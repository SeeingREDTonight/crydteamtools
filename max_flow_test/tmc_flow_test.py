# TMC Flow Test — Adaptive max-volumetric-flow detection for extruders
# credits:
#   Steven (Fragmon) — Crydteam
#   YouTube: https://www.youtube.com/@crydteamprinting
#
# License: GPLv3

import logging
import math
import os
import time
import json


def _pstdev(values):
    """Population standard deviation. Local implementation to avoid
    importing Python's stdlib `statistics` module — Kalico ships an
    extras/statistics.py plugin that shadows it when this file is
    imported from klipper/klippy/extras."""
    n = len(values)
    if n == 0:
        return 0.0
    mean = sum(values) / n
    variance = sum((x - mean) ** 2 for x in values) / n
    return math.sqrt(variance)


def _sanitize_csv(name):
    """Make a thermistor object name CSV-column-safe: ASCII alnum
    plus underscore. 'temperature_sensor chamber' -> 'temperature_sensor_chamber'.
    """
    out = []
    for ch in name:
        if ch.isalnum() or ch == '_':
            out.append(ch)
        else:
            out.append('_')
    return ''.join(out) or 'sensor'


SAMPLE_INTERVAL = 0.02    # 50 Hz polling — stall clicks last ~10-30 ms,
                          # at 20 Hz most of them fell between samples
THERMAL_EVERY_N = 12      # thermal snapshot every Nth SG sample (~4 Hz)
STARTUP_TRIM_FRACTION = 0.10  # drop first 10 % of each run's SG samples
                              # (spin-up transient reads near 0 on SG2
                              # and poisons sg_min / dip detection)
TAIL_TRIM_FRACTION = 0.05     # drop last 5 % of each run: during the
                              # deceleration ramp SG is equally invalid
                              # and produced ~2 false dip/peak samples
                              # per run (CSV 2026-08-09_09-10-42)
MIN_HOTEND_TEMP = 180.0
MODULE_NAME = "TMC Flow Test"
MODULE_VERSION = "1.2.2"
SG_MIN_INFORMATIVE = 50   # below this SG value, readings are noise


# ─── Driver-specific trigger profiles ──────────────────────────────────
# Each TMC chip family has its own SG noise characteristics, scale, and
# CoolStep behaviour. Rather than embedding "if self.sg2_driver / else"
# branches throughout the trigger code, all tunable thresholds are
# centralized here in a profile per driver. To adjust sensitivity for
# one driver family without affecting the others, edit only its profile.
#
# The base profile contains the values that have been validated against
# real TMC5160 test runs (CSVs from 11-37, 12-40, 12-56, 13-20). Other
# drivers START as identical copies and can be tuned independently from
# real data.

class TriggerProfile:
    """Threshold values used by the slip-detection triggers.

    Inherit and override only the fields that need tuning for a
    specific driver family. Values here are the validated TMC5160
    defaults — DO NOT change them without re-running the TMC5160
    regression suite.
    """

    # ─── _check_cv_spike thresholds ────────────────────────────────
    # Pattern (a) "high-variance trip" — last CV above this absolute
    # value triggers immediately, regardless of baseline.
    CV_HIGH_VARIANCE = 10.0

    # Pattern (b) "ratio jump from immediate baseline"
    CV_JUMP_RATIO_COARSE = 2.5
    CV_JUMP_RATIO_BISECT = 2.0
    CV_JUMP_MIN_COARSE = 5.0
    CV_JUMP_MIN_BISECT = 4.0

    # Pattern (c) "rising trend" — three consecutive ≥1.3× steps
    CV_RISING_RATIO = 1.3
    CV_RISING_MIN_PRIOR = 2.5
    CV_RISING_MIN_LAST_COARSE = 5.0
    CV_RISING_MIN_LAST_BISECT = 4.0

    # Pattern (d) "vs coarse-phase baseline" (bisection only)
    CV_VS_COARSE_RATIO = 2.0
    CV_VS_COARSE_MIN = 5.0

    # Single-step low-baseline jump check
    CV_LOWBASE_RATIO = 1.5

    # ─── _check_iqr_spread thresholds ──────────────────────────────
    # Pattern (a) "ratio vs immediate prior steps"
    IQR_RATIO_COARSE = 3.0
    IQR_RATIO_BISECT = 1.7
    IQR_RATIO_MIN_ABS = 12   # don't fire on tiny absolute IQRs
    # Pattern (a) CV cross-check: when True, pattern (a) only fires
    # if CV ALSO confirms elevated noise relative to prior steps.
    # Helps drivers (e.g. SG4 with active CoolStep) where CS regulation
    # widens within-step IQR without affecting run-to-run reproducibility
    # (CV stays low). TMC5160 was validated with this OFF — leave it
    # off for SG2 unless real data shows false positives there too.
    IQR_RATIO_REQUIRE_CV = False
    IQR_RATIO_CV_FLOOR_COARSE = 3.0   # absolute CV floor, coarse phase
    IQR_RATIO_CV_FLOOR_BISECT = 2.5   # absolute CV floor, bisection

    # Pattern (b) "vs coarse-phase median"
    IQR_VS_COARSE_RATIO = 2.5
    IQR_VS_COARSE_MIN_ABS = 18
    IQR_VS_COARSE_REQUIRE_CV = True   # CV cross-check on by default
    IQR_VS_COARSE_CV_FLOOR = 3.0
    IQR_VS_COARSE_CV_RATIO = 1.5

    # Pattern (c) "absolute" — trigger if IQR ≥ this in bisection
    IQR_ABSOLUTE_TRIGGER = 25

    # ─── _check_sg_max_spike thresholds ────────────────────────────
    SG_MAX_RATIO_TO_MEDIAN = 3.0
    SG_MAX_RATIO_TO_COARSE = 1.3
    SG_MAX_ABS_GAP = 200            # absolute gap floor (sg2_driver)
    SG_MAX_ABS_GAP_SG4 = 150        # absolute gap floor (sg4 chips)
    SG_MAX_BIG_RATIO = 4.0          # alt path — extreme ratio
    SG_MAX_BIG_GAP = 300            # alt path — extreme gap

    # ─── _check_sg_max_spike COARSE-phase thresholds ──────────────
    # The coarse-phase path uses a separate, "gap-jump" criterion.
    # It fires when (sg_max - sg_median) makes a big jump above the
    # baseline gap from earlier coarse steps. This catches stick-slip
    # stalls where the medians stay compact (so IQR/CV/median triggers
    # all miss it) but sg_max suddenly records repeated decoupling
    # spikes inside the runs. Conservative defaults to avoid early
    # noise-driven false positives.
    COARSE_GAP_JUMP_RATIO = 2.0      # current gap >= 2x prior baseline gap
    COARSE_GAP_JUMP_ABS_FLOOR = 350  # current gap must clear this absolute floor
    COARSE_GAP_JUMP_PREV_FRACTION = 0.7  # prev step's gap must also be ≥ 0.7×floor
    # Inverse direction (rare case where SG rises with load)
    SG_MIN_RATIO_TO_MEDIAN = 3.0
    SG_MIN_ABS_GAP = 80

    # ─── _check_dip_collapse thresholds (SG2 only) ─────────────────
    # A stick-slip stall drives SG_RESULT toward 0 for a few samples —
    # far too few to move median/IQR/CV. After the spin-up trim,
    # samples below DIP_FRACTION × run-median count as collapse dips.
    # The trigger fires when the step's dip RATE clears both an
    # absolute floor and a multiple of the prior-step baseline.
    DIP_FRACTION = 0.25             # dip = sample < this × run median
    DIP_MIN_COUNT = 3               # min dips in the step to consider
    DIP_RATE_MIN = 0.01             # ≥ 1 % of the step's samples
    DIP_BASE_RATIO = 3.0            # and ≥ 3× prior-step baseline rate

    # ─── _check_peak_unload thresholds (SG2 only) ──────────────────
    # The EXTRUDER counterpart of the dip: when the gear grinds
    # through the filament, the motor is UNLOADED and SG spikes UP
    # toward its no-load value while the median stays smooth
    # (validated: CSV 2026-08-09_08-37-15, full grind at flow=50 with
    # sg_max 1.6× median but compact median/IQR/CV). Samples above
    # (1 + PEAK_FRACTION) × run-median count as unload peaks —
    # counting samples instead of using sg_max makes a single stray
    # spike harmless.
    # Thresholds validated against two real runs: full grind at
    # flow=50 (08-37-15, ~1 % estimated peak rate) must fire, while
    # 4 stray samples in 1241 during clean operation at flow=40
    # (09-10-42, decel artefacts) must NOT.
    PEAK_FRACTION = 0.45            # peak = sample > 1.45 × run median
    PEAK_MIN_COUNT = 8              # min peaks in the step to consider
    PEAK_RATE_MIN = 0.006           # ≥ 0.6 % of the step's samples
    PEAK_BASE_RATIO = 3.0           # and ≥ 3× prior-step baseline rate

    # ─── _check_thermal_runaway thresholds (all drivers) ───────────
    # Continuous grinding keeps the motor load HIGH (the gear digs into
    # stationary filament), so every SG trigger above can stay quiet.
    # The give-away is thermal: the melt zone stops absorbing the
    # commanded flow, so the hotend climbs above target even though the
    # heater was struggling one step earlier.
    # Validated against 8 real runs: clean ones peak at 0.7-4.9 °C above
    # target, the grinding run reached 17.9 °C.
    TEMP_RUNAWAY_ABS = 5.0          # °C above target, absolute floor
    TEMP_RUNAWAY_JUMP = 4.0         # AND °C above the earlier steps'
                                    # own offset (protects setups that
                                    # simply run warm)
    # Heater pinned at this duty for a whole step = the hotend is at its
    # thermal limit; higher flows are not sustainable in a real print.
    PWM_SATURATED = 0.95

    # ─── _check_heater_collapse thresholds (all drivers) ──────────
    # Melting more plastic takes more heat, so heater duty must RISE
    # with flow. If the commanded flow goes up but duty falls sharply,
    # the filament has stopped arriving (stall, grind or clog). This
    # also catches a motor that loses sync at move start and then sits
    # still in a spinning field — SG keeps tracking the commanded step
    # rate there, so every SG trigger stays blind, and a strong heater
    # PID can hold the overshoot under TEMP_RUNAWAY_ABS.
    # Observed: TMC5160 run (2026-11-12) with a stall at move start:
    # duty 0.378 → 0.268 (−29 %) at the stall step, vs +5..10 % per
    # healthy step; repeat run −31 %.
    # PLA strip (2026-11-12 19:39) happened mid-step: −11 % / −0.040 at
    # the strip step. Biggest drop across three clean runs: −5 % /
    # −0.019. So 10 % / 0.03 catches the strip step with ~2x margin.
    PWM_COLLAPSE_FRAC = 0.10        # duty < (1 − this) × lower-flow step
    PWM_COLLAPSE_MIN_ABS = 0.03     # AND an absolute drop of ≥ 3 % duty

    # ─── _check_run_outlier thresholds ─────────────────────────────
    OUTLIER_MAD_RATIO = 4.0         # deviation ≥ this × MAD
    OUTLIER_MIN_REL = 0.08          # AND ≥ 8 % of median

    # ─── _is_borderline thresholds ─────────────────────────────────
    BORDER_CV_LOW = 4.0
    BORDER_CV_HIGH = 7.0
    BORDER_CV_RATIO = 1.5
    BORDER_IQR_LOW = 15
    BORDER_IQR_HIGH = 25
    BORDER_IQR_RATIO = 1.7
    BORDER_IQR_CV_FLOOR = 2.0       # CV must be at least this for IQR-only borderline

    # ─── SG-level / step-jump basics ───────────────────────────────
    # Used by _sg_min_informative and _sg_jump_threshold.
    # Note: sg2_driver returns -1 (no gating) for SG_MIN_INFORMATIVE,
    # while SG4 returns the SG_MIN_INFORMATIVE module constant.
    SG_JUMP_THRESHOLD = 5           # SG2: small jumps significant

    # ─── Plateau trigger threshold ─────────────────────────────────
    # The plateau trigger fires when the cumulative SG-load over 2
    # steps falls short of expected_2step × this fraction. Lower
    # values (e.g. 0.3) mean the trigger is more permissive — only
    # really flat sections fire. Higher values (e.g. 0.6) mean any
    # mild deceleration of the SG decline fires.
    #
    # TMC5160 baseline of 0.5 was tuned with real Sherpa data. TMC2240
    # SG2 has a steeper SG-vs-load saturation curve (each delta is
    # half the previous), so 0.5 fires very early on the natural
    # saturation slope. Lower this fraction for TMC2240.
    PLATEAU_RATIO = 0.5             # cumulative-load < this × expected_2step

    # ─── Plateau saturation-skip threshold ─────────────────────────
    # The plateau trigger compares the current step's SG-delta to the
    # average of the prior 3 deltas. If any prior step had SG values
    # close to the 1023 saturation ceiling, those step's delta is
    # artificially inflated (clipping artefact) — and the plateau
    # trigger then thinks "load should have dropped MUCH more than
    # it did" and falsely fires on naturally saturating curves.
    #
    # Skip plateau evaluation entirely if any prior SG median exceeds
    # this threshold. The validated TMC5160 11-37-30 baseline has
    # max prior SG ~536 (well below 700). A test starting at flow=10
    # with SGT≥15 typically sees SG=850-1023 on the first step —
    # that's the case we need to skip.
    PLATEAU_SATURATION_SKIP = 700   # skip plateau if prior SG > this

    # ─── Warmup-skip threshold ─────────────────────────────────────
    # Per-step, the first repetition often shows different SG behaviour
    # than the others (motor transitions from cold-stop to extruding,
    # filament path settles, etc.). If run 1's SG average deviates from
    # the rest by MORE than this fraction, it's flagged as warmup and
    # excluded from the median/IQR/CV stats. The TMC5160 baseline
    # (10 %) was tuned with real Sherpa Mini data; TMC2240 needs a
    # more aggressive threshold because its run-to-run drift is
    # systematically higher (typically 3-6 %) — without dropping
    # warmup, the run-outlier trigger fires repeatedly on these
    # warmup runs and aborts the test far below real slip.
    WARMUP_DRIFT_THRESHOLD = 0.10   # 10 % default (TMC5160 baseline)

    # ─── Auto-SGT-tuning ranges (SG2 drivers) ──────────────────────
    # Used by the optional pre-test SGT auto-tuning phase. Values
    # apply to drivers using the signed `sgt` field (TMC5160, TMC2130,
    # TMC2240 in SG2 mode). TMC2209 (SG4) uses different field/range
    # and is not auto-tuned in this version.
    #
    # Probe strategy: extrude AT the user's START flow (the same flow
    # the test will begin coarse sweeping from). Calibrating SGT for
    # exactly this load gives the most representative SG-vs-load
    # baseline — which matters because SGT effects are flow-dependent
    # and a probe-derived SGT only matches the test if the probe flow
    # matches the start flow.
    #
    # Targets are evaluated at the start-flow load:
    #   - SGT_LOW_TARGET_MIN: SG should be at least this for adequate
    #                         dynamic range. Below → bump SGT up.
    #   - SGT_LOW_TARGET_MAX: SG should NOT exceed this. Slightly
    #                         below the 1023 saturation ceiling so a
    #                         truly-saturated reading is detected.
    SGT_LOW_TARGET_MIN = 600        # SG at start_flow ≥ this
    SGT_LOW_TARGET_MAX = 1022       # SG at start_flow ≤ this (1023=saturated)
    SGT_RANGE_MIN = -64             # absolute SGT lower bound (SG2)
    SGT_RANGE_MAX = 63              # absolute SGT upper bound (SG2)
    SGT_AUTOTUNE_PROBE_DURATION = 5.0  # seconds per probe (matches test)
    SGT_AUTOTUNE_PROBE_REPEATS = 5     # repetitions per probe (statistical)
    SGT_AUTOTUNE_MAX_ITERATIONS = 5    # safety cap


class TMC5160Profile(TriggerProfile):
    """TMC5160 / TMC2130 / TMC2660 — SG2 + SpreadCycle.

    These are the validated production values. Do not change without
    running the regression suite against the historical CSVs.
    """
    pass  # uses base defaults (= validated 5160 values)


class TMC2240Profile(TriggerProfile):
    """TMC2240 — SG2 + SpreadCycle (high-torque path).

    The plugin runs TMC2240 in the same SG2/SpreadCycle path as the
    TMC5160, because the SG4/StealthChop path delivers ~50 % less peak
    torque (confirmed empirically: same Sherpa Mini reaches 110 mm³/s
    on TMC5160 SG2 but only 60 mm³/s on TMC2240 SG4).

    The TMC2240 SG2 hardware path produces nearly identical SG signal
    characteristics to the TMC5160 — but with one systematic
    difference: the first repetition of every measurement step shows
    a 3-6 % SG drift compared to runs 2-5. This warmup pattern is
    consistent across all flow steps (confirmed with CSV 21-02-45).
    The default TMC5160 warmup-skip threshold (10 %) doesn't catch
    these drifts in the coarse phase but does catch them in
    bisection / verify, where the run-outlier trigger then
    erroneously interprets the warmup run as slip — producing a
    max-flow result that's roughly half the real hardware capability.

    Fix: lower the warmup-skip threshold so run 1 is reliably
    excluded throughout the test. All other thresholds inherit from
    the validated TMC5160 baseline.

    Second fix: lower the plateau-trigger threshold. The TMC2240 SG2
    saturation curve is steeper than the TMC5160 — observed ratios
    of 0.23 between consecutive deltas during natural saturation. The
    TMC5160 default of 0.5 fires far too early on this curve. Setting
    PLATEAU_RATIO to 0.2 makes the trigger fire only when the trend
    truly flattens (which is when slip actually begins).
    """
    WARMUP_DRIFT_THRESHOLD = 0.04   # 4 % — catches the systematic drift
    PLATEAU_RATIO = 0.2             # steeper saturation than TMC5160

    # Third fix (validated against CSV 2026-08-08_10-03-53): audible
    # stalls at flow=70 produced a real CV jump (3.4 % = 3.3× the
    # 1.0 % baseline) that stayed under the TMC5160 floor of 5.0.
    # TMC2240 SG2 runs are tighter (clean CV ~0.7-1.3 %), so a lower
    # absolute floor is safe for the coarse phase — but paired with a
    # stricter ratio: replaying CSV 2026-05-02_19-37-17 showed a 2.9×
    # jump to CV 3.0 % at flow=80 on a run whose trend only flattened
    # at 95, and firing there would truncate the result.
    CV_JUMP_MIN_COARSE = 3.0
    CV_JUMP_RATIO_COARSE = 3.0


class TMC2209Profile(TriggerProfile):
    """TMC2209 — SG4 in StealthChop or SpreadCycle.

    Empirical findings on TMC2209 (Sherpa Mini, validated 2026-04-30):

    Unlike the SG2 drivers (TMC5160 / TMC2240), TMC2209 SG4_RESULT
    behaves very differently and the SG-median trigger logic that
    works on SG2 does NOT apply:

    1. SG_RESULT scale is 0-510 (not 0-1023). Trinamic spec says
       "higher = lower load" same as SG2, but in practice SG_RESULT
       on TMC2209 only produces meaningful values above a hardware-
       specific minimum velocity — below that it sticks in a "bias
       region" of 6-22 regardless of actual load.

    2. SG-median MAGNITUDE is unreliable as a slip indicator. On
       validated TMC2209 hardware, SG actually INCREASED with load
       (from 6 at low flow to 86 at near-stall), the opposite of
       what Trinamic's spec describes. This is why Klipper's
       sensorless homing on TMC2209 only uses the binary DIAG-pin
       trigger (SG < 2*SGTHRS), not the SG_RESULT magnitude itself.

    3. SG run-to-run CV (variance) IS a reliable slip indicator.
       At validated stall point (flow=140, motor physically slipped):
         CV jumped from 2.4% (no slip) to 23.3% (full stall).
       That's a 10x increase — easily detectable.

    Strategy: trigger primarily on CV-spike and IQR-widening. Disable
    SG-median-based triggers (snap-back, plateau, max-spike) since
    they assume the SG2 "smooth load curve" model that doesn't apply.
    """
    # ─── SG-median triggers DISABLED ──────────────────────────────
    # Set thresholds extremely high so they never fire. SG2 logic
    # (snap-back / plateau / max-spike) doesn't apply to SG4.
    SG_JUMP_THRESHOLD = 99999       # disable sudden-jump trigger
    SG_MAX_SPIKE_RATIO = 99.0       # disable max-spike (decoupling) trigger
    SG_MAX_SPIKE_RATIO_BISECT = 99.0
    COARSE_GAP_JUMP_RATIO = 99.0     # disable coarse gap-jump for SG4
    COARSE_GAP_JUMP_ABS_FLOOR = 99999
    COARSE_GAP_JUMP_PREV_FRACTION = 99.0
    PLATEAU_RATIO = 0.0             # 0 → never fire plateau
    PLATEAU_SATURATION_SKIP = 0     # n/a since plateau disabled

    # ─── CV-based triggers TIGHTENED ──────────────────────────────
    # Empirically: CV ~2-7 % during clean operation, jumps to 20-30 %
    # at slip. Set CV_HIGH_VARIANCE high enough to ignore noise but
    # low enough to catch the slip jump.
    CV_HIGH_VARIANCE = 12.0         # absolute CV trigger (vs 5.0 on SG2)
    CV_JUMP_RATIO_BISECT = 3.0      # 3x baseline jump in bisection
    CV_JUMP_RATIO_COARSE = 4.0      # 4x for coarse (noise-tolerant)
    CV_RISING_BISECT = 8.0          # rising CV trend trigger
    CV_RISING_COARSE = 12.0
    CV_VS_COARSE_BISECT = 4.0       # 4x coarse-baseline CV
    CV_VS_COARSE_COARSE = 6.0
    CV_LOWBASE_RATIO = 5.0          # only trigger if base CV > 1 %
    CV_FLOOR = 1.0

    # ─── IQR-based triggers TIGHTENED ─────────────────────────────
    # IQR also widens dramatically at slip. On the 0-510 scale the
    # absolute IQR threshold needs to be lower than for SG2.
    IQR_RATIO_COARSE = 5.0
    IQR_RATIO_BISECT = 4.0
    IQR_VS_COARSE_BISECT = 6.0
    IQR_ABSOLUTE_TRIGGER = 50       # raw SG units (vs 25 on SG2)
    REQUIRE_CV_CROSS_CHECK = True   # IQR alone needs CV confirmation

    # ─── Border / sanity ──────────────────────────────────────────
    BORDER_CV_LOW = 6.0             # CV between 6-12% = borderline
    BORDER_CV_HIGH = 12.0
    BORDER_IQR_LOW = 30
    BORDER_IQR_HIGH = 50

    # ─── Dip-collapse / peak-unload triggers DISABLED ─────────────
    # SG4 sticks in a low "bias region" at low velocity and its
    # magnitude doesn't follow the SG2 collapse/unload model — low or
    # high samples are normal, not stalls.
    DIP_MIN_COUNT = 999999
    PEAK_MIN_COUNT = 999999

    # ─── Outlier detection (single-run slip in 5 reps) ────────────
    OUTLIER_MAD_RATIO = 4.0
    OUTLIER_MIN_REL = 0.15          # relaxed for SG4 noise (vs 0.08 on SG2)

    # ─── Warmup-skip ──────────────────────────────────────────────
    WARMUP_DRIFT_THRESHOLD = 0.15   # SG4 has higher first-run drift


def get_trigger_profile(driver_type):
    """Return the appropriate TriggerProfile subclass for a given driver."""
    if driver_type == 'tmc2240':
        return TMC2240Profile
    if driver_type == 'tmc2209':
        return TMC2209Profile
    # Default: TMC5160-style SG2 profile (covers tmc5160, tmc2130, tmc2660,
    # and any unknown driver — safest fallback).
    return TMC5160Profile


class _HardStop(Exception):
    """Raised after a hard filament failure (heater collapse or BDWidth
    feed): the filament is stripped, stalled or clogged, so bisecting
    further would only grind it. The run has already been finalised
    and reported when this is raised."""


class TMCFlowTest:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.reactor = self.printer.get_reactor()
        self.gcode = self.printer.lookup_object('gcode')

        self.stepper_name = config.get('extruder_stepper', 'extruder')
        self.filament_diameter = config.getfloat(
            'filament_diameter', 1.75, above=0.)
        self.melt_zone_length = config.getfloat(
            'melt_zone_length', 42.0, above=0.)
        self.min_hotend_temp = config.getfloat(
            'min_hotend_temp', MIN_HOTEND_TEMP, above=0.)

        # Fixed low reference flow (mm³/s) measured BEFORE the coarse
        # sweep whenever START is well above it. The two reference
        # steps give every baseline-driven trigger (CV, dips, peaks,
        # IQR, gap) clean slip-free data to compare against and anchor
        # the SGT auto-tune — without them a high START would leave
        # the early sweep steps trigger-blind or with slip-polluted
        # baselines. 0 disables the reference phase.
        self.reference_flow = config.getfloat(
            'reference_flow', 15.0, minval=0.)

        # Part-cooling fan speed during the test (0-100 %). Set in
        # config or override via `FAN_SPEED=` on TMC_FLOW_FIND_MAX.
        # Once the test starts the value is locked — printer-cooling
        # fans actually have non-trivial impact on max flow (heater
        # has to do ~10-20 % more work to maintain temperature when
        # fan is at 100 % vs off), so consistency matters more than
        # flexibility here.
        self.test_fan_speed = config.getfloat(
            'test_fan_speed', 0.0, minval=0.0, maxval=100.0)
        # Optional: explicitly name the fan to control. Default is
        # the printer's primary part-cooling fan via M106.
        self.fan_object_name = config.get('fan_object_name', '').strip()

        # Optional list of additional thermistors to log alongside the
        # extruder's. Useful for chamber, heatbreak, ambient, etc.
        # Comma-separated list of object names like:
        #   extra_thermistors: chamber, heater_generic heatbreak
        # Each name must be the printer-config-section header (see
        # `[temperature_sensor name]` and similar). Names that don't
        # resolve are silently skipped at test start.
        extras_raw = config.get('extra_thermistors', '').strip()
        if extras_raw:
            self.extra_thermistors = [
                n.strip() for n in extras_raw.split(',') if n.strip()]
        else:
            self.extra_thermistors = []

        # Optional chamber / filter fans to run for the whole test
        # (air filtering only — they do not change the measurement).
        # Comma-separated names; the type is detected at run start:
        #   [fan_generic X]      → SET_FAN_SPEED FAN=X SPEED=1
        #   [temperature_fan X]  → target set to its min_temp so it runs
        # Previous speed / target is restored when the test ends.
        cf_raw = config.get('chamber_fans', '').strip()
        self.chamber_fans = [n.strip() for n in cf_raw.split(',')
                             if n.strip()]

        config_dir = os.path.expanduser('~/printer_data/config')
        if not os.path.isdir(config_dir):
            config_dir = os.path.expanduser('~')
        default_dir = os.path.join(config_dir, 'Flowtest')
        # expanduser: the settings template documents a '~/...' path, and
        # without this a user-supplied value would create a literal '~'
        # directory next to klippy instead of writing into the home dir.
        self.output_dir = os.path.expanduser(
            config.get('output_dir', default_dir))

        self.filament_area = math.pi * (self.filament_diameter / 2.0) ** 2

        # Driver detection state
        self.tmc = None
        self.driver_type = None   # 'tmc2240', 'tmc2209', 'tmc5160', etc.
        self.is_2240 = False
        self.is_2209 = False
        self.is_5160 = False
        # SG2 family (TMC5160, TMC2130, TMC2660) — StallGuard2 in SpreadCycle,
        # SG_RESULT in DRV_STATUS bits 0-9, threshold field 'sgt' (signed).
        self.sg2_driver = False
        self.sg4_available = False  # SG4_RESULT register (TMC2240)
        # Trigger threshold profile — replaced in _lookup_tmc with the
        # driver-appropriate subclass. Default = TMC5160 (validated).
        self.profile = TMC5160Profile

        # Sample buffers
        self.samples_sg = []
        self.samples_time = []
        self.sampling_active = False
        self.sample_timer = None
        self.sample_start_time = 0.0

        # Main command
        self.gcode.register_command(
            'TMC_FLOW_FIND_MAX', self.cmd_TMC_FLOW_FIND_MAX,
            desc='Find max volumetric flow rate via StallGuard slip '
                 'detection (Coarse → Bisection → Verification)')
        self.gcode.register_command(
            'TMC_FLOW_STATUS', self.cmd_TMC_FLOW_STATUS,
            desc='Show current TMC StallGuard diagnostic values')
        self.gcode.register_command(
            'TMC_FLOW_GUI', self.cmd_TMC_FLOW_GUI,
            desc='Write the beginner-friendly control panel (HTML) with '
                 'live driver/config values into the output directory')
        self.gcode.register_command(
            'TMC_FLOW_SGT_CALIBRATE', self.cmd_TMC_FLOW_SGT_CALIBRATE,
            desc='Auto-tune the StallGuard SGT threshold at a low '
                 'reference flow and stage the result for SAVE_CONFIG '
                 '(persisted in printer.cfg)')
        self.gcode.register_command(
            'TMC_FLOW_TEST_SG_VARIANTS',
            self.cmd_TMC_FLOW_TEST_SG_VARIANTS,
            desc='TMC2209-only diagnostic: probe SG_RESULT in both '
                 'StealthChop and SpreadCycle to determine empirically '
                 'which chopper mode produces a usable signal on this '
                 'specific hardware')

    # ─── TMC driver lookup ──────────────────────────────────────────

    def _lookup_tmc(self):
        if self.tmc is not None:
            return
        candidates = ['tmc2240', 'tmc5160', 'tmc2209', 'tmc2226',
                      'tmc2130', 'tmc2208', 'tmc2660']
        for drv in candidates:
            obj_name = '%s %s' % (drv, self.stepper_name)
            tmc = self.printer.lookup_object(obj_name, None)
            if tmc is not None:
                self.tmc = tmc
                self.driver_type = drv
                self.is_2240 = (drv == 'tmc2240')
                self.is_2209 = (drv == 'tmc2209')
                self.is_5160 = (drv == 'tmc5160')
                # TMC2240 is treated as an SG2 driver (SpreadCycle path).
                # Its SG4/StealthChop path delivers significantly less
                # torque at high speeds — confirmed empirically: the same
                # Sherpa Mini extruder reaches ~110 mm³/s on TMC5160 SG2
                # but only ~60 mm³/s on TMC2240 SG4. Since this plugin
                # exists to find MAX flow, we always use the high-torque
                # SG2 path on TMC2240.
                self.sg2_driver = drv in ('tmc5160', 'tmc2130',
                                          'tmc2660', 'tmc2240')
                # Pick the trigger threshold profile for this driver
                # family. Each profile is independently tunable; the
                # TMC5160 profile is the validated baseline.
                self.profile = get_trigger_profile(drv)
                # Check for SG4_RESULT register (TMC2240 only)
                try:
                    self.sg4_available = (
                        self.is_2240
                        and 'SG4_RESULT' in self.tmc.mcu_tmc.name_to_reg)
                except AttributeError:
                    self.sg4_available = False
                logging.info(
                    "tmc_flow_test: using %s for stepper '%s' "
                    "(SG4=%s, is_2209=%s, sg2=%s, profile=%s)",
                    drv, self.stepper_name, self.sg4_available,
                    self.is_2209, self.sg2_driver,
                    self.profile.__name__)
                return
        raise self.gcode.error(
            "tmc_flow_test: no TMC driver found for stepper '%s'."
            % self.stepper_name)

    # ─── Driver-specific field accessors ────────────────────────────

    def _get_sg_threshold_field_name(self):
        """Return the correct SG-threshold field name for this driver.

        TMC5160 / TMC2130 / TMC2660 / TMC2240 (SG2 path) use 'sgt'.
        TMC2209 uses 'sgthrs'.
        """
        if self.sg2_driver:
            # TMC5160 / TMC2130 / TMC2660 / TMC2240 — StallGuard2 threshold
            # (signed -64..63)
            return 'sgt'
        # TMC2209, TMC2226, TMC2208 use 'sgthrs'
        return 'sgthrs'

    def _get_sg_label(self):
        """Human-readable label for the SG signal.

        TMC2240 in SG2 mode uses SG_RESULT in DRV_STATUS just like the
        TMC5160 — the SG4_RESULT register is unused on the SG2 path.
        """
        return 'SG_RESULT'

    # ─── TMC config validation ──────────────────────────────────────

    def _check_tmc_config(self):
        """Verify TMC driver is configured correctly for the SG-based test.

        Returns (problems, infos).
        """
        problems = []
        infos = []
        if self.tmc is None:
            return ([('tmc', None, 'No TMC driver found')], infos)

        def get(name):
            try:
                return self.tmc.fields.get_field(name)
            except (KeyError, AttributeError):
                return None

        tpwmthrs = get('tpwmthrs')
        tcoolthrs = get('tcoolthrs')
        semin = get('semin')
        en_pwm_mode = get('en_pwm_mode')
        en_spread_cycle = get('en_spreadcycle')
        sg_thrs_field = self._get_sg_threshold_field_name()
        sg_thrs_val = get(sg_thrs_field)

        # ─── Chopper-mode check (driver-specific) ───
        # SG4 family (TMC2209 only now): StallGuard4 needs StealthChop ON
        # SG2 family (TMC5160, TMC2130, TMC2660, TMC2240):
        #   StallGuard2 needs SpreadCycle
        if self.is_2209:
            stealthchop_active = True
            stealthchop_indicator = None
            if en_spread_cycle is not None:
                stealthchop_active = (en_spread_cycle == 0)
                stealthchop_indicator = (
                    'en_spreadcycle', en_spread_cycle, 'should be 0')

            if not stealthchop_active and stealthchop_indicator:
                problems.append(
                    (stealthchop_indicator[0], stealthchop_indicator[1],
                     'StealthChop is not active. StallGuard4 needs '
                     'StealthChop ON.\n'
                     'Add to your [%s extruder] section:\n'
                     '  stealthchop_threshold: 999999' % self.driver_type))
            elif (tpwmthrs is not None and tpwmthrs > 0
                  and tpwmthrs < 0x10000):
                # Mid-range tpwmthrs would switch to SpreadCycle at higher
                # speeds.
                problems.append(
                    ('tpwmthrs', tpwmthrs,
                     'StealthChop only active below a velocity threshold '
                     '(tpwmthrs=%d). At higher flows the driver switches '
                     'to SpreadCycle and breaks StallGuard4.\n'
                     'Set:\n'
                     '  stealthchop_threshold: 999999' % tpwmthrs))
        elif self.sg2_driver:
            # TMC5160 / TMC2130 / TMC2660 / TMC2240 — SG2 requires SpreadCycle.
            #
            # Klipper's tmc.py writes tpwmthrs = 0xFFFFF (= 1048575) as the
            # default when stealthchop_threshold is NOT set in the config —
            # together with en_pwm_mode=0 this is the "SpreadCycle at all
            # speeds" state we actually want for SG2. So:
            #
            #   • en_pwm_mode is the authoritative signal: 0 = SpreadCycle.
            #   • Only flag tpwmthrs when en_pwm_mode is 1 AND tpwmthrs is in
            #     the mid-range (0 < tpwmthrs < 0xFFFFF), which would make
            #     the driver switch chopper mode based on velocity.
            #   • tpwmthrs == 0xFFFFF means "never enter StealthChop" and is
            #     the Klipper default → not a problem.
            stealthchop_default_tpwmthrs = 0xFFFFF  # 1048575
            if en_pwm_mode is not None and en_pwm_mode == 1:
                # Build TMC2240-specific advice if applicable
                if self.is_2240:
                    extra = (
                        '\nNote: this plugin uses TMC2240 in '
                        'SG2/SpreadCycle mode for max torque (~50 %% '
                        'higher peak flow than SG4/StealthChop). If you '
                        'use sensorless homing, set the chopper mode '
                        'temporarily for the test only.')
                else:
                    extra = ''
                if (tpwmthrs is not None and 0 < tpwmthrs
                        < stealthchop_default_tpwmthrs):
                    problems.append(
                        ('tpwmthrs', tpwmthrs,
                         '%s StallGuard2 needs SpreadCycle at all speeds, '
                         'but tpwmthrs=%d enables StealthChop below that '
                         'velocity.\n'
                         'Remove stealthchop_threshold from your [%s '
                         'extruder] section (or set it to 0).%s'
                         % (self.driver_type.upper(), tpwmthrs,
                            self.driver_type, extra)))
                else:
                    problems.append(
                        ('en_pwm_mode', en_pwm_mode,
                         'StealthChop is active. %s StallGuard2 only works '
                         'in SpreadCycle mode.\n'
                         'Remove stealthchop_threshold from your [%s '
                         'extruder] section (do not set it to 0 — remove '
                         'the line).%s'
                         % (self.driver_type.upper(), self.driver_type,
                            extra)))
            # else: en_pwm_mode == 0 → SpreadCycle is active. tpwmthrs is
            # irrelevant in this case (the driver never enters StealthChop),
            # whether it's 0 or 0xFFFFF.

        # ─── tcoolthrs check (StallGuard gate) ───
        if tcoolthrs == 0:
            problems.append(
                ('tcoolthrs', tcoolthrs,
                 'StallGuard reading disabled (tcoolthrs=0).\n'
                 'Add to your [%s extruder] section:\n'
                 '  coolstep_threshold: 0.5' % self.driver_type))

        # ─── SG threshold check ───
        # Only enforced for TMC2209 (SG4 path). For SG2 drivers
        # (TMC5160/2130/2660/2240), sgt is signed (-64..63) and
        # SG_RESULT can be read regardless of sgt; we don't depend on
        # the hardware stop trigger.
        if self.is_2209 and sg_thrs_val == 0:
            problems.append(
                (sg_thrs_field, sg_thrs_val,
                 'SGTHRS is 0. StallGuard trigger inactive.\n'
                 'Add to your [tmc2209 extruder] section:\n'
                 '  driver_SGTHRS: 100'))

        # ─── CoolStep informational check ───
        # The plugin runs SG-only triggers. CoolStep can still be enabled
        # in the user's config (it doesn't break the test) but it may
        # dampen SG signal during high-load steps and produce a less
        # conservative max-flow result.
        if semin is not None and semin > 0:
            infos.append(
                "Note: CoolStep is active (driver_SEMIN=%d). The test "
                "will run with your configured SEMIN, but CoolStep may "
                "dampen SG signal during high-load steps. For the most "
                "conservative max-flow result, set driver_SEMIN: 0 "
                "temporarily." % semin)
        else:
            infos.append(
                "CoolStep is disabled (driver_SEMIN=0). Motor runs at "
                "constant IRUN — recommended for clean SG signal.")

        # Driver info
        if self.is_2209:
            infos.append(
                "Driver: TMC2209 detected (uses SG_RESULT, sgthrs).")
        elif self.is_2240:
            infos.append(
                "Driver: TMC2240 detected — using SG2/SpreadCycle path "
                "(SG_RESULT via DRV_STATUS, sgt threshold) for max "
                "torque. SG4/StealthChop path is intentionally bypassed "
                "since it delivers ~50%% less peak flow.")
        elif self.sg2_driver:
            infos.append(
                "Driver: %s detected (StallGuard2 in SpreadCycle, "
                "SG_RESULT via DRV_STATUS, sgt threshold)."
                % self.driver_type.upper())

        return (problems, infos)

    # ─── SG sampling ─────────────────────────────────────────────────

    def _read_sg(self):
        """Read StallGuard value directly from the driver register.

        TMC2209: SG_RESULT register (dedicated SG4 register)
        TMC5160 / TMC2130 / TMC2660 / TMC2240 (SG2 path):
            DRV_STATUS bits 0-9 (sg_result field)
        Other drivers: fallback via get_status()

        Note: TMC2240 has a separate SG4_RESULT register which we
        deliberately ignore — we run TMC2240 in SpreadCycle/SG2 mode
        for max torque, and SG2 lives in DRV_STATUS just like TMC5160.
        """
        # TMC2209 has its own SG_RESULT register (SG4)
        if self.is_2209:
            try:
                reg_val = self.tmc.mcu_tmc.get_register('SG_RESULT')
                return reg_val & 0x3FF
            except Exception as e:
                logging.debug(
                    "tmc_flow_test: SG_RESULT read failed: %s", e)
                return None

        # SG2 family (TMC5160 / TMC2130 / TMC2660 / TMC2240):
        # SG_RESULT is in DRV_STATUS bits 0-9
        if self.sg2_driver:
            try:
                reg_val = self.tmc.mcu_tmc.get_register('DRV_STATUS')
                return reg_val & 0x3FF
            except Exception as e:
                logging.debug(
                    "tmc_flow_test: DRV_STATUS read failed: %s", e)
                return None

        # Fallback for any other unknown driver
        try:
            drv = self.tmc.get_status(self.reactor.monotonic())
            if 'drv_status' in drv and isinstance(drv['drv_status'], dict):
                return drv['drv_status'].get('sg_result')
            return drv.get('sg_result')
        except Exception:
            return None

    def _start_sampling(self):
        self.samples_sg = []
        self.samples_time = []
        # Thermal samples are taken at the same cadence as SG, so we
        # only capture data DURING active extrusion. The first thermal
        # sample at t=0 establishes the baseline; later samples show
        # the heater's response to the load.
        self.samples_thermal = []
        self.sample_start_time = self.reactor.monotonic()
        self.sampling_active = True
        # Cadence for thermal sampling — SG sampling is 50 Hz which is
        # overkill for thermal. Sample thermal every Nth SG sample.
        self._thermal_sample_counter = 0
        self.sample_timer = self.reactor.register_timer(
            self._sample_callback, self.reactor.NOW)

    def _stop_sampling(self):
        self.sampling_active = False
        if self.sample_timer is not None:
            self.reactor.unregister_timer(self.sample_timer)
            self.sample_timer = None

    def _sample_callback(self, eventtime):
        if not self.sampling_active:
            return self.reactor.NEVER
        sg = self._read_sg()
        rel_t = eventtime - self.sample_start_time
        # For drivers we read directly from registers, accept any non-None
        # value (including 0). For fallback drivers via get_status, we
        # filter 0/None as those usually mean "not yet polled".
        direct_read = self.is_2209 or self.sg2_driver
        if sg is not None:
            if direct_read or sg > 0:
                self.samples_sg.append(sg)
                self.samples_time.append(rel_t)
        # Thermal sampling at lower cadence (~4 Hz) — temp + PWM don't
        # change fast enough to need the full SG sample rate.
        self._thermal_sample_counter += 1
        if self._thermal_sample_counter >= THERMAL_EVERY_N:
            self._thermal_sample_counter = 0
            self.samples_thermal.append(self._get_thermal_snapshot())
        return eventtime + SAMPLE_INTERVAL

    # ─── Statistics ─────────────────────────────────────────────────

    @staticmethod
    def _stats(samples):
        """Median + IQR + basic stats. Returns None if empty."""
        if not samples:
            return None
        sorted_s = sorted(samples)
        n = len(sorted_s)

        def percentile(p):
            if n == 1:
                return sorted_s[0]
            k = (n - 1) * p / 100.0
            f = int(k)
            c = min(f + 1, n - 1)
            return sorted_s[f] + (k - f) * (sorted_s[c] - sorted_s[f])

        return {
            'min': sorted_s[0], 'max': sorted_s[-1],
            'avg': sum(sorted_s) / n,
            'median': percentile(50),
            'p25': percentile(25),
            'p75': percentile(75),
            'std': _pstdev(sorted_s) if n > 1 else 0.0,
            'n': n,
        }

    # ─── CSV / HTML output ──────────────────────────────────────────

    def _write_csv(self, path, results, meta):
        with open(path, 'w') as f:
            f.write("# TMC Flow Test v%s results\n" % MODULE_VERSION)
            f.write("# Plugin by Steven (Fragmon) — Crydteam\n")
            f.write("# YouTube: https://www.youtube.com/@crydteamprinting\n")
            tmc_settings = None
            for k, v in meta.items():
                if k == 'tmc_settings':
                    tmc_settings = v
                    continue
                f.write("# %s: %s\n" % (k, v))
            # TMC settings as a separate, readable block at the end
            # of the comment header — useful as a paper trail.
            if tmc_settings:
                f.write("#\n# TMC driver settings at test start:\n")
                for label, value, raw in tmc_settings:
                    f.write("#   %s = %s  [%s]\n" % (label, value, raw))

            # Determine extra-thermistor names that appear in any
            # result, so we can emit a stable column order.
            extra_names = []
            seen = set()
            for r in results:
                th = r.get('thermal') or {}
                for name in (th.get('extras') or {}).keys():
                    if name not in seen:
                        seen.add(name)
                        extra_names.append(name)

            # CSV header
            base_header = ("phase,flow_mm3s,sg_median,sg_p25,sg_p75,sg_avg,"
                           "sg_min,sg_max,sg_n,n_repeats,sg_run_cv_pct,"
                           "sg_dips,sg_dip_rate_pct,"
                           "sg_peaks,sg_peak_rate_pct,"
                           "run_sg_avgs,"
                           "temp_target,temp_start,temp_end,temp_min,"
                           "temp_avg,temp_drop,"
                           "pwm_min,pwm_max,pwm_avg,"
                           "tmc_otpw,tmc_ot")
            extra_header = "".join(
                ",extra_%s_min,extra_%s_avg,extra_%s_max"
                % (_sanitize_csv(n), _sanitize_csv(n), _sanitize_csv(n))
                for n in extra_names)
            data_header = (",t_s,idle_pwm,melt_duty,"
                           "melt_duty_per_flow,melt_fit_pct,"
                           "bdw_cmd_mm,bdw_meas_mm,bdw_ratio,"
                           "enc_gap_mm")
            f.write(base_header + extra_header + data_header + "\n")

            # Heat per mm³ vs the fitted low-flow line, in % (only once
            # the run has finished and the fit is in meta).
            def fit_pct(r):
                fa = meta.get('melt_fit_intercept')
                fb = meta.get('melt_fit_slope')
                pwm_r = (r.get('thermal') or {}).get('pwm_avg')
                if fa is None or not fb or pwm_r is None or not r['flow']:
                    return None
                p0 = meta.get('idle_pwm_start')
                cur = _idle_at(r.get('t_s'))
                drift = (cur - p0) if (p0 is not None
                                       and cur is not None) else 0.0
                return (pwm_r - drift - fa) / (fb * r['flow']) * 100

            # Idle duty at a step's time: linear between the start and
            # end baselines when both exist (corrects chamber drift).
            def _idle_at(t):
                p0 = meta.get('idle_pwm_start')
                p1 = meta.get('idle_pwm_end')
                t1 = meta.get('idle_t_end_s')
                if p0 is None:
                    return p1
                if p1 is None or not t1 or t is None:
                    return p0
                frac = min(max(t / t1, 0.0), 1.0)
                return p0 + (p1 - p0) * frac

            for r in results:
                sg = r.get('sg') or {}
                rc = r.get('run_consistency') or {}
                run_sg = r.get('run_sg_avgs') or []
                phase = r.get('phase', 'coarse')
                th = r.get('thermal') or {}

                def fmt(d, key):
                    v = d.get(key, '')
                    if isinstance(v, float):
                        return "%.1f" % v
                    return str(v)

                def fmt_t(d, key, dec=1):
                    v = d.get(key)
                    if v is None:
                        return ''
                    return "%.*f" % (dec, v)

                def fmt_pwm(d, key):
                    v = d.get(key)
                    if v is None:
                        return ''
                    return "%.3f" % v

                def fmt_flag(d, key):
                    v = d.get(key)
                    if v is None:
                        return ''
                    return "%d" % v

                dip = r.get('dip') or {}
                peak = r.get('peak') or {}
                base_row = ("%s,%.2f,%s,%s,%s,%s,%s,%s,%s,"
                            "%d,%s,%s,%s,%s,%s,%s,"
                            "%s,%s,%s,%s,%s,%s,"
                            "%s,%s,%s,"
                            "%s,%s") % (
                    phase,
                    r['flow'],
                    fmt(sg, 'median'), fmt(sg, 'p25'), fmt(sg, 'p75'),
                    fmt(sg, 'avg'), sg.get('min', ''), sg.get('max', ''),
                    sg.get('n', 0),
                    len(run_sg),
                    "%.1f" % rc.get('sg_cv', 0) if rc else '',
                    str(dip.get('count', '')),
                    "%.2f" % (dip['rate'] * 100) if dip else '',
                    str(peak.get('count', '')),
                    "%.2f" % (peak['rate'] * 100) if peak else '',
                    '|'.join("%.1f" % v for v in run_sg),
                    fmt_t(th, 'temp_target'),
                    fmt_t(th, 'temp_start'),
                    fmt_t(th, 'temp_end'),
                    fmt_t(th, 'temp_min'),
                    fmt_t(th, 'temp_avg'),
                    fmt_t(th, 'temp_drop', 2),
                    fmt_pwm(th, 'pwm_min'),
                    fmt_pwm(th, 'pwm_max'),
                    fmt_pwm(th, 'pwm_avg'),
                    fmt_flag(th, 'tmc_otpw_any'),
                    fmt_flag(th, 'tmc_ot_any'),
                )
                # Append extra-thermistor columns in the same order
                # as the header.
                extras = (th.get('extras') or {})
                extra_cols = []
                for name in extra_names:
                    e = extras.get(name) or {}
                    extra_cols.append(fmt_t(e, 'min'))
                    extra_cols.append(fmt_t(e, 'avg'))
                    extra_cols.append(fmt_t(e, 'max'))
                idle = _idle_at(r.get('t_s'))
                pwm_a = th.get('pwm_avg')
                melt = (pwm_a - idle) if (idle is not None
                                          and pwm_a is not None) else None
                b = r.get('bdw') or {}

                def f_or(v, fmt_s):
                    return (fmt_s % v) if v is not None else ''
                data_cols = [
                    f_or(r.get('t_s'), '%.1f'),
                    f_or(idle, '%.4f'),
                    f_or(melt, '%.4f'),
                    f_or(melt / r['flow'] if (melt is not None
                                              and r['flow']) else None,
                         '%.6f'),
                    f_or(fit_pct(r), '%.1f'),
                    f_or(b.get('cmd_mm'), '%.1f'),
                    f_or(b.get('meas_mm'), '%.1f'),
                    f_or(b.get('ratio'), '%.4f'),
                    f_or((r.get('enc') or {}).get('gap_mm'), '%.1f'),
                ]
                f.write(base_row
                        + ("," + ",".join(extra_cols)
                           if extra_cols else "")
                        + "," + ",".join(data_cols)
                        + "\n")

    def _write_html(self, path, results, meta, limit_reason,
                    final_result=None):
        """Compact HTML report with chart.

        final_result (optional dict): when present, renders the prominent
        result panel with the optimal flow value. Expected keys:
            max_safe (float)         -- mm³/s, the optimal value
            verify_cv (float)        -- run-to-run CV in percent
            quality (str)            -- e.g. 'good (stable)'
            stop_reason (str|None)   -- last trigger reason (optional)
            trigger_events (list)    -- chronological trigger / borderline
                                         / verify-fail events for the
                                         decision-trail panel
            baseline_stats (dict)    -- "healthy" CV/IQR ranges from the
                                         coarse phase (used by the
                                         decision-trail and chart bands)
        """
        # ─── Build chart datasets (sorted by flow rate) ─────────────
        # Test phases run NOT in flow order: bisection revisits lower
        # flows after a coarse trigger, verification re-tests the safe
        # value. Plotting in time order makes the line zigzag and
        # confuses interpretation. We therefore:
        #   1. Deduplicate by flow — when the same flow appears multiple
        #      times (e.g. verify after bisect at the same flow), the
        #      LATEST entry wins (it's the most-confirmed measurement).
        #   2. Sort by flow ascending.
        # The Test Details table below keeps chronological order so the
        # decision flow is still inspectable.
        results_by_flow = {}
        for r in results:
            results_by_flow[r['flow']] = r  # later entry overwrites
        results_for_charts = sorted(results_by_flow.values(),
                                     key=lambda r: r['flow'])

        flows = [r['flow'] for r in results_for_charts]
        phases = [r.get('phase', 'coarse') for r in results_for_charts]
        sg_label = self._get_sg_label()
        sg_median = [r['sg']['median'] if r['sg'] else None
                     for r in results_for_charts]
        sg_p25 = [r['sg']['p25'] if r['sg'] else None
                  for r in results_for_charts]
        sg_p75 = [r['sg']['p75'] if r['sg'] else None
                  for r in results_for_charts]
        sg_avg = [r['sg']['avg'] if r['sg'] else None
                  for r in results_for_charts]

        # Thermal data extraction — all values may be None per step.
        # The chart's JS code checks if any non-null exists before
        # rendering the chart at all.
        def th(r, key):
            t = r.get('thermal') or {}
            return t.get(key)
        temp_actual = [th(r, 'temp_avg') for r in results_for_charts]
        temp_target = [th(r, 'temp_target') for r in results_for_charts]
        temp_min    = [th(r, 'temp_min') for r in results_for_charts]
        temp_drop   = [th(r, 'temp_drop') for r in results_for_charts]
        pwm_avg     = [th(r, 'pwm_avg') for r in results_for_charts]
        pwm_max     = [th(r, 'pwm_max') for r in results_for_charts]
        tmc_otpw    = [th(r, 'tmc_otpw_any') for r in results_for_charts]
        tmc_ot      = [th(r, 'tmc_ot_any') for r in results_for_charts]

        # Extra thermistors: collect a stable column-order across all
        # steps, then per-step avg per sensor as parallel arrays.
        extras_seen = []
        extras_seen_set = set()
        for r in results_for_charts:
            t = r.get('thermal') or {}
            for n in (t.get('extras') or {}).keys():
                if n not in extras_seen_set:
                    extras_seen_set.add(n)
                    extras_seen.append(n)
        extras_per_step = {}
        for name in extras_seen:
            vals = []
            for r in results_for_charts:
                t = r.get('thermal') or {}
                extras = t.get('extras') or {}
                e = extras.get(name) or {}
                vals.append(e.get('avg'))
            extras_per_step[name] = vals

        # Intra-run drift load (Phase 3 thermal indicator)
        def ir(r, key):
            v = r.get('intra_run') or {}
            return v.get(key)
        intra_drift = [ir(r, 'mean_drift_load') for r in results_for_charts]

        # CV per step (run-to-run variance). One of the strongest slip
        # signals. Used by the variance bar chart in the new layout.
        cv_data = []
        for r in results_for_charts:
            rc = r.get('run_consistency') or {}
            cv_data.append(rc.get('sg_cv'))

        # Linear feed speed and residence time per step. Computed from
        # melt_zone_length and filament_area. Both purely informational
        # — no triggers depend on them.
        linear_speeds = []
        residence_times = []
        for r in results_for_charts:
            flow = r['flow']
            if self.filament_area > 0:
                ls = flow / self.filament_area  # mm/s
                linear_speeds.append(ls)
                if ls > 0 and self.melt_zone_length > 0:
                    residence_times.append(self.melt_zone_length / ls)
                else:
                    residence_times.append(None)
            else:
                linear_speeds.append(None)
                residence_times.append(None)

        if final_result is not None:
            max_safe = final_result.get('max_safe')
            verify_cv = final_result.get('verify_cv', 0.0)
            quality = final_result.get('quality', '')
            stop_reason = final_result.get('stop_reason')
            trigger_events = final_result.get('trigger_events') or []
            baseline_stats = final_result.get('baseline_stats')

            stop_html = ''
            if stop_reason:
                stop_html = (
                    '<div class="stop-line">Stop trigger that ended the '
                    'search: <em>%s</em></div>' % stop_reason)

            summary_html = (
                '<div class="summary final">'
                '<h2>Maximum Safe Volumetric Flow</h2>'
                '<div class="big-number">%.1f<span class="unit">'
                'mm³/s</span></div>'
                '<p class="result-explainer">'
                'This is the highest extrusion speed (in cubic millimeters '
                'of plastic per second) where your motor still '
                'reliably grips the filament. Beyond this point, the '
                'extruder gear starts to slip on the filament — '
                'producing under-extrusion, layer issues or print '
                'failures.'
                '</p>'
                '<div class="quality-line">Verification quality: '
                '<strong>%s</strong> '
                '<span class="quality-tooltip" title="How consistent the '
                'measurement was across the 5 verify repetitions. Lower '
                'CV = more consistent = more trust in the result.">'
                '(CV = %.1f%%) ⓘ</span></div>'
                '<div class="recommendations">'
                '<h3>Use these values in your slicer</h3>'
                '<div class="rec-table">'
                '<span class="rec-label">Conservative (80%%)</span>'
                '<span class="rec-value">%.1f mm³/s</span>'
                '<span class="rec-note">recommended for slicer</span>'
                '<span class="rec-label">Aggressive (90%%)</span>'
                '<span class="rec-value">%.1f mm³/s</span>'
                '<span class="rec-note">only with margin</span>'
                '</div>'
                '<p class="rec-explainer">'
                'The 80%% value is what most people should use as their '
                'slicer\'s "max volumetric speed" — it leaves headroom '
                'for filament variability, temperature fluctuations and '
                'long prints. The 90%% value is more aggressive and '
                'should only be used after long-term print testing.'
                '</p>'
                '</div>'
                '%s'
                '</div>'
                % (max_safe, quality, verify_cv,
                   max_safe * 0.8, max_safe * 0.9, stop_html))
        elif limit_reason:
            summary_html = (
                '<div class="summary"><h2>Result</h2>'
                '<p>Stop reason: <strong>%s</strong></p></div>'
                % limit_reason)
            trigger_events = []
            baseline_stats = None
            max_safe = None
        else:
            summary_html = (
                '<div class="summary"><h2>Result</h2>'
                '<p>Test completed without trigger.</p></div>')
            trigger_events = []
            baseline_stats = None
            max_safe = None

        # ─── Decision-trail panel: WHY the result is what it is ───
        decision_trail_html = ''
        if trigger_events or baseline_stats:
            kind_labels = {
                'trigger': ('🚨', 'TRIGGER',
                            'Slip detected — phase ended at this flow'),
                'borderline-retest': ('🟡', 'BORDERLINE',
                                      'In gray zone — re-measured for '
                                      'confirmation'),
                'borderline-confirmed': ('🚨', 'BORDERLINE CONFIRMED',
                                         'Re-test also borderline — '
                                         'treated as trigger'),
                'verify-borderline': ('🟡', 'VERIFY BORDERLINE',
                                      'Verify in gray zone — re-measured'),
                'verify-fail': ('⚠️', 'VERIFY FAILED',
                                'Verify trigged at this flow — fell back '
                                'into bisection'),
            }
            event_rows = []
            for ev in trigger_events:
                icon, label, _desc = kind_labels.get(
                    ev['kind'], ('•', ev['kind'].upper(), ''))
                metric_bits = []
                if ev.get('cv') is not None:
                    metric_bits.append('CV %.1f%%' % ev['cv'])
                if ev.get('iqr') is not None:
                    metric_bits.append('IQR %.0f' % ev['iqr'])
                if ev.get('sg_median') is not None:
                    metric_bits.append('%s %.0f'
                                       % (sg_label, ev['sg_median']))
                metrics_str = ' · '.join(metric_bits) if metric_bits else '—'
                event_rows.append(
                    '<tr class="ev-%s">'
                    '<td class="ev-icon">%s</td>'
                    '<td class="ev-phase">%s</td>'
                    '<td class="ev-flow">%.1f mm³/s</td>'
                    '<td class="ev-label">%s</td>'
                    '<td class="ev-metrics">%s</td>'
                    '<td class="ev-reason">%s</td>'
                    '</tr>'
                    % (ev['kind'], icon,
                       ev['phase'], ev['flow'], label,
                       metrics_str, ev['reason']))

            # Baseline summary row
            baseline_block = ''
            if baseline_stats:
                baseline_block = (
                    '<div class="baseline-stats">'
                    '<h3>Healthy ranges (from %d coarse-phase steps)</h3>'
                    '<p class="baseline-explainer">'
                    'During the coarse sweep before slip started, the '
                    'plugin saw consistent measurements — CV and IQR '
                    'stayed in the ranges shown below. Once any later '
                    'step jumped <strong>far outside</strong> these '
                    '"healthy" ranges (e.g. CV doubling, IQR widening '
                    'a lot), that was the plugin\'s signal that slip '
                    'had started. Think of it like measuring a '
                    'person\'s normal heart rate at rest — if it '
                    'suddenly doubles, something\'s changed.'
                    '</p>'
                    '<table class="baseline-table">'
                    '<thead><tr><th>Metric</th><th>Min</th><th>Median</th>'
                    '<th>Mean</th><th>Max</th>'
                    '<th>Trigger floor (bisection)</th></tr></thead>'
                    '<tbody>'
                    '<tr><td>Run-to-run CV (%%)</td>'
                    '<td>%.1f</td><td>%.1f</td><td>%.1f</td><td>%.1f</td>'
                    '<td>≥ 4 %% AND ≥ 2.0× median, OR ≥ 7 %% absolute</td>'
                    '</tr>'
                    '<tr><td>IQR (P75 − P25, raw SG units)</td>'
                    '<td>%.0f</td><td>%.0f</td><td>%.0f</td><td>%.0f</td>'
                    '<td>≥ 18 AND ≥ 2.5× median, OR ≥ 25 absolute</td>'
                    '</tr>'
                    '</tbody></table>'
                    '<p class="baseline-legend"><em>'
                    'CV = how much the 5 repetitions of one step varied. '
                    'IQR = how spread out the SG samples were within '
                    'one step. Both are explained in the glossary above.'
                    '</em></p>'
                    '</div>'
                    % (baseline_stats['n_steps'],
                       baseline_stats['cv_min'],
                       baseline_stats['cv_median'],
                       baseline_stats['cv_mean'],
                       baseline_stats['cv_max'],
                       baseline_stats['iqr_min'],
                       baseline_stats['iqr_median'],
                       baseline_stats['iqr_mean'],
                       baseline_stats['iqr_max']))

            events_block = ''
            if event_rows:
                events_block = (
                    '<h3>Decision events (chronological)</h3>'
                    '<p class="trail-explainer">Each event below is a '
                    'point where the plugin detected something unusual '
                    'and changed its decision. The bottom-most '
                    'event-row is what produced the final result.</p>'
                    '<div class="trail-table-wrap">'
                    '<table class="trail-table">'
                    '<thead><tr>'
                    '<th></th><th>Phase</th><th>Flow</th><th>Verdict</th>'
                    '<th>Metrics</th><th>Reason</th>'
                    '</tr></thead><tbody>%s</tbody></table>'
                    '</div>'
                    % ''.join(event_rows))

            decision_trail_html = (
                '<details class="decision-trail" open>'
                '<summary>Why this value? — Decision trail</summary>'
                '<p class="trail-intro">'
                'This section shows <strong>exactly why</strong> the '
                'plugin chose the final value. The "healthy ranges" '
                'table tells you what normal CV/IQR looks like during '
                'the slip-free portion of the test — the plugin '
                'treats anything well outside those ranges as a slip '
                'signal. The events table below lists every '
                'decision moment in chronological order; the '
                'bottom-most row is what produced the final result.'
                '</p>'
                '%s%s'
                '</details>' % (baseline_block, events_block))

        # Build meta info, but skip the tmc_settings list (rendered
        # separately below as its own collapsible section).
        meta_html = ''.join(
            '<div><strong>%s:</strong> %s</div>' % (k, v)
            for k, v in meta.items() if k != 'tmc_settings')

        # Build TMC settings block (collapsible <details> for compactness)
        tmc_settings = meta.get('tmc_settings') or []
        if tmc_settings:
            tmc_rows = ''.join(
                '<tr><td style="text-align:left">%s</td>'
                '<td style="text-align:right;font-family:monospace">%s</td>'
                '<td style="text-align:left;color:#888;font-family:monospace">'
                '%s</td></tr>'
                % (label, value, raw)
                for label, value, raw in tmc_settings)
            tmc_settings_html = (
                '<details class="tmc-settings"><summary>'
                'TMC driver settings at test start (%d fields) — '
                'click to expand</summary>'
                '<table class="tmc-settings-table"><thead><tr>'
                '<th style="text-align:left">Setting</th>'
                '<th style="text-align:right">Value</th>'
                '<th style="text-align:left">Field</th>'
                '</tr></thead><tbody>%s</tbody></table>'
                '</details>'
                % (len(tmc_settings), tmc_rows))
        else:
            tmc_settings_html = ''

        # Data table
        rows = []
        for r in results:
            sg = r.get('sg') or {}
            rc = r.get('run_consistency') or {}
            cv_str = ''
            cv_class = ''
            if rc and 'sg_cv' in rc:
                cv = rc['sg_cv']
                cv_str = "%.1f%%" % cv
                if cv > 25:
                    cv_class = ' style="background:#ffcdd2"'
                elif cv > 10:
                    cv_class = ' style="background:#fff9c4"'

            def fmt(d, key, fs="%.1f"):
                v = d.get(key)
                if v is None or v == '':
                    return '-'
                return fs % v

            th = r.get('thermal') or {}
            dip = r.get('dip') or {}
            peak = r.get('peak') or {}
            ex = self._step_extras(r, meta)

            def opt(v, fs):
                return (fs % v) if v is not None else '-'
            fail_style = (' style="background:rgba(248,81,73,0.18)"'
                          if r.get('hard_fail') else '')
            rows.append(
                "<tr%s><td>%s</td><td>%.1f</td><td><b>%s</b></td>"
                "<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"
                "<td>%d</td><td%s>%s</td><td>%s</td><td>%s</td>"
                "<td>%s</td><td>%s</td><td>%s</td>"
                "<td>%s</td><td>%s</td><td>%s</td><td>%s</td>"
                "<td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
                    fail_style,
                    r.get('phase', ''), r['flow'], fmt(sg, 'median'),
                    fmt(sg, 'p25'), fmt(sg, 'p75'), fmt(sg, 'avg'),
                    opt(sg.get('min'), '%s'), opt(sg.get('max'), '%s'),
                    sg.get('n', 0), cv_class, cv_str or '-',
                    opt(dip.get('count'), '%d'), opt(peak.get('count'), '%d'),
                    opt(th.get('temp_avg'), '%.1f'),
                    opt(th.get('pwm_avg'), '%.3f'),
                    opt(r.get('t_s'), '%.0f'),
                    opt(ex['idle'], '%.3f'),
                    opt(ex['melt'], '%.3f'),
                    opt(ex['melt_per_flow'], '%.5f'),
                    opt(ex['fit_pct'], '%.0f%%'),
                    opt(ex['bdw_cmd'], '%.0f'),
                    opt(ex['bdw_meas'], '%.0f'),
                    opt(ex['bdw_pct'], '%.0f%%'),
                    opt((r.get('enc') or {}).get('gap_mm'), '%.1f')))

        table_header = (
            "<th>Phase</th><th>Flow (mm³/s)</th><th>%s median</th>"
            "<th>%s P25</th><th>%s P75</th><th>%s avg</th>"
            "<th>min</th><th>max</th>"
            "<th>n</th><th>Inter-run CV</th><th>Dips</th><th>Peaks</th>"
            "<th>Temp avg (°C)</th><th>Heater duty</th><th>t (s)</th>"
            "<th>Idle duty</th><th>Melt duty</th><th>Melt / mm³</th>"
            "<th>vs low-flow line</th>"
            "<th>BDW cmd (mm)</th><th>BDW meas (mm)</th><th>BDW feed</th>"
            "<th>Encoder gap (mm)</th>"
            % (sg_label, sg_label, sg_label, sg_label))

        table = ("<table><thead><tr>" + table_header
                 + "</tr></thead><tbody>"
                 + "".join(rows) + "</tbody></table>")

        # ─── NEW LAYOUT: Build hero, insights and timeline blocks ────

        # Hero block: big number, slicer recommendations, status pill
        if max_safe is not None:
            verify_cv_for_hero = (final_result.get('verify_cv', 0.0)
                                  if final_result else 0.0)
            verify_repeats_for_hero = (final_result.get('verify_repeats', 0)
                                       if final_result else 0)
            quality_for_hero = (final_result.get('quality', '')
                                if final_result else '')
            status_class = 'success'
            status_text = 'Verified — CV %.1f%% over %d runs' % (
                verify_cv_for_hero, verify_repeats_for_hero)
            if str(quality_for_hero).startswith('stopped'):
                status_text = ('Stopped on filament failure — last step '
                               'before the filament stopped')
            if (quality_for_hero in ('borderline', 'fragile')
                    or str(quality_for_hero).startswith('stopped')):
                status_class = 'warning'
            hero_html = (
                '<div class="hero">'
                '<div class="hero-grid">'
                '<div>'
                '<div class="hero-label">Maximum Safe Flow</div>'
                '<div class="hero-value">%.0f<span class="unit">mm³/s'
                '</span></div>'
                '<div class="hero-status %s"><span>%s</span></div>'
                '</div>'
                '<div class="slicer-recos">'
                '<div class="reco-card conservative">'
                '<div class="reco-percent">80%%</div>'
                '<div class="reco-value">%.0f</div>'
                '<div class="reco-label">conservative<br>everyday use</div>'
                '</div>'
                '<div class="reco-card aggressive">'
                '<div class="reco-percent">90%%</div>'
                '<div class="reco-value">%.0f</div>'
                '<div class="reco-label">aggressive<br>after validation</div>'
                '</div>'
                '</div>'
                '</div>'
                '</div>'
                % (max_safe, status_class, status_text,
                   max_safe * 0.8, max_safe * 0.9))
        else:
            stop_text = limit_reason or 'Test ended without verified result'
            hero_html = (
                '<div class="hero">'
                '<div class="hero-grid">'
                '<div>'
                '<div class="hero-label">Test Result</div>'
                '<div class="hero-value" style="font-size:48px">No verified value</div>'
                '<div class="hero-status warning"><span>%s</span></div>'
                '</div>'
                '<div></div>'
                '</div>'
                '</div>'
                % stop_text)

        # Melt limit (print quality) and fit warning under the hero
        if meta.get('melt_limit'):
            hero_html += (
                '<div class="summary" style="margin-top:12px">'
                '<strong>Melt limit (print quality):</strong> %s'
                '%s</div>'
                % (meta['melt_limit'],
                   ('<br><em>Note: %s</em>' % meta['melt_fit_warning'])
                   if meta.get('melt_fit_warning') else ''))
        if meta.get('slip_onset'):
            hero_html += (
                '<div class="summary" style="margin-top:12px">'
                '<strong>Slip start (BDWidth feed):</strong> %s</div>'
                % meta['slip_onset'])

        # Insights cards: 4 quick-glance summaries
        # 1) Result Quality
        if max_safe is not None:
            verify_cv_ic = (final_result.get('verify_cv', 0.0)
                            if final_result else 0.0)
            quality_ic = (final_result.get('quality', '')
                          if final_result else '')
            verify_n_ic = (final_result.get('verify_repeats', 0)
                           if final_result else 0)
            quality_card_class = 'success'
            quality_value = 'Verified'
            if quality_ic in ('borderline', 'fragile'):
                quality_card_class = 'warning'
                quality_value = quality_ic.title()
            elif quality_ic == '':
                quality_value = 'Done'
            quality_detail = ('CV %.1f%% across %d verify runs'
                              % (verify_cv_ic, verify_n_ic))
            if str(quality_ic).startswith('stopped'):
                quality_card_class = 'warning'
                quality_value = 'Stopped'
                quality_detail = ('filament failure — last step before the stop, '
                                  'not verified')
            quality_card = (
                '<div class="insight-card %s">'
                '<div class="insight-icon">✓</div>'
                '<div class="insight-title">Result Quality</div>'
                '<div class="insight-value">%s</div>'
                '<div class="insight-detail">%s</div></div>'
                % (quality_card_class, quality_value, quality_detail))
        else:
            quality_card = (
                '<div class="insight-card warning">'
                '<div class="insight-icon">!</div>'
                '<div class="insight-title">Result Quality</div>'
                '<div class="insight-value">Inconclusive</div>'
                '<div class="insight-detail">No verified value</div></div>')

        # 2) First Trigger
        first_trig_event = None
        for ev in (trigger_events or []):
            if ev.get('kind') == 'trigger':
                first_trig_event = ev
                break
        if first_trig_event:
            trig_metrics = []
            if first_trig_event.get('cv') is not None:
                trig_metrics.append('CV %.1f%%' % first_trig_event['cv'])
            if first_trig_event.get('iqr') is not None:
                trig_metrics.append('IQR %.0f' % first_trig_event['iqr'])
            trig_metric_str = ' · '.join(trig_metrics) if trig_metrics else '—'
            trigger_card = (
                '<div class="insight-card info">'
                '<div class="insight-icon">⊘</div>'
                '<div class="insight-title">First Trigger</div>'
                '<div class="insight-value">%.0f mm³/s</div>'
                '<div class="insight-detail">%s</div></div>'
                % (first_trig_event['flow'], trig_metric_str))
        else:
            trigger_card = (
                '<div class="insight-card info">'
                '<div class="insight-icon">⊘</div>'
                '<div class="insight-title">First Trigger</div>'
                '<div class="insight-value">— </div>'
                '<div class="insight-detail">no trigger fired</div></div>')

        # 3) Thermal Watch
        max_pwm = max((p for p in pwm_max if p is not None), default=None)
        max_drop = max((d for d in temp_drop if d is not None), default=None)
        any_otpw = any(tmc_otpw)
        any_ot = any(tmc_ot)
        # Cold extrusion onset for the insight card. The user-supplied
        # hint takes precedence over auto-detection (visual observation
        # is more reliable than the heater-PWM heuristic).
        cold_extrusion_hint_val = meta.get('cold_extrusion_hint', 0)
        if cold_extrusion_hint_val > 0:
            cold_onset_for_card = cold_extrusion_hint_val
            cold_onset_source = 'user'
        else:
            cold_onset_for_card = self._detect_cold_extrusion_onset(results)
            cold_onset_source = 'auto'

        # Compute peak stress score for the insight card
        score_pairs = self._compute_thermal_stress_per_step(results)
        score_values = [s for _, s in score_pairs if s is not None]
        peak_stress = max(score_values) if score_values else None

        if max_pwm is None and max_drop is None:
            thermal_card = (
                '<div class="insight-card info">'
                '<div class="insight-icon">∅</div>'
                '<div class="insight-title">Thermal Watch</div>'
                '<div class="insight-value">No data</div>'
                '<div class="insight-detail">heater readings unavailable'
                '</div></div>')
        else:
            # The stress-score-based assessment is the PRIMARY indicator.
            # Specific PWM/drop/driver flags add secondary detail.
            if peak_stress is not None:
                if peak_stress < 30:
                    therm_class = 'success'
                    therm_label = 'Stable'
                elif peak_stress < 60:
                    therm_class = 'warning'
                    therm_label = 'Moderate stress'
                else:
                    therm_class = 'danger'
                    therm_label = 'High thermal stress'
            else:
                therm_class = 'success'
                therm_label = 'Stable'

            therm_detail_bits = []
            if peak_stress is not None:
                therm_detail_bits.append('peak stress %.0f/100' % peak_stress)
            if cold_onset_for_card is not None:
                source_tag = ' (user)' if cold_onset_source == 'user' else ''
                therm_detail_bits.append(
                    'onset @ %.0f mm³/s%s' % (cold_onset_for_card, source_tag))
            if max_pwm is not None and max_pwm >= 0.85:
                therm_detail_bits.append('PWM peak %.0f%%' % (max_pwm * 100))
            if max_drop is not None and max_drop >= 1.5:
                therm_detail_bits.append('drop %.1f °C' % max_drop)
            # Recovery deficit: did the heater catch up between reps?
            # If not, that's an even stronger thermal-saturation signal
            # than the drop during a single extrusion.
            recovery_deficits = [
                (r.get('thermal') or {}).get('recovery_deficit_max')
                for r in results]
            recovery_deficits = [d for d in recovery_deficits if d is not None]
            max_recovery = max(recovery_deficits) if recovery_deficits else None
            if max_recovery is not None and max_recovery >= 1.0:
                therm_detail_bits.append(
                    'recovery gap %.1f °C' % max_recovery)
                if max_recovery >= 3.0 and therm_class == 'success':
                    therm_class = 'warning'
                    therm_label = 'Heater not recovering'
            # Override label if driver thermal flags were raised
            if any_ot:
                therm_class = 'danger'
                therm_label = 'Driver overheated'
            elif any_otpw and therm_class == 'success':
                therm_class = 'warning'
                therm_label = 'Driver hot (≥120 °C)'
            therm_detail = (' · '.join(therm_detail_bits)
                            if therm_detail_bits else 'no concerns')
            thermal_card = (
                '<div class="insight-card %s">'
                '<div class="insight-icon">⚠</div>'
                '<div class="insight-title">Thermal Watch</div>'
                '<div class="insight-value">%s</div>'
                '<div class="insight-detail">%s</div></div>'
                % (therm_class, therm_label, therm_detail))

        # 4) Driver Config — concise summary
        driver_label_short = self.driver_type.upper()
        sg_variant = 'SG2' if self.sg2_driver else 'SG4'
        chopper_mode = 'SpreadCycle'  # default for our supported configs
        try:
            if self.tmc is not None:
                if self.is_2209:
                    en_sc = self.tmc.fields.get_field('en_spreadcycle')
                    chopper_mode = ('SpreadCycle' if en_sc == 1
                                    else 'StealthChop')
                else:
                    en_pwm = self.tmc.fields.get_field('en_pwm_mode')
                    chopper_mode = ('StealthChop' if en_pwm == 1
                                    else 'SpreadCycle')
        except (KeyError, AttributeError):
            pass
        # Try to grab a few key values from tmc_settings for the detail line
        tmc_settings_dict = {raw: value for label, value, raw
                              in (meta.get('tmc_settings') or [])}
        detail_bits = []
        if 'sgt' in tmc_settings_dict:
            detail_bits.append('SGT=%s' % tmc_settings_dict['sgt'])
        if 'sgthrs' in tmc_settings_dict:
            detail_bits.append('SGTHRS=%s' % tmc_settings_dict['sgthrs'])
        if 'irun' in tmc_settings_dict:
            detail_bits.append('IRUN=%s' % tmc_settings_dict['irun'])
        config_detail = ', '.join(detail_bits) if detail_bits else 'see details'
        config_card = (
            '<div class="insight-card info">'
            '<div class="insight-icon">⚙</div>'
            '<div class="insight-title">Driver Config</div>'
            '<div class="insight-value">%s / %s</div>'
            '<div class="insight-detail">%s · %s</div></div>'
            % (sg_variant, chopper_mode, driver_label_short, config_detail))

        insights_html = (
            '<div class="insights">%s%s%s%s</div>'
            % (quality_card, trigger_card, thermal_card, config_card))

        # Decision Timeline — friendly summary of phase transitions
        timeline_items = []
        # Group results by phase to produce one entry per coarse-run /
        # bisect-step / verify event.
        coarse_results = [r for r in results if r.get('phase') == 'coarse']
        bisect_results = [r for r in results if r.get('phase') == 'bisect']
        verify_results = [r for r in results if r.get('phase') == 'verify']

        # First coarse summary
        if coarse_results:
            non_trigger_coarse = []
            triggered_coarse = None
            for cr in coarse_results:
                # detect trigger event for this flow+phase
                triggered_here = False
                for ev in (trigger_events or []):
                    if (ev['phase'] == 'coarse'
                            and abs(ev['flow'] - cr['flow']) < 0.01
                            and ev.get('kind') in ('trigger',
                                                   'borderline-confirmed')):
                        triggered_here = True
                        triggered_coarse = (cr, ev)
                        break
                if not triggered_here:
                    non_trigger_coarse.append(cr)
            if non_trigger_coarse:
                f_low = non_trigger_coarse[0]['flow']
                f_high = non_trigger_coarse[-1]['flow']
                first_sg = (non_trigger_coarse[0]['sg'].get('median')
                            if non_trigger_coarse[0].get('sg') else None)
                last_sg = (non_trigger_coarse[-1]['sg'].get('median')
                           if non_trigger_coarse[-1].get('sg') else None)
                last_cv = (non_trigger_coarse[-1].get('run_consistency')
                           or {}).get('sg_cv')
                detail_parts = ['%d coarse steps' % len(non_trigger_coarse)]
                if first_sg is not None and last_sg is not None:
                    detail_parts.append(
                        '%s went from %.0f to %.0f' % (
                            sg_label, first_sg, last_sg))
                if last_cv is not None:
                    detail_parts.append('final CV %.1f%%' % last_cv)
                timeline_items.append(
                    '<div class="timeline-item">'
                    '<div class="timeline-phase coarse">Coarse sweep</div>'
                    '<div class="timeline-content">'
                    '<div class="timeline-flow">%.0f → %.0f mm³/s</div>'
                    '<div class="timeline-detail">%s</div>'
                    '</div></div>'
                    % (f_low, f_high, '. '.join(detail_parts) + '.'))
            if triggered_coarse:
                cr, ev = triggered_coarse
                timeline_items.append(
                    '<div class="timeline-item">'
                    '<div class="timeline-phase trigger">Trigger</div>'
                    '<div class="timeline-content">'
                    '<div class="timeline-flow">%.0f mm³/s — slip detected</div>'
                    '<div class="timeline-detail">%s</div>'
                    '</div></div>'
                    % (cr['flow'], ev.get('reason', 'trigger fired')))

        # Bisection items
        if bisect_results:
            for br in bisect_results:
                cv_b = (br.get('run_consistency') or {}).get('sg_cv')
                # Was this bisect step a borderline-confirmed trigger?
                trig_here = None
                for ev in (trigger_events or []):
                    if (ev['phase'] == 'bisect'
                            and abs(ev['flow'] - br['flow']) < 0.01):
                        trig_here = ev
                        break
                if trig_here:
                    label_text = ('triggered (re-test confirmed)'
                                  if trig_here.get('kind')
                                     == 'borderline-confirmed'
                                  else trig_here.get('kind', 'trigger'))
                    timeline_items.append(
                        '<div class="timeline-item">'
                        '<div class="timeline-phase trigger">Bisect</div>'
                        '<div class="timeline-content">'
                        '<div class="timeline-flow">%.0f mm³/s — %s</div>'
                        '<div class="timeline-detail">%s</div>'
                        '</div></div>'
                        % (br['flow'], label_text,
                           trig_here.get('reason', '')))
                else:
                    cv_str = ('CV %.1f%%' % cv_b) if cv_b is not None else '—'
                    timeline_items.append(
                        '<div class="timeline-item">'
                        '<div class="timeline-phase bisect">Bisect</div>'
                        '<div class="timeline-content">'
                        '<div class="timeline-flow">%.0f mm³/s — passed</div>'
                        '<div class="timeline-detail">%s</div>'
                        '</div></div>'
                        % (br['flow'], cv_str))

        # Verify items
        if verify_results:
            last_verify = verify_results[-1]
            cv_v = (last_verify.get('run_consistency') or {}).get('sg_cv')
            v_failed = False
            for ev in (trigger_events or []):
                if ev.get('kind') == 'verify-fail':
                    v_failed = True
                    break
            if v_failed:
                timeline_items.append(
                    '<div class="timeline-item">'
                    '<div class="timeline-phase trigger">Verify</div>'
                    '<div class="timeline-content">'
                    '<div class="timeline-flow">%.0f mm³/s — verify failed</div>'
                    '<div class="timeline-detail">Re-bisected to find a stable point.</div>'
                    '</div></div>'
                    % last_verify['flow'])
            else:
                cv_str = ('CV %.1f%%' % cv_v) if cv_v is not None else '—'
                n_repeats = len(last_verify.get('run_sg_avgs') or [])
                timeline_items.append(
                    '<div class="timeline-item">'
                    '<div class="timeline-phase verify">Verify</div>'
                    '<div class="timeline-content">'
                    '<div class="timeline-flow">%.0f mm³/s — confirmed</div>'
                    '<div class="timeline-detail">%d repetitions, %s. '
                    'Result is verified.</div>'
                    '</div></div>'
                    % (last_verify['flow'], n_repeats, cv_str))

        if timeline_items:
            timeline_html = (
                '<div class="section">'
                '<div class="section-header">'
                '<div class="section-title">Decision Timeline</div>'
                '<div class="section-meta">why this value</div>'
                '</div>'
                '<div class="timeline">%s</div>'
                '</div>'
                % ''.join(timeline_items))
        else:
            timeline_html = ''

        # Driver label for header
        driver_label_full = '%s / %s' % (self.driver_type, self.stepper_name)

        # Decide CV trigger threshold for the variance bar chart
        cv_trigger_threshold = self.profile.CV_HIGH_VARIANCE

        html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TMC Flow Test — %(timestamp)s</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;700&family=Manrope:wght@400;500;600;700;800&display=swap" rel="stylesheet">
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0"></script>
<script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation@3.0.1"></script>
<style>
:root {
    --bg-primary: #0e1117; --bg-secondary: #161b22;
    --bg-elevated: #1c2128; --bg-card: #21262d;
    --bg-card-hover: #30363d; --border: #30363d;
    --border-bright: #484f58; --text-primary: #e6edf3;
    --text-secondary: #8b949e; --text-tertiary: #6e7681;
    --accent: #58a6ff; --accent-dim: #1f6feb;
    --success: #3fb950; --success-bg: rgba(63, 185, 80, 0.12);
    --warning: #d29922; --warning-bg: rgba(210, 153, 34, 0.12);
    --danger: #f85149; --danger-bg: rgba(248, 81, 73, 0.12);
    --info: #58a6ff; --info-bg: rgba(88, 166, 255, 0.12);
    --font-mono: 'JetBrains Mono', monospace;
    --font-sans: 'Manrope', system-ui, sans-serif;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: var(--font-sans); background: var(--bg-primary);
    color: var(--text-primary); line-height: 1.5; min-height: 100vh;
    background-image:
      radial-gradient(at 20%% 0%%, rgba(31, 111, 235, 0.08) 0px, transparent 50%%),
      radial-gradient(at 80%% 0%%, rgba(88, 166, 255, 0.04) 0px, transparent 50%%);
    background-attachment: fixed;
}
.container { max-width: 1280px; margin: 0 auto; padding: 0 24px; }
header {
    border-bottom: 1px solid var(--border); padding: 20px 0; margin-bottom: 32px;
    background: rgba(14, 17, 23, 0.85); backdrop-filter: blur(8px);
    position: sticky; top: 0; z-index: 100;
}
header .container {
    display: flex; align-items: center; justify-content: space-between; gap: 24px;
}
.brand { display: flex; align-items: center; gap: 12px; }
.brand-logo {
    width: 36px; height: 36px; border-radius: 8px;
    background: linear-gradient(135deg, var(--accent) 0%%, var(--accent-dim) 100%%);
    display: grid; place-items: center; font-weight: 800; color: white;
    font-family: var(--font-mono); font-size: 14px; letter-spacing: -0.05em;
}
.brand-text { line-height: 1.2; }
.brand-name { font-weight: 700; font-size: 15px; letter-spacing: -0.01em; }
.brand-subtitle { color: var(--text-tertiary); font-size: 12px; font-weight: 500; }
.header-meta { display: flex; gap: 16px; font-size: 12px;
               color: var(--text-tertiary); font-family: var(--font-mono); }
.header-meta-item { display: flex; align-items: center; gap: 6px; }
.header-meta-item .dot { width: 6px; height: 6px; border-radius: 50%%; background: var(--success); }
.hero {
    background: linear-gradient(135deg, rgba(63, 185, 80, 0.08) 0%%,
                rgba(31, 111, 235, 0.05) 100%%);
    border: 1px solid var(--border); border-radius: 16px;
    padding: 40px 48px; margin-bottom: 32px; position: relative; overflow: hidden;
}
.hero::before {
    content: ''; position: absolute; top: 0; right: 0; width: 240px; height: 240px;
    background: radial-gradient(circle, rgba(63, 185, 80, 0.15) 0%%, transparent 70%%);
    pointer-events: none;
}
.hero-grid { display: grid; grid-template-columns: 1fr 320px; gap: 48px;
             align-items: center; position: relative; }
.hero-label { font-size: 12px; text-transform: uppercase; letter-spacing: 0.15em;
              color: var(--text-tertiary); font-weight: 600; margin-bottom: 8px; }
.hero-value { font-family: var(--font-mono); font-size: 88px; font-weight: 700;
              line-height: 1; letter-spacing: -0.04em; color: var(--text-primary);
              display: flex; align-items: baseline; gap: 12px; margin-bottom: 16px; }
.hero-value .unit { font-size: 24px; color: var(--text-secondary);
                    font-weight: 500; letter-spacing: 0; }
.hero-status {
    display: inline-flex; align-items: center; gap: 8px;
    padding: 6px 12px; background: var(--success-bg);
    border: 1px solid rgba(63, 185, 80, 0.3); border-radius: 999px;
    font-size: 12px; font-weight: 600; color: var(--success);
    text-transform: uppercase; letter-spacing: 0.05em;
}
.hero-status.warning { background: var(--warning-bg);
                       border-color: rgba(210, 153, 34, 0.3); color: var(--warning); }
.hero-status.warning::before { background: var(--warning); box-shadow: 0 0 8px var(--warning); }
.hero-status::before { content: ''; width: 8px; height: 8px; border-radius: 50%%;
                       background: var(--success); box-shadow: 0 0 8px var(--success); }
.slicer-recos { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; }
.reco-card { background: var(--bg-card); border: 1px solid var(--border);
             border-radius: 12px; padding: 16px; text-align: center;
             transition: border-color 0.2s; }
.reco-card:hover { border-color: var(--border-bright); }
.reco-percent { font-size: 11px; color: var(--text-tertiary);
                text-transform: uppercase; letter-spacing: 0.1em; font-weight: 600; }
.reco-value { font-family: var(--font-mono); font-size: 28px; font-weight: 700;
              color: var(--text-primary); line-height: 1.1; margin: 4px 0; }
.reco-label { font-size: 11px; color: var(--text-secondary); }
.reco-card.conservative .reco-percent { color: var(--success); }
.reco-card.aggressive .reco-percent { color: var(--warning); }
.insights { display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px;
            margin-bottom: 32px; }
.insight-card { background: var(--bg-card); border: 1px solid var(--border);
                border-radius: 10px; padding: 16px; transition: border-color 0.2s; }
.insight-card:hover { border-color: var(--border-bright); }
.insight-icon { display: inline-flex; align-items: center; justify-content: center;
                width: 28px; height: 28px; border-radius: 6px; margin-bottom: 10px;
                font-size: 14px; }
.insight-card.success .insight-icon { background: var(--success-bg); color: var(--success); }
.insight-card.warning .insight-icon { background: var(--warning-bg); color: var(--warning); }
.insight-card.info .insight-icon { background: var(--info-bg); color: var(--info); }
.insight-card.danger .insight-icon { background: var(--danger-bg); color: var(--danger); }
.insight-title { font-size: 11px; text-transform: uppercase; letter-spacing: 0.1em;
                 color: var(--text-tertiary); margin-bottom: 4px; font-weight: 600; }
.insight-value { font-family: var(--font-mono); font-size: 20px; font-weight: 700;
                 color: var(--text-primary); line-height: 1.2; margin-bottom: 4px; }
.insight-detail { font-size: 12px; color: var(--text-secondary); }
.section { margin-bottom: 32px; }
.section-header { display: flex; align-items: baseline; justify-content: space-between;
                  margin-bottom: 16px; border-bottom: 1px solid var(--border);
                  padding-bottom: 8px; }
.section-title { font-size: 18px; font-weight: 700; color: var(--text-primary);
                 letter-spacing: -0.01em; }
.section-meta { font-size: 12px; color: var(--text-tertiary); font-family: var(--font-mono); }
.setup-notes-section .section-meta { font-family: inherit; }
.notes-export-btn {
    background: var(--bg-card); color: var(--text-primary);
    border: 1px solid var(--border); border-radius: 6px;
    padding: 6px 12px; font-family: inherit; font-size: 12px;
    cursor: pointer; transition: all 0.15s;
}
.notes-export-btn:hover {
    border-color: var(--accent); color: var(--accent);
}
.setup-notes-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 12px;
    background: var(--bg-card); border: 1px solid var(--border);
    border-radius: 12px; padding: 20px;
}
.setup-field { display: flex; flex-direction: column; gap: 4px; }
.setup-field-wide { grid-column: 1 / -1; }
.setup-field label {
    font-size: 11px; font-weight: 600;
    color: var(--text-tertiary); letter-spacing: 0.04em;
    text-transform: uppercase;
}
.setup-field input, .setup-field textarea {
    background: var(--bg-page); color: var(--text-primary);
    border: 1px solid var(--border); border-radius: 6px;
    padding: 8px 10px; font-family: inherit; font-size: 13px;
    transition: border-color 0.15s;
    width: 100%%; box-sizing: border-box;
}
.setup-field input:focus, .setup-field textarea:focus {
    outline: none; border-color: var(--accent);
}
.setup-field textarea { resize: vertical; min-height: 60px; }
.setup-notes-help {
    margin-top: 12px; font-size: 12px; color: var(--text-tertiary);
    line-height: 1.5;
}
.tabs { display: flex; gap: 4px; margin-bottom: 16px; border-bottom: 1px solid var(--border); }
.tab { padding: 8px 16px; background: transparent; border: none;
       color: var(--text-secondary); font-family: inherit; font-size: 13px;
       font-weight: 600; cursor: pointer; border-bottom: 2px solid transparent;
       margin-bottom: -1px; transition: all 0.2s; }
.tab:hover { color: var(--text-primary); }
.tab.active { color: var(--accent); border-bottom-color: var(--accent); }
.tab-panel { display: none; }
.tab-panel.active { display: block; }
.chart-card { background: var(--bg-card); border: 1px solid var(--border);
              border-radius: 12px; padding: 24px; margin-bottom: 16px; }
.chart-explainer { font-size: 13px; color: var(--text-secondary);
                   margin-bottom: 20px; line-height: 1.6;
                   padding-left: 12px; border-left: 3px solid var(--accent-dim); }
.chart-explainer strong { color: var(--text-primary); }
canvas { max-height: 400px; }
.timeline { background: var(--bg-card); border: 1px solid var(--border);
            border-radius: 12px; padding: 24px; }
.timeline-item { display: grid; grid-template-columns: 100px 1fr; gap: 16px;
                 padding: 14px 0; border-bottom: 1px solid var(--border); }
.timeline-item:last-child { border-bottom: none; }
.timeline-phase { font-family: var(--font-mono); font-size: 11px; text-transform: uppercase;
                  letter-spacing: 0.1em; color: var(--text-tertiary); font-weight: 600;
                  padding-top: 2px; }
.timeline-phase.coarse { color: var(--info); }
.timeline-phase.bisect { color: var(--warning); }
.timeline-phase.verify { color: var(--success); }
.timeline-phase.trigger { color: var(--danger); }
.timeline-content { font-size: 13px; }
.timeline-flow { font-family: var(--font-mono); font-weight: 700;
                 color: var(--text-primary); margin-bottom: 4px; }
.timeline-detail { color: var(--text-secondary); }
details { background: var(--bg-secondary); border: 1px solid var(--border);
          border-radius: 10px; margin-bottom: 8px; }
details summary { padding: 14px 20px; cursor: pointer; font-weight: 600;
                  color: var(--text-primary); font-size: 14px; display: flex;
                  align-items: center; justify-content: space-between;
                  list-style: none; user-select: none; }
details summary::-webkit-details-marker { display: none; }
details summary::after { content: '↓'; font-family: var(--font-mono);
                         color: var(--text-tertiary); transition: transform 0.2s;
                         font-size: 14px; }
details[open] summary::after { transform: rotate(180deg); color: var(--accent); }
details > div { padding: 0 20px 20px; color: var(--text-secondary);
                font-size: 13px; line-height: 1.7; }
details table { font-size: 12px; }
table { width: 100%%; border-collapse: collapse; font-family: var(--font-mono);
        font-size: 12px; }
table th { text-align: left; padding: 10px 8px; border-bottom: 1px solid var(--border-bright);
           color: var(--text-tertiary); font-weight: 600; text-transform: uppercase;
           letter-spacing: 0.05em; font-size: 11px; background: var(--bg-elevated);
           position: sticky; top: 0; }
table td { padding: 8px; border-bottom: 1px solid var(--border); color: var(--text-secondary); }
table tr:hover td { background: var(--bg-card-hover); color: var(--text-primary); }
.row-coarse td:first-child { color: var(--info); font-weight: 600; }
.row-bisect td:first-child { color: var(--warning); font-weight: 600; }
.row-verify td:first-child { color: var(--success); font-weight: 600; }
footer { border-top: 1px solid var(--border); padding: 24px 0; margin-top: 48px;
         color: var(--text-tertiary); font-size: 12px; }
footer .container { display: flex; justify-content: space-between;
                    align-items: center; flex-wrap: wrap; gap: 16px; }
footer a { color: var(--accent); text-decoration: none; }
footer a:hover { text-decoration: underline; }
.footer-credits { display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
.footer-credits .author { color: var(--text-secondary); font-weight: 600; }
@media (max-width: 900px) {
    .hero { padding: 24px; }
    .hero-grid { grid-template-columns: 1fr; gap: 24px; }
    .hero-value { font-size: 64px; }
    .insights { grid-template-columns: repeat(2, 1fr); }
    header .container { flex-direction: column; align-items: flex-start; }
    .header-meta { flex-wrap: wrap; }
}
@media (max-width: 600px) {
    .container { padding: 0 16px; }
    .insights { grid-template-columns: 1fr; }
    .timeline-item { grid-template-columns: 1fr; gap: 4px; }
    .tabs { overflow-x: auto; }
}
</style>
</head>
<body>

<header>
  <div class="container">
    <div class="brand">
      <div class="brand-logo">TMC</div>
      <div class="brand-text">
        <div class="brand-name">TMC Flow Test</div>
        <div class="brand-subtitle">Adaptive max-volumetric-flow detection</div>
      </div>
    </div>
    <div class="header-meta">
      <div class="header-meta-item"><span class="dot"></span><span>v%(version)s</span></div>
      <div class="header-meta-item"><span>%(timestamp)s</span></div>
      <div class="header-meta-item"><span>%(driver_label)s</span></div>
    </div>
  </div>
</header>

<div class="container">

%(hero_html)s

%(insights_html)s

<div class="section setup-notes-section">
  <div class="section-header">
    <div class="section-title">Setup Notes</div>
    <div class="section-meta">
      <button id="exportNotesBtn" class="notes-export-btn"
              title="Save these fields into the HTML so reopening keeps them"
              onclick="saveSetupNotes()">💾 Bake into HTML</button>
    </div>
  </div>
  <div class="setup-notes-grid">
    <div class="setup-field">
      <label>Hotend</label>
      <input type="text" data-key="hotend"
             placeholder="e.g. Rapido HF, Volcano…">
    </div>
    <div class="setup-field">
      <label>Extruder</label>
      <input type="text" data-key="extruder"
             placeholder="e.g. Sherpa Mini, Orbiter v2…">
    </div>
    <div class="setup-field">
      <label>Nozzle</label>
      <input type="text" data-key="nozzle"
             placeholder="e.g. Bondtech CHT 0.4 brass…">
    </div>
    <div class="setup-field">
      <label>Filament</label>
      <input type="text" data-key="filament"
             placeholder="e.g. Polymaker PolyLite PLA, black, 1.75…">
    </div>
    <div class="setup-field">
      <label>Filament manufacturer</label>
      <input type="text" data-key="manufacturer"
             placeholder="e.g. Polymaker, Sunlu, Prusament…">
    </div>
    <div class="setup-field setup-field-wide">
      <label>Notes (free text)</label>
      <textarea rows="3" data-key="notes"
                placeholder="Anything else worth recording for later: chamber temp, filament age, drying schedule…"></textarea>
    </div>
  </div>
  <div class="setup-notes-help">
    Fill in any fields you want to record alongside this report. Click
    <strong>💾 Bake into HTML</strong> to write the values into this file's
    HTML — that way the report can be saved/shared and the values stay with it.
    Otherwise the fields persist only in your browser.
  </div>
</div>

<div class="section">
  <div class="section-header">
    <div class="section-title">Visual Analysis</div>
    <div class="section-meta">%(measurement_count)s measurements</div>
  </div>
  <div class="tabs">
    <button class="tab active" data-target="sg" onclick="showTab(event, 'sg')">StallGuard signal</button>
    <button class="tab" data-target="thermal" onclick="showTab(event, 'thermal')">Thermal profile</button>
    <button class="tab" data-target="variance" onclick="showTab(event, 'variance')">Run-to-run variance</button>
    <button class="tab" data-target="melt" onclick="showTab(event, 'melt')">Melt &amp; feed</button>
    <button class="tab" data-target="compare" onclick="showTab(event, 'compare')">Compare runs</button>
  </div>
  <div id="tab-sg" class="tab-panel active">
    <div class="chart-card">
      <div class="chart-explainer">
        The <strong>median SG line</strong> shows motor load at each tested
        flow. The <strong>shaded band</strong> covers the P25–P75 sample range.
        <strong>Lower line = more load.</strong> Sudden jumps, spikes upward,
        or wide bands indicate slip is starting. Phase boundaries marked as
        vertical lines.
      </div>
      <canvas id="sgChart"></canvas>
    </div>
  </div>
  <div id="tab-thermal" class="tab-panel">
    <div class="chart-card" id="thermalChartCard" style="display:none;">
      <div class="chart-explainer">
        Tracks <strong>hotend temperature</strong> (red),
        <strong>heater PWM</strong> (orange) and <strong>residence time</strong>
        in melt zone (green) per step. <strong>Two thermal limits to watch:</strong>
        <em>heater power</em> (PWM saturates at 1.0 + temperature falls)
        or <em>filament speed</em> (residence below your hotend's class minimum).
        Reference lines: V6 ≥1.5 s, Volcano ≥0.6 s, CHT/HF ≥0.3 s.
      </div>
      <canvas id="thermalChart"></canvas>
    </div>
    <div class="chart-card" id="thermalChartEmpty">
      <p style="color: var(--text-tertiary); text-align: center; padding: 40px;">
        No thermal data captured for this run.
      </p>
    </div>
  </div>
  <div id="tab-variance" class="tab-panel">
    <div class="chart-card">
      <div class="chart-explainer">
        <strong>Run-to-run CV</strong> measures how much the repeated runs at
        each flow varied from each other. Low CV (&lt;3 %%) means measurements
        are repeatable. When CV spikes upward, the motor is producing
        inconsistent SG readings between runs — the classic slip signature.
        Bars are coloured: green (clean), yellow (borderline), red (above CV
        trigger threshold).
      </div>
      <canvas id="cvChart"></canvas>
    </div>
  </div>
  <div id="tab-compare" class="tab-panel">
    <div class="chart-card">
      <div class="chart-explainer">
        <strong>Compare runs.</strong> Load CSVs from other runs of this
        test (e.g. the same filament at other temperatures) to overlay them
        on this one. Everything is computed in your browser from the CSVs —
        the melt limit uses the same low-flow line fit as the plugin, so
        older CSVs work too. The CSV that produced this report, and files
        already loaded, are rejected. Solid line = this report, dashed =
        loaded runs.
      </div>
      <input type="file" id="cmpFile" accept=".csv,text/csv" multiple>
      <div id="cmpMsg" style="margin-top:8px; color: var(--text-tertiary);"></div>
      <table style="margin-top:12px">
        <thead><tr><th>Run</th><th>Material</th><th>Target temp</th><th>Timestamp</th>
          <th>Melt limit (mm³/s)</th><th>Result (mm³/s)</th><th>Status</th>
          <th></th></tr></thead>
        <tbody id="cmpTable"></tbody>
      </table>
    </div>
    <div class="chart-card">
      <div class="chart-explainer"><strong>Heat vs low-flow line</strong> —
        the melt limit is where a run drops below the dashed line.</div>
      <canvas id="cmpMeltChart"></canvas>
    </div>
    <div class="chart-card">
      <div class="chart-explainer"><strong>BDWidth feed</strong> — filament
        the motion sensor saw vs commanded.</div>
      <canvas id="cmpBdwChart"></canvas>
    </div>
    <div class="chart-card">
      <div class="chart-explainer"><strong>StallGuard median</strong> per
        flow step.</div>
      <canvas id="cmpSgChart"></canvas>
    </div>
  </div>
  <div id="tab-melt" class="tab-panel">
    <div class="chart-card" id="meltChartCard">
      <div class="chart-explainer">
        <strong>Heat vs the low-flow line</strong> (blue): heater duty above
        idle compared with the straight line fitted through the healthy
        low-flow steps. Flat near 100 %% = the plastic fully melts; the
        first drop below the dashed 85 %% line is the <em>melt limit</em>
        (print quality). <strong>BDWidth feed</strong> (orange): filament
        the motion sensor saw entering the extruder vs commanded. A sharp
        drop means slip or a stall; the red dashed line is the feed-check
        limit. Sensor accuracy depends on the filament surface.
      </div>
      <canvas id="meltChart"></canvas>
    </div>
    <div class="chart-card" id="meltChartEmpty" style="display:none;">
      <p style="color: var(--text-tertiary); text-align: center; padding: 40px;">
        No melt or BDWidth data for this run (needs IDLE &gt; 0 and/or a
        bdwidth sensor).
      </p>
    </div>
  </div>
</div>

%(timeline_html)s

<div class="section">
  <div class="section-header">
    <div class="section-title">Test Details</div>
    <div class="section-meta">click to expand</div>
  </div>
  <details>
    <summary>Full data table</summary>
    <div>%(data_table)s</div>
  </details>
  <details>
    <summary>Test configuration</summary>
    <div>%(meta_html)s</div>
  </details>
  <details>
    <summary>TMC driver settings at test start</summary>
    <div>%(tmc_settings_html)s</div>
  </details>
  %(decision_trail_html)s
</div>

<div class="section">
  <div class="section-header">
    <div class="section-title">Reference</div>
    <div class="section-meta">explanations</div>
  </div>

  <details>
    <summary>What does "max volumetric flow" actually mean?</summary>
    <div><p>It's the maximum amount of plastic your extruder can push per
    second while still maintaining grip on the filament. Above this value
    the motor either slips on the gear teeth or extrudes filament that
    wasn't fully melted. Both lead to under-extrusion in real prints.
    Slicers use this value to cap volumetric speed even when you set high
    print speeds.</p></div>
  </details>

  <details>
    <summary>StallGuard, SG_RESULT, CV — what's measured?</summary>
    <div><p><strong>StallGuard</strong> is a TMC driver feature that
    measures motor torque without external sensors. It works by analyzing
    the back-EMF of the coils. <strong>SG_RESULT</strong> is its output —
    on TMC5160 it ranges 0–1023 (lower = more load), on TMC2209 it's
    0–510. The plugin samples this value many times per second.</p>
    <p><strong>Run-to-run CV</strong> is the coefficient of variation
    between repeated measurements at the same flow. Low CV (&lt;3 %%)
    means measurements are repeatable. High CV means something is
    inconsistent — usually slip starting.</p></div>
  </details>

  <details>
    <summary>Why the 80%% / 90%% recommendations?</summary>
    <div><p>The "Maximum Safe" value is what the plugin found just below
    where slip starts. The 80%% recommendation gives you a 20%% safety
    buffer — recommended for everyday use because filament thickness
    varies, hotend temperature drifts, and longer prints stress the system
    more. The 90%% value is closer to the limit and should only be used
    after you've validated it with several long real-world prints.</p></div>
  </details>

  <details>
    <summary>Hotend temperature & heater PWM</summary>
    <div><p>The thermal chart logs the hotend's behavior during each step.
    <strong>Temperature drops</strong> at high flow happen when the heater
    can't melt filament fast enough. <strong>Heater PWM</strong> reaching
    1.0 means Klipper is asking for full power; if it stays there for
    several seconds AND temperature falls, the system is near thermal
    saturation. <strong>Driver thermal flags</strong> (OTPW ≥120 °C,
    OT ≥150 °C) appear as triangles when the TMC chip itself gets hot.
    None of this triggers the test or affects slip detection — purely
    informational.</p></div>
  </details>

  <details>
    <summary>Residence time in the melt zone</summary>
    <div><p>The number of seconds each piece of filament spends inside
    the heated melt zone. Computed as
    <code>melt_zone_length / linear_feed_speed</code>. As flow increases,
    residence time shrinks — even with plenty of heater power, very short
    residence times mean the filament core stays cold.</p>
    <p>Each hotend design has a <strong>minimum residence time</strong>
    needed for clean melting. The reference lines in the chart show these
    minimums, ordered from weakest to strongest:</p>
    <ul>
    <li><strong>V6-class</strong> — minimum ~1.5 s residence</li>
    <li><strong>Volcano-class</strong> — minimum ~0.6 s residence</li>
    <li><strong>CHT / HF-class</strong> — minimum ~0.3 s residence</li>
    <li><strong>Goliath / Bondtech-CHT</strong> — &lt; 0.3 s ok</li>
    </ul>
    <p>The tooltip shows the <strong>minimum required hotend class</strong>
    at each flow — anything stronger also works. Example: 0.8 s residence
    needs at minimum a Volcano, but a CHT or Goliath would also handle it
    fine. If your hotend is below the minimum class shown at the slip point,
    the bottleneck was filament speed (not heater wattage).</p></div>
  </details>

  <details>
    <summary>Coloured zones &amp; thermal stress score</summary>
    <div><p>The chart backgrounds are tinted in three zones based on a
    heuristic <strong>thermal stress score</strong> calculated per step.
    The score combines five signals (each weighted differently):</p>
    <ul>
    <li><strong>Heater PWM level</strong> — how hard the heater is working
    (0–30 points). PWM 30%% adds nothing, 95%%+ adds full 30.</li>
    <li><strong>Temperature drop from target</strong> — how much actual
    temp falls below target (0–25 points). 0 °C drops adds nothing,
    5 °C+ adds full 25.</li>
    <li><strong>PWM rising trend</strong> — heater taking on more work
    relative to early-test baseline (0–15 points).</li>
    <li><strong>PWM peak hits</strong> — saturation events even briefly
    (0–10 points).</li>
    <li><strong>Intra-run SG drift</strong> — load drifting upward
    DURING individual 5-second runs (0–20 points). This is a strong
    cold-extrusion signature: as filament heats up too slowly, motor
    load increases continuously through a single run rather than
    staying steady. Pure motor slip jumps abruptly; cold extrusion
    drifts gradually.</li>
    </ul>
    <p>Score interpretation:
    <strong style="color:#3fb950">0–30 = green / stable</strong>,
    <strong style="color:#d29922">30–60 = yellow / moderate stress</strong>,
    <strong style="color:#f85149">60–100 = red / high stress (cold extrusion likely)</strong>.</p>
    <p>The yellow zone in the charts marks where the score first crosses
    30. <strong>Important:</strong> This is a heuristic combining heater
    behavior with within-run load trends. It does NOT trigger or modify
    the slip-detection logic. Hover any thermal-chart point to see the
    per-step score AND the intra-run drift value (positive %% = load
    increased during the run, the cold-extrusion fingerprint).</p>
    <p><strong>About the data:</strong> All thermal values (PWM, temperature)
    are sampled at 4 Hz <em>during active extrusion only</em> — never
    during the cooldown gaps between reps. This avoids the otherwise
    confusing pattern where PWM drops to baseline and temperature
    overshoots while the heater is idle.</p></div>
  </details>
</div>

</div>

<footer>
  <div class="container">
    <div class="footer-credits">
      <span>Plugin by</span>
      <span class="author">Steven (Fragmon) — Crydteam</span>
      <span>·</span>
      <a href="https://www.youtube.com/@crydteamprinting">YouTube</a>
      <span>·</span>
      <a href="https://github.com/Fragmon/klipper_max_flow_test">GitHub</a>
    </div>
    <div>Generated %(timestamp)s · v%(version)s · GPL-3.0</div>
  </div>
</footer>

<script>
function showTab(event, name) {
    document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
    document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
    event.currentTarget.classList.add('active');
    document.getElementById('tab-' + name).classList.add('active');
    // Charts built while their tab was hidden place their annotation
    // lines against a zero-size chart area; redraw them once visible.
    if (name === 'compare' && window._cmpRender) { window._cmpRender(); }
    if (typeof Chart !== 'undefined') {
        document.querySelectorAll('#tab-' + name + ' canvas').forEach(c => {
            const ch = Chart.getChart(c);
            if (ch) { ch.resize(); ch.update(); }
        });
    }
}

// ─── Setup-notes persistence ────────────────────────────────────────
// Two layers of persistence:
// (1) Live-edits in localStorage so the fields survive a page reload.
// (2) "Bake into HTML" button writes the values into the document
//     itself, so saving the HTML preserves the values.
const NOTES_KEY = 'tmc_flow_notes_' + (document.title || 'report');
function loadSetupNotes() {
    document.querySelectorAll('.setup-field [data-key]').forEach(el => {
        // (a) prefer a value baked into the HTML attribute
        const baked = el.getAttribute('data-baked');
        if (baked !== null && baked !== '') {
            el.value = baked;
            return;
        }
        // (b) otherwise restore from localStorage
        try {
            const saved = localStorage.getItem(NOTES_KEY + ':' + el.dataset.key);
            if (saved !== null) el.value = saved;
        } catch (e) {}
    });
    // Persist on every change
    document.querySelectorAll('.setup-field [data-key]').forEach(el => {
        el.addEventListener('input', () => {
            try {
                localStorage.setItem(
                    NOTES_KEY + ':' + el.dataset.key, el.value);
            } catch (e) {}
        });
    });
}
function saveSetupNotes() {
    // Bake current values into the document's HTML so a Save-As keeps
    // them. We write the value as `data-baked="..."` on each input.
    document.querySelectorAll('.setup-field [data-key]').forEach(el => {
        el.setAttribute('data-baked', el.value);
        // Also set the `value` attribute for inputs / textContent for
        // textarea so the serialized HTML reflects the live value.
        if (el.tagName.toLowerCase() === 'textarea') {
            el.textContent = el.value;
        } else {
            el.setAttribute('value', el.value);
        }
    });
    const btn = document.getElementById('exportNotesBtn');
    const orig = btn.textContent;
    btn.textContent = '✓ Baked. Now File → Save Page As…';
    setTimeout(() => { btn.textContent = orig; }, 4000);
}
window.addEventListener('DOMContentLoaded', loadSetupNotes);

const flows = %(flows)s;
const phases = %(phases)s;
const sgMedian = %(sg_median)s;
const sgP25 = %(sg_p25)s;
const sgP75 = %(sg_p75)s;
const sgAvg = %(sg_avg)s;
const cvData = %(cv_data)s;
const tempActual = %(temp_actual)s;
const tempTarget = %(temp_target)s;
const tempMin = %(temp_min)s;
const tempDrop = %(temp_drop)s;
const pwmAvg = %(pwm_avg)s;
const pwmMax = %(pwm_max)s;
const tmcOtpw = %(tmc_otpw)s;
const tmcOt = %(tmc_ot)s;
const linearSpeeds = %(linear_speeds)s;
const residenceTimes = %(residence_times)s;
const chartEvents = %(chart_events)s;
const triggerFlow = %(trigger_flow)s;
const verifyFlow = %(verify_flow)s;
const cvTriggerThreshold = %(cv_trigger)s;
const coldExtrusionHint = %(cold_extrusion_hint)s;
const coldOnsetFlow = %(cold_onset_flow)s;
const stressScores = %(stress_scores)s;
const intraDrift = %(intra_drift)s;
const extrasNames = %(extras_names)s;
const extrasData = %(extras_data)s;
const meltFitPct = %(melt_fit_pct)s;
const bdwPct = %(bdw_pct)s;
const meltPerFlow = %(melt_per_flow)s;
const bdwMinPct = %(bdw_min_pct)s;
const meltKneePct = %(melt_knee_pct)s;
const ownCsv = %(own_csv)s;

Chart.defaults.color = '#8b949e';
Chart.defaults.borderColor = '#30363d';
Chart.defaults.font.family = "'Manrope', system-ui, sans-serif";

// ─── Helper: pair y-values with their x-value (flow) ───────────────
// Without this, Chart.js would treat the flow array as categorical
// labels and space them evenly even when the actual flow steps are
// uneven (coarse=5, bisect=2, verify=1, etc.). Pairing as {x,y}
// objects lets us use a true linear axis and the points sit at their
// real flow values.
function asPoints(yArr) {
    const out = [];
    for (let i = 0; i < flows.length; i++) {
        if (yArr[i] !== null && yArr[i] !== undefined) {
            out.push({ x: flows[i], y: yArr[i] });
        }
    }
    return out;
}

// Build phase-boundary annotations for SG chart
const phaseAnnotations = {};
let lastPhase = null;
for (let i = 0; i < phases.length; i++) {
    if (phases[i] !== lastPhase && lastPhase !== null) {
        phaseAnnotations['phase_' + i] = {
            type: 'line',
            xMin: flows[i] - 0.5, xMax: flows[i] - 0.5,
            borderColor: 'rgba(110, 118, 129, 0.4)',
            borderWidth: 1, borderDash: [2, 4],
        };
    }
    lastPhase = phases[i];
}

// Build trigger event annotations
const eventAnnotations = {};
chartEvents.forEach((ev, i) => {
    eventAnnotations['ev_' + i] = {
        type: 'point',
        xValue: ev.flow,
        yValue: ev.value,
        backgroundColor: 'rgba(248, 81, 73, 0.85)',
        borderColor: '#f85149',
        radius: 6, borderWidth: 2,
    };
});

if (triggerFlow !== null) {
    eventAnnotations['triggerLine'] = {
        type: 'line',
        xMin: triggerFlow, xMax: triggerFlow,
        borderColor: '#f85149', borderWidth: 2, borderDash: [4,4],
        label: { display: true, content: 'first trigger',
                 position: 'start', backgroundColor: '#f85149',
                 color: 'white', font: {size: 10}, padding: 4 },
    };
}
if (verifyFlow !== null) {
    eventAnnotations['verifyLine'] = {
        type: 'line',
        xMin: verifyFlow, xMax: verifyFlow,
        borderColor: '#3fb950', borderWidth: 2,
        label: { display: true, content: 'verified',
                 position: 'end', backgroundColor: '#3fb950',
                 color: 'white', font: {size: 10}, padding: 4 },
    };
}
// Note: the cold-extrusion hint line is now drawn by buildZoneAnnotations()
// as part of the colored-zone system — that way the hint and the
// auto-detection share the same visual marker.

const sgAnnotations = Object.assign({}, phaseAnnotations, eventAnnotations);

// ─── Color zones (visual hint, no test impact) ───────────────────
// Three coloured background bands divide the flow range:
//   green:  start → cold-extrusion onset (or slip if no cold)
//   yellow: cold-extrusion onset → slip trigger
//   red:    slip trigger → max
// Each chart that wants the zones merges these annotations into its
// own annotations object.
function buildZoneAnnotations() {
    const z = {};
    // Determine zone boundaries based on which detections fired.
    const flowMin = %(flow_min_for_zones)s;
    const flowMax = %(flow_max_for_zones)s;
    // The user-supplied COLD_EXTRUSION_HINT takes precedence over the
    // auto-detection — visual observation of curling/sputtering is more
    // reliable than the heater-PWM heuristic for many setups (especially
    // strong heaters on weaker hotends, where curling happens long
    // before any heater-side stress is visible).
    const coldFlow = (coldExtrusionHint > 0) ? coldExtrusionHint
                                              : coldOnsetFlow;
    const slipFlow = triggerFlow;

    // Three states are possible:
    //  1. cold + slip: green / yellow / red
    //  2. only slip:   green / red (no yellow zone)
    //  3. only cold:   green / yellow (no red zone)
    //  4. neither:     all green
    let greenEnd, yellowEnd;
    if (coldFlow !== null && slipFlow !== null) {
        if (coldFlow < slipFlow) {
            greenEnd = coldFlow;
            yellowEnd = slipFlow;
        } else {
            // Cold detected after slip — unusual, treat as just slip
            greenEnd = slipFlow;
            yellowEnd = null;
        }
    } else if (slipFlow !== null) {
        greenEnd = slipFlow;
        yellowEnd = null;
    } else if (coldFlow !== null) {
        greenEnd = coldFlow;
        yellowEnd = flowMax;
    } else {
        greenEnd = flowMax;
        yellowEnd = null;
    }

    // Green zone (safe)
    z['zoneGreen'] = {
        type: 'box',
        xMin: flowMin, xMax: greenEnd,
        backgroundColor: 'rgba(63, 185, 80, 0.06)',
        borderWidth: 0,
        drawTime: 'beforeDatasetsDraw',
    };
    // Yellow zone (cold-extrusion suspected)
    if (yellowEnd !== null && yellowEnd > greenEnd) {
        z['zoneYellow'] = {
            type: 'box',
            xMin: greenEnd, xMax: yellowEnd,
            backgroundColor: 'rgba(210, 153, 34, 0.10)',
            borderWidth: 0,
            drawTime: 'beforeDatasetsDraw',
        };
    }
    // Red zone (slip / stall)
    const redStart = (yellowEnd !== null) ? yellowEnd : greenEnd;
    if (redStart < flowMax && slipFlow !== null) {
        z['zoneRed'] = {
            type: 'box',
            xMin: redStart, xMax: flowMax,
            backgroundColor: 'rgba(248, 81, 73, 0.08)',
            borderWidth: 0,
            drawTime: 'beforeDatasetsDraw',
        };
    }
    // Cold-onset marker line. Label depends on source:
    //   - if the user supplied COLD_EXTRUSION_HINT, it's based on
    //     visual observation
    //   - otherwise auto-detection from heater-PWM saturation
    if (coldFlow !== null) {
        const labelText = (coldExtrusionHint > 0)
            ? '❄ cold-ext. (user)'
            : '❄ cold-ext. (auto)';
        z['coldOnsetLine'] = {
            type: 'line',
            xMin: coldFlow, xMax: coldFlow,
            borderColor: 'rgba(210, 153, 34, 0.7)',
            borderWidth: 2, borderDash: [5, 3],
            label: { display: true,
                     content: labelText,
                     position: 'start',
                     backgroundColor: 'rgba(210, 153, 34, 0.85)',
                     color: 'white', font: {size: 10}, padding: 3 },
        };
    }
    return z;
}

const zoneAnnotations = buildZoneAnnotations();
Object.assign(sgAnnotations, zoneAnnotations);

// SG chart
new Chart(document.getElementById('sgChart'), {
    type: 'line',
    data: { datasets: [
        { label: 'P75', data: asPoints(sgP75),
          borderColor: 'rgba(88, 166, 255, 0.4)',
          backgroundColor: 'rgba(88, 166, 255, 0.08)',
          fill: '+1', borderDash: [3,3], pointRadius: 2 },
        { label: 'median', data: asPoints(sgMedian),
          borderColor: '#58a6ff', borderWidth: 3, pointRadius: 5, fill: false },
        { label: 'P25', data: asPoints(sgP25),
          borderColor: 'rgba(88, 166, 255, 0.4)',
          fill: false, borderDash: [3,3], pointRadius: 2 },
        { label: 'avg', data: asPoints(sgAvg),
          borderColor: 'rgba(139, 148, 158, 0.6)',
          borderDash: [6,3], fill: false, borderWidth: 1, pointRadius: 0 },
    ]},
    options: {
        responsive: true, interaction: { mode: 'index', intersect: false },
        scales: {
            x: { type: 'linear',
                 title: { display: true, text: 'Flow Rate (mm³/s)' } },
            y: { title: { display: true, text: '%(sg_label)s' } },
        },
        plugins: {
            legend: { position: 'top' },
            annotation: { annotations: sgAnnotations },
        },
    },
});

// Thermal chart — only if any thermal data exists
const hasThermalData = tempActual.some(v => v !== null);
if (hasThermalData) {
    document.getElementById('thermalChartCard').style.display = 'block';
    document.getElementById('thermalChartEmpty').style.display = 'none';

    const otpwPoints = [];
    const otPoints = [];
    for (let i = 0; i < flows.length; i++) {
        if (tmcOtpw[i]) otpwPoints.push({x: flows[i], y: 0.05});
        if (tmcOt[i])   otPoints.push({x: flows[i], y: 0.05});
    }

    const validTemps = tempActual.concat(tempMin)
        .filter(v => v !== null && v !== undefined);
    const validTargets = tempTarget.filter(v => v !== null && v !== undefined);
    let tempYMin = null, tempYMax = null;
    if (validTemps.length > 0) {
        tempYMin = Math.min(...validTemps) - 2;
        tempYMax = Math.max(...validTemps, ...validTargets) + 2;
    }

    // Build extra-thermistor datasets dynamically (e.g. chamber,
    // heatbreak). Each gets a unique colour from a small palette so
    // up to ~6 extras read clearly. Toggle via legend click.
    const extraColours = [
        '#56d364', // green
        '#79c0ff', // light blue
        '#d2a8ff', // purple
        '#ffa657', // orange
        '#a5a5a5', // grey
        '#ff7b72', // salmon
    ];
    const extraDatasets = [];
    for (let i = 0; i < extrasNames.length; i++) {
        const name = extrasNames[i];
        const data = extrasData[name] || [];
        // Skip sensors that produced no usable data
        if (!data.some(v => v !== null && v !== undefined)) continue;
        // Update tempYMin/tempYMax to cover the extra sensor range
        const valid = data.filter(v => v !== null && v !== undefined);
        if (valid.length > 0) {
            const lo = Math.min(...valid);
            const hi = Math.max(...valid);
            if (tempYMin === null || lo - 2 < tempYMin) tempYMin = lo - 2;
            if (tempYMax === null || hi + 2 > tempYMax) tempYMax = hi + 2;
        }
        extraDatasets.push({
            label: name,
            data: asPoints(data),
            borderColor: extraColours[i %% extraColours.length],
            borderWidth: 2,
            borderDash: [5, 4],
            fill: false,
            pointRadius: 3,
            yAxisID: 'yTemp',
        });
    }

    new Chart(document.getElementById('thermalChart'), {
        type: 'line',
        data: { datasets: [
            { label: 'temp target', data: asPoints(tempTarget),
              borderColor: 'rgba(88, 166, 255, 0.6)',
              borderDash: [4, 3], borderWidth: 1.5,
              fill: false, pointRadius: 0, yAxisID: 'yTemp' },
            { label: 'temp actual', data: asPoints(tempActual),
              borderColor: '#f85149', borderWidth: 2.5, pointRadius: 4,
              fill: false, yAxisID: 'yTemp' },
            { label: 'temp min', data: asPoints(tempMin),
              borderColor: 'rgba(255, 152, 0, 0.6)',
              borderDash: [2, 4], borderWidth: 1, fill: false,
              pointRadius: 2, yAxisID: 'yTemp' },
            ...extraDatasets,
            { label: 'heater PWM (avg)', data: asPoints(pwmAvg),
              borderColor: '#d29922',
              backgroundColor: 'rgba(210, 153, 34, 0.18)',
              fill: 'origin', borderWidth: 2, pointRadius: 3, yAxisID: 'yPwm' },
            { label: 'heater PWM (max)', data: asPoints(pwmMax),
              borderColor: 'rgba(230, 81, 0, 0.7)',
              borderDash: [3, 3], borderWidth: 1, pointRadius: 0,
              fill: false, yAxisID: 'yPwm' },
            { label: 'residence (s)', data: asPoints(residenceTimes),
              borderColor: '#3fb950', borderWidth: 2.5, pointRadius: 4,
              fill: false, yAxisID: 'yResidence' },
            { label: 'driver OTPW (≥120 °C)', data: otpwPoints,
              showLine: false, pointStyle: 'triangle', pointRadius: 8,
              backgroundColor: 'rgba(255, 152, 0, 0.9)',
              borderColor: '#e65100', yAxisID: 'yPwm' },
            { label: 'driver OT (≥150 °C)', data: otPoints,
              showLine: false, pointStyle: 'triangle', pointRadius: 10,
              backgroundColor: 'rgba(244, 67, 54, 0.9)',
              borderColor: '#b71c1c', yAxisID: 'yPwm' },
        ]},
        options: {
            responsive: true, interaction: { mode: 'index', intersect: false },
            scales: {
                x: { type: 'linear',
                     title: { display: true, text: 'Flow Rate (mm³/s)' } },
                yTemp: {
                    type: 'linear', position: 'left',
                    title: { display: true, text: 'Temperature (°C)' },
                    min: tempYMin, max: tempYMax,
                },
                yPwm: {
                    type: 'linear', position: 'right',
                    title: { display: true, text: 'PWM (0-1)' },
                    min: 0, max: 1.0,
                    grid: { drawOnChartArea: false },
                },
                yResidence: {
                    type: 'linear', position: 'right',
                    title: { display: true, text: 'Residence (s)' },
                    min: 0, grid: { drawOnChartArea: false },
                    afterFit: function(s) { s.width = 50; },
                },
            },
            plugins: {
                legend: { position: 'top' },
                tooltip: {
                    callbacks: {
                        afterBody: function(items) {
                            // Resolve which flow step we're on. With
                            // {x,y} data we can't trust dataIndex (it
                            // depends on the dataset's null-skipping),
                            // so look up by the actual x value (flow).
                            const flow = items[0].parsed.x;
                            const idx = flows.indexOf(flow);
                            if (idx < 0) return [];
                            const lines = [];
                            const drop = tempDrop[idx];
                            if (drop !== null && drop !== undefined && drop > 0.1) {
                                lines.push('Drop from target: ' + drop.toFixed(1) + ' °C');
                            }
                            const ls = linearSpeeds[idx];
                            const rt = residenceTimes[idx];
                            if (ls !== null && ls !== undefined) {
                                lines.push('Linear feed: ' + ls.toFixed(1) + ' mm/s');
                            }
                            if (rt !== null && rt !== undefined) {
                                let cls = '';
                                // Minimum hotend class needed to maintain
                                // this residence time — anything stronger
                                // (e.g. CHT > Volcano > V6) also works.
                                if (rt >= 1.5) cls = 'V6 or better';
                                else if (rt >= 0.6) cls = 'Volcano or better (CHT/Goliath OK)';
                                else if (rt >= 0.3) cls = 'CHT/HF or better';
                                else cls = 'Goliath / Bondtech-CHT only';
                                lines.push('Min. hotend: ' + cls);
                            }
                            // Thermal stress score (heuristic 0-100)
                            const ss = stressScores[idx];
                            if (ss !== null && ss !== undefined) {
                                let stressLabel = '';
                                if (ss < 30) stressLabel = 'low';
                                else if (ss < 60) stressLabel = 'moderate';
                                else stressLabel = 'high';
                                lines.push('Thermal stress: ' + ss.toFixed(0) +
                                           '/100 (' + stressLabel + ')');
                            }
                            // Intra-run SG drift (cold extrusion fingerprint)
                            const drift = intraDrift[idx];
                            if (drift !== null && drift !== undefined) {
                                let driftLabel = '';
                                if (Math.abs(drift) < 0.5) driftLabel = 'stable';
                                else if (drift > 3) driftLabel = 'load growing — cold suspected';
                                else if (drift > 1) driftLabel = 'load growing slightly';
                                else if (drift < -1) driftLabel = 'load decreasing';
                                else driftLabel = 'minimal change';
                                lines.push('Intra-run drift: ' +
                                           (drift > 0 ? '+' : '') +
                                           drift.toFixed(1) + '%% (' +
                                           driftLabel + ')');
                            }
                            return lines;
                        }
                    }
                },
                annotation: {
                    annotations: Object.assign({
                        v6Residence: {
                            type: 'line', yMin: 1.5, yMax: 1.5, yScaleID: 'yResidence',
                            borderColor: 'rgba(63, 185, 80, 0.4)', borderDash: [4, 4],
                            label: { display: true, content: 'V6 min 1.5s', position: 'start',
                                     backgroundColor: 'rgba(63, 185, 80, 0.7)',
                                     color: 'white', font: {size: 9}, padding: 2 },
                        },
                        volcanoResidence: {
                            type: 'line', yMin: 0.6, yMax: 0.6, yScaleID: 'yResidence',
                            borderColor: 'rgba(210, 153, 34, 0.4)', borderDash: [4, 4],
                            label: { display: true, content: 'Volcano min 0.6s', position: 'start',
                                     backgroundColor: 'rgba(210, 153, 34, 0.7)',
                                     color: 'white', font: {size: 9}, padding: 2 },
                        },
                        chtResidence: {
                            type: 'line', yMin: 0.3, yMax: 0.3, yScaleID: 'yResidence',
                            borderColor: 'rgba(248, 81, 73, 0.4)', borderDash: [4, 4],
                            label: { display: true, content: 'CHT min 0.3s', position: 'start',
                                     backgroundColor: 'rgba(248, 81, 73, 0.7)',
                                     color: 'white', font: {size: 9}, padding: 2 },
                        },
                    }, zoneAnnotations),
                },
            },
        },
    });
}

// CV chart — bars colored by severity
// Built as {x, y} pairs on a linear axis so bars sit at their real
// flow values, not in evenly-spaced category slots.
const cvBars = [];
const cvColors = [];
for (let i = 0; i < flows.length; i++) {
    const v = cvData[i];
    if (v === null || v === undefined) continue;
    cvBars.push({ x: flows[i], y: v });
    cvColors.push(
        v >= cvTriggerThreshold ? '#f85149'
        : v >= cvTriggerThreshold * 0.6 ? '#d29922'
        : '#3fb950');
}
new Chart(document.getElementById('cvChart'), {
    type: 'bar',
    data: { datasets: [
        { label: 'run-to-run CV (%%)', data: cvBars,
          backgroundColor: cvColors, borderWidth: 0 },
    ]},
    options: {
        responsive: true,
        scales: {
            x: { type: 'linear',
                 title: { display: true, text: 'Flow Rate (mm³/s)' } },
            y: { title: { display: true, text: 'CV (%%)' }, min: 0 },
        },
        plugins: {
            legend: { display: false },
            annotation: {
                annotations: Object.assign({
                    cvTrigger: {
                        type: 'line', yMin: cvTriggerThreshold, yMax: cvTriggerThreshold,
                        borderColor: 'rgba(248, 81, 73, 0.5)', borderDash: [4, 4],
                        label: { display: true,
                                 content: 'CV trigger ' + cvTriggerThreshold + '%%',
                                 position: 'end', backgroundColor: '#f85149',
                                 color: 'white', font: {size: 10}, padding: 3 },
                    },
                }, zoneAnnotations),
            },
        },
    },
});

(function () {
    const meltPts = flows.map((f, i) => (meltFitPct[i] === null ? null
        : {x: f, y: meltFitPct[i]})).filter(p => p !== null);
    const bdwPts = flows.map((f, i) => (bdwPct[i] === null ? null
        : {x: f, y: bdwPct[i]})).filter(p => p !== null);
    if (!meltPts.length && !bdwPts.length) {
        document.getElementById('meltChartCard').style.display = 'none';
        document.getElementById('meltChartEmpty').style.display = '';
        return;
    }
    const ann = {};
    if (meltPts.length) {
        ann.knee = { type: 'line', scaleID: 'y', value: meltKneePct,
            borderColor: 'rgba(56, 139, 253, 0.6)', borderDash: [6, 4],
            label: { display: true, content: 'melt limit ' + meltKneePct + '%%',
                     position: 'start', backgroundColor: '#388bfd',
                     color: 'white', font: {size: 10}, padding: 3 } };
    }
    if (bdwPts.length && bdwMinPct > 0) {
        ann.feed = { type: 'line', scaleID: 'y', value: bdwMinPct,
            borderColor: 'rgba(248, 81, 73, 0.6)', borderDash: [6, 4],
            label: { display: true, content: 'feed check ' + bdwMinPct + '%%',
                     position: 'end', backgroundColor: '#f85149',
                     color: 'white', font: {size: 10}, padding: 3 } };
    }
    new Chart(document.getElementById('meltChart'), {
        type: 'line',
        data: { datasets: [
            { label: 'Heat vs low-flow line (%%)', data: meltPts,
              borderColor: '#388bfd', backgroundColor: '#388bfd',
              tension: 0, pointRadius: 3 },
            { label: 'BDWidth feed (%%)', data: bdwPts,
              borderColor: '#d29922', backgroundColor: '#d29922',
              tension: 0, pointRadius: 3 },
        ]},
        options: {
            responsive: true,
            scales: {
                x: { type: 'linear',
                     title: { display: true, text: 'Flow Rate (mm³/s)' } },
                y: { title: { display: true, text: '%% of expected' },
                     suggestedMin: 0, suggestedMax: 115 },
            },
            plugins: { annotation: { annotations: ann } },
        },
    });
})();

// ─── Compare runs: load other CSVs of this test and overlay them ───
(function () {
    const PALETTE = ['#388bfd', '#f0883e', '#3fb950', '#a371f7',
                     '#f85149', '#d29922', '#39c5cf', '#db61a2'];
    const runs = [];          // {name, meta, rows, label, color}
    const charts = {};
    const msg = document.getElementById('cmpMsg');

    function parseCsv(text) {
        const meta = {};
        const lines = text.split(/\\r?\\n/);
        let header = null;
        const rows = [];
        for (const line of lines) {
            if (!line.trim()) continue;
            if (line.startsWith('#')) {
                const m = line.match(/^#\\s*([A-Za-z0-9_]+):\\s*(.*)$/);
                if (m) meta[m[1]] = m[2].trim();
                continue;
            }
            const cells = line.split(',');
            if (!header) { header = cells; continue; }
            const r = {};
            header.forEach((h, i) => { r[h] = cells[i]; });
            rows.push(r);
        }
        if (!header || header.indexOf('flow_mm3s') < 0
                || header.indexOf('pwm_avg') < 0) {
            throw new Error('not a TMC flow test CSV (no flow_mm3s / '
                            + 'pwm_avg columns)');
        }
        return {meta: meta, rows: rows};
    }
    const num = v => (v === undefined || v === '' || v === null)
        ? null : (isNaN(+v) ? null : +v);

    // Same line-fit melt analysis as the plugin (works on older CSVs
    // too): duty = a + b·flow through steps 2–5, drift-corrected with
    // idle_pwm when present; %% of that line per step.
    function analyse(run) {
        const rows = run.rows.filter(r => num(r.flow_mm3s) && num(r.pwm_avg) !== null);
        const idle0 = num(run.meta.idle_pwm_start);
        const pts = rows.map(r => {
            let pwm = num(r.pwm_avg);
            const ip = num(r.idle_pwm);
            if (idle0 !== null && ip !== null) pwm -= (ip - idle0);
            return {flow: num(r.flow_mm3s), pwm: pwm, row: r};
        });
        const fit = pts.slice(1, 5);
        run.fitPct = [];
        run.knee = null; run.lastOk = null;
        if (fit.length >= 3 && pts.length >= fit.length + 2) {
            const n = fit.length;
            const sx = fit.reduce((s, p) => s + p.flow, 0);
            const sy = fit.reduce((s, p) => s + p.pwm, 0);
            const sxx = fit.reduce((s, p) => s + p.flow * p.flow, 0);
            const sxy = fit.reduce((s, p) => s + p.flow * p.pwm, 0);
            const den = n * sxx - sx * sx;
            if (den > 0) {
                const b = (n * sxy - sx * sy) / den;
                const a = (sy - b * sx) / n;
                if (b > 0) {
                    run.fitPct = pts.map(p => ({x: p.flow,
                        y: (p.pwm - a) / (b * p.flow) * 100}));
                    const fitMax = Math.max(...fit.map(p => p.flow));
                    const later = run.fitPct.filter(p => p.x > fitMax + 0.001)
                        .sort((u, v) => u.x - v.x);
                    const k = later.find(p => p.y < meltKneePct);
                    if (k) {
                        run.knee = k.x;
                        const ok = run.fitPct.filter(p => p.x < k.x - 0.001
                                                      && p.y >= meltKneePct);
                        run.lastOk = ok.length ? Math.max(...ok.map(p => p.x)) : null;
                    }
                }
            }
        }
        run.bdw = rows.filter(r => num(r.bdw_ratio) !== null)
            .map(r => ({x: num(r.flow_mm3s), y: num(r.bdw_ratio) * 100}));
        run.sg = rows.filter(r => num(r.sg_median) !== null)
            .map(r => ({x: num(r.flow_mm3s), y: num(r.sg_median)}))
            .sort((u, v) => u.x - v.x);
        run.fitPct.sort((u, v) => u.x - v.x);
        run.bdw.sort((u, v) => u.x - v.x);
        const tt = (run.meta.hotend_temp || '').match(/target\\s+([\\d.]+)/);
        const t0 = rows.length ? num(rows[0].temp_target) : null;
        run.temp = tt ? +tt[1] : t0;
        run.label = (run.temp !== null ? run.temp.toFixed(0) + ' °C' : '? °C')
            + ' · ' + (run.meta.timestamp || run.name);
    }

    function matOf(run) {
        const m = (run.meta.material || '').trim();
        return (m && m !== 'not specified') ? m : '';
    }
    function dataset(run, key) {
        return { label: run.label, data: run[key], borderColor: run.color,
                 backgroundColor: run.color, tension: 0, pointRadius: 3,
                 borderDash: run.own ? [] : [5, 3] };
    }
    function draw(id, key, yTitle, ann) {
        const ds = runs.filter(r => r[key] && r[key].length).map(r => dataset(r, key));
        const card = document.getElementById(id).closest('.chart-card');
        card.style.display = ds.length ? '' : 'none';
        if (charts[id]) charts[id].destroy();
        if (!ds.length) return;
        charts[id] = new Chart(document.getElementById(id), {
            type: 'line', data: { datasets: ds },
            options: { responsive: true,
                scales: { x: { type: 'linear', title: { display: true,
                                  text: 'Flow Rate (mm³/s)' } },
                          y: { title: { display: true, text: yTitle } } },
                plugins: { annotation: { annotations: ann || {} } } } });
    }
    function render() {
        draw('cmpMeltChart', 'fitPct', '%% of low-flow line',
             { knee: { type: 'line', scaleID: 'y', value: meltKneePct,
                       borderColor: 'rgba(56,139,253,0.5)', borderDash: [6, 4] } });
        draw('cmpBdwChart', 'bdw', 'BDWidth feed (%%)',
             bdwMinPct > 0 ? { feed: { type: 'line', scaleID: 'y',
                       value: bdwMinPct, borderColor: 'rgba(248,81,73,0.5)',
                       borderDash: [6, 4] } } : {});
        draw('cmpSgChart', 'sg', 'SG median');
        const tb = document.getElementById('cmpTable');
        tb.innerHTML = runs.map(r =>
            '<tr><td><span style="display:inline-block;width:10px;height:10px;'
            + 'border-radius:50%%;background:' + r.color + '"></span> '
            + (r.own ? '<b>this report</b>' : r.name) + '</td>'
            + '<td>' + (matOf(r) || '-') + '</td>'
            + '<td>' + (r.temp !== null ? r.temp.toFixed(0) + ' °C' : '-') + '</td>'
            + '<td>' + (r.meta.timestamp || '-') + '</td>'
            + '<td>' + (r.lastOk !== null ? r.lastOk.toFixed(1) : (r.fitPct.length ? 'above ' + Math.max(...r.fitPct.map(p => p.x)).toFixed(0) : '-')) + '</td>'
            + '<td>' + (r.meta.result_max_flow || '-') + '</td>'
            + '<td>' + (r.meta.result_quality || '-') + '</td>'
            + '<td>' + (r.own ? '' : '<button class="cmp-rm" data-name="'
                + r.name.replace(/"/g, '') + '">remove</button>') + '</td></tr>'
        ).join('');
        tb.querySelectorAll('.cmp-rm').forEach(b => b.onclick = () => {
            const i = runs.findIndex(r => !r.own && r.name === b.dataset.name);
            if (i >= 0) runs.splice(i, 1);
            render();
        });
    }
    function add(name, text, own) {
        const p = parseCsv(text);
        const run = { name: name, meta: p.meta, rows: p.rows, own: own,
                      text: text.trim(),
                      color: PALETTE[runs.length %% PALETTE.length] };
        analyse(run);
        runs.push(run);
    }

    try { add('this report', ownCsv, true); }
    catch (e) { msg.textContent = 'Could not read this report\\'s own data: ' + e.message; }

    document.getElementById('cmpFile').addEventListener('change', ev => {
        const files = Array.from(ev.target.files || []);
        const notes = [];
        let pending = files.length;
        if (!pending) return;
        files.forEach(f => {
            const rd = new FileReader();
            rd.onload = () => {
                const text = String(rd.result);
                try {
                    const p = parseCsv(text);
                    const own = runs.find(r => r.own);
                    const ts = p.meta.timestamp;
                    if (own && (text.trim() === own.text
                            || (ts && ts === own.meta.timestamp))) {
                        notes.push(f.name + ': this is the run that produced '
                                   + 'this report — skipped.');
                    } else if (runs.some(r => !r.own && (r.text === text.trim()
                            || (ts && r.meta.timestamp === ts)))) {
                        notes.push(f.name + ': already loaded — skipped.');
                    } else {
                        add(f.name, text, false);
                        const run = runs[runs.length - 1];
                        if (own && matOf(run) && matOf(own)
                                && matOf(run) !== matOf(own)) {
                            notes.push(f.name + ': different material ('
                                       + matOf(run) + ' vs ' + matOf(own)
                                       + ').');
                        } else if (own && (!matOf(run) || !matOf(own))) {
                            notes.push(f.name + ': material not recorded '
                                       + 'in both runs — check it is the '
                                       + 'same filament.');
                        }
                        if (own && run.meta.driver && own.meta.driver
                                && run.meta.driver !== own.meta.driver) {
                            notes.push(f.name + ': different driver/extruder ('
                                       + run.meta.driver + ').');
                        }
                    }
                } catch (e) {
                    notes.push(f.name + ': ' + e.message);
                }
                if (--pending === 0) {
                    msg.textContent = notes.join('  ');
                    ev.target.value = '';
                    render();
                }
            };
            rd.readAsText(f);
        });
    });
    window._cmpRender = render;
})();
</script>
</body>
</html>"""
        # Build trigger-event annotations for the chart. We attach
        # each event to the LAST results-array index whose flow matches
        # (since re-tests replace earlier entries, this is well-defined
        # except for verify-fail where the failed verify still sits in
        # results — in that case we attach to the verify entry directly).
        chart_events = []
        if final_result is not None:
            for ev in (final_result.get('trigger_events') or []):
                # Find the matching results index (search in reverse so
                # we hit the last/most-recent entry for this flow+phase).
                idx = None
                for i in range(len(results) - 1, -1, -1):
                    rr = results[i]
                    if (rr.get('phase') == ev['phase']
                            and abs(rr.get('flow', 0) - ev['flow'])
                            < 0.01):
                        idx = i
                        break
                if idx is None:
                    continue
                chart_events.append({
                    'idx': idx,
                    'flow': ev['flow'],
                    'phase': ev['phase'],
                    'kind': ev['kind'],
                    'reason': ev['reason'],
                    'cv': ev.get('cv'),
                    'iqr': ev.get('iqr'),
                })

        # Healthy-range bands for the chart (drawn as horizontal
        # background regions if baseline_stats is available).
        healthy_iqr = None
        if final_result is not None:
            bs = final_result.get('baseline_stats')
            if bs and 'iqr_min' in bs and 'iqr_max' in bs:
                healthy_iqr = {
                    'min': bs['iqr_min'],
                    'max': max(bs['iqr_max'],
                               bs['iqr_median'] * 1.5),
                    'median': bs['iqr_median'],
                }

        # First trigger flow + verify flow for chart annotations
        first_trigger_flow = None
        for ev in (trigger_events or []):
            if ev.get('kind') == 'trigger':
                first_trigger_flow = ev['flow']
                break
        verify_flow = None
        if max_safe is not None:
            verify_flow = max_safe

        # Cold-extrusion onset detection (purely visual, no test impact).
        # Returns the flow at which heater saturation + temperature drop
        # were first observed together. May be None.
        cold_onset_flow = self._detect_cold_extrusion_onset(results)

        # Continuous thermal stress score per step (0-100). Used to colour
        # the charts and to drive the cold-extrusion onset detection above.
        # See _compute_thermal_stress_per_step for the scoring breakdown.
        stress_scores_pairs = self._compute_thermal_stress_per_step(results)
        stress_scores = [s for _, s in stress_scores_pairs]
        max_stress = max((s for s in stress_scores if s is not None),
                         default=None)
        max_stress_flow = None
        if max_stress is not None:
            for flow, score in stress_scores_pairs:
                if score == max_stress:
                    max_stress_flow = flow
                    break

        # The chart's coloured background zones use these boundaries.
        # min_flow / max_flow span the X-axis of the data we have.
        if results:
            min_flow_for_zones = min(r['flow'] for r in results)
            max_flow_for_zones = max(r['flow'] for r in results)
        else:
            min_flow_for_zones = 0
            max_flow_for_zones = 100

        rendered = html % {
            'timestamp': meta.get('timestamp', '-'),
            'version': MODULE_VERSION,
            'driver_label': driver_label_full,
            'measurement_count': len(results),
            'meta_html': meta_html,
            'tmc_settings_html': tmc_settings_html,
            'decision_trail_html': decision_trail_html,
            'hero_html': hero_html,
            'insights_html': insights_html,
            'timeline_html': timeline_html,
            'sg_label': sg_label,
            'data_table': table,
            'flows': json.dumps(flows),
            'phases': json.dumps(phases),
            'sg_median': json.dumps(sg_median),
            'sg_p25': json.dumps(sg_p25),
            'sg_p75': json.dumps(sg_p75),
            'sg_avg': json.dumps(sg_avg),
            'cv_data': json.dumps(cv_data),
            'melt_fit_pct': json.dumps(
                [self._step_extras(r, meta)['fit_pct']
                 for r in results_for_charts]),
            'bdw_pct': json.dumps(
                [self._step_extras(r, meta)['bdw_pct']
                 for r in results_for_charts]),
            'melt_per_flow': json.dumps(
                [self._step_extras(r, meta)['melt_per_flow']
                 for r in results_for_charts]),
            'bdw_min_pct': json.dumps(
                (getattr(self, '_bdw_min_ratio', 0.0) or 0.0) * 100),
            'melt_knee_pct': json.dumps(self.MELT_KNEE_FRAC * 100),
            'own_csv': json.dumps(self._csv_text(results, meta)).replace(
                '</', '<\\/'),
            'cv_trigger': json.dumps(cv_trigger_threshold),
            'chart_events': json.dumps(chart_events),
            'trigger_flow': json.dumps(first_trigger_flow),
            'verify_flow': json.dumps(verify_flow),
            'cold_onset_flow': json.dumps(cold_onset_flow),
            'stress_scores': json.dumps(stress_scores),
            'flow_min_for_zones': json.dumps(min_flow_for_zones),
            'flow_max_for_zones': json.dumps(max_flow_for_zones),
            'cold_extrusion_hint': json.dumps(
                meta.get('cold_extrusion_hint', 0)),
            # Thermal data — see _aggregate_thermal_samples
            'temp_actual': json.dumps(temp_actual),
            'temp_target': json.dumps(temp_target),
            'temp_min': json.dumps(temp_min),
            'temp_drop': json.dumps(temp_drop),
            'pwm_avg': json.dumps(pwm_avg),
            'pwm_max': json.dumps(pwm_max),
            'tmc_otpw': json.dumps(tmc_otpw),
            'tmc_ot': json.dumps(tmc_ot),
            'intra_drift': json.dumps(intra_drift),
            # Filament-speed / residence-time helpers
            'linear_speeds': json.dumps(linear_speeds),
            'residence_times': json.dumps(residence_times),
            # Extra thermistors (chamber, heatbreak, etc.). Each one
            # becomes a separate dotted line on the thermal chart.
            'extras_names': json.dumps(extras_seen),
            'extras_data': json.dumps(extras_per_step),
        }
        with open(path, 'w') as f:
            f.write(rendered)

    def _csv_text(self, results, meta):
        """This run's CSV as a string (embedded in the HTML report so the
        Compare tab can treat it exactly like a loaded CSV)."""
        import tempfile
        fd, tmp = tempfile.mkstemp(suffix='.csv')
        os.close(fd)
        try:
            self._write_csv(tmp, results, meta)
            with open(tmp) as f:
                return f.read()
        except Exception:
            return ''
        finally:
            try:
                os.remove(tmp)
            except OSError:
                pass

    def _save_report(self, results, meta, timestamp, limit_reason,
                     no_html, gcmd=None, announce=True,
                     final_result=None):
        """Save CSV and HTML report. Safe to call repeatedly."""
        if not results:
            return
        if final_result:
            # Recorded in the CSV header too, so the HTML comparison can
            # show other runs' results from their CSV alone.
            if final_result.get('max_safe') is not None:
                meta['result_max_flow'] = '%.1f' % final_result['max_safe']
            if final_result.get('quality'):
                meta['result_quality'] = final_result['quality']
        try:
            os.makedirs(self.output_dir, exist_ok=True)
            csv_path = os.path.join(
                self.output_dir, 'tmc_flow_%s.csv' % timestamp)
            self._write_csv(csv_path, results, meta)
            if announce and gcmd is not None:
                gcmd.respond_info("CSV saved: %s" % csv_path)
            if not no_html:
                html_path = os.path.join(
                    self.output_dir,
                    'tmc_flow_%s.html' % timestamp)
                self._write_html(html_path, results, meta, limit_reason,
                                 final_result=final_result)
                if announce and gcmd is not None:
                    gcmd.respond_info("HTML saved: %s" % html_path)
        except Exception as e:
            if gcmd is not None:
                gcmd.respond_info(
                    "Warning: report write failed: %s" % e)
            logging.exception("tmc_flow_test: report write failed")

    # ─── Single flow step measurement ───────────────────────────────

    # ─── Thermal data collection ───────────────────────────────────
    # Read hotend temperature + heater PWM + TMC driver thermal flags.
    # These values are LOGGED ONLY (CSV columns) — no triggers, no
    # auto-detection, no behavioural change. They give the user
    # context to spot thermal issues during analysis (e.g. a slip
    # trigger that coincided with a temp drop is suspicious).

    def _get_thermal_snapshot(self):
        """Return a dict with current hotend + driver thermal state.
        All keys may be None if the relevant subsystem isn't available
        — treated as 'no data' downstream.
        """
        snap = {
            'temp_actual': None,    # current hotend temp (°C)
            'temp_target': None,    # target hotend temp (°C)
            'heater_pwm': None,     # current heater PWM, 0.0-1.0
            'tmc_otpw': None,       # TMC over-temperature warning flag (≥120 °C)
            'tmc_ot': None,         # TMC over-temperature error flag (≥150 °C)
            'extras': {},           # name -> temperature for extra thermistors
        }

        # Hotend temp + PWM via Klipper extruder/heater objects
        try:
            extruder = self.printer.lookup_object('extruder')
            heater = extruder.get_heater()
            cur_t, target_t = heater.get_temp(self.reactor.monotonic())
            snap['temp_actual'] = float(cur_t)
            snap['temp_target'] = float(target_t)
            # last_pwm_value is always present; for both PID and MPC
            # control methods Klipper exposes the actual fired PWM
            try:
                pwm = heater.last_pwm_value
                if pwm is not None:
                    snap['heater_pwm'] = float(pwm)
            except AttributeError:
                pass
        except Exception:
            pass

        # TMC driver thermal warning flags via DRV_STATUS
        if self.tmc is not None:
            for field, key in (('otpw', 'tmc_otpw'), ('ot', 'tmc_ot')):
                try:
                    val = self.tmc.fields.get_field(field)
                    if val is not None:
                        snap[key] = int(val)
                except (KeyError, AttributeError):
                    pass

        # Extra thermistors (chamber, heatbreak, ambient, etc.).
        # Each one is looked up by its full Klipper object name.
        # Resolution failures are silent so a printer.cfg change that
        # removes a thermistor doesn't break the test mid-run.
        if self.extra_thermistors:
            now = self.reactor.monotonic()
            for name in self.extra_thermistors:
                try:
                    obj = self.printer.lookup_object(name)
                    # temperature_sensor: has .last_temp; generic
                    # heaters: get_temp(eventtime) returns (cur, target).
                    if hasattr(obj, 'get_temp'):
                        result = obj.get_temp(now)
                        # heater-style returns tuple (cur, target)
                        cur = (result[0] if isinstance(result, tuple)
                               else result)
                        snap['extras'][name] = float(cur)
                    elif hasattr(obj, 'last_temp'):
                        snap['extras'][name] = float(obj.last_temp)
                except Exception:
                    pass  # not a temp sensor or doesn't exist

        return snap

    @staticmethod
    def _aggregate_thermal_samples(samples):
        """Reduce a list of per-rep thermal snapshots to summary stats.
        Each input is a dict from _get_thermal_snapshot. Returns aggregate
        across all reps with min/max/avg where meaningful.
        """
        agg = {
            'temp_start': None, 'temp_end': None,
            'temp_min': None, 'temp_max': None, 'temp_avg': None,
            'temp_target': None, 'temp_drop': None,
            'pwm_min': None, 'pwm_max': None, 'pwm_avg': None,
            'tmc_otpw_any': 0, 'tmc_ot_any': 0,
            'extras': {},
        }
        if not samples:
            return agg

        # Temperature aggregation
        temps = [s['temp_actual'] for s in samples
                 if s.get('temp_actual') is not None]
        if temps:
            agg['temp_start'] = temps[0]
            agg['temp_end'] = temps[-1]
            agg['temp_min'] = min(temps)
            agg['temp_max'] = max(temps)
            agg['temp_avg'] = sum(temps) / len(temps)

        targets = [s['temp_target'] for s in samples
                   if s.get('temp_target') is not None]
        if targets:
            agg['temp_target'] = targets[-1]  # last target wins
            if temps:
                agg['temp_drop'] = max(targets[-1] - min(temps), 0.0)

        # PWM aggregation
        pwms = [s['heater_pwm'] for s in samples
                if s.get('heater_pwm') is not None]
        if pwms:
            agg['pwm_min'] = min(pwms)
            agg['pwm_max'] = max(pwms)
            agg['pwm_avg'] = sum(pwms) / len(pwms)

        # Driver thermal warning flags — flagged if ANY sample saw it
        agg['tmc_otpw_any'] = int(any(
            s.get('tmc_otpw') for s in samples
            if s.get('tmc_otpw') is not None))
        agg['tmc_ot_any'] = int(any(
            s.get('tmc_ot') for s in samples
            if s.get('tmc_ot') is not None))

        # Extra thermistors: aggregate min/max/avg per sensor name
        extras_agg = {}
        # Collect all extra-sensor names that appeared in any sample
        all_extras = set()
        for s in samples:
            extras = s.get('extras') or {}
            all_extras.update(extras.keys())
        for name in all_extras:
            vals = [s['extras'][name] for s in samples
                    if s.get('extras')
                    and name in s['extras']
                    and s['extras'][name] is not None]
            if vals:
                extras_agg[name] = {
                    'min': min(vals),
                    'max': max(vals),
                    'avg': sum(vals) / len(vals),
                }
        agg['extras'] = extras_agg

        return agg

    @staticmethod
    def _compute_thermal_stress_per_step(results):
        """Compute a 0-100 thermal stress score for each step.

        The score combines five signals to estimate how close the
        hotend system is to a cold-extrusion regime:

        - PWM absolute level (0-30 points): how hard the heater is working
        - Temperature drop (0-25 points): how much actual temp falls below target
        - PWM rising trend vs baseline (0-15 points): heater taking on more work
        - PWM peak proximity (0-10 points): hits at saturation
        - Intra-run SG drift (0-20 points): SG drifting load-ward DURING a run,
          a strong cold-extrusion fingerprint (see _compute_intra_run_trend)

        This is a HEURISTIC. It does not replace direct visual inspection
        of the strand at the nozzle, and it does not trigger anything in
        the test. Its only purpose is to colour-code the chart so the
        user gets a visual hint of WHERE in the flow range thermal stress
        starts to build.

        Returns a list of (flow, score) tuples — one per step, in
        results order. Score is None if no thermal data was available
        for that step.
        """
        # Build PWM baseline from the first 3 coarse steps
        baseline_pwms = []
        for r in results:
            if r.get('phase', 'coarse') != 'coarse':
                continue
            t = r.get('thermal') or {}
            p = t.get('pwm_avg')
            if p is not None:
                baseline_pwms.append(p)
                if len(baseline_pwms) >= 3:
                    break
        baseline = (sum(baseline_pwms) / len(baseline_pwms)
                    if baseline_pwms else None)

        scores = []
        for r in results:
            flow = r['flow']
            t = r.get('thermal') or {}
            ir = r.get('intra_run') or {}
            pwm_avg = t.get('pwm_avg')
            pwm_max = t.get('pwm_max')
            drop = t.get('temp_drop')
            mean_drift_load = ir.get('mean_drift_load') if ir else None
            consistent_growing = ir.get('consistent_growing_load') \
                if ir else 0

            # If we have NEITHER thermal nor intra-run, mark as no data
            if pwm_avg is None and drop is None and mean_drift_load is None:
                scores.append((flow, None))
                continue

            score = 0.0

            # 1) PWM absolute level (0-30)
            if pwm_avg is not None:
                # PWM=0.30 ≈ 0 stress, PWM=0.95 ≈ ~30
                score += min(30, max(0, (pwm_avg - 0.3) * 46))

            # 2) Temperature drop (0-25)
            if drop is not None:
                # 0 °C = 0, 5 °C = 25
                score += min(25, max(0, drop * 5))

            # 3) PWM rising trend (0-15)
            if pwm_avg is not None and baseline is not None:
                delta = pwm_avg - baseline
                # +0.20 PWM = +15 stress
                score += min(15, max(0, delta * 75))

            # 4) PWM peak saturation (0-10)
            if pwm_max is not None:
                if pwm_max >= 0.99:
                    score += 10
                elif pwm_max >= 0.85:
                    score += 5

            # 5) Intra-run SG drift toward higher load (0-20)
            # Positive drift_load = SG drifted toward higher-load
            # direction during the run = filament getting harder to push
            # = classic cold-extrusion fingerprint
            if mean_drift_load is not None:
                if mean_drift_load >= 5.0:
                    score += 20  # very strong gradient
                elif mean_drift_load >= 3.0:
                    score += 15
                elif mean_drift_load >= 1.5:
                    score += 10
                elif mean_drift_load >= 0.5:
                    score += 5
                # Bonus if multiple runs all show growing load
                # (consistency check — rules out a single fluke)
                if consistent_growing >= 3:
                    score += 5

            scores.append((flow, min(100, round(score, 1))))
        return scores

    @staticmethod
    def _detect_cold_extrusion_onset(results):
        """Find the flow where thermal stress first crosses the
        'yellow' threshold (>= 30). This is purely informational and
        does not affect any test logic."""
        for flow, score in TMCFlowTest._compute_thermal_stress_per_step(
                results):
            if score is not None and score >= 30:
                return flow
        return None

    @staticmethod
    def _compute_intra_run_trend(run_sg, sg2_driver=True):
        """Compute time-trend metrics from one run's SG samples.

        Within a single 5-second run, samples are recorded in time
        order (~50 Hz). If SG drifts during the run it tells us
        something the per-run-median can't: a steady slope is the
        signature of cold extrusion (load grows as filament cools
        further), while a sudden jump near the end is motor slip
        (mechanical limit hit).

        Returns dict with:
            slope          — linear regression slope of SG vs sample index
            early_avg      — mean of first third of samples
            late_avg       — mean of last third
            drift_pct      — (late_avg - early_avg) / median * 100
            drift_pct_load — same, oriented so positive = load increasing
                             (negative SG drift on SG2 drivers, positive on SG4)
            n_samples      — count after warmup-skip
        Returns None if too few samples.
        """
        n_raw = len(run_sg)
        if n_raw < 20:  # need at least ~1 second of samples
            return None
        # Skip first 10 % as warmup (motor speed-up, sensor settling)
        skip = max(2, n_raw // 10)
        samples = run_sg[skip:]
        n = len(samples)
        if n < 10:
            return None

        # Linear regression slope (SG units per sample-index)
        mean_x = (n - 1) / 2.0
        mean_y = sum(samples) / n
        num = sum((i - mean_x) * (s - mean_y)
                  for i, s in enumerate(samples))
        den = sum((i - mean_x) ** 2 for i in range(n))
        slope = (num / den) if den > 0 else 0.0

        # Early vs late thirds
        third = max(1, n // 3)
        early = samples[:third]
        late = samples[-third:]
        early_avg = sum(early) / len(early)
        late_avg = sum(late) / len(late)
        drift = late_avg - early_avg
        median_approx = mean_y if mean_y > 1 else 1.0
        drift_pct = (drift / median_approx * 100.0)

        # Orient drift so POSITIVE = load INCREASING during the run
        # (which is the cold-extrusion signature).
        # SG2 drivers (TMC5160): SG goes DOWN as load increases
        # SG4 drivers (TMC2209/2240): SG goes UP as load increases
        if sg2_driver:
            drift_pct_load = -drift_pct  # invert
        else:
            drift_pct_load = drift_pct

        return {
            'slope': slope,
            'early_avg': early_avg,
            'late_avg': late_avg,
            'drift_pct': drift_pct,
            'drift_pct_load': drift_pct_load,
            'n_samples': n,
        }

    @staticmethod
    def _aggregate_intra_run_trends(per_run_trends):
        """Reduce per-run trend dicts to a step-level aggregate.

        Returns aggregate dict or None if no valid per-run data.
        Keys:
            mean_drift_load     — average load-direction drift across runs
                                   (positive = load increased during runs)
            max_drift_load      — single largest load-increasing drift
            consistent_negative — count of runs where load clearly increased
                                   (drift_load > 1 %), suggests cold-extrusion
            n_valid_runs        — number of runs that had usable data
        """
        valid = [t for t in per_run_trends if t is not None]
        if not valid:
            return None
        load_drifts = [t['drift_pct_load'] for t in valid]
        return {
            'mean_drift_load': sum(load_drifts) / len(load_drifts),
            'max_drift_load': max(load_drifts),
            'consistent_growing_load':
                sum(1 for d in load_drifts if d > 1.0),
            'n_valid_runs': len(valid),
        }

    # ─── Data-collection helpers (no triggers) ─────────────────────

    def _measure_idle_baseline(self, gcmd, settle_s, idle_s, label):
        """Hold temperature without extruding, then average the heater
        duty. This is the zero-flow point for the heat-per-mm³ columns.
        Returns {'pwm', 'temp', 't_mid'} or None."""
        gcmd.respond_info(
            "\n>>> Idle heater baseline (%s): settling %.0f s, then "
            "measuring %.0f s — no extrusion <<<" % (label, settle_s, idle_s))
        self.gcode.run_script_from_command("M400")
        if settle_s > 0:
            self.gcode.run_script_from_command(
                "G4 P%d" % int(settle_s * 1000))
        t_start = self.reactor.monotonic()
        pwms, temps = [], []
        n = max(1, int(idle_s / 0.5))
        for _ in range(n):
            snap = self._get_thermal_snapshot()
            if snap.get('heater_pwm') is not None:
                pwms.append(snap['heater_pwm'])
            if snap.get('temp_actual') is not None:
                temps.append(snap['temp_actual'])
            self.gcode.run_script_from_command("G4 P500")
        if not pwms:
            gcmd.respond_info("  Idle baseline: no heater PWM available.")
            return None
        res = {'pwm': sum(pwms) / len(pwms),
               'temp': (sum(temps) / len(temps)) if temps else 0.0,
               't_mid': (t_start + self.reactor.monotonic()) / 2.0}
        gcmd.respond_info(
            "  Idle heater duty (%s): %.1f %% at %.1f°C (%d samples)"
            % (label, res['pwm'] * 100, res['temp'], len(pwms)))
        return res

    def _resolve_bdwidth(self, gcmd, param):
        """Find a bdwidth sensor (markniu/bdwidth). Reads its internal
        motion accumulator, which is not part of its public status, so
        anything missing simply disables the motion columns."""
        if param.lower() in ('none', 'off', '0'):
            return None
        obj = None
        name = None
        for full_name, o in self.printer.lookup_objects():
            if not full_name.startswith('bdwidth'):
                continue
            short = full_name.split(None, 1)[-1]
            if param.lower() == 'auto' or param in (short, full_name):
                obj, name = o, full_name
                break
        if obj is None:
            if param.lower() != 'auto':
                gcmd.respond_info("BDWidth '%s' not found — motion "
                                  "columns disabled." % param)
            return None
        if not hasattr(obj, 'actual_total_move'):
            gcmd.respond_info("BDWidth '%s' has no motion accumulator "
                              "(module version?) — motion columns "
                              "disabled." % name)
            return None
        mode = str(getattr(obj, 'is_active', '?'))
        info = {'obj': obj, 'name': name,
                'short': name.split(None, 1)[-1],
                'orig_mode': mode, 'mode': mode,
                'counts_per_mm': float(getattr(obj, 'linear_motion', 42.8)),
                'sample_time': float(getattr(obj, 'sample_time', 1.0))}
        gcmd.respond_info("BDWidth motion data: using '%s' (mode=%s)."
                          % (name, mode))
        return info

    def _bdwidth_prepare(self, gcmd, bdw):
        """Switch the sensor to motion-only for the test: width
        compensation would change the extrusion multiplier (so the
        delivered flow differs from the commanded one), and with the
        sensor disabled its motion counter does not update at all."""
        # bdwidth's own runout check prints two console lines per poll
        # while the filament is stalled, and could pause the printer
        # mid-test. The flow test has its own feed check, so push the
        # runout distance out of reach for the duration.
        obj = bdw['obj']
        if hasattr(obj, 'runout_delay_length'):
            bdw['orig_runout_delay'] = obj.runout_delay_length
            obj.runout_delay_length = 1.0e9
            try:
                obj._update_filament_runout_pos()
            except Exception:
                pass
        if bdw['orig_mode'].strip() == 'motion':
            return
        self.gcode.run_script_from_command(
            "SET_BDWIDTH NAME=%s COMMAND=ENABLE_MOTION\nM221 S100"
            % bdw['short'])
        bdw['mode'] = str(getattr(bdw['obj'], 'is_active', 'motion'))
        gcmd.respond_info(
            "BDWidth '%s' switched to motion-only for the test (was %s); "
            "extrusion multiplier reset to 100 %%. Restored at the end."
            % (bdw['name'], bdw['orig_mode']))

    def _bdwidth_restore(self, gcmd, bdw):
        """Put the sensor back in the mode it had before the test."""
        obj = bdw['obj']
        if 'orig_runout_delay' in bdw:
            obj.runout_delay_length = bdw.pop('orig_runout_delay')
            try:
                obj._update_filament_runout_pos()
            except Exception:
                pass
        orig = bdw['orig_mode'].strip()
        if orig == bdw['mode'].strip():
            return
        if 'disable' in orig:
            cmd = 'DISABLE'
        elif orig == 'width':
            cmd = 'ENABLE_WIDTH'
        elif orig == 'motion':
            cmd = 'ENABLE_MOTION'
        else:
            cmd = 'ENABLE'           # 'all'
        try:
            self.gcode.run_script_from_command(
                "SET_BDWIDTH NAME=%s COMMAND=%s" % (bdw['short'], cmd))
            bdw['mode'] = orig
            gcmd.respond_info("BDWidth '%s' restored to mode=%s."
                              % (bdw['name'], orig))
        except Exception as e:
            gcmd.respond_info("Could not restore BDWidth mode (%s) — "
                              "run SET_BDWIDTH NAME=%s COMMAND=%s"
                              % (e, bdw['short'], cmd))

    # ─── MPC filament constants (Kalico) ─────────────────────────
    # Kalico's MPC feed-forward needs the filament's density and heat
    # capacity: MPC_SET HEATER=extruder FILAMENT_DENSITY=… 
    # FILAMENT_HEAT_CAPACITY=…  If the config has the _SET_MPC_MATERIAL
    # macro from Kalico's MPC docs, that macro (and its table) is used;
    # otherwise the same table from those docs is built in.
    MPC_MATERIALS = {
        'PLA': (1.25, 2.20), 'PETG': (1.27, 2.20), 'PC+ABS': (1.15, 2.20),
        'ABS': (1.06, 2.40), 'ASA': (1.07, 2.10), 'PA6': (1.12, 2.50),
        'PA': (1.15, 2.50), 'PC': (1.20, 1.90), 'TPU': (1.21, 2.00),
        'TPU-90A': (1.15, 2.00), 'TPU-95A': (1.22, 2.00),
        'ABS-CF': (1.11, 2.40), 'ASA-CF': (1.11, 2.10),
        'PA6-CF': (1.19, 2.50), 'PC+ABS-CF': (1.22, 2.20),
        'PC+CF': (1.36, 1.90), 'PLA-CF': (1.29, 2.20),
        'PETG-CF': (1.30, 2.20),
    }

    def _mpc_control(self):
        try:
            heater = self.printer.lookup_object('extruder').get_heater()
        except Exception:
            return None
        ctrl = getattr(heater, 'control', None)
        if ctrl is None or not all(hasattr(ctrl, a) for a in (
                'const_filament_density', 'const_filament_heat_capacity')):
            return None
        return ctrl

    def _mpc_apply_material(self, gcmd, material):
        ctrl = self._mpc_control()
        if ctrl is None:
            gcmd.respond_info(
                "MATERIAL=%s recorded. The hotend is not under MPC "
                "control, so no heater model was changed." % material)
            return None
        orig = (float(ctrl.const_filament_density),
                float(ctrl.const_filament_heat_capacity))
        state = {'orig': orig, 'applied': False}
        macro = self.printer.lookup_object(
            'gcode_macro _SET_MPC_MATERIAL', None)
        try:
            if macro is not None:
                self.gcode.run_script_from_command(
                    "_SET_MPC_MATERIAL MATERIAL=%s HEATER=extruder"
                    % material)
                source = 'your _SET_MPC_MATERIAL macro'
            elif material in self.MPC_MATERIALS:
                d, c = self.MPC_MATERIALS[material]
                self.gcode.run_script_from_command(
                    "MPC_SET HEATER=extruder FILAMENT_DENSITY=%.2f "
                    "FILAMENT_HEAT_CAPACITY=%.2f" % (d, c))
                source = 'built-in Kalico table'
            else:
                gcmd.respond_info(
                    "MATERIAL=%s is not in the MPC material table — MPC "
                    "keeps its current filament constants (%.2f / %.2f)."
                    % (material, orig[0], orig[1]))
                return state
        except Exception as e:
            gcmd.respond_info("Could not set MPC for %s (%s) — MPC "
                              "unchanged." % (material, e))
            return state
        state.update({
            'applied': True, 'source': source,
            'density': float(ctrl.const_filament_density),
            'heat_capacity': float(ctrl.const_filament_heat_capacity)})
        # The macro reports the change itself; only the built-in path
        # needs a line of its own.
        if macro is None:
            gcmd.respond_info(
                "MPC set for %s: density %.2f, heat capacity %.2f "
                "(restored to %.2f / %.2f after the test)."
                % (material, state['density'], state['heat_capacity'],
                   orig[0], orig[1]))
        return state

    # ─── Chamber / filter fans (config: chamber_fans) ─────────────
    def _chamber_fans_start(self, gcmd):
        """Switch on the configured chamber fans; return what is needed
        to put them back. Unknown names and other fan types are
        reported and skipped."""
        state = []
        if not self.chamber_fans:
            return state
        now = self.reactor.monotonic()
        on = []
        for name in self.chamber_fans:
            parts = name.split(None, 1)
            short = parts[-1]
            kinds = (('fan_generic', 'temperature_fan') if len(parts) == 1
                     else (parts[0],))
            obj = kind = None
            for k in kinds:
                o = self.printer.lookup_object('%s %s' % (k, short), None)
                if o is not None:
                    obj, kind = o, k
                    break
            if obj is None:
                gcmd.respond_info(
                    "chamber_fans: '%s' is not a [fan_generic] or "
                    "[temperature_fan] — skipped." % name)
                continue
            try:
                st = obj.get_status(now)
                if kind == 'fan_generic':
                    prev = float(st.get('speed', 0.0))
                    self.gcode.run_script_from_command(
                        "SET_FAN_SPEED FAN=%s SPEED=1" % short)
                else:
                    prev = float(st.get('target', 0.0))
                    low = float(getattr(obj, 'min_temp', 0.0))
                    self.gcode.run_script_from_command(
                        "SET_TEMPERATURE_FAN_TARGET TEMPERATURE_FAN=%s "
                        "TARGET=%.1f" % (short, low))
                state.append({'name': short, 'kind': kind, 'prev': prev})
                on.append(short)
            except Exception as e:
                gcmd.respond_info("chamber_fans: could not switch on "
                                  "'%s' (%s) — skipped." % (name, e))
        if on:
            gcmd.respond_info("Chamber fans on for the test: %s "
                              "(restored at the end)." % ", ".join(on))
        return state

    def _chamber_fans_restore(self, gcmd, state):
        while state:
            f = state.pop()
            try:
                if f['kind'] == 'fan_generic':
                    self.gcode.run_script_from_command(
                        "SET_FAN_SPEED FAN=%s SPEED=%.3f"
                        % (f['name'], f['prev']))
                else:
                    self.gcode.run_script_from_command(
                        "SET_TEMPERATURE_FAN_TARGET TEMPERATURE_FAN=%s "
                        "TARGET=%.1f" % (f['name'], f['prev']))
            except Exception as e:
                gcmd.respond_info("Could not restore fan '%s' (%s)."
                                  % (f['name'], e))

    def _mpc_restore(self, gcmd, state):
        if not state.get('applied'):
            return
        d, c = state['orig']
        try:
            self.gcode.run_script_from_command(
                "MPC_SET HEATER=extruder FILAMENT_DENSITY=%.4f "
                "FILAMENT_HEAT_CAPACITY=%.4f" % (d, c))
            state['applied'] = False
        except Exception as e:
            gcmd.respond_info("Could not restore MPC filament constants "
                              "(%s) — run MPC_SET HEATER=extruder "
                              "FILAMENT_DENSITY=%.2f "
                              "FILAMENT_HEAT_CAPACITY=%.2f"
                              % (e, d, c))

    # ─── Generic filament motion sensor (encoder) ─────────────────
    # Klipper's [filament_motion_sensor] (BTT SFS and similar) only
    # pulses every few mm and keeps no counter, but it records where the
    # extruder was at the last pulse: filament_runout_pos =
    # extruder position at the last pulse + detection_length. After a
    # step, "extruder position now − that" is how far the extruder
    # pushed with no filament movement seen.

    def _resolve_encoder(self, gcmd, param, gap_param):
        if param.lower() in ('none', 'off', '0'):
            return None
        for full_name, o in self.printer.lookup_objects():
            if not full_name.startswith('filament_motion_sensor'):
                continue
            short = full_name.split(None, 1)[-1]
            if param.lower() != 'auto' and param not in (short, full_name):
                continue
            # This test always drives [extruder]; in auto mode ignore
            # sensors bound to other extruders.
            if (param.lower() == 'auto'
                    and getattr(o, 'extruder_name', 'extruder')
                    not in (None, 'extruder')):
                continue
            if not all(hasattr(o, a) for a in (
                    'filament_runout_pos', 'detection_length',
                    '_get_extruder_pos')):
                gcmd.respond_info("Encoder '%s' does not expose its "
                                  "position (Klipper version?) — skipped."
                                  % full_name)
                return None
            det = float(o.detection_length)
            limit = gap_param if gap_param > 0 else max(2.0 * det, 10.0)
            helper = getattr(o, 'runout_helper', None)
            gcmd.respond_info("Filament motion data: using encoder '%s' "
                              "(stop when the extruder pushes %.1f mm "
                              "without a pulse)." % (full_name, limit))
            return {'obj': o, 'name': full_name, 'short': short,
                    'detection_length': det, 'gap_limit': limit,
                    'orig_enabled': getattr(helper, 'sensor_enabled',
                                            None)}
        if param.lower() != 'auto':
            gcmd.respond_info("Encoder '%s' not found — skipped." % param)
        return None

    def _encoder_prepare(self, gcmd, enc):
        """Disable the sensor's runout action (pause / runout_gcode)
        for the test; position tracking keeps running."""
        if enc['orig_enabled']:
            try:
                self.gcode.run_script_from_command(
                    "SET_FILAMENT_SENSOR SENSOR=%s ENABLE=0" % enc['short'])
            except Exception:
                pass

    def _encoder_restore(self, gcmd, enc):
        if enc['orig_enabled']:
            try:
                self.gcode.run_script_from_command(
                    "SET_FILAMENT_SENSOR SENSOR=%s ENABLE=1" % enc['short'])
            except Exception as e:
                gcmd.respond_info("Could not re-enable %s (%s) — run "
                                  "SET_FILAMENT_SENSOR SENSOR=%s ENABLE=1"
                                  % (enc['name'], e, enc['short']))

    def _encoder_step(self, enc):
        o = enc['obj']
        try:
            self.gcode.run_script_from_command("M400")
            now_pos = float(o._get_extruder_pos())
            last_pulse = float(o.filament_runout_pos) - float(
                o.detection_length)
            gap = max(now_pos - last_pulse, 0.0)
        except Exception:
            gap = None
        return {'gap_mm': gap, 'limit_mm': enc['gap_limit']}

    def _check_encoder_feed(self, results):
        if not results:
            return None
        e = results[-1].get('enc') or {}
        gap = e.get('gap_mm')
        limit = e.get('limit_mm')
        if gap is None or not limit or gap <= limit:
            return None
        return ("filament motion sensor saw no movement for the last "
                "%.1f mm of extrusion (limit %.1f mm) — the filament "
                "stopped moving (stall, slip or clog)" % (gap, limit))

    def _bdwidth_total(self, bdw):
        """Wait until the sensor has polled past all motion so far,
        then return its accumulated count."""
        wait_ms = int((bdw['sample_time'] * 1.5 + 0.3) * 1000)
        self.gcode.run_script_from_command("M400\nG4 P%d" % wait_ms)
        try:
            return float(bdw['obj'].actual_total_move)
        except Exception:
            return None

    @staticmethod
    def _bdwidth_step(bdw, before, after, cmd_mm):
        out = {'cmd_mm': cmd_mm, 'meas_mm': None, 'ratio': None}
        if before is None or after is None or bdw['counts_per_mm'] <= 0:
            return out
        meas = (after - before) / bdw['counts_per_mm']
        out['meas_mm'] = meas
        out['ratio'] = meas / cmd_mm if cmd_mm > 0 else None
        return out

    def _measure_step(self, gcmd, target_flow, step_duration, repeat,
                      skip_warmup=True, rep_check=None, quiet=False):
        """Run a single flow measurement (multiple repetitions, aggregate).

        rep_check(reps_done, rep_len_mm, mm_per_s) is called after every
        rep but the last; a non-empty return ends the step early (the
        filament clearly stopped — no point grinding the other reps).
        quiet=True skips the summary line (the caller prints its own)."""
        early_stop = None
        mm_per_sec = target_flow / self.filament_area
        feed_rate = mm_per_sec * 60.0
        extrude_length = mm_per_sec * step_duration

        per_run_sg = []
        run_sg_avgs = []
        per_run_trends = []  # intra-run drift metrics
        # Thermal samples now taken DURING extrusion (via the sample
        # timer) and recovery samples taken BETWEEN reps. Old design
        # took snapshots before+after each rep, but those mixed
        # extrusion-stress with recovery values from the cooldown
        # gap, giving misleading averages.
        active_thermal_samples = []   # samples during active extrusion
        recovery_thermal_samples = []  # samples between reps (idle)

        for rep in range(repeat):
            # Capture ONE recovery snapshot just before extrusion starts
            # (gives us the heater's idle/recovery state for diagnostics)
            recovery_thermal_samples.append(self._get_thermal_snapshot())

            self._start_sampling()
            try:
                self.gcode.run_script_from_command(
                    "G1 E%.4f F%.1f\nM400" % (extrude_length, feed_rate))
            finally:
                self._stop_sampling()

            # Pull thermal samples that were captured DURING the
            # extrusion (4 Hz cadence) — these are the only ones that
            # represent actual heater stress under load.
            active_thermal_samples.extend(self.samples_thermal)

            run_sg_raw = list(self.samples_sg)
            # Trim the spin-up AND spin-down transients: while the
            # motor accelerates from / decelerates to standstill, SG2
            # readings are invalid in both directions — they poison
            # sg_min/sg_max and the dip/peak counters. (The intra-run
            # trend code below does its own skip on the raw list.)
            n_raw = len(run_sg_raw)
            if n_raw >= 20:
                skip = max(2, int(n_raw * STARTUP_TRIM_FRACTION))
                tail = max(1, int(n_raw * TAIL_TRIM_FRACTION))
                run_sg = run_sg_raw[skip:n_raw - tail]
            else:
                run_sg = run_sg_raw
            per_run_sg.append(run_sg)
            if run_sg:
                run_sg_avgs.append(sum(run_sg) / len(run_sg))

            # Compute intra-run trend metrics (Phase 3: detect graduel
            # load increase within a single run — the signature of cold
            # extrusion as opposed to abrupt motor slip).
            per_run_trends.append(
                self._compute_intra_run_trend(
                    run_sg_raw, sg2_driver=self.sg2_driver))

            if rep < repeat - 1:
                if rep_check is not None:
                    why = rep_check(rep + 1, extrude_length, mm_per_sec)
                    if why:
                        early_stop = {'reps': rep + 1, 'of': repeat,
                                      'reason': why}
                        break
                self.gcode.run_script_from_command("G4 P300")

        # Warmup-skip: if first run deviates more than the profile's
        # WARMUP_DRIFT_THRESHOLD from the rest, exclude it. The
        # threshold is driver-specific because run-to-run drift
        # varies significantly between TMC chip families:
        #   TMC5160: ~10 % default (validated against historical CSVs)
        #   TMC2240: 4 %  (the chip shows systematic 3-6 % drift on run 1)
        warmup_dropped = False
        included_indices = list(range(len(run_sg_avgs)))
        if skip_warmup and len(run_sg_avgs) >= 3:
            run1 = run_sg_avgs[0]
            rest_mean = sum(run_sg_avgs[1:]) / len(run_sg_avgs[1:])
            if rest_mean > 0:
                deviation = abs(run1 - rest_mean) / rest_mean
                if deviation > self.profile.WARMUP_DRIFT_THRESHOLD:
                    warmup_dropped = True
                    included_indices = list(range(1, len(run_sg_avgs)))

        # Aggregate from included runs
        agg_sg = []
        for idx in included_indices:
            agg_sg.extend(per_run_sg[idx])

        sg_stats = self._stats(agg_sg)

        # Collapse-dip and unload-peak counting (SG2 only): samples
        # far BELOW the run median are brief motor-stall events;
        # samples far ABOVE it are unload events — the extruder gear
        # grinding through the filament removes the load and SG snaps
        # toward its no-load value. Thanks to the spin-up trim above,
        # a clean run has essentially none of either.
        dip_agg = None
        peak_agg = None
        if self.sg2_driver:
            dip_count = 0
            peak_count = 0
            agg_n = 0
            for idx in included_indices:
                samples = per_run_sg[idx]
                if len(samples) < 10:
                    continue
                run_median = sorted(samples)[len(samples) // 2]
                dip_thresh = run_median * self.profile.DIP_FRACTION
                peak_thresh = run_median * (
                    1.0 + self.profile.PEAK_FRACTION)
                dip_count += sum(1 for s in samples if s < dip_thresh)
                peak_count += sum(1 for s in samples
                                  if s > peak_thresh)
                agg_n += len(samples)
            if agg_n > 0:
                dip_agg = {'count': dip_count, 'n': agg_n,
                           'rate': dip_count / agg_n}
                peak_agg = {'count': peak_count, 'n': agg_n,
                            'rate': peak_count / agg_n}

        included_sg_avgs = [run_sg_avgs[i] for i in included_indices
                            if i < len(run_sg_avgs)]
        run_consistency = None
        if len(included_sg_avgs) > 1:
            sg_run_std = _pstdev(included_sg_avgs)
            sg_run_mean = sum(included_sg_avgs) / len(included_sg_avgs)
            sg_cv = (sg_run_std / sg_run_mean * 100.0
                     if sg_run_mean > 0 else 0)
            run_consistency = {
                'sg_run_std': sg_run_std,
                'sg_cv': sg_cv,
                'warmup_dropped': warmup_dropped,
            }

        # Aggregate thermal samples — only the ACTIVE extrusion ones
        # define pwm_avg/max/min and the temp_drop. Recovery samples
        # are kept for diagnostics (did the heater recover between reps?)
        # Done before the summary line so it can report heater
        # saturation.
        thermal_agg = self._aggregate_thermal_samples(active_thermal_samples)

        # Compact summary line
        sg_med_str = "%.0f" % sg_stats['median'] if sg_stats else 'n/a'
        cv_str = ("%.1f%%" % run_consistency['sg_cv']
                  if run_consistency else 'n/a')
        warmup_str = ' [run 1 excluded as warmup]' if warmup_dropped else ''
        dip_str = (' | dips = %d' % dip_agg['count']
                   if dip_agg and dip_agg['count'] > 0 else '')
        if peak_agg and peak_agg['count'] > 0:
            dip_str += ' | peaks = %d' % peak_agg['count']
        # Heater pinned = the hotend, not the motor, is now the limit.
        pwm_avg = (thermal_agg or {}).get('pwm_avg')
        if pwm_avg is not None and pwm_avg >= self.profile.PWM_SATURATED:
            dip_str += ' | heater at 100%% (thermal limit reached)'
        if not quiet:
            gcmd.respond_info(
                "  %.1f mm³/s | SG median = %s | "
                "run-to-run CV = %s%s%s"
                % (target_flow, sg_med_str, cv_str, dip_str, warmup_str))

        # Augment with recovery diagnostics (mainly: did the heater
        # come back to target between reps?)
        if recovery_thermal_samples and thermal_agg:
            recovery_temps = [s['temp_actual']
                              for s in recovery_thermal_samples
                              if s.get('temp_actual') is not None]
            recovery_targets = [s['temp_target']
                                for s in recovery_thermal_samples
                                if s.get('temp_target') is not None]
            if recovery_temps and recovery_targets:
                # How far below target was the heater at the start of
                # each rep (averaged)? Indicates whether cooldown was
                # enough, or heater is struggling between extrusions.
                deltas = [(t - a) for a, t in
                          zip(recovery_temps, recovery_targets)]
                thermal_agg['recovery_deficit_avg'] = (
                    sum(deltas) / len(deltas))
                thermal_agg['recovery_deficit_max'] = max(deltas)

        return {
            'flow': target_flow,
            'sg': sg_stats,
            'run_consistency': run_consistency,
            'run_sg_avgs': run_sg_avgs,
            'warmup_dropped': warmup_dropped,
            'dip': dip_agg,
            'peak': peak_agg,
            'thermal': thermal_agg,
            'intra_run': self._aggregate_intra_run_trends(per_run_trends),
            'early_stop': early_stop,
            'summary': {'sg': sg_med_str, 'cv': cv_str,
                        'dips': (dip_agg or {}).get('count'),
                        'peaks': (peak_agg or {}).get('count'),
                        'warmup': warmup_dropped},
        }

    # ─── Trigger detection — SG-only mode ──────────────────────────

    def _sg_min_informative(self):
        """Lowest SG value still considered informative.

        For SG4 drivers (TMC2240/2209) low SG values mean 'almost no
        load' and are noisy; we ignore them. For SG2 drivers (TMC5160 et
        al.) low SG means HIGH load — exactly the regime we care about
        — so don't gate on it.
        """
        if self.sg2_driver:
            return -1  # never reject on this gate
        return SG_MIN_INFORMATIVE

    def _sg_jump_threshold(self):
        """Minimum raw |SG delta| considered a real (non-noise) jump.

        Uses the driver-specific profile value, which captures the
        relationship between SG scale and jump magnitudes for that
        chip family.
        """
        return self.profile.SG_JUMP_THRESHOLD

    def _check_triggers_sg(self, results):
        """SG-based slip detection.

        SG can either RISE or FALL with increasing motor load depending on
        driver chip, chopper mode, motor, and CoolStep state. Rather than
        hard-coding a direction per driver, we let the data decide:
        average sign of the recent step-to-step deltas tells us which way
        SG is trending. A slip then shows up as a sharp move in the
        OPPOSITE direction (motor decouples from load → SG snaps back to
        its 'no-load' value).

        Two triggers, both fire only when stepping UP:
          1. SG reload jump: SG moves against the established trend by
             more than the typical step magnitude (and >15 raw units)
          2. Plateau over 2 steps: the trend stalls — cumulative
             trend-direction movement is less than half the typical
        """
        # Direction-independent triggers (dip / CV / IQR) get a usable
        # baseline after 2 prior steps; the trend triggers below guard
        # themselves with len(sg_deltas) >= 3 (needs 5 results). The
        # old flat len<5 gate left the first 4 coarse steps entirely
        # trigger-blind — a stall on step 4 went unnoticed.
        # Heater-collapse needs only ONE lower-flow step to compare
        # against, and it does not read StallGuard, so it runs before
        # the history gate below.
        # Both of these mean the filament physically stopped. They are
        # flagged HARD so the run ends instead of bisecting on chewed
        # filament (see _run_find_max: check()).
        self._hard_failure = None
        heater_reason = self._check_heater_collapse(results)
        if heater_reason:
            self._hard_failure = heater_reason
            return heater_reason

        # Feed check: needs no history at all, so it also covers a
        # stall on the very first step.
        feed_reason = self._check_bdwidth_feed(results)
        if feed_reason:
            self._hard_failure = feed_reason
            return feed_reason
        enc_reason = self._check_encoder_feed(results)
        if enc_reason:
            self._hard_failure = enc_reason
            return enc_reason

        if not results or len(results) < 3:
            return None

        # Thermal trigger first: it is the only one that still works
        # when the extruder grinds continuously (see the docstring
        # there), and it must not be skipped by the SG-specific guards
        # below.
        thermal_reason = self._check_thermal_runaway(results)
        if thermal_reason:
            return thermal_reason

        r = results[-1]
        sg_stats = r['sg']
        if sg_stats is None or sg_stats['n'] == 0:
            return None
        sg_med = sg_stats['median']
        target_flow = r['flow']
        sg_label = self._get_sg_label()

        if sg_med <= self._sg_min_informative():
            return None

        prev_flow = results[-2].get('flow', 0)
        going_up = target_flow >= prev_flow - 0.001
        # SG trend triggers (snap-back, abnormal jump) are direction-
        # sensitive — they require a monotonic up-sweep to compare
        # against. Bisection probes can step DOWN, so we skip the
        # trend triggers in that case but still run the
        # direction-independent CV-spike and IQR-spread checks.

        # ─── Trend-direction triggers (snap-back / over-jump / plateau)
        # ─── only valid during monotonic up-sweep (coarse phase). In
        # ─── bisection, the trend baseline is meaningless because we
        # ─── jump around the slip point.
        if going_up:
            # Trend/plateau logic works on the sweep steps only —
            # reference steps are excluded because the flow jump from
            # the reference block to the sweep start would poison the
            # delta baseline (the trigger checks never run DURING the
            # reference phase, so tr[-1] is always the current step).
            tr = [r for r in results
                  if r.get('phase', 'coarse') != 'ref']
            # Raw step-to-step SG deltas from recent history.
            sg_deltas = []
            for j in range(max(1, len(tr) - 5), len(tr) - 1):
                rj, rj_prev = tr[j], tr[j-1]
                if rj.get('sg') and rj_prev.get('sg'):
                    sg_deltas.append(
                        rj['sg']['median'] - rj_prev['sg']['median'])

            if len(sg_deltas) >= 3:
                # Use median (not mean) of prior deltas as the
                # baseline. Median is robust against outliers like
                # the saturation-region first delta — e.g. starting
                # at flow=10 with a high SGT can produce a -418
                # delta that poisons mean-based predictions for
                # several steps and makes the natural saturation
                # curve look like a "plateau" trigger.
                sorted_deltas = sorted(sg_deltas)
                expected_delta = sorted_deltas[len(sorted_deltas) // 2]
                actual_delta = sg_med - tr[-2]['sg']['median']

                # Trend direction (+1 = SG rises with load, -1 = SG falls
                # with load, 0 = flat / inconclusive). Use 1 raw unit as
                # the noise-floor so a tiny non-zero average doesn't lock
                # us into a spurious sign.
                if expected_delta > 1.0:
                    trend_sign = +1
                elif expected_delta < -1.0:
                    trend_sign = -1
                else:
                    trend_sign = 0

                # In trend-direction terms, slip can manifest two ways:
                #   • "Over-jump"   — load_signal jumps in the same
                #                      direction as the trend, but >2×
                #                      the typical step. This is what
                #                      TMC2240/2209 (SG4) often show
                #                      empirically: SG keeps rising but
                #                      suddenly MUCH faster than before.
                #   • "Snap-back"   — load_signal moves OPPOSITE to the
                #                      trend (motor decouples, SG returns
                #                      toward its no-load value). This is
                #                      the textbook SG2 stall signature
                #                      seen on TMC5160 etc.
                # Both should fire; we test each independently.
                if trend_sign != 0:
                    expected_load = expected_delta * trend_sign     # > 0
                    actual_load = actual_delta * trend_sign         # any sign

                    # Snap-back (slip via decoupling)
                    if (actual_load < -expected_load
                            and abs(actual_delta)
                            > self._sg_jump_threshold()):
                        if trend_sign > 0:
                            return ("%s reload jump: %+.0f (expected to "
                                    "keep rising ~+%.0f) — slip detected"
                                    % (sg_label, actual_delta,
                                       expected_load))
                        return ("%s reload jump: %+.0f (expected to keep "
                                "falling ~-%.0f) — slip detected"
                                % (sg_label, actual_delta, expected_load))

                    # Over-jump (slip via abnormal acceleration)
                    if (actual_load > expected_load * 2.0
                            and abs(actual_delta)
                            > self._sg_jump_threshold()):
                        if trend_sign > 0:
                            return ("%s abnormal jump: %+.0f vs expected "
                                    "~+%.0f (%.1fx larger) — slip detected"
                                    % (sg_label, actual_delta,
                                       expected_load,
                                       actual_load / expected_load))
                        return ("%s abnormal drop: %+.0f vs expected "
                                "~-%.0f (%.1fx larger) — slip detected"
                                % (sg_label, actual_delta, expected_load,
                                   actual_load / expected_load))

                    # Trigger 2: plateau over 2 steps, only meaningful
                    # when the trend is sizeable. The threshold is
                    # driver-specific (PLATEAU_RATIO in the profile)
                    # because the steepness of the SG-vs-load
                    # saturation curve varies between chip families.
                    #
                    # Restricted to the coarse phase only: in
                    # bisection / verify, flow probes can step in
                    # either direction (e.g. 45 → 42 → 44 → 43), and
                    # the trend baseline assumption breaks down. The
                    # `going_up` check above doesn't catch this
                    # because individual bisect probes are still
                    # locally upward — but the prior-step deltas in
                    # bisection mode reflect bisection bracket logic,
                    # not the smooth coarse sweep this trigger
                    # assumes.
                    #
                    # Also skip plateau if any prior coarse step had
                    # SG values in the saturation region (close to
                    # 1023). Saturation clipping artificially inflates
                    # the prior delta — the SG-vs-load curve looks
                    # MUCH steeper than reality, so the trigger then
                    # interprets normal saturation behaviour later in
                    # the sweep as a "stalled trend".
                    last_phase = tr[-1].get('phase', 'coarse')
                    prior_sg_max = max(
                        tr[j]['sg']['median']
                        for j in range(max(0, len(tr) - 5),
                                       len(tr) - 1)
                        if tr[j].get('sg'))
                    saturation_skip = (
                        prior_sg_max
                        > self.profile.PLATEAU_SATURATION_SKIP)
                    if (expected_load > 5
                            and last_phase == 'coarse'
                            and not saturation_skip):
                        # Single-step plateau check (NEW): if the most
                        # recent step's SG-delta is essentially flat
                        # compared to the recent baseline, fire even
                        # without needing a 2-step cumulative window.
                        # An abrupt plateau (e.g. typical -25 deltas
                        # then suddenly +0.5) is a strong slip signal
                        # and shouldn't be smoothed over 2 steps.
                        # Threshold is half the PLATEAU_RATIO — must
                        # be a much clearer single-step plateau to
                        # fire than a borderline 2-step accumulation.
                        single_threshold_factor = (
                            self.profile.PLATEAU_RATIO / 2.0)
                        single_threshold = (expected_load
                                            * single_threshold_factor)
                        if actual_load < single_threshold:
                            direction = ("rising" if trend_sign > 0
                                         else "falling")
                            return ("%s single-step plateau: %s trend "
                                    "flattened abruptly (only %.0f vs "
                                    "expected %.0f load-units in this "
                                    "step) — strong sign that flow "
                                    "increase no longer adds motor load"
                                    % (sg_label, direction,
                                       actual_load, expected_load))

                        # Cumulative 2-step check (existing logic)
                        prev_actual = (tr[-2]['sg']['median']
                                       - tr[-3]['sg']['median'])
                        cumulative_load = (actual_load
                                           + prev_actual * trend_sign)
                        expected_2step = expected_load * 2
                        if cumulative_load < (expected_2step
                                              * self.profile.PLATEAU_RATIO):
                            direction = ("rising" if trend_sign > 0
                                         else "falling")
                            return ("%s plateau over 2 steps: %s trend "
                                    "stalled (only %.0f vs typical %.0f "
                                    "load-units) — flow no longer "
                                    "increasing motor load (slip starting)"
                                    % (sg_label, direction,
                                       cumulative_load, expected_2step))

        # ─── Direction-independent triggers (always evaluated) ───
        # CV-spike and IQR-spread fire on the per-step variance signal
        # itself, not on a trend baseline, so they make sense for both
        # coarse and bisection probes.

        # Trigger 0: collapse dips (SG2 only). Brief stick-slip stalls
        # drive SG toward 0 for a few samples — the distribution stats
        # below absorb them completely, so check the dip counter first.
        dip_reason = self._check_dip_collapse(results, sg_label)
        if dip_reason:
            return dip_reason

        # Trigger 0b: unload peaks (SG2 only) — the extruder-grinding
        # mirror image of the collapse dip.
        peak_reason = self._check_peak_unload(results, sg_label)
        if peak_reason:
            return peak_reason

        # Trigger 3: run-to-run CV spike. Intermittent slip shows up as
        # a sudden burst of variance between repeats at the same flow,
        # even when the median doesn't move much. Fire when the latest
        # CV is both elevated in absolute terms AND much higher than the
        # baseline of recent steps.
        cv_reason = self._check_cv_spike(results, sg_label)
        if cv_reason:
            return cv_reason

        # The remaining distribution triggers keep the original 5-step
        # maturity requirement: early coarse steps often sit in the SG
        # saturation region with spiky max/IQR, and replaying the
        # validated CSVs showed false fires (e.g. 5160 gap-jump at
        # flow=40 on a run whose real limit was ~114) when these run
        # on a 2-step baseline.
        if len(results) < 5:
            return None

        # Trigger 4: IQR/spread anomaly. Catches "quiet" stalls where
        # the median looks normal but the distribution widens (brief
        # intermittent stall absorbed by median over 5 repeats).
        iqr_reason = self._check_iqr_spread(results, sg_label)
        if iqr_reason:
            return iqr_reason

        # Trigger 5: SG max spike — bisection / verify only. Catches
        # decoupling events where one of the repeats briefly snaps SG
        # to its no-load value, but the median absorbs it.
        spike_reason = self._check_sg_max_spike(results, sg_label)
        if spike_reason:
            return spike_reason

        # Trigger 6: per-run outlier — bisection / verify only. Catches
        # the case where 1 of N repeats has a clearly different SG/CS
        # profile from the others (intermittent slip in just one run).
        outlier_reason = self._check_run_outlier(results, sg_label)
        if outlier_reason:
            return outlier_reason

        return None

    MELT_KNEE_FRAC = 0.85   # heat per mm³ below this × baseline = knee

    def _step_extras(self, r, meta):
        """Per-step data-collection values for reports: idle duty at the
        step's time, melt duty (above idle), melt duty per mm³, % of the
        fitted low-flow line (once the fit is in meta) and BDWidth."""
        pwm = (r.get('thermal') or {}).get('pwm_avg')
        idle = self._idle_pwm_at(meta, r.get('t_s'))
        melt = (pwm - idle) if (pwm is not None and idle is not None) else None
        flow = r.get('flow') or 0.0
        fit_pct = None
        fa = meta.get('melt_fit_intercept')
        fb = meta.get('melt_fit_slope')
        if fa is not None and fb and pwm is not None and flow:
            p0 = meta.get('idle_pwm_start')
            drift = (idle - p0) if (p0 is not None
                                     and idle is not None) else 0.0
            fit_pct = (pwm - drift - fa) / (fb * flow) * 100
        b = r.get('bdw') or {}
        ratio = b.get('ratio')
        return {
            'idle': idle,
            'melt': melt,
            'melt_per_flow': (melt / flow) if (melt is not None
                                               and flow) else None,
            'fit_pct': fit_pct,
            'bdw_cmd': b.get('cmd_mm'),
            'bdw_meas': b.get('meas_mm'),
            'bdw_pct': (ratio * 100) if ratio is not None else None,
        }

    @staticmethod
    def _idle_pwm_at(meta, t):
        p0 = meta.get('idle_pwm_start')
        p1 = meta.get('idle_pwm_end')
        t1 = meta.get('idle_t_end_s')
        if p0 is None:
            return p1
        if p1 is None or not t1 or t is None:
            return p0
        frac = min(max(t / t1, 0.0), 1.0)
        return p0 + (p1 - p0) * frac

    MELT_FIT_STEPS = 4      # steps 2..5 in run order define the line
    MELT_OFFSET_WARN = 0.05 # |fit intercept − measured idle| above this

    def _melt_fit(self, results, meta):
        """Straight line heater-duty = a + b·flow through the healthy
        low-flow steps (steps 2–5 in run order; step 1 reads low while
        things warm up). Duty is drift-corrected with the start/end idle
        baselines. Fitting the line instead of dividing by the measured
        idle absorbs the constant per-rep offset that otherwise makes
        heat-per-mm³ slide down from the first step on."""
        p0 = meta.get('idle_pwm_start')
        pts = []
        for r in results:
            pwm = (r.get('thermal') or {}).get('pwm_avg')
            if pwm is None or not r.get('flow') or r.get('hard_fail'):
                continue
            drift = 0.0
            if p0 is not None:
                cur_idle = self._idle_pwm_at(meta, r.get('t_s'))
                if cur_idle is not None:
                    drift = cur_idle - p0
            pts.append((r['flow'], pwm - drift))
        fit_pts = pts[1:1 + self.MELT_FIT_STEPS]
        if len(fit_pts) < 3 or len(pts) < len(fit_pts) + 2:
            return None
        n = len(fit_pts)
        sx = sum(f for f, _ in fit_pts)
        sy = sum(v for _, v in fit_pts)
        sxx = sum(f * f for f, _ in fit_pts)
        sxy = sum(f * v for f, v in fit_pts)
        den = n * sxx - sx * sx
        if den <= 0:
            return None
        b = (n * sxy - sx * sy) / den
        a = (sy - b * sx) / n
        return {'a': a, 'b': b, 'pts': pts,
                'fit_max_flow': max(f for f, _ in fit_pts)}

    def _melt_knee(self, results, meta):
        """Quality limit from the heater: each step's duty above the
        fitted intercept vs what the low-flow line predicts. The first
        flow where it drops below MELT_KNEE_FRAC is the melt limit.
        Report only — never stops the test."""
        fit = self._melt_fit(results, meta)
        if fit is None:
            return None
        a, b = fit['a'], fit['b']
        frac = self.MELT_KNEE_FRAC
        out = {'intercept': a, 'slope': b, 'frac': frac, 'warn': None}
        if b <= 0:
            out['warn'] = ('low-flow heater duty does not rise with flow '
                           '(slope %.5f) — melt limit not available' % b)
            out['knee'] = None
            out['max_tested'] = max(f for f, _ in fit['pts'])
            return out
        idle0 = meta.get('idle_pwm_start')
        if idle0 is not None and abs(a - idle0) > self.MELT_OFFSET_WARN:
            out['warn'] = ('fit intercept %.3f is %.3f away from the '
                           'measured idle %.3f — treat the melt limit '
                           'with caution' % (a, a - idle0, idle0))
        later = sorted((f, (v - a) / (b * f)) for f, v in fit['pts']
                       if f > fit['fit_max_flow'] + 0.001)
        knee = next(((f, rt) for f, rt in later if rt < frac), None)
        if knee is None:
            out['knee'] = None
            out['max_tested'] = max(f for f, _ in fit['pts'])
            return out
        ok = [f for f, v in fit['pts']
              if f < knee[0] - 0.001 and (v - a) / (b * f) >= frac]
        out.update({'knee': knee[0], 'knee_ratio': knee[1],
                    'last_ok': max(ok) if ok else None})
        return out

    def _check_bdwidth_feed(self, results):
        """Filament did not move: the BDWidth motion sensor (before the
        extruder) saw far less filament than the step commanded. Catches
        a motor that loses sync and sits still, a stripped/slipping
        drive and a hard clog, whatever StallGuard reads.

        The sensor's scale and per-step noise are loose (±5 %, one
        −25 % dropout seen in 10 steps), so the threshold is coarse and
        the ratio is scale-corrected by the median of the earlier steps
        (clamped to 0.8–1.2; 1.0 until two earlier steps exist).
        """
        min_ratio = getattr(self, '_bdw_min_ratio', 0.0) or 0.0
        if min_ratio <= 0 or not results:
            return None
        cur = results[-1].get('bdw') or {}
        ratio = cur.get('ratio')
        if ratio is None:
            return None
        prior = sorted(r['bdw']['ratio'] for r in results[:-1]
                       if (r.get('bdw') or {}).get('ratio') is not None)
        if len(prior) >= 2:
            n = len(prior)
            scale = (prior[n // 2] if n % 2
                     else (prior[n // 2 - 1] + prior[n // 2]) / 2.0)
            scale = min(max(scale, 0.8), 1.2)
        else:
            scale = 1.0
        rel = ratio / scale
        if rel >= min_ratio:
            return None
        return ("BDWidth saw %.0f %% of the expected filament feed "
                "(%.1f of %.1f mm, sensor scale %.2f) — below the %.0f %% "
                "limit, so the filament stopped moving (stall, slip or "
                "clog)" % (rel * 100, cur.get('meas_mm') or 0.0,
                           cur.get('cmd_mm') or 0.0, scale,
                           min_ratio * 100))

    SLIP_ONSET_FRAC = 0.90   # feed below this × low-flow reading = slip

    def _slip_onset(self, results):
        """Where partial slip starts, from BDWidth (report only).

        Feed ratios are scaled by the median of the first three steps
        (the sensor's own scale for this filament, clamped 0.8–1.2).
        Slip start = the lowest flow from which every higher step reads
        below SLIP_ONSET_FRAC, so a single noisy dropout lower down is
        ignored."""
        steps = [(r['flow'], r['bdw']['ratio']) for r in results
                 if (r.get('bdw') or {}).get('ratio') is not None]
        if len(steps) < 3:
            return None
        first = sorted(v for _, v in steps[:3])
        scale = min(max(first[1], 0.8), 1.2)
        frac = self.SLIP_ONSET_FRAC
        by_flow = sorted((f, v / scale) for f, v in steps)
        onset = None
        for i, (f, rel) in enumerate(by_flow):
            if all(rr < frac for _, rr in by_flow[i:]):
                onset = (f, rel)
                break
        out = {'frac': frac, 'scale': scale,
               'max_tested': by_flow[-1][0]}
        if onset is None:
            out.update({'onset': None})
            return out
        clean = [f for f, rr in by_flow
                 if f < onset[0] - 0.001 and rr >= frac]
        out.update({'onset': onset[0], 'onset_rel': onset[1],
                    'last_clean': max(clean) if clean else None})
        return out

    def _live_melt_expect(self, results):
        """Expected heater duty for the latest step from the low-flow
        line fitted so far (steps 2–5 in run order; no drift correction
        yet — the end idle comes later). None while still fitting."""
        pts = [(r['flow'], r['thermal']['pwm_avg']) for r in results
               if (r.get('thermal') or {}).get('pwm_avg') is not None
               and r.get('flow') and not r.get('hard_fail')]
        if len(pts) <= 1 + self.MELT_FIT_STEPS:
            return None
        fit = pts[1:1 + self.MELT_FIT_STEPS]
        n = len(fit)
        sx = sum(f for f, _ in fit)
        sy = sum(v for _, v in fit)
        sxx = sum(f * f for f, _ in fit)
        sxy = sum(f * v for f, v in fit)
        den = n * sxx - sx * sx
        if den <= 0:
            return None
        b = (n * sxy - sx * sy) / den
        a = (sy - b * sx) / n
        return {'a': a, 'b': b}

    def _check_heater_collapse(self, results):
        """Detect that the hotend stopped melting the commanded flow.

        Compares this step's average heater duty (active extrusion
        samples only) with the step at the highest flow BELOW this one.
        That reference is always a step that passed: once a step fails,
        every later probe sits below it. Valid in coarse, bisection and
        verify alike.

        Fires when duty fell by more than PWM_COLLAPSE_FRAC and by at
        least PWM_COLLAPSE_MIN_ABS even though more flow was commanded.
        """
        if not results or len(results) < 2:
            return None
        p = self.profile
        cur = results[-1]
        flow = cur.get('flow')
        pwm = (cur.get('thermal') or {}).get('pwm_avg')
        if flow is None or pwm is None:
            return None
        ref = None
        for r in results[:-1]:
            rf = r.get('flow')
            rp = (r.get('thermal') or {}).get('pwm_avg')
            if rf is None or rp is None or rf >= flow - 0.001:
                continue
            # highest lower flow wins; on a tie the most recent one
            if ref is None or rf >= ref[0]:
                ref = (rf, rp)
        if ref is None:
            return None
        ref_flow, ref_pwm = ref
        drop = ref_pwm - pwm
        if drop < p.PWM_COLLAPSE_MIN_ABS:
            return None
        if pwm >= ref_pwm * (1.0 - p.PWM_COLLAPSE_FRAC):
            return None
        return ("heater duty fell to %.0f %% at %.1f mm³/s from %.0f %% "
                "at %.1f mm³/s — more flow was commanded but less heat "
                "is being absorbed, so filament stopped arriving "
                "(stall, grind or clog)"
                % (pwm * 100, flow, ref_pwm * 100, ref_flow))

    def _check_thermal_runaway(self, results):
        """Detect that the extruder stopped delivering the commanded flow.

        Every trigger above reads the StallGuard signal. Continuous
        grinding defeats all of them: the hobbed gear chews into
        filament that no longer advances, which keeps the motor load
        high and steady — medians stay smooth, no dips, no unload
        peaks, CV unchanged.

        What does change is heat. Melting X mm³/s absorbs a matching
        amount of power; if the material stops arriving, that heat sink
        disappears and the hotend overshoots its target even though the
        heater was maxed out one step earlier. Validated on a TMC2240
        run that ground at 80 mm³/s: +2 °C at 80, +7.5 at 90, +17.9 at
        100, while eight clean historical runs never exceeded +4.9.

        Both an absolute floor and a jump over the run's own baseline
        must be cleared, so a hotend that simply runs a few degrees warm
        does not trip it.
        """
        p = self.profile
        th = (results[-1].get('thermal') or {})
        target = th.get('temp_target')
        avg = th.get('temp_avg')
        if not target or avg is None or target <= 0:
            return None
        offset = avg - target
        if offset < p.TEMP_RUNAWAY_ABS:
            return None
        prior = []
        for r in results[:-1]:
            t = r.get('thermal') or {}
            if t.get('temp_target') and t.get('temp_avg') is not None:
                prior.append(t['temp_avg'] - t['temp_target'])
        if len(prior) < 2:
            return None
        s = sorted(prior)
        n = len(s)
        base = s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0
        if offset - base < p.TEMP_RUNAWAY_JUMP:
            return None
        return ("hotend ran %.1f °C above target (%.1f vs %.1f °C) while "
                "the earlier steps held %+.1f °C — the melt zone is no "
                "longer absorbing the commanded flow, so the extruder is "
                "not delivering it (grinding)"
                % (offset, avg, target, base))

    def _check_dip_collapse(self, results, sg_label):
        """Detect brief SG collapse dips — the stick-slip stall signature.

        On SG2 drivers a lost-step event snaps SG_RESULT toward 0 for
        ~1-2 samples per click. Median, IQR and CV barely move (a few
        dips among hundreds of samples), so those triggers stay quiet
        even while the motor audibly rattles. Count samples below
        DIP_FRACTION × run-median (spin-up transient already trimmed
        in _measure_step) and fire when the step's dip rate clears the
        absolute floor AND a multiple of the prior-step baseline.
        """
        if not self.sg2_driver:
            return None
        p = self.profile
        dip = results[-1].get('dip')
        if not dip or dip['n'] == 0:
            return None
        rate = dip['rate']
        if dip['count'] < p.DIP_MIN_COUNT or rate < p.DIP_RATE_MIN:
            return None
        prior_rates = sorted(
            r['dip']['rate'] for r in results[:-1]
            if r.get('dip') and r['dip']['n'] > 0)
        baseline = (prior_rates[len(prior_rates) // 2]
                    if prior_rates else 0.0)
        if rate >= max(p.DIP_RATE_MIN, baseline * p.DIP_BASE_RATIO):
            return ("%s collapse dips: %d of %d samples below %.0f%% of "
                    "the run median (%.1f%% of samples vs baseline "
                    "%.1f%%) — brief stall events (stick-slip)"
                    % (sg_label, dip['count'], dip['n'],
                       p.DIP_FRACTION * 100, rate * 100, baseline * 100))
        return None

    def _check_peak_unload(self, results, sg_label):
        """Detect unload peaks — the extruder-grinding signature.

        When the drive gear grinds through the filament the motor is
        briefly UNLOADED and SG spikes toward its no-load value while
        median/IQR/CV stay compact (the mirror image of the collapse
        dip). Counting samples above (1+PEAK_FRACTION) × run-median
        makes a single stray sg_max outlier harmless — the trigger
        needs a sustained peak rate above the prior-step baseline.
        """
        if not self.sg2_driver:
            return None
        p = self.profile
        peak = results[-1].get('peak')
        if not peak or peak['n'] == 0:
            return None
        rate = peak['rate']
        if peak['count'] < p.PEAK_MIN_COUNT or rate < p.PEAK_RATE_MIN:
            return None
        prior_rates = sorted(
            r['peak']['rate'] for r in results[:-1]
            if r.get('peak') and r['peak']['n'] > 0)
        baseline = (prior_rates[len(prior_rates) // 2]
                    if prior_rates else 0.0)
        if rate >= max(p.PEAK_RATE_MIN, baseline * p.PEAK_BASE_RATIO):
            return ("%s unload peaks: %d of %d samples above %.0f%% of "
                    "the run median (%.2f%% of samples vs baseline "
                    "%.2f%%) — motor repeatedly unloaded (gear "
                    "grinding through the filament)"
                    % (sg_label, peak['count'], peak['n'],
                       (1.0 + p.PEAK_FRACTION) * 100,
                       rate * 100, baseline * 100))
        return None

    def _check_cv_spike(self, results, sg_label):
        """Detect a run-to-run CV spike — intermittent slip signature.

        A clean motor produces tight repeats (CV usually <3% on stable
        SG2 setups, may go up to 5% on noisier SG4). When the motor
        starts intermittently slipping, individual repeats diverge
        sharply, even if the median across all samples looks 'normal'.

        Triggers (any one fires):
          (a) "high-variance trip" — last CV >= 10% absolute (regardless
              of baseline). Catches sudden chaotic slip.
          (b) "low-baseline jump"  — last CV >= ratio * avg of prior
              N steps AND last CV >= absolute floor. Ratio and floor
              are tighter in coarse, looser in bisection (where we're
              already near slip).
          (c) "rising CV trend"   — CV grew across 2 consecutive
              steps (each >=1.3x prior).
          (d) "bisection absolute" — in bisection / verify, ANY CV
              >= 2x the *coarse-phase median CV* (or >= 5% absolute,
              whichever is higher). Catches the case where
              intermittent slip already started during coarse and
              "polluted" the immediate baseline, making patterns (b)
              and (c) under-react.

        Baseline strategy: when we're in bisection/verify, the most
        meaningful baseline is the *coarse phase* — it represents the
        clean, slip-free motor behaviour. Using only the immediate 3
        prior steps as baseline lets earlier elevated bisection CVs
        raise the threshold and mask later anomalies.
        """
        last = results[-1]
        last_rc = last.get('run_consistency') or {}
        if 'sg_cv' not in last_rc:
            return None
        last_cv = last_rc['sg_cv']

        last_phase = last.get('phase', 'coarse')
        in_bisection = last_phase in ('bisect', 'verify')

        # Build two baselines:
        #   - "immediate": last 3 steps regardless of phase (sensitive
        #     to recent trend)
        #   - "coarse": median of CVs from the COARSE phase only
        #     (slip-free reference, only available once we have data)
        prior_cvs = []
        for r in results[-4:-1]:
            rc = r.get('run_consistency') or {}
            if 'sg_cv' in rc:
                prior_cvs.append(rc['sg_cv'])

        coarse_cvs = []
        for r in results[:-1]:
            if r.get('phase', 'coarse') not in ('coarse', 'ref'):
                continue
            rc = r.get('run_consistency') or {}
            if 'sg_cv' in rc:
                coarse_cvs.append(rc['sg_cv'])
        # Median of coarse CVs, excluding the LAST coarse step (which
        # may be the elevated step that triggered bisection).
        coarse_baseline_cv = None
        if len(coarse_cvs) >= 4:
            sorted_cv = sorted(coarse_cvs[:-1])
            n = len(sorted_cv)
            coarse_baseline_cv = (sorted_cv[n // 2] if n % 2
                                  else (sorted_cv[n // 2 - 1]
                                        + sorted_cv[n // 2]) / 2)

        # Pattern (a): high-variance trip — works without baseline
        if last_cv >= self.profile.CV_HIGH_VARIANCE:
            if not prior_cvs:
                return ("run-to-run %s CV spiked to %.1f%% — repeats "
                        "diverging, intermittent slip"
                        % (sg_label, last_cv))
            avg_prior_cv = sum(prior_cvs) / len(prior_cvs)
            if last_cv >= self.profile.CV_LOWBASE_RATIO * avg_prior_cv:
                return ("run-to-run %s CV spiked to %.1f%% (baseline "
                        "~%.1f%% over previous 3 steps) — repeats "
                        "diverging, intermittent slip"
                        % (sg_label, last_cv, avg_prior_cv))

        if len(prior_cvs) < 3:
            return None
        avg_prior_cv = sum(prior_cvs) / len(prior_cvs)

        # Pattern (b): low-baseline jump vs immediate prior
        b_ratio = (self.profile.CV_JUMP_RATIO_BISECT if in_bisection
                   else self.profile.CV_JUMP_RATIO_COARSE)
        b_min_cv = (self.profile.CV_JUMP_MIN_BISECT if in_bisection
                    else self.profile.CV_JUMP_MIN_COARSE)
        if last_cv >= b_min_cv and last_cv >= b_ratio * avg_prior_cv:
            return ("run-to-run %s CV jumped to %.1f%% (%.1fx baseline "
                    "of %.1f%% over previous 3 steps) — repeats "
                    "diverging, intermittent slip"
                    % (sg_label, last_cv,
                       last_cv / avg_prior_cv, avg_prior_cv))

        # Pattern (c): rising CV trend across 2 consecutive steps.
        c_min_cv = (self.profile.CV_RISING_MIN_LAST_BISECT if in_bisection
                    else self.profile.CV_RISING_MIN_LAST_COARSE)
        c_ratio = self.profile.CV_RISING_RATIO
        if (last_cv >= c_min_cv
                and len(prior_cvs) >= 2
                and prior_cvs[-1] >= c_ratio * prior_cvs[-2]
                and last_cv >= c_ratio * prior_cvs[-1]
                and prior_cvs[-1] >= self.profile.CV_RISING_MIN_PRIOR):
            return ("run-to-run %s CV rising across 3 steps "
                    "(%.1f%% → %.1f%% → %.1f%%) — gradual slip onset"
                    % (sg_label, prior_cvs[-2], prior_cvs[-1], last_cv))

        # Pattern (d): bisection absolute vs coarse baseline.
        # In bisection/verify, compare against the slip-free coarse
        # phase median. This catches cases where one elevated bisect
        # step has already pulled the immediate baseline up, masking
        # subsequent elevated values.
        if (in_bisection
                and coarse_baseline_cv is not None
                and last_cv >= self.profile.CV_VS_COARSE_MIN
                and last_cv >=
                    self.profile.CV_VS_COARSE_RATIO * coarse_baseline_cv):
            return ("run-to-run %s CV %.1f%% in bisection vs coarse-"
                    "phase baseline %.1f%% (%.1fx) — slip signature "
                    "above clean-extrusion noise floor"
                    % (sg_label, last_cv, coarse_baseline_cv,
                       last_cv / coarse_baseline_cv))

        return None

    def _check_iqr_spread(self, results, sg_label):
        """Detect an IQR/spread anomaly — quiet stall signature.

        Sometimes a motor stalls briefly (start or end of move) and
        recovers, leaving a median that looks normal but a wider
        distribution. The IQR (P75-P25) widens when this happens.

        Four triggers, in order:
          (a) "ratio vs immediate" — current IQR >= ratio * avg of
              prior 3 IQRs. Ratio is 3.0 in coarse, 1.7 in bisection.
          (b) "ratio vs coarse baseline" (bisection only) — current
              IQR >= 2.0 * the median IQR from the COARSE phase.
              Catches cases where earlier elevated bisect IQRs have
              raised the immediate baseline and masked subsequent
              anomalies.
          (c) "absolute" — current IQR >= 25 raw units in bisection /
              verify. A motor near slip should have median, P25 and
              P75 within a few units of each other; an IQR of 25+
              indicates intermittent stall regardless of context.
          (d) "cumulative growth" (coarse only) — IQR has roughly
              doubled relative to the early-test baseline (median of
              first 3 coarse IQRs) AND now exceeds an absolute floor.
              Catches GRADUAL widening that single-step ratios miss.
              Strong pre-slip indicator: stick-slip events build up
              over multiple steps before the abrupt fail.
        """
        last = results[-1]
        sg_stats = last.get('sg') or {}
        if 'p25' not in sg_stats or 'p75' not in sg_stats:
            return None
        last_iqr = sg_stats['p75'] - sg_stats['p25']

        prior_iqrs = []
        for r in results[-4:-1]:
            sg_p = r.get('sg') or {}
            if 'p25' in sg_p and 'p75' in sg_p:
                prior_iqrs.append(sg_p['p75'] - sg_p['p25'])

        last_phase = last.get('phase', 'coarse')
        in_bisection = last_phase in ('bisect', 'verify')

        # (c) absolute — fires unconditionally on any genuinely large
        # IQR during bisection / verify. Cheap and reliable.
        if in_bisection and last_iqr >= self.profile.IQR_ABSOLUTE_TRIGGER:
            return ("%s spread widened: IQR %.0f raw units in %s — "
                    "samples no longer cluster, intermittent stall"
                    % (sg_label, last_iqr, last_phase))

        # (d) cumulative growth — coarse phase only. Catches gradually
        # widening IQR over many steps (typical slow stick-slip onset
        # that single-step ratio triggers miss). Requires at least 5
        # coarse steps of history so the early baseline is meaningful.
        if not in_bisection:
            coarse_iqrs = []
            for r in results[:-1]:
                if r.get('phase', 'coarse') not in ('coarse', 'ref'):
                    continue
                sg_p = r.get('sg') or {}
                if 'p25' in sg_p and 'p75' in sg_p:
                    coarse_iqrs.append(sg_p['p75'] - sg_p['p25'])
            # Need enough history for a stable baseline + clear margin
            # before the current step.
            if len(coarse_iqrs) >= 5:
                # Early baseline: median of first 3 IQRs (after the
                # very first which can be warmup-noisy)
                early_iqrs = sorted(coarse_iqrs[:4])
                # take the 2nd and 3rd smallest as a robust early baseline
                early_baseline = (early_iqrs[1] + early_iqrs[2]) / 2.0
                # Threshold: IQR has ≥ 2× early baseline AND exceeds
                # an absolute floor (so very-small early IQRs of 2-3
                # don't trigger on noise alone). Floor is set
                # deliberately above IQR_RATIO_MIN_ABS so this trigger
                # is independent of the per-step ratio path.
                growth_floor_abs = max(
                    self.profile.IQR_RATIO_MIN_ABS, 12)
                if (early_baseline >= 0.5
                        and last_iqr >= growth_floor_abs
                        and last_iqr >= 2.0 * early_baseline):
                    return ("%s spread grew over sweep: IQR climbed "
                            "from baseline ~%.0f (early steps) to %.0f "
                            "(%.1fx) — gradual pre-slip widening, motor "
                            "is stuttering more often as flow rises"
                            % (sg_label, early_baseline, last_iqr,
                               last_iqr / early_baseline))

        if len(prior_iqrs) < 3:
            return None
        avg_prior_iqr = sum(prior_iqrs) / len(prior_iqrs)
        if avg_prior_iqr < 1:
            return None

        # (a) ratio vs immediate prior steps.
        # If the profile requires CV cross-check (e.g. SG4 drivers with
        # active CoolStep), confirm with the run-to-run CV — without
        # elevated CV, a wide IQR is likely a CoolStep-current-transition
        # effect, not slip. Genuine slip raises both within-step spread
        # AND run-to-run variance.
        ratio_threshold = (self.profile.IQR_RATIO_BISECT if in_bisection
                           else self.profile.IQR_RATIO_COARSE)
        if (last_iqr >= self.profile.IQR_RATIO_MIN_ABS
                and last_iqr >= ratio_threshold * avg_prior_iqr):
            cv_confirms = True
            cv_note = ""
            if self.profile.IQR_RATIO_REQUIRE_CV:
                last_rc = last.get('run_consistency') or {}
                last_cv_a = last_rc.get('sg_cv')
                cv_floor = (self.profile.IQR_RATIO_CV_FLOOR_BISECT
                            if in_bisection
                            else self.profile.IQR_RATIO_CV_FLOOR_COARSE)
                # Build prior-CV baseline same way as prior-IQR baseline
                prior_cvs_a = []
                for r in results[-4:-1]:
                    rc = r.get('run_consistency') or {}
                    if 'sg_cv' in rc:
                        prior_cvs_a.append(rc['sg_cv'])
                avg_prior_cv = (sum(prior_cvs_a) / len(prior_cvs_a)
                                if prior_cvs_a else 0.0)
                # CV must be elevated above absolute floor AND above
                # baseline (or baseline must be tiny). 1.4× ratio is
                # less strict than the value-based jump triggers in
                # _check_cv_spike — IQR pattern (a) only needs CV to
                # *not contradict* the spread anomaly.
                cv_confirms = (last_cv_a is not None
                               and last_cv_a >= cv_floor
                               and (avg_prior_cv < 0.5
                                    or last_cv_a >= 1.4 * avg_prior_cv))
                if cv_confirms and last_cv_a is not None:
                    cv_note = (" (CV %.1f%% confirms, prior ~%.1f%%)"
                               % (last_cv_a, avg_prior_cv))
            if cv_confirms:
                return ("%s spread widened: IQR %.0f vs baseline ~%.0f "
                        "(%.1fx) — quiet/intermittent stall in samples"
                        "%s"
                        % (sg_label, last_iqr, avg_prior_iqr,
                           last_iqr / avg_prior_iqr, cv_note))

        # (b) ratio vs coarse baseline (bisection only).
        # Cross-check with CV: a pure spread anomaly with NORMAL CV is
        # often just statistical fluctuation (e.g. a single outlying
        # sample stretched the IQR but the run-to-run averages are
        # tight). Genuine slip widens BOTH spread AND CV. So we only
        # fire when CV also confirms elevated noise relative to coarse.
        if in_bisection:
            coarse_iqrs = []
            coarse_cvs = []
            for r in results[:-1]:
                if r.get('phase', 'coarse') not in ('coarse', 'ref'):
                    continue
                sg_p = r.get('sg') or {}
                if 'p25' in sg_p and 'p75' in sg_p:
                    coarse_iqrs.append(sg_p['p75'] - sg_p['p25'])
                rc = r.get('run_consistency') or {}
                if 'sg_cv' in rc:
                    coarse_cvs.append(rc['sg_cv'])
            if len(coarse_iqrs) >= 4 and len(coarse_cvs) >= 4:
                # Median of all but the last coarse step (which may be
                # the elevated step that triggered bisection).
                sorted_iqr = sorted(coarse_iqrs[:-1])
                n = len(sorted_iqr)
                median_coarse_iqr = (sorted_iqr[n // 2] if n % 2
                                     else (sorted_iqr[n // 2 - 1]
                                           + sorted_iqr[n // 2]) / 2)
                sorted_cv = sorted(coarse_cvs[:-1])
                m = len(sorted_cv)
                median_coarse_cv = (sorted_cv[m // 2] if m % 2
                                    else (sorted_cv[m // 2 - 1]
                                          + sorted_cv[m // 2]) / 2)
                last_rc = last.get('run_consistency') or {}
                last_cv = last_rc.get('sg_cv', 0.0)
                # CV cross-check: must be elevated vs coarse AND above
                # an absolute floor. If CV is below baseline-ish, it's
                # not slip — it's a wide-distribution fluke.
                if self.profile.IQR_VS_COARSE_REQUIRE_CV:
                    cv_confirms = (
                        last_cv >= self.profile.IQR_VS_COARSE_CV_FLOOR
                        and median_coarse_cv >= 0.5
                        and last_cv >=
                            (self.profile.IQR_VS_COARSE_CV_RATIO
                             * median_coarse_cv))
                else:
                    cv_confirms = True
                if (median_coarse_iqr >= 1
                        and last_iqr >= self.profile.IQR_VS_COARSE_MIN_ABS
                        and last_iqr >=
                            (self.profile.IQR_VS_COARSE_RATIO
                             * median_coarse_iqr)
                        and cv_confirms):
                    return ("%s spread IQR %.0f in %s vs coarse-phase "
                            "median IQR ~%.0f (%.1fx); CV %.1f%% also "
                            "elevated vs coarse %.1f%% — slip widening"
                            % (sg_label, last_iqr, last_phase,
                               median_coarse_iqr,
                               last_iqr / median_coarse_iqr,
                               last_cv, median_coarse_cv))

        return None

    def _check_sg_max_spike(self, results, sg_label):
        """Detect a no-load decoupling spike inside a single step.

        When the motor briefly decouples from the load during one of
        the repeats, SG snaps toward its no-load value. The 5-second
        per-run average + median over 5 runs absorbs the brief spike,
        but sg_max / sg_min records it.

        Direction is data-driven, not driver-hardcoded: across all
        three driver families (TMC5160 SG2, TMC2240 SG4, TMC2209 SG4)
        SG normally FALLS as load increases. So decoupling = SG snaps
        UP toward the no-load value, which means sg_max is the right
        signal regardless of driver. (We still keep a sg_min check as
        fallback in case a particular setup runs the other direction.)

        Fires in coarse / bisection / verify. Coarse uses tighter
        thresholds because natural sg_max variation is wider there.
        """
        last = results[-1]
        phase = last.get('phase', 'coarse')
        # Coarse path uses elevated thresholds (see below). The
        # legacy "bisect/verify only" gate has been removed because a
        # type of stick-slip stall produces compact medians/IQRs but
        # very large sg_max recovery spikes — invisible to the
        # standard triggers but obvious in sg_max if you look at it.
        sg = last.get('sg') or {}
        if 'max' not in sg or 'min' not in sg or 'median' not in sg:
            return None
        last_max = sg['max']
        last_min = sg['min']
        last_med = sg['median']

        # Build baseline of "typical" sg_max and sg_min from the
        # coarse phase. We need both because the no-load direction
        # depends on the actual SG-vs-load slope of THIS hardware.
        coarse_maxes = []
        coarse_mins = []
        coarse_meds = []
        for r in results[:-1]:
            if r.get('phase', 'coarse') not in ('coarse', 'ref'):
                continue
            sg_p = r.get('sg') or {}
            if 'max' in sg_p:
                coarse_maxes.append(sg_p['max'])
            if 'min' in sg_p:
                coarse_mins.append(sg_p['min'])
            if 'median' in sg_p:
                coarse_meds.append(sg_p['median'])
        if len(coarse_maxes) < 4 or len(coarse_meds) < 4:
            # Coarse gap-jump path can work with fewer steps because
            # the criterion is much stricter (gap_ratio AND abs_floor).
            # Need at least 3 prior steps to have a meaningful baseline.
            if (phase == 'coarse'
                    and len(coarse_maxes) >= 3
                    and len(coarse_meds) >= 3):
                pass  # continue with coarse-only check below
            else:
                return None

        # Determine SG-vs-load trend direction: compare first vs last
        # coarse-phase median (excluding the last coarse step which
        # likely triggered). If SG dropped from start to end, normal
        # = falling-with-load → decoupling = max spikes UP. Otherwise
        # = rising-with-load → decoupling = min drops DOWN.
        first_med = coarse_meds[0]
        # Use median of last 2-3 (excluding very last, which triggered)
        if len(coarse_meds) >= 3:
            late_meds = coarse_meds[-3:-1] if len(coarse_meds) > 3 else coarse_meds[-2:-1]
            late_med = sum(late_meds) / len(late_meds)
        else:
            late_med = coarse_meds[-1]
        sg_falls_with_load = late_med < first_med

        def _median(xs):
            s = sorted(xs); n = len(s)
            return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2

        if sg_falls_with_load:
            # Normal case for all common driver setups. Decoupling →
            # SG snaps UP. Compare sg_max to typical coarse sg_max.
            median_coarse_max = _median(coarse_maxes[:-1])
            # Use a separate, stricter "max-gap-jump" path tuned to
            # the coarse phase: we want to fire when the gap (max-med)
            # makes a big absolute jump above the prior coarse-baseline
            # gap. This is the signature seen in stick-slip stalls
            # where the medians stay compact but sg_max suddenly
            # records repeated decoupling spikes.
            prior_gaps = [m - md for m, md in zip(coarse_maxes[:-1],
                                                  coarse_meds[:-1])
                          if m > md]
            baseline_gap = _median(prior_gaps) if prior_gaps else 0
            current_gap = last_max - last_med
            # Pick threshold based on phase
            if phase == 'coarse':
                # Conservative: gap must be >=2x prior coarse baseline
                # AND the absolute floor must be cleared. Avoids early
                # noise-driven false positives.
                gap_ratio_threshold = (
                    self.profile.COARSE_GAP_JUMP_RATIO)
                gap_abs_floor = (
                    self.profile.COARSE_GAP_JUMP_ABS_FLOOR)
            else:
                # Bisect/verify: legacy a/b paths below
                gap_ratio_threshold = 0
                gap_abs_floor = 0

            absolute_floor = (self.profile.SG_MAX_ABS_GAP
                              if self.sg2_driver
                              else self.profile.SG_MAX_ABS_GAP_SG4)
            ratio_vs_coarse = (last_max
                               / max(median_coarse_max, 1))
            ratio_vs_med = last_max / max(last_med, 1)
            big_gap = (last_max - last_med) >= self.profile.SG_MAX_BIG_GAP

            # COARSE: gap-jump path — most sensitive to the
            #         stick-slip-with-compact-median signature.
            # Requires TWO consecutive steps with elevated gaps to
            # avoid false-positives from a single-sample sg_max
            # outlier (one stray spike in 564 samples = 0.18 % of
            # the signal which is just noise, not real slip).
            if phase == 'coarse':
                # Compute previous step's gap if it exists
                prev_gap = 0
                if len(coarse_maxes) >= 2 and len(coarse_meds) >= 2:
                    # coarse_maxes[-1] is the second-most-recent step
                    prev_gap = coarse_maxes[-1] - coarse_meds[-1]
                # Both prev and current must be ratio*baseline above
                # the historical baseline AND clear the absolute floor.
                prev_ratio = prev_gap / max(baseline_gap, 1)
                cur_ratio = current_gap / max(baseline_gap, 1)
                if (baseline_gap > 0
                        and cur_ratio >= gap_ratio_threshold
                        and prev_ratio >= gap_ratio_threshold
                        and current_gap >= gap_abs_floor
                        and prev_gap >= gap_abs_floor *
                            self.profile.COARSE_GAP_JUMP_PREV_FRACTION):
                    return ("%s max-median gap %.0f in this step is "
                            "%.1fx the prior coarse baseline (%.0f), "
                            "with the previous step (gap %.0f) also "
                            "elevated — persistent stall-recovery "
                            "spikes inside the runs (motor decoupling)"
                            % (sg_label, current_gap, cur_ratio,
                               baseline_gap, prev_gap))
                # Coarse stops here (don't fire on the bisect/verify
                # paths below)
                return None

            # BISECT/VERIFY: legacy two-path detection
            fires_a = (ratio_vs_med >= self.profile.SG_MAX_RATIO_TO_MEDIAN
                       and ratio_vs_coarse >=
                            self.profile.SG_MAX_RATIO_TO_COARSE
                       and last_max - last_med >= absolute_floor)
            fires_b = (ratio_vs_med >= self.profile.SG_MAX_BIG_RATIO
                       and big_gap)
            if fires_a or fires_b:
                return ("%s max %.0f in this step is %.1fx the "
                        "median (%.0f) and %.1fx the typical coarse "
                        "max (%.0f) — at least one repeat showed "
                        "motor decoupling spike (no-load value)"
                        % (sg_label, last_max, ratio_vs_med, last_med,
                           ratio_vs_coarse, median_coarse_max))
        else:
            # Rare case: SG rises with load. Decoupling → SG snaps
            # DOWN. Compare sg_min to typical coarse sg_min.
            if not coarse_mins:
                return None
            median_coarse_min = _median(coarse_mins[:-1])
            if (last_min > 0
                    and last_min <=
                        last_med / self.profile.SG_MIN_RATIO_TO_MEDIAN
                    and last_med - last_min >=
                        self.profile.SG_MIN_ABS_GAP):
                return ("%s min %.0f is %.1fx below the median (%.0f)"
                        " — at least one repeat dropped to no-load "
                        "value (decoupling)"
                        % (sg_label, last_min,
                           last_med / max(last_min, 1), last_med))

        return None

    def _check_run_outlier(self, results, sg_label):
        """Detect when one of the per-run averages is an outlier.

        Each measurement step runs N repeats (default 5). We have the
        list of per-run SG averages. If the run-to-run variance is
        dominated by a SINGLE outlying run (one repeat much off, the
        others tight), the median absorbs it but it's clear evidence
        of intermittent slip.

        Detection:
          - At least 4 valid runs
          - Find the run whose SG_avg deviates most from the median of
            the other runs
          - If that deviation is > 4x the median absolute deviation
            (MAD) of the rest AND >= 8 % of the median, it's an outlier

        If warmup was dropped during stats aggregation (run 1 deviated
        too much from the rest), this check operates only on the
        warmup-included runs (runs 2..N) — otherwise it would
        re-flag the warmup run as an outlier, defeating the purpose
        of the warmup-skip.

        Only fires in bisect or verify (where we care about
        reproducibility, not coarse exploration).
        """
        last = results[-1]
        if last.get('phase') not in ('bisect', 'verify'):
            return None
        all_run_sg = last.get('run_sg_avgs') or []
        warmup_dropped = last.get('warmup_dropped', False)

        # If warmup was excluded from stats, also exclude it from
        # outlier analysis. The displayed run number compensates for
        # the offset so reports refer to the original run index.
        if warmup_dropped and len(all_run_sg) >= 1:
            run_sg = all_run_sg[1:]
            run_idx_offset = 1
        else:
            run_sg = all_run_sg
            run_idx_offset = 0

        if len(run_sg) < 4:
            return None

        # Find the run furthest from the median of the rest
        n = len(run_sg)
        outlier_idx = -1
        outlier_dev = 0.0
        median_rest = 0.0
        mad_rest = 0.0
        for i in range(n):
            others = [run_sg[j] for j in range(n) if j != i]
            sorted_o = sorted(others)
            m = len(sorted_o)
            med_o = (sorted_o[m // 2] if m % 2
                     else (sorted_o[m // 2 - 1] + sorted_o[m // 2]) / 2)
            devs = [abs(x - med_o) for x in others]
            sorted_d = sorted(devs)
            mad_o = (sorted_d[m // 2] if m % 2
                     else (sorted_d[m // 2 - 1] + sorted_d[m // 2]) / 2)
            dev_i = abs(run_sg[i] - med_o)
            if dev_i > outlier_dev:
                outlier_dev = dev_i
                outlier_idx = i
                median_rest = med_o
                mad_rest = mad_o

        if outlier_idx < 0 or median_rest <= 0:
            return None

        # Outlier criterion: deviation >= 4x MAD AND >= 8% of median
        # (avoids tripping on tiny absolute differences when MAD is
        # near zero).
        if mad_rest < 0.5:
            mad_rest = 0.5  # noise floor
        if not (outlier_dev >= self.profile.OUTLIER_MAD_RATIO * mad_rest
                and outlier_dev >=
                    self.profile.OUTLIER_MIN_REL * median_rest):
            return None

        # SG-only outlier detection. Index reporting uses the
        # original run number (1-based), accounting for the warmup
        # offset so the user can identify which CSV row corresponds.
        original_run_num = outlier_idx + 1 + run_idx_offset
        return ("run %d is an outlier (SG_avg %.1f vs %.1f for the "
                "other repeats, %.1fx MAD) — at least one of the "
                "repeats stalled while the others ran clean"
                % (original_run_num, run_sg[outlier_idx],
                   median_rest, outlier_dev / mad_rest))

    def _is_borderline(self, results):
        """Detect if the latest bisection step sits in a 'gray zone'
        — not clearly safe but not clearly stalled either.

        Returns a string explaining why it's borderline, or None if
        the result is conclusive.

        Used by the bisection loop to decide whether to re-test the
        same flow before classifying it as safe or as a trigger.

        Heuristic: look at CV and IQR vs the coarse-phase baselines.
        Confidence ranges:
          - CV in coarse: typically 1-3% on a stable SG2 setup
          - CV >= 7% in bisection AND >= 3x coarse baseline: clearly slip
          - CV < 4% in bisection: clearly safe
          - CV in 4-6.9% range, OR 2-3x coarse baseline: BORDERLINE
          - IQR >= 25 in bisection: clearly slip (already trigger (c))
          - IQR < 15: clearly safe
          - IQR 15-24: BORDERLINE
        """
        if not results:
            return None
        last = results[-1]
        if last.get('phase') not in ('bisect', 'verify'):
            return None
        last_rc = last.get('run_consistency') or {}
        last_cv = last_rc.get('sg_cv')
        sg_stats = last.get('sg') or {}
        if 'p25' not in sg_stats or 'p75' not in sg_stats:
            return None
        last_iqr = sg_stats['p75'] - sg_stats['p25']

        # Build coarse-phase reference
        coarse_cvs = []
        coarse_iqrs = []
        for r in results[:-1]:
            if r.get('phase', 'coarse') not in ('coarse', 'ref'):
                continue
            rc = r.get('run_consistency') or {}
            sg_p = r.get('sg') or {}
            if 'sg_cv' in rc:
                coarse_cvs.append(rc['sg_cv'])
            if 'p25' in sg_p and 'p75' in sg_p:
                coarse_iqrs.append(sg_p['p75'] - sg_p['p25'])
        if len(coarse_cvs) < 4:
            return None  # not enough baseline for comparison

        # Use median-of-coarse (excluding last coarse step which may be
        # the elevated one that triggered bisection)
        sorted_cv = sorted(coarse_cvs[:-1])
        n = len(sorted_cv)
        coarse_med_cv = (sorted_cv[n // 2] if n % 2
                         else (sorted_cv[n // 2 - 1]
                               + sorted_cv[n // 2]) / 2)
        sorted_iqr = sorted(coarse_iqrs[:-1])
        n = len(sorted_iqr)
        coarse_med_iqr = (sorted_iqr[n // 2] if n % 2
                          else (sorted_iqr[n // 2 - 1]
                                + sorted_iqr[n // 2]) / 2)

        # CV gray zone: elevated but not screaming
        cv_borderline = (last_cv is not None
                         and self.profile.BORDER_CV_LOW <= last_cv
                                                       < self.profile.BORDER_CV_HIGH
                         and (coarse_med_cv < 1
                              or last_cv >=
                                  self.profile.BORDER_CV_RATIO
                                  * coarse_med_cv))

        # IQR gray zone: wide but not extreme. Cross-check with CV — a
        # widened IQR with a very LOW CV is usually just one outlier
        # sample stretching the percentiles, not real intermittent
        # slip. Genuine slip widens both spread and run-to-run noise.
        # We require CV at least BORDER_IQR_CV_FLOOR (more than typical
        # clean noise) for IQR-only borderlines to count.
        iqr_borderline = (self.profile.BORDER_IQR_LOW <= last_iqr
                                                      < self.profile.BORDER_IQR_HIGH
                          and (coarse_med_iqr < 1
                               or last_iqr >=
                                   self.profile.BORDER_IQR_RATIO
                                   * coarse_med_iqr)
                          and last_cv is not None
                          and last_cv >= self.profile.BORDER_IQR_CV_FLOOR)

        if cv_borderline and iqr_borderline:
            return ("CV %.1f%% (vs coarse ~%.1f%%) AND IQR %.0f "
                    "(vs coarse ~%.0f) both in borderline range"
                    % (last_cv, coarse_med_cv,
                       last_iqr, coarse_med_iqr))
        if cv_borderline:
            return ("CV %.1f%% in borderline range (coarse baseline "
                    "~%.1f%%, clear trigger >=7%%)"
                    % (last_cv, coarse_med_cv))
        if iqr_borderline:
            return ("IQR %.0f in borderline range (coarse baseline "
                    "~%.0f, clear trigger >=25)"
                    % (last_iqr, coarse_med_iqr))
        return None

    def _coarse_baseline_stats(self, results):
        """Compute "healthy" CV and IQR ranges from the coarse phase.

        Used by the HTML report to show the user what was considered
        normal during the slip-free portion of the test, so trigger
        events can be interpreted in context.

        Returns dict with keys:
            cv_min, cv_median, cv_max, cv_mean
            iqr_min, iqr_median, iqr_max, iqr_mean
            n_steps
        Returns None if there's not enough coarse data.
        """
        coarse_cvs = []
        coarse_iqrs = []
        for r in results:
            if r.get('phase', 'coarse') not in ('coarse', 'ref'):
                continue
            rc = r.get('run_consistency') or {}
            sg_p = r.get('sg') or {}
            if 'sg_cv' in rc:
                coarse_cvs.append(rc['sg_cv'])
            if 'p25' in sg_p and 'p75' in sg_p:
                coarse_iqrs.append(sg_p['p75'] - sg_p['p25'])
        # Drop the LAST coarse step (likely the one that triggered).
        # Only do so if we have enough data; otherwise keep them all.
        if len(coarse_cvs) >= 4:
            coarse_cvs = coarse_cvs[:-1]
        if len(coarse_iqrs) >= 4:
            coarse_iqrs = coarse_iqrs[:-1]

        if not coarse_cvs or not coarse_iqrs:
            return None

        def _median(xs):
            s = sorted(xs)
            n = len(s)
            return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2

        return {
            'cv_min': min(coarse_cvs),
            'cv_max': max(coarse_cvs),
            'cv_median': _median(coarse_cvs),
            'cv_mean': sum(coarse_cvs) / len(coarse_cvs),
            'iqr_min': min(coarse_iqrs),
            'iqr_max': max(coarse_iqrs),
            'iqr_median': _median(coarse_iqrs),
            'iqr_mean': sum(coarse_iqrs) / len(coarse_iqrs),
            'n_steps': len(coarse_cvs),
        }

    # ─── Helper: rotation_distance lookup ───────────────────────────

    def _get_rotation_distance(self, extruder):
        try:
            rd = extruder.extruder_stepper.stepper.get_rotation_distance()
            return rd[0] if isinstance(rd, tuple) else rd
        except AttributeError:
            pass
        try:
            rd = extruder.stepper.get_rotation_distance()
            return rd[0] if isinstance(rd, tuple) else rd
        except AttributeError:
            pass
        try:
            cfg = self.printer.lookup_object('configfile')
            settings = cfg.get_status(
                self.reactor.monotonic())['settings']
            return float(settings[self.stepper_name]['rotation_distance'])
        except Exception:
            return None

    # ─── TMC_FLOW_STATUS — diagnostic ───────────────────────────────

    def cmd_TMC_FLOW_STATUS(self, gcmd):
        """Show current TMC SG values and configuration status."""
        self._lookup_tmc()
        activate = gcmd.get_int('ACTIVATE', 1, minval=0, maxval=1)

        if activate:
            gcmd.respond_info("Activating motor (1 mm extrusion)...")
            self.gcode.run_script_from_command(
                "M83\nG1 E1 F60\nM400\nG92 E0")

        gcmd.respond_info(
            "===== TMC %s Status (stepper '%s') ====="
            % (self.driver_type or 'unknown', self.stepper_name))

        sg = self._read_sg()
        sg_label = self._get_sg_label()
        gcmd.respond_info(
            "%s: %s (range 0-510, lower = more load)"
            % (sg_label, str(sg) if sg is not None else 'n/a'))

        # Run config check
        problems, infos = self._check_tmc_config()
        if problems:
            gcmd.respond_info("Configuration issues found:")
            for fname, val, _desc in problems:
                gcmd.respond_info("  ⚠ %s = %s" % (fname, val))
        else:
            gcmd.respond_info("✓ Configuration looks good.")
        for info in infos:
            gcmd.respond_info(info)

    def _snapshot_tmc_settings(self):
        """Read TMC register state most relevant to this test.

        Captures only the fields that affect StallGuard / CoolStep
        behaviour during the flow test, so the report has a short and
        focused paper trail of the configuration that produced the
        result. Missing fields are skipped silently.
        """
        snapshot = []
        if self.tmc is None:
            return snapshot

        def get(name):
            try:
                return self.tmc.fields.get_field(name)
            except (KeyError, AttributeError):
                return None

        def add(label, name):
            v = get(name)
            if v is None:
                return
            snapshot.append((label, str(v), name))

        # Identity
        snapshot.append(("Driver", str(self.driver_type or '?'),
                         'driver_type'))
        snapshot.append(("Stepper", str(self.stepper_name or '?'),
                         'stepper_name'))

        # Chopper mode (decides which StallGuard engine works)
        add("en_pwm_mode (1=StealthChop)", 'en_pwm_mode')
        add("en_spreadcycle (1=SpreadCycle, TMC2209)", 'en_spreadcycle')

        # Velocity gates that govern StallGuard / CoolStep activity
        add("TPWMTHRS (StealthChop velocity threshold)", 'tpwmthrs')
        add("TCOOLTHRS (CoolStep / StallGuard velocity threshold)",
            'tcoolthrs')
        add("THIGH (full-step velocity threshold)", 'thigh')

        # StallGuard sensitivity
        add("SGTHRS (TMC2209 StallGuard4 threshold)", 'sgthrs')
        add("SGT (StallGuard2 threshold, signed -64..+63)", 'sgt')
        add("sg4_thrs (StallGuard4 threshold, TMC2240)", 'sg4_thrs')
        add("sg4_filt_en (SG4 filter enabled)", 'sg4_filt_en')
        add("SFILT (SG2 filter enabled)", 'sfilt')

        # CoolStep configuration
        add("SEMIN (CoolStep lower threshold = SEMIN*32)", 'semin')
        add("SEMAX (CoolStep upper threshold offset)", 'semax')
        add("SEUP (current-up step size)", 'seup')
        add("SEDN (current-down step size)", 'sedn')
        add("SEIMIN (CoolStep min current; 0=1/2 IRUN, 1=1/4)",
            'seimin')

        # Run-time current state
        add("IRUN (run current scale 0..31)", 'irun')

        return snapshot

    # ─── Auto-SGT tuning ────────────────────────────────────────────

    def _can_autotune_sgt(self):
        """Return True if the current driver supports automatic SGT
        tuning. Currently only signed-SGT (SG2) drivers are supported:
        TMC5160, TMC2130, TMC2240. TMC2209 (SG4) and TMC2660 are
        skipped — they have different field semantics.
        """
        if self.tmc is None:
            return False
        if self.is_2209:
            return False  # SG4 with SGTHRS — not yet supported
        # Verify the sgt field actually exists on this driver
        try:
            self.tmc.fields.get_field('sgt')
            return True
        except (KeyError, AttributeError):
            return False

    def _set_sgt(self, value):
        """Set the SGT field via SET_TMC_FIELD. Returns True on success."""
        try:
            self.gcode.run_script_from_command(
                "SET_TMC_FIELD STEPPER=%s FIELD=sgt VALUE=%d"
                % (self.stepper_name, int(value)))
            return True
        except Exception as e:
            logging.exception("tmc_flow_test: failed to set sgt: %s" % e)
            return False

    def _probe_sg_at_flow(self, gcmd, flow, duration, repeat,
                           skip_warmup=True):
        """Run a multi-repeat extrusion at the given flow and return
        SG statistics (median, p25, p75, cv). Used for the auto-SGT
        tuning probe. Wraps _measure_step so we get the full statistical
        treatment (warmup-skip, median, IQR, run-to-run CV) rather than
        a noisy single-shot average.

        skip_warmup: if True (default), the first repetition is
        excluded if it deviates from the rest. Set to False for
        diagnostic measurements where you want the unfiltered raw
        signal — useful when the SG signal itself is degenerate (e.g.
        TMC2209 SG4 producing only a few discrete values).

        Returns dict with 'median', 'cv', 'min', 'max' or None on failure.
        """
        # _measure_step does its own gcode run + sampling + warmup-skip
        result = self._measure_step(gcmd, flow, duration, repeat,
                                     skip_warmup=skip_warmup)
        if not result or not result.get('sg'):
            return None
        sg = result['sg']
        rc = result.get('run_consistency') or {}
        return {
            'median': sg.get('median'),
            'p25': sg.get('p25'),
            'p75': sg.get('p75'),
            'min': sg.get('min'),
            'max': sg.get('max'),
            'cv': rc.get('sg_cv', 0),
            'warmup_dropped': result.get('warmup_dropped', False),
        }

    def _autotune_sgt(self, gcmd, start_flow):
        """Pre-test SGT auto-tuning probe.

        Runs short low-flow probes and adjusts SGT until SG_RESULT
        sits in the profile's target range (SGT_LOW_TARGET_MIN .. MAX).
        This guarantees the main test starts with adequate dynamic
        range, which is the single biggest factor for accurate
        slip detection.

        Probe flow strategy: probes AT the user's START flow.
        Reason: SGT effects on the SG-vs-load curve are flow-dependent.
        If we probe at flow=5 but the test actually starts at flow=10,
        the SGT we tune for flow=5 will produce a less-than-target SG
        value at flow=10 — leaving the early test steps with too
        little dynamic range, which makes the natural SG-vs-load
        saturation curve look like a "plateau" trigger and produces
        false-positive slip detections in the coarse phase.
        Probing at start_flow makes the calibration directly
        relevant to where the test begins.

        SGT semantics (signed -64..+63):
          higher SGT  →  less sensitive  →  higher SG_RESULT values
          lower SGT   →  more sensitive  →  lower SG_RESULT values

        On exit, the SGT register holds the tuned value. Returns
        (original_sgt, final_sgt) so the caller can either restore the
        original (recommended for safety) or keep the tuned value
        for the test run. The console output recommends the tuned
        value for permanent inclusion in printer.cfg.
        """
        if not self._can_autotune_sgt():
            gcmd.respond_info(
                "Auto-SGT skipped: driver doesn't support sgt field "
                "tuning (TMC2209/SG4 not yet supported).")
            return None, None

        try:
            original_sgt = self.tmc.fields.get_field('sgt')
        except (KeyError, AttributeError):
            gcmd.respond_info(
                "Auto-SGT skipped: could not read current sgt value.")
            return None, None

        sgt_min = self.profile.SGT_RANGE_MIN
        sgt_max = self.profile.SGT_RANGE_MAX
        target_min = self.profile.SGT_LOW_TARGET_MIN
        target_max = self.profile.SGT_LOW_TARGET_MAX
        # Probe AT the user's START flow, not a fixed lower value.
        # The reason: SGT calibration is flow-dependent. If we probe
        # at flow=5 but the test starts at flow=10, the SGT we tuned
        # for flow=5 will produce a smaller SG range at flow=10 than
        # expected — leaving us in or near the target band ONLY at
        # flow=5, and below it from flow=10 onwards. That triggers
        # plateau false-positives later (the curve looks flat because
        # SG has already collapsed to the lower end of the scale).
        # Probing at start_flow guarantees the scale is calibrated
        # for the load the test actually begins with.
        probe_flow = start_flow
        probe_duration = self.profile.SGT_AUTOTUNE_PROBE_DURATION
        probe_repeats = self.profile.SGT_AUTOTUNE_PROBE_REPEATS
        max_iter = self.profile.SGT_AUTOTUNE_MAX_ITERATIONS

        gcmd.respond_info(
            "===== Auto-SGT tuning =====\n"
            "Probing SG_RESULT at %.1f mm³/s (%d reps × %.1f s)\n"
            "Target range at low load: SG_RESULT %d-%d (median)\n"
            "Current SGT: %d (range %d..%d)"
            % (probe_flow, probe_repeats, probe_duration,
               target_min, target_max,
               original_sgt, sgt_min, sgt_max))

        cur_sgt = original_sgt
        history = []  # [(sgt, sg_median, cv), ...]

        for iteration in range(max_iter):
            stats = self._probe_sg_at_flow(gcmd, probe_flow,
                                            probe_duration, probe_repeats)
            if stats is None or stats.get('median') is None:
                gcmd.respond_info(
                    "  iteration %d: SGT=%d → no SG samples (skipping)"
                    % (iteration + 1, cur_sgt))
                break

            sg_med = stats['median']
            sg_cv = stats.get('cv', 0)
            sg_max = stats.get('max', 0)
            history.append((cur_sgt, sg_med, sg_cv))

            # In-range check: median must be in target range AND no
            # samples saturated. Saturation at probe-flow guarantees
            # saturation at the test's actual start_flow (which is
            # higher than probe), so we treat any 1023 samples as
            # "SGT too high" even if the median looks acceptable.
            saturated = sg_max >= 1023
            in_range = (target_min <= sg_med <= target_max
                        and not saturated)

            if in_range:
                status = "✓ in range"
            elif saturated and target_min <= sg_med <= target_max:
                status = "saturated samples — SGT too high"
            else:
                status = "needs adjustment"
            sat_note = " [max=1023 saturated!]" if saturated else ""
            gcmd.respond_info(
                "  iteration %d: SGT=%d → SG median=%.0f, CV=%.1f%%, "
                "max=%.0f (%s)%s"
                % (iteration + 1, cur_sgt, sg_med, sg_cv, sg_max,
                   status, sat_note))

            if in_range:
                break

            # Decide adjustment direction. Saturation forces "lower SGT"
            # regardless of median (the median is misleading if upper
            # samples are clipped to 1023).
            if saturated:
                # Always lower SGT when samples saturated
                step = 5 if sg_med > 1000 else 3
                new_sgt = max(sgt_min, cur_sgt - step)
            elif sg_med < target_min:
                # Too sensitive — raise SGT
                step = 5 if sg_med < target_min / 2 else 3
                new_sgt = min(sgt_max, cur_sgt + step)
            else:
                # Above target_max but not saturated — lower SGT
                step = 3
                new_sgt = max(sgt_min, cur_sgt - step)

            if new_sgt == cur_sgt:
                gcmd.respond_info(
                    "  Reached SGT bound (%d) without entering target "
                    "range — stopping." % cur_sgt)
                break

            cur_sgt = new_sgt
            if not self._set_sgt(cur_sgt):
                gcmd.respond_info(
                    "  Failed to write SGT=%d — aborting auto-tune."
                    % cur_sgt)
                cur_sgt = original_sgt
                self._set_sgt(original_sgt)
                break

        else:
            gcmd.respond_info(
                "  Reached max iterations (%d) without entering target "
                "range." % max_iter)

        # Final report. The loop above uses (median in range AND not
        # saturated) as the success condition — and breaks out of the
        # loop on the first successful iteration. So if the loop
        # exited on success, the LAST history entry is the in-range
        # one. We use the same in-range definition here.
        last_sg = history[-1][1] if history else None
        final_sgt = cur_sgt
        in_target = (last_sg is not None
                     and target_min <= last_sg <= target_max)
        if in_target:
            if final_sgt != original_sgt:
                gcmd.respond_info(
                    "Auto-SGT: tuned to SGT=%d (was %d).\n"
                    "→ For a permanent fix, add to your "
                    "[%s extruder] section:\n"
                    "    driver_SGT: %d\n"
                    "Continuing test with tuned SGT..."
                    % (final_sgt, original_sgt,
                       self.driver_type, final_sgt))
            else:
                gcmd.respond_info(
                    "Auto-SGT: current SGT=%d already optimal — "
                    "no change needed."
                    % original_sgt)
        else:
            gcmd.respond_info(
                "Auto-SGT: could not reach target range. "
                "Continuing with SGT=%d (last SG median=%s).\n"
                "  → Consider manually tuning SGT and re-running."
                % (final_sgt,
                   "%.0f" % last_sg if last_sg is not None else "n/a"))

        return original_sgt, final_sgt

    def _save_sgt_to_config(self, gcmd, sgt_value):
        """Stage driver_SGT for Klipper's SAVE_CONFIG mechanism.

        The value lands in printer.cfg's autosave block when the user
        runs SAVE_CONFIG (which also restarts Klipper) — from then on
        it survives every restart, overriding any driver_SGT in the
        main config section.
        """
        section = "%s %s" % (self.driver_type, self.stepper_name)
        try:
            configfile = self.printer.lookup_object('configfile')
            configfile.set(section, 'driver_SGT', int(sgt_value))
            gcmd.respond_info(
                "SGT %d staged for [%s].\n"
                ">>> Run SAVE_CONFIG to write it into printer.cfg "
                "(Klipper restarts). The value is then permanent."
                % (int(sgt_value), section))
            return True
        except Exception as e:
            gcmd.respond_info(
                "Could not stage SGT for SAVE_CONFIG (%s). Add it "
                "manually to [%s]:\n    driver_SGT: %d"
                % (e, section, int(sgt_value)))
            return False

    def _restore_sgt_if_needed(self, gcmd, original_sgt, final_sgt,
                                keep_sgt, announce=True,
                                save_sgt=False):
        """Restore SGT to its original value unless keep_sgt/save_sgt
        is set or the auto-tune was a no-op. Safe to call multiple
        times."""
        if (original_sgt is None
                or final_sgt is None
                or final_sgt == original_sgt):
            return  # nothing to do
        if save_sgt:
            if announce:
                gcmd.respond_info(
                    "Auto-SGT: SAVE_SGT=1 — keeping tuned SGT=%d "
                    "active and staging it for SAVE_CONFIG."
                    % final_sgt)
            self._save_sgt_to_config(gcmd, final_sgt)
            return
        if keep_sgt:
            if announce:
                gcmd.respond_info(
                    "Auto-SGT: KEEP_SGT=1 — leaving SGT=%d active "
                    "until next FIRMWARE_RESTART." % final_sgt)
            return
        if self._set_sgt(original_sgt):
            if announce:
                gcmd.respond_info(
                    "Auto-SGT: restored original SGT=%d "
                    "(test ran with tuned SGT=%d). To make the "
                    "tuned value permanent, add to your config:\n"
                    "    driver_SGT: %d"
                    % (original_sgt, final_sgt, final_sgt))
        else:
            if announce:
                gcmd.respond_info(
                    "Auto-SGT: WARNING failed to restore original "
                    "SGT=%d. Run FIRMWARE_RESTART to reset."
                    % original_sgt)

    def cmd_TMC_FLOW_GUI(self, gcmd):
        """Write the interactive control panel (tmc_flow_gui.html) into
        the output directory, with the live config baked in.

        Everything here must survive a cold printer and a missing /
        unsupported driver — the whole point of the page is to TELL the
        user what is wrong, so it must never fail to be written.
        """
        src = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           'tmc_flow_gui.html')
        if not os.path.isfile(src):
            raise gcmd.error(
                "tmc_flow_test: GUI template not found at %s — pull the "
                "latest crydteamtools repo." % src)
        try:
            with open(src, encoding='utf-8') as f:
                html = f.read()
        except Exception as e:
            # Anything that is not a gcode.error would put Klipper into
            # shutdown ("Internal error on command"), so convert it.
            raise gcmd.error(
                "tmc_flow_test: cannot read the GUI template %s: %s"
                % (src, e))
        if '/*CFG*/null' not in html:
            # A silent no-op replace would ship a page claiming "no live
            # printer data" while this command reports success.
            raise gcmd.error(
                "tmc_flow_test: GUI template at %s has no /*CFG*/null "
                "placeholder — the file was modified or reformatted." % src)

        # Driver lookup is best-effort: _lookup_tmc raises when no
        # StallGuard-capable driver is present, which the page reports.
        try:
            self._lookup_tmc()
        except Exception:
            pass

        def field(name):
            if self.tmc is None:
                return None
            try:
                return self.tmc.fields.get_field(name)
            except (KeyError, AttributeError):
                return None

        problems = []
        try:
            probs, _infos = self._check_tmc_config()
            for fname, val, desc in probs:
                # First line of the description is the human summary.
                headline = str(desc).split('\n')[0]
                problems.append("%s (%s = %s)" % (headline, fname, val))
        except Exception as e:
            logging.debug("tmc_flow_test: GUI config check failed: %s", e)

        hotend_temp = None
        try:
            extruder = self.printer.lookup_object('extruder', None)
            if extruder is not None:
                hotend_temp = round(
                    extruder.get_heater().get_temp(
                        self.reactor.monotonic())[0], 1)
        except Exception:
            pass

        cfg = {
            'version': MODULE_VERSION,
            'driver': self.driver_type,
            'stepper': self.stepper_name,
            # Driver-derived values are only meaningful once the lookup
            # succeeded — before that the class defaults would report a
            # plausible but wrong threshold field.
            'sgt': (field(self._get_sg_threshold_field_name())
                    if self.tmc is not None else None),
            'can_autotune_sgt': bool(self._can_autotune_sgt()),
            'problems': problems,
            'filament_diameter': self.filament_diameter,
            'melt_zone_length': self.melt_zone_length,
            'min_hotend_temp': self.min_hotend_temp,
            'reference_flow': self.reference_flow,
            'test_fan_speed': self.test_fan_speed,
            'hotend_temp': hotend_temp,
            'output_dir': self.output_dir,
        }
        # json.dumps escapes for a JS string, not for HTML: a value
        # containing "</script>" would end the script block.
        cfg_json = json.dumps(cfg).replace('</', '<\\/')
        html = html.replace('/*CFG*/null', cfg_json)

        try:
            os.makedirs(self.output_dir, exist_ok=True)
        except Exception as e:
            raise gcmd.error(
                "tmc_flow_test: cannot create output dir %s: %s"
                % (self.output_dir, e))
        dst = os.path.join(self.output_dir, 'tmc_flow_gui.html')
        # Write via a temp file so an interrupted write can't leave a
        # half-written control panel behind. os.makedirs succeeds on an
        # existing but unwritable directory, so the real permission /
        # ENOSPC / read-only-remount failures surface HERE — and an
        # uncaught OSError would shut the printer down.
        tmp = dst + '.tmp'
        try:
            with open(tmp, 'w', encoding='utf-8') as f:
                f.write(html)
            os.replace(tmp, dst)
        except Exception as e:
            try:
                os.remove(tmp)
            except Exception:
                pass
            raise gcmd.error(
                "tmc_flow_test: cannot write the control panel to %s: %s"
                % (dst, e))
        gcmd.respond_info(
            "Control panel written to:\n  %s\n"
            "Open it via your web UI's file browser (Flowtest folder) or "
            "any browser. Re-run TMC_FLOW_GUI after config or driver "
            "changes to refresh the baked-in values." % dst)

    def cmd_TMC_FLOW_SGT_CALIBRATE(self, gcmd):
        """Standalone SGT tuning with persistence.

        Runs the same probe-based auto-tune the main test uses,
        anchored at a low reference flow, then stages the result as
        driver_SGT for SAVE_CONFIG. SAVE=0 keeps the value active for
        this session only (like KEEP_SGT on the main test).
        """
        default_flow = (self.reference_flow
                        if self.reference_flow > 0 else 15.0)
        flow = gcmd.get_float('FLOW', default_flow, above=0.)
        save = gcmd.get_int('SAVE', 1, minval=0, maxval=1)

        if not self._can_autotune_sgt():
            raise gcmd.error(
                "SGT auto-tune not supported on this driver "
                "(needs the signed sgt field: TMC5160/2130/2240).")

        # Hotend must be at printing temperature — this extrudes.
        extruder = self.printer.lookup_object('extruder')
        heater = extruder.get_heater()
        cur_temp, _ = heater.get_temp(self.reactor.monotonic())
        target_temp = heater.target_temp
        if cur_temp < self.min_hotend_temp:
            raise gcmd.error(
                "Hotend too cold: %.1f°C (min %.1f°C). Heat it to "
                "printing temperature first."
                % (cur_temp, self.min_hotend_temp))
        if target_temp > 0 and cur_temp < target_temp - 5.0:
            raise gcmd.error(
                "Hotend not at target: %.1f°C (target %.1f°C)."
                % (cur_temp, target_temp))

        gcmd.respond_info(
            "──── SGT calibration ────\n"
            "  driver: %s on %s | probe flow: %.1f mm³/s [FLOW]\n"
            "  result is %s"
            % (self.driver_type, self.stepper_name, flow,
               "staged for SAVE_CONFIG [SAVE=1]" if save
               else "session-only [SAVE=0]"))

        self.gcode.run_script_from_command("M83\nG92 E0")
        original_sgt, final_sgt = self._autotune_sgt(gcmd, flow)
        if original_sgt is None or final_sgt is None:
            gcmd.respond_info(
                "SGT calibration did not produce a result — see the "
                "messages above.")
            return

        if final_sgt == original_sgt:
            gcmd.respond_info(
                "SGT %d already in the target range — nothing to "
                "change." % original_sgt)
            return

        if save:
            self._save_sgt_to_config(gcmd, final_sgt)
        else:
            gcmd.respond_info(
                "Tuned SGT=%d stays active until FIRMWARE_RESTART "
                "(was %d). To persist it, re-run with SAVE=1 or add:\n"
                "    driver_SGT: %d" % (final_sgt, original_sgt,
                                        final_sgt))

    # ─── TMC_FLOW_TEST_SG_VARIANTS ─────────────────────────────────
    # Diagnostic for TMC2209 only. Trinamic states StallGuard4 is
    # "intended for StealthChop mode, only" but doesn't categorically
    # state it returns nothing in SpreadCycle. Different production
    # batches and motor combinations give wildly different behaviour
    # in SpreadCycle. This command runs the same probe in BOTH chopper
    # modes and tells the user empirically what their hardware does.
    #
    # Why this matters: TMC2209 in SpreadCycle delivers ~50 % more
    # peak torque than in StealthChop, so if a particular hardware
    # combo happens to produce a usable SG signal in SpreadCycle, the
    # user could get significantly higher max-flow results.

    def _save_tmc2209_chopper_state(self):
        """Snapshot TMC2209 chopper-mode-relevant fields for restore.

        TMC2209 uses different field names than TMC5160:
          - TMC2209 GCONF Bit 2 = 'en_spreadcycle' (1=SpreadCycle, 0=StealthChop)
          - TMC5160 GCONF Bit 2 = 'en_pwm_mode'    (1=StealthChop, 0=SpreadCycle)
        Note the inverted semantics! This method is TMC2209-only.

        Returns dict with 'en_spreadcycle', 'tpwmthrs', 'sgthrs',
        'tcoolthrs' (or None values where read failed).
        """
        snapshot = {}
        for field in ('en_spreadcycle', 'tpwmthrs',
                       'sgthrs', 'tcoolthrs'):
            try:
                snapshot[field] = self.tmc.fields.get_field(field)
            except (KeyError, AttributeError):
                snapshot[field] = None
        return snapshot

    def _set_tmc_field_safe(self, field, value):
        """SET_TMC_FIELD wrapper. Returns True on success."""
        try:
            self.gcode.run_script_from_command(
                "SET_TMC_FIELD STEPPER=%s FIELD=%s VALUE=%d"
                % (self.stepper_name, field, int(value)))
            return True
        except Exception as e:
            logging.exception(
                "tmc_flow_test: failed to set %s: %s" % (field, e))
            return False

    def _restore_tmc2209_chopper(self, snapshot, gcmd):
        """Restore the chopper state captured by snapshot. Best-effort:
        logs failures but tries to set every field that has a value."""
        ok = True
        for field in ('en_spreadcycle', 'tpwmthrs',
                       'sgthrs', 'tcoolthrs'):
            if snapshot.get(field) is not None:
                if not self._set_tmc_field_safe(field, snapshot[field]):
                    ok = False
        if ok:
            gcmd.respond_info(
                "Original chopper config restored.")
        else:
            gcmd.respond_info(
                "WARNING: chopper restore had failures. Run "
                "FIRMWARE_RESTART to fully reset.")
        return ok

    def _evaluate_sg4_quality(self, low_stats, high_stats):
        """Decide if a (low-flow, high-flow) pair of probes shows a
        usable SG4 signal for slip detection.

        IMPORTANT: SG4 on TMC2209 behaves differently from SG2. The
        original SG2 logic (load drops SG_RESULT proportionally) does
        NOT apply. Empirically, SG_RESULT on TMC2209:
          - Sticks in a "bias region" (6-22) at low velocity
          - Often INCREASES (not decreases) with load
          - But CV (run-to-run variance) reliably EXPLODES at slip

        We declare the signal USABLE if ANY of the following holds —
        each one alone is enough evidence that downstream slip
        detection has something to work with:

          (A) CV-spike: high_cv ≥ 10 % AND CV ratio ≥ 3×
              The CV-based triggers will fire at real slip.

          (B) Range separation: P25-P75 boxes of the two probes do
              NOT overlap. This means even single noisy samples can
              be classified as "low-flow region" vs "high-flow region"
              — the strongest possible evidence the signal carries
              load information.

          (C) Magnitude change: |delta| ≥ 30 SG-units. Down from the
              previous 50 — anything ≥ 30 is well above typical
              measurement noise (~3–5 units stddev) and gives
              magnitude-based triggers headroom.

          (D) Ratio change: high/low ratio ≥ 2× (or ≤ 0.5× inverse).
              Catches cases where absolute delta is small but the
              relative change is large (e.g. low=10, high=25:
              delta=15 but ratio=2.5×).

        We accept either direction of SG-median change — some TMC2209
        hardware shows inverted behaviour (SG up at load). The
        plugin's CV/IQR-based triggers work regardless of SG direction.

        Returns dict with 'usable' (bool), 'reason' (str), 'delta'
        (signed median difference high − low), 'cv_ratio' (high_cv /
        low_cv).
        """
        if not low_stats or not high_stats:
            return {'usable': False,
                    'reason': 'one or both probes returned no samples',
                    'delta': None, 'cv_ratio': None}

        low_med = low_stats.get('median', 0)
        high_med = high_stats.get('median', 0)
        low_p25 = low_stats.get('p25', low_med)
        low_p75 = low_stats.get('p75', low_med)
        high_p25 = high_stats.get('p25', high_med)
        high_p75 = high_stats.get('p75', high_med)
        low_cv = low_stats.get('cv', 100)
        high_cv = high_stats.get('cv', 100)
        delta = high_med - low_med
        # Avoid division by zero for very small CV
        cv_ratio = high_cv / max(low_cv, 0.5)
        # Median ratio (always in [1, ∞) by symmetric definition)
        if low_med > 0 and high_med > 0:
            med_ratio = max(high_med / low_med, low_med / high_med)
        else:
            med_ratio = 1.0

        # Both probes returning zero means no signal at all
        if low_med == 0 and high_med == 0:
            return {'usable': False,
                    'reason': 'SG_RESULT stuck at 0 in both probes — '
                              'driver returns no SG signal',
                    'delta': delta, 'cv_ratio': cv_ratio}

        # Path A: CV-spike (primary slip signature on TMC2209)
        if cv_ratio >= 3.0 and high_cv >= 10.0:
            return {'usable': True,
                    'reason': ('CV-spike detected: CV jumped from '
                               '%.1f%% (low) to %.1f%% (high), %.1fx '
                               'increase. CV-based slip triggers '
                               'will catch real slip events.'
                               % (low_cv, high_cv, cv_ratio)),
                    'delta': delta, 'cv_ratio': cv_ratio}

        # Path B: P25-P75 ranges of the two probes don't overlap.
        # This is the strongest possible evidence of load sensitivity:
        # even individual samples are classifiable as "low" vs "high".
        ranges_separate = (low_p75 < high_p25) or (high_p75 < low_p25)
        if ranges_separate:
            return {'usable': True,
                    'reason': ('clean range separation: low-flow P25-P75 '
                               '= %d-%d, high-flow P25-P75 = %d-%d. The '
                               'two distributions don\'t overlap — '
                               'strong load-correlated signal.'
                               % (low_p25, low_p75, high_p25, high_p75)),
                    'delta': delta, 'cv_ratio': cv_ratio}

        # Path C: large absolute magnitude change
        if abs(delta) >= 30:
            direction = ('drops' if delta < 0 else 'rises')
            return {'usable': True,
                    'reason': ('SG-median %s by %d units between low '
                               'and high flow (well above noise floor). '
                               'Signal is responsive to load (note: '
                               'TMC2209 SG direction can be inverted '
                               'vs. SG2 spec — this is normal).'
                               % (direction, abs(delta))),
                    'delta': delta, 'cv_ratio': cv_ratio}

        # Path D: relative ratio change (catches small-magnitude /
        # high-relative-change cases)
        if med_ratio >= 2.0:
            return {'usable': True,
                    'reason': ('SG-median changed by %.1fx between low '
                               'and high flow (%d → %d). Relative '
                               'change is large enough for magnitude-'
                               'based triggers to work.'
                               % (med_ratio, low_med, high_med)),
                    'delta': delta, 'cv_ratio': cv_ratio}

        # Neither CV nor SG-median moved enough — signal is dead
        if low_cv > 25 or high_cv > 25:
            return {'usable': False,
                    'reason': ('high noise without load-correlated '
                               'change: CV %.1f%% / %.1f%%, SG delta '
                               'only %d, ratio %.1fx. Pure noise, '
                               'not signal.'
                               % (low_cv, high_cv, delta, med_ratio)),
                    'delta': delta, 'cv_ratio': cv_ratio}
        return {'usable': False,
                'reason': ('signal too flat: SG delta = %+d, ratio '
                           '%.1fx, CV ratio %.1fx, ranges %d-%d vs '
                           '%d-%d (overlap). Need any of: |delta|≥30, '
                           'ratio≥2×, range separation, or CV-spike. '
                           'Try higher HIGH_FLOW or higher run_current.'
                           % (delta, med_ratio, cv_ratio,
                              low_p25, low_p75, high_p25, high_p75)),
                'delta': delta, 'cv_ratio': cv_ratio}

    def cmd_TMC_FLOW_TEST_SG_VARIANTS(self, gcmd):
        """Empirically probe whether TMC2209 SG4 produces usable
        readings in SpreadCycle mode (in addition to the documented
        StealthChop mode).

        Test method:
          1. Save current chopper state
          2. Force StealthChop, probe SG at LOW_FLOW and HIGH_FLOW
          3. Force SpreadCycle, probe SG at LOW_FLOW and HIGH_FLOW
          4. Compare: usable signal needs (a) non-zero median, (b)
             load-sensitive Δmedian ≥ 30, (c) run-to-run CV < 25 %
          5. Restore original chopper state
          6. Print verdict and recommendation
        """
        self._lookup_tmc()
        if self.driver_type != 'tmc2209':
            gcmd.respond_info(
                "TMC_FLOW_TEST_SG_VARIANTS is only meaningful for "
                "TMC2209 (current driver: %s). For TMC5160/TMC2240 "
                "the plugin already runs in the optimal SG2 + "
                "SpreadCycle mode." % self.driver_type)
            return

        low_flow = gcmd.get_float('LOW_FLOW', 5.0,
                                   minval=1.0, maxval=30.0)
        high_flow = gcmd.get_float('HIGH_FLOW', 20.0,
                                    minval=5.0, maxval=150.0)
        duration = gcmd.get_float('DURATION', 5.0,
                                   minval=2.0, maxval=15.0)
        repeat = gcmd.get_int('REPEAT', 5, minval=3, maxval=10)
        sgthrs = gcmd.get_int('SGTHRS', 100, minval=0, maxval=255)

        if high_flow <= low_flow:
            raise gcmd.error(
                "HIGH_FLOW (%.1f) must be greater than LOW_FLOW (%.1f)"
                % (high_flow, low_flow))

        # Hotend-temp safety check — same pattern as the main test
        extruder = self.printer.lookup_object('extruder')
        heater = extruder.get_heater()
        cur_temp, _ = heater.get_temp(self.reactor.monotonic())
        target_temp = heater.target_temp
        if cur_temp < self.min_hotend_temp:
            raise gcmd.error(
                "Hotend too cold: %.1f °C < %.1f °C minimum"
                % (cur_temp, self.min_hotend_temp))

        gcmd.respond_info(
            "===== TMC2209 SG_RESULT chopper-mode comparison =====\n"
            "This test will run %d × %.1f s probes at flow=%.1f mm³/s\n"
            "and flow=%.1f mm³/s in BOTH StealthChop and SpreadCycle.\n"
            "It empirically determines whether SG_RESULT works in\n"
            "SpreadCycle on YOUR specific hardware.\n"
            "Expected runtime: ~%.0f minutes.\n"
            "Hotend: %.1f °C / target %.1f °C\n"
            "------------------------------------------------"
            % (repeat, duration, low_flow, high_flow,
               (4 * repeat * (duration + 0.5)) / 60.0,
               cur_temp, target_temp))

        # Snapshot current chopper config
        snapshot = self._save_tmc2209_chopper_state()
        gcmd.respond_info(
            "Saved chopper state: en_spreadcycle=%s tpwmthrs=%s "
            "sgthrs=%s tcoolthrs=%s"
            % (snapshot['en_spreadcycle'], snapshot['tpwmthrs'],
               snapshot['sgthrs'], snapshot['tcoolthrs']))

        # Set extruder relative + zero
        self.gcode.run_script_from_command("M83\nG92 E0\nM400")
        # Ensure SGTHRS is a sane value for the test (irrelevant for
        # SG_RESULT reads themselves, but needed if user has it at 0)
        self._set_tmc_field_safe('sgthrs', sgthrs)
        # TCOOLTHRS must be > 0 for SG to be active (Trinamic spec).
        # Use 0xFFFFF (= 1048575, max 20-bit value) so SG is active
        # at any motor velocity.
        self._set_tmc_field_safe('tcoolthrs', 0xFFFFF)
        # Sync to ensure the fields are committed before we proceed
        self.gcode.run_script_from_command("M400")

        results = {}
        try:
            # ─── Phase A: StealthChop ──────────────────────────────
            # TMC2209 GCONF Bit 2 'en_spreadcycle':
            #   0 = StealthChop (Trinamic default for TMC2209)
            #   1 = SpreadCycle
            # TPWMTHRS = 0xFFFFF means "never auto-switch to SpreadCycle
            # at higher velocity" — we want pure StealthChop here.
            gcmd.respond_info(">>> Phase A: StealthChop (documented)")
            self._set_tmc_field_safe('en_spreadcycle', 0)
            self._set_tmc_field_safe('tpwmthrs', 0xFFFFF)
            # Settle delay + sync. Re-assert M83 + G92 E0 because some
            # SET_TMC_FIELD operations can disturb the extruder state
            # (or Klipper may reset it as part of the field write).
            self.gcode.run_script_from_command(
                "M400\nG4 P1500\nM83\nG92 E0\nM400")

            gcmd.respond_info(
                "  StealthChop: low-flow probe at %.1f mm³/s..." % low_flow)
            sc_low = self._probe_sg_at_flow(
                gcmd, low_flow, duration, repeat, skip_warmup=False)
            # Reset extruder state between probes too — the rep-loop
            # in _measure_step relies on M83 still being active.
            self.gcode.run_script_from_command(
                "M400\nG4 P2000\nM83\nG92 E0\nM400")
            gcmd.respond_info(
                "  StealthChop: high-flow probe at %.1f mm³/s..." % high_flow)
            sc_high = self._probe_sg_at_flow(
                gcmd, high_flow, duration, repeat, skip_warmup=False)

            results['stealthchop'] = {
                'low': sc_low, 'high': sc_high,
                'eval': self._evaluate_sg4_quality(sc_low, sc_high),
            }
            # Inter-phase rest
            self.gcode.run_script_from_command(
                "M400\nG4 P3000\nM83\nG92 E0\nM400")

            # ─── Phase B: SpreadCycle ──────────────────────────────
            # Force SpreadCycle by setting en_spreadcycle=1.
            # TPWMTHRS irrelevant when en_spreadcycle=1 since GCONF
            # already forces SpreadCycle.
            gcmd.respond_info(">>> Phase B: SpreadCycle (experimental)")
            self._set_tmc_field_safe('en_spreadcycle', 1)
            self._set_tmc_field_safe('tpwmthrs', 0)
            self.gcode.run_script_from_command(
                "M400\nG4 P1500\nM83\nG92 E0\nM400")

            gcmd.respond_info(
                "  SpreadCycle: low-flow probe at %.1f mm³/s..." % low_flow)
            sp_low = self._probe_sg_at_flow(
                gcmd, low_flow, duration, repeat, skip_warmup=False)
            self.gcode.run_script_from_command(
                "M400\nG4 P2000\nM83\nG92 E0\nM400")
            gcmd.respond_info(
                "  SpreadCycle: high-flow probe at %.1f mm³/s..." % high_flow)
            sp_high = self._probe_sg_at_flow(
                gcmd, high_flow, duration, repeat, skip_warmup=False)

            results['spreadcycle'] = {
                'low': sp_low, 'high': sp_high,
                'eval': self._evaluate_sg4_quality(sp_low, sp_high),
            }
        finally:
            # ─── Always restore original chopper state ─────────────
            self._restore_tmc2209_chopper(snapshot, gcmd)

        # ─── Report ────────────────────────────────────────────────
        def fmt_probe(p):
            if not p or p.get('median') is None:
                return "no samples"
            return ("median=%.0f, p25=%.0f, p75=%.0f, min=%.0f, max=%.0f, "
                    "CV=%.1f%%"
                    % (p.get('median', 0), p.get('p25', 0),
                       p.get('p75', 0), p.get('min', 0),
                       p.get('max', 0), p.get('cv', 0)))

        gcmd.respond_info(
            "================================================\n"
            "===== Results =====")

        for mode_name, mode_label in [('stealthchop', 'StealthChop'),
                                      ('spreadcycle', 'SpreadCycle')]:
            data = results.get(mode_name, {})
            ev = data.get('eval', {})
            usable = ev.get('usable', False)
            verdict = "✓ USABLE" if usable else "✗ NOT USABLE"
            gcmd.respond_info(
                "%s: %s\n"
                "  flow=%.1f mm³/s: %s\n"
                "  flow=%.1f mm³/s: %s\n"
                "  → %s"
                % (mode_label, verdict, low_flow,
                   fmt_probe(data.get('low')),
                   high_flow,
                   fmt_probe(data.get('high')),
                   ev.get('reason', 'no evaluation')))

        # ─── Recommendation ────────────────────────────────────────
        sc_eval = results['stealthchop']['eval']
        sp_eval = results['spreadcycle']['eval']
        sc_ok = sc_eval.get('usable', False)
        sp_ok = sp_eval.get('usable', False)

        gcmd.respond_info("===== Recommendation =====")
        if sc_ok and sp_ok:
            sc_cvr = sc_eval.get('cv_ratio') or 0
            sp_cvr = sp_eval.get('cv_ratio') or 0
            # Both work — pick the one with cleaner slip signature
            if sp_cvr > sc_cvr * 1.2:
                gcmd.respond_info(
                    "✓ Both StealthChop AND SpreadCycle produce a\n"
                    "  usable slip signal on YOUR hardware.\n"
                    "  SpreadCycle has stronger CV-spike at slip\n"
                    "  (CV ratio %.1fx vs %.1fx in StealthChop).\n"
                    "  → SpreadCycle gives full motor torque AND a\n"
                    "    detectable slip signal — use it for max flow:\n"
                    "      stealthchop_threshold: 0\n"
                    "      driver_SGTHRS: %d\n"
                    "  → Run TMC_FLOW_FIND_MAX MAX=150 START=30\n"
                    "  → CAUTION: SG4 in SpreadCycle is unsupported by\n"
                    "    Trinamic spec. Validate with several long\n"
                    "    prints before trusting the slicer value."
                    % (sp_cvr, sc_cvr, sgthrs))
            else:
                gcmd.respond_info(
                    "✓ Both StealthChop AND SpreadCycle produce a\n"
                    "  usable slip signal on YOUR hardware.\n"
                    "  StealthChop has stronger or similar CV-spike\n"
                    "  (CV ratio %.1fx vs %.1fx).\n"
                    "  → Use the documented StealthChop mode:\n"
                    "      stealthchop_threshold: 999999\n"
                    "      driver_SGTHRS: %d\n"
                    "  → Run TMC_FLOW_FIND_MAX MAX=150 START=30"
                    % (sc_cvr, sp_cvr, sgthrs))
        elif sc_ok:
            gcmd.respond_info(
                "✓ StealthChop produces a usable slip signal on\n"
                "  YOUR hardware (CV ratio %.1fx between low/high flow).\n"
                "✗ SpreadCycle did not show a clean slip signature.\n"
                "  → Use the documented StealthChop config:\n"
                "      stealthchop_threshold: 999999\n"
                "      driver_SGTHRS: %d\n"
                "  → Run TMC_FLOW_FIND_MAX MAX=150 START=30\n"
                "    (START=30 because TMC2209 SG4 has a low-velocity\n"
                "    bias region — start above it for clean readings)"
                % (sc_eval.get('cv_ratio') or 0, sgthrs))
        elif sp_ok:
            gcmd.respond_info(
                "✓ SpreadCycle produces a usable slip signal on YOUR\n"
                "  hardware (CV ratio %.1fx between low/high flow).\n"
                "✗ StealthChop did not show a clean slip signature.\n"
                "  → This is unusual but the test confirms SpreadCycle\n"
                "    works for slip detection here.\n"
                "  → Force SpreadCycle in your config:\n"
                "      stealthchop_threshold: 0\n"
                "      driver_SGTHRS: %d\n"
                "  → Run TMC_FLOW_FIND_MAX MAX=150 START=30\n"
                "  → CAUTION: SG4 in SpreadCycle is unsupported per\n"
                "    Trinamic spec. Validate before relying on result."
                % (sp_eval.get('cv_ratio') or 0, sgthrs))
        else:
            gcmd.respond_info(
                "✗ Neither chopper mode produced a usable slip signal\n"
                "  on this hardware combination.\n"
                "  Diagnostics to try:\n"
                "  1. Was the motor actually slipping at HIGH_FLOW?\n"
                "     If not, increase HIGH_FLOW (try 150) so the\n"
                "     test reaches actual stall.\n"
                "  2. Increase run_current (e.g. 1.0 A) for more\n"
                "     load differentiation.\n"
                "  3. Check with DUMP_TMC STEPPER=extruder during a\n"
                "     long G1 E... command — does SG_RESULT change?\n"
                "  4. If SG_RESULT stays at 0 or one fixed value:\n"
                "     hardware-side SG4 problem (some TMC2209 clones\n"
                "     have non-functional SG4). Try a different\n"
                "     TMC2209 board, or switch to TMC2240.")

    # ─── Fan + extra-thermistor helpers ─────────────────────────────

    def _capture_current_fan_speed(self):
        """Return current part-cooling fan speed as 0-100 percent.
        Returns 0.0 if no fan is found or status can't be read.
        """
        fan_obj = self._lookup_fan_object()
        if fan_obj is None:
            return 0.0
        try:
            now = self.reactor.monotonic()
            status = fan_obj.get_status(now)
            speed_frac = status.get('speed', 0.0)
            return float(speed_frac) * 100.0
        except Exception:
            return 0.0

    def _set_part_cooling_fan(self, percent, gcmd, quiet=False):
        """Set part-cooling fan to a given 0-100 percent.
        Uses M106 by default (drives the printer's primary fan).
        Falls back to no-op if no fan is configured.
        """
        # Clamp into valid M106 range
        percent = max(0.0, min(100.0, float(percent)))
        m106_value = int(round(percent / 100.0 * 255.0))
        try:
            if m106_value <= 0:
                self.gcode.run_script_from_command("M107")
            else:
                self.gcode.run_script_from_command(
                    "M106 S%d" % m106_value)
        except Exception as e:
            if not quiet:
                gcmd.respond_info(
                    "Could not set fan speed to %.0f %% (%s) — "
                    "continuing without fan control." % (percent, e))

    def _lookup_fan_object(self):
        """Best-effort look-up of the printer's part-cooling fan object.

        Tries:
        1. The user-named object via [tmc_flow_test] fan_object_name.
        2. 'fan' (Klipper's built-in [fan] section).
        3. Any 'fan_generic <name>' that mentions 'part' in its name.
        Returns None if nothing usable was found.
        """
        if self.fan_object_name:
            try:
                return self.printer.lookup_object(self.fan_object_name)
            except Exception:
                pass
        # Default: M106 controls the [fan] section (Klipper's
        # built-in primary part-cooling fan).
        try:
            return self.printer.lookup_object('fan')
        except Exception:
            pass
        # Heuristic fallback: any fan_generic with 'part' in the name
        try:
            objs = self.printer.lookup_objects('fan_generic')
            for name, obj in objs:
                if 'part' in name.lower():
                    return obj
        except Exception:
            pass
        return None

    def _resolve_extra_thermistors(self, gcmd):
        """Verify configured extra thermistors actually exist.
        Logs a notice for each one resolved, and a warning for
        each one that doesn't exist.
        """
        if not self.extra_thermistors:
            return
        ok = []
        missing = []
        for name in self.extra_thermistors:
            try:
                obj = self.printer.lookup_object(name)
                if hasattr(obj, 'get_temp') or hasattr(obj, 'last_temp'):
                    ok.append(name)
                else:
                    missing.append((name, "no temperature interface"))
            except Exception as e:
                missing.append((name, str(e)))
        if ok:
            gcmd.respond_info(
                "Extra thermistors logged this run: %s"
                % ", ".join(ok))
        for name, why in missing:
            gcmd.respond_info(
                "Extra thermistor '%s' not available (%s) — skipping."
                % (name, why))

    # ─── Main commands ──────────────────────────────────────────────

    def cmd_TMC_FLOW_FIND_MAX(self, gcmd):
        """Run the StallGuard-based flow test."""
        self._lookup_tmc()
        self._mpc_state = None
        self._fan_state = []
        try:
            self._run_find_max(gcmd)
        except _HardStop:
            pass  # already reported and cleaned up
        finally:
            if getattr(self, '_fan_state', None):
                self._chamber_fans_restore(gcmd, self._fan_state)
            # MATERIAL= changed MPC: never leave it changed, even if the
            # run aborted before its own cleanup.
            st = getattr(self, '_mpc_state', None)
            if st is not None and st.get('applied'):
                self._mpc_restore(gcmd, st)

    def _run_find_max(self, gcmd):
        """Common bracket-bisection algorithm using SG triggers."""
        self._lookup_tmc()

        start_flow = gcmd.get_float('START', 10.0, above=0.)
        max_flow = gcmd.get_float('MAX', 80.0, above=start_flow)
        coarse_step = gcmd.get_float('COARSE_STEP', 10.0, above=0.)
        min_step = gcmd.get_float('MIN_STEP', 1.0, above=0.)
        step_duration = gcmd.get_float('DURATION', 5.0, above=0.5)
        repeat = gcmd.get_int('REPEAT', 5, minval=1, maxval=10)
        verify_repeats = gcmd.get_int(
            'VERIFY_REPEATS', 5, minval=1, maxval=10)
        purge = gcmd.get_float('PURGE', 0.0, minval=0.)
        cooldown = gcmd.get_float('COOLDOWN', 15.0, minval=0., maxval=300.)
        max_bisect = gcmd.get_int(
            'MAX_BISECT_STEPS', 6, minval=2, maxval=15)
        no_html = gcmd.get_int('NO_HTML', 0, minval=0, maxval=1)
        skip_check = gcmd.get_int(
            'SKIP_TMC_CHECK', 0, minval=0, maxval=1)
        auto_sgt = gcmd.get_int('AUTO_SGT', 1, minval=0, maxval=1)
        keep_sgt = gcmd.get_int('KEEP_SGT', 0, minval=0, maxval=1)
        save_sgt = gcmd.get_int('SAVE_SGT', 0, minval=0, maxval=1)
        cold_extrusion_hint = gcmd.get_float(
            'COLD_EXTRUSION_HINT', 0.0, minval=0., maxval=200.)
        ref_flow = gcmd.get_float('REF_FLOW', self.reference_flow,
                                  minval=0.)
        # Data-collection additions (recorded only, no triggers):
        #   TEMP=<°C>     heat the hotend and wait (0 = use current target)
        #   SETTLE=<s>    hold temperature this long before each idle
        #                 baseline so the PID output is steady
        #   IDLE=<s>      idle heater-duty baseline length, measured at
        #                 the start and again at the end (0 = off)
        #   BDWIDTH=<name|auto|none>  bdwidth sensor for motion data
        test_temp = gcmd.get_float('TEMP', 0.0, minval=0., maxval=400.)
        settle_s = gcmd.get_float('SETTLE', 20.0, minval=0., maxval=900.)
        idle_s = gcmd.get_float('IDLE', 60.0, minval=0., maxval=600.)
        bdw_param = gcmd.get('BDWIDTH', 'auto').strip()
        #   BDW_MIN_RATIO=<0..1>  fail a step when BDWidth saw less than
        #                 this fraction of the expected filament feed
        #                 (scale-corrected); 0 = record only, no trigger
        self._bdw_min_ratio = gcmd.get_float(
            'BDW_MIN_RATIO', 0.6, minval=0., maxval=0.95)
        #   ENCODER=<name|auto|none>  Klipper [filament_motion_sensor]
        #                 (e.g. BTT SFS) used when no bdwidth is found
        #   ENCODER_GAP=<mm>  fail a step when the extruder advanced this
        #                 far past the sensor's last pulse; 0 = auto
        #                 (2 × the sensor's detection_length, min 10 mm)
        enc_param = gcmd.get('ENCODER', 'auto').strip()
        #   MATERIAL=<PLA|PETG|ABS|…>  recorded in the CSV/HTML (the
        #                 Compare tab checks it) and, on an MPC hotend,
        #                 used to set MPC's filament density / heat
        #                 capacity for the test. Blank = unchanged.
        material = gcmd.get('MATERIAL', '').strip().upper()
        enc_gap = gcmd.get_float('ENCODER_GAP', 0.0, minval=0.)

        # Part-cooling fan speed for the test. Defaults to the value
        # in [tmc_flow_test] config (or 0 if not set). Override per
        # invocation with FAN_SPEED=<percent>.
        fan_speed_pct = gcmd.get_float(
            'FAN_SPEED', self.test_fan_speed, minval=0.0, maxval=100.0)

        if min_step >= coarse_step:
            raise gcmd.error(
                "MIN_STEP (%.1f) must be smaller than COARSE_STEP (%.1f)"
                % (min_step, coarse_step))

        # ─── TMC config check ───
        if not skip_check:
            problems, infos = self._check_tmc_config()
            if problems:
                msg = "\n=== TMC Configuration Issue(s) ===\n"
                for fname, val, desc in problems:
                    msg += "Problem: %s (current value: %s)\n%s\n\n" % (
                        fname, val, desc)
                msg += ("After updating printer.cfg, run:\n"
                        "  FIRMWARE_RESTART\n"
                        "  TMC_FLOW_FIND_MAX\n\n"
                        "To skip this check (advanced):\n"
                        "  TMC_FLOW_FIND_MAX SKIP_TMC_CHECK=1")
                raise gcmd.error(msg)
            gcmd.respond_info("TMC configuration check passed.")
            for info in infos:
                gcmd.respond_info(info)

        # ─── Material → MPC filament constants (MATERIAL=) ───
        mpc_state = None
        if material:
            mpc_state = self._mpc_apply_material(gcmd, material)
            self._mpc_state = mpc_state

        # ─── Hotend heat-up (TEMP=) ───
        if test_temp > 0:
            if test_temp < self.min_hotend_temp:
                raise gcmd.error(
                    "TEMP=%.0f is below min_hotend_temp (%.0f)."
                    % (test_temp, self.min_hotend_temp))
            gcmd.respond_info(
                "Heating hotend to %.1f°C and waiting..." % test_temp)
            # Not M109: macro sets (e.g. Demon) often override it with a
            # loose window (target-7), which then fails the "within 5°C"
            # check below. Wait to within 1°C with stock commands.
            self.gcode.run_script_from_command(
                "SET_HEATER_TEMPERATURE HEATER=extruder TARGET=%.1f\n"
                "TEMPERATURE_WAIT SENSOR=extruder MINIMUM=%.1f "
                "MAXIMUM=%.1f"
                % (test_temp, test_temp - 1.0, test_temp + 5.0))

        # ─── Hotend temp check ───
        extruder = self.printer.lookup_object('extruder')
        heater = extruder.get_heater()
        cur_temp, _ = heater.get_temp(self.reactor.monotonic())
        target_temp = heater.target_temp
        if cur_temp < self.min_hotend_temp:
            raise gcmd.error(
                "Hotend too cold: %.1f°C (min %.1f°C)."
                % (cur_temp, self.min_hotend_temp))
        if target_temp > 0 and cur_temp < target_temp - 5.0:
            raise gcmd.error(
                "Hotend not at target: %.1f°C (target %.1f°C). "
                "Wait for M109 or run M109 S%d before testing."
                % (cur_temp, target_temp, int(target_temp)))

        # ─── Auto-SGT tuning ───
        # Optional pre-test step that probes SG_RESULT at low flow
        # and adjusts SGT to land in the profile's healthy target
        # range. Set AUTO_SGT=0 to skip. By default the original
        # SGT is RESTORED after the test so printer.cfg behavior
        # doesn't silently change; pass KEEP_SGT=1 to leave the
        # tuned value active until the next FIRMWARE_RESTART.
        #
        # IMPORTANT: _measure_step issues `G1 E<value>` commands which
        # are interpreted in the current extruder mode. The main test
        # (further below) sets `M83\nG92 E0` (relative + reset)
        # before the coarse loop, but auto-SGT runs BEFORE that block,
        # so we set the same mode here. Without this, only the first
        # repetition extrudes anything — subsequent reps see "G1 E5"
        # as "go TO position 5" and the extruder is already there.
        # Reference phase plan: two fixed low-flow steps before the
        # sweep, only when START is clearly above them (otherwise the
        # sweep itself starts low and refs would just cost time).
        ref_flows = []
        if ref_flow > 0:
            ref2 = ref_flow + 5.0
            if start_flow >= ref2 + 5.0:
                ref_flows = [ref_flow, ref2]

        original_sgt = None
        final_sgt = None
        if auto_sgt and self._can_autotune_sgt():
            self.gcode.run_script_from_command("M83\nG92 E0")
            # Anchor the SGT calibration at the reference flow when a
            # reference phase runs — probing at a high START would
            # calibrate against near-slip load.
            sgt_probe_flow = ref_flows[0] if ref_flows else start_flow
            original_sgt, final_sgt = self._autotune_sgt(
                gcmd, sgt_probe_flow)

        rotation_distance = self._get_rotation_distance(extruder)
        if rotation_distance is None:
            raise gcmd.error("Could not determine rotation_distance")

        timestamp = time.strftime('%Y-%m-%d_%H-%M-%S')
        sg_label = self._get_sg_label()
        meta = {
            'timestamp': timestamp,
            'material': material or 'not specified',
            'driver': '%s on %s (%s register, sample rate %.0f Hz)' % (
                self.driver_type, self.stepper_name, sg_label,
                1.0 / SAMPLE_INTERVAL),
            'algorithm': 'ADAPTIVE BISECTION (coarse=%.1f, min=%.1f, '
                         'repeats=%d)' % (coarse_step, min_step, repeat),
            'hotend_temp': '%.1f °C (target %.1f °C)' % (
                cur_temp, target_temp),
            'filament_diameter': '%.2f mm' % self.filament_diameter,
            'melt_zone_length': '%.1f mm' % self.melt_zone_length,
            'rotation_distance': '%.4f mm' % rotation_distance,
            'flow_range': '%.1f → %.1f mm³/s (adaptive)' % (
                start_flow, max_flow),
            'reference_flow': ('%.1f + %.1f mm³/s'
                               % (ref_flows[0], ref_flows[1])
                               if ref_flows else 'skipped'),
            'step_duration': '%.1f s' % step_duration,
            'cold_extrusion_hint': cold_extrusion_hint,
            'tmc_settings': self._snapshot_tmc_settings(),
        }

        # ─── Banner ───
        gcmd.respond_info(
            "===== TMC Flow Test v%s =====\n"
            "Range: %.0f → %.0f mm³/s | Coarse: %.0f | Min: %.0f mm³/s\n"
            "Each measurement: %d reps × %.1f s (median over all samples)\n"
            "Cool-down between phases: %.0f s | Hotend: %.1f°C "
            "(target %.1f°C)\n"
            "Sample rate: %.0f Hz | Driver: %s | Register: %s\n"
            "------------------------------------------------\n"
            "Algorithm:\n"
            "  Phase 1 (Coarse): increase flow in big steps until trigger\n"
            "  Phase 2 (Bisection): narrow to ±%.0f mm³/s by halving\n"
            "  Phase 3 (Verification): confirm with %d reps\n"
            "================================================"
            % (MODULE_VERSION,
               start_flow, max_flow, coarse_step, min_step,
               repeat, step_duration, cooldown,
               cur_temp, target_temp,
               1.0 / SAMPLE_INTERVAL, self.driver_type, sg_label,
               min_step, verify_repeats))

        if cold_extrusion_hint > 0:
            gcmd.respond_info(
                "Cold-extrusion hint: %.1f mm³/s — will be marked as "
                "a vertical line in the HTML report. The plugin's slip "
                "detection ignores this value; it's purely informational."
                % cold_extrusion_hint)

        # ─── Resolve extra thermistors and store original fan speed ───
        # Both happen here so they're locked in BEFORE any extrusion
        # starts. Once the test runs, neither value should change.
        self._resolve_extra_thermistors(gcmd)
        original_fan_pct = self._capture_current_fan_speed()
        self._set_part_cooling_fan(fan_speed_pct, gcmd)
        if fan_speed_pct > 0:
            gcmd.respond_info(
                "Part-cooling fan locked at %.0f %% for the test "
                "(was %.0f %% — restored at end). Fan speed materially "
                "affects max flow; consistency matters more than "
                "flexibility during a max-flow test."
                % (fan_speed_pct, original_fan_pct))

        # ─── BDWidth motion sensor (data only) ───
        bdw = self._resolve_bdwidth(gcmd, bdw_param)
        if bdw is not None:
            self._bdwidth_prepare(gcmd, bdw)
            meta['bdwidth'] = ('%s (test mode=%s, was %s, %.1f counts/mm, '
                               'sample %.1f s)' % (
                                   bdw['name'], bdw['mode'],
                                   bdw['orig_mode'], bdw['counts_per_mm'],
                                   bdw['sample_time']))
        else:
            meta['bdwidth'] = 'not used'
        # Generic encoder (Klipper filament_motion_sensor) only when no
        # bdwidth: bdwidth measures the feed ratio, an encoder can only
        # tell "the filament stopped".
        enc = None
        if bdw is None:
            enc = self._resolve_encoder(gcmd, enc_param, enc_gap)
            if enc is not None:
                self._encoder_prepare(gcmd, enc)
                meta['motion_sensor'] = (
                    '%s (filament_motion_sensor, detection_length %.1f mm,'
                    ' stop when > %.1f mm without a pulse)'
                    % (enc['name'], enc['detection_length'],
                       enc['gap_limit']))
        if bdw is None and enc is None:
            meta['motion_sensor'] = 'none'
            gcmd.respond_info(
                "No filament motion sensor (bdwidth or "
                "filament_motion_sensor) found — the feed check is "
                "skipped. Heater collapse still catches stalls and "
                "strips.")

        # ─── Chamber / filter fans (config: chamber_fans) ───
        # Before the idle baseline, so both idle readings see them.
        self._fan_state = self._chamber_fans_start(gcmd)
        if self._fan_state:
            meta['chamber_fans'] = ', '.join(
                f['name'] for f in self._fan_state)

        # ─── Idle heater baseline (start) ───
        # Fan is already locked, so the idle losses match the test.
        if idle_s > 0:
            ib = self._measure_idle_baseline(gcmd, settle_s, idle_s,
                                             'start')
            if ib:
                meta['idle_pwm_start'] = round(ib['pwm'], 4)
                meta['idle_temp_start'] = round(ib['temp'], 2)
        t0 = self.reactor.monotonic()
        if mpc_state is not None and mpc_state.get('applied'):
            meta['mpc_filament'] = (
                'density %.2f g/cm³, heat capacity %.2f J/g/K (%s; '
                'restored to %.2f / %.2f after the test)'
                % (mpc_state['density'], mpc_state['heat_capacity'],
                   mpc_state['source'], mpc_state['orig'][0],
                   mpc_state['orig'][1]))
        meta['idle_settle_s'] = settle_s
        meta['idle_measure_s'] = idle_s

        if purge > 0:
            self.gcode.run_script_from_command(
                "M83\nG1 E%.2f F300\nG92 E0" % purge)
        self.gcode.run_script_from_command("M83\nG92 E0")

        results = []
        trigger_events_ref = []   # same list object as trigger_events

        def bdw_scale():
            """Sensor scale for this filament: median feed ratio of the
            steps so far (clamped 0.8–1.2; 1.0 until two exist)."""
            prior = sorted(r['bdw']['ratio'] for r in results
                           if (r.get('bdw') or {}).get('ratio')
                           is not None)
            if len(prior) < 2:
                return 1.0
            n = len(prior)
            m = (prior[n // 2] if n % 2
                 else (prior[n // 2 - 1] + prior[n // 2]) / 2.0)
            return min(max(m, 0.8), 1.2)

        def measure_annotated(flow, reps):
            """_measure_step plus the data-only extras: time since the
            start baseline (for idle-drift correction), filament motion
            (BDWidth / encoder), a per-rep stop when the filament has
            clearly stopped, and one combined console line."""
            bdw_before = (self._bdwidth_total(bdw)
                          if bdw is not None else None)
            min_ratio = getattr(self, '_bdw_min_ratio', 0.0) or 0.0

            def rep_check(done, rep_len, mm_s):
                # No wait here: the sensor may lag by up to one sample
                # period, so that much is taken off the expected length
                # (errs towards NOT stopping).
                if bdw is not None and bdw_before is not None \
                        and min_ratio > 0:
                    try:
                        now = float(bdw['obj'].actual_total_move)
                    except Exception:
                        return None
                    cmd = rep_len * done
                    exp = max(cmd - bdw['sample_time'] * mm_s, cmd * 0.5)
                    meas = (now - bdw_before) / bdw['counts_per_mm']
                    rel = meas / exp / bdw_scale()
                    if rel < min_ratio:
                        return ("BDWidth saw %.0f %% of the expected feed "
                                "after rep %d (%.0f of %.0f mm)"
                                % (rel * 100, done, meas, cmd))
                if enc is not None:
                    e = self._encoder_step(enc)
                    gap = e.get('gap_mm')
                    if gap is not None and gap > enc['gap_limit']:
                        return ("encoder saw no movement for %.0f mm "
                                "(limit %.1f mm) after rep %d"
                                % (gap, enc['gap_limit'], done))
                return None

            check_fn = (rep_check if (bdw is not None or enc is not None)
                        else None)
            r = self._measure_step(
                gcmd, flow, step_duration, reps, rep_check=check_fn,
                quiet=True)
            r['t_s'] = self.reactor.monotonic() - t0
            es = r.get('early_stop')
            reps_done = es['reps'] if es else reps
            parts = []
            sm = r.get('summary') or {}
            sg_txt = "SG %s (CV %s" % (sm.get('sg', 'n/a'),
                                       sm.get('cv', 'n/a'))
            sg_txt += ", dips %d" % (sm.get('dips') or 0)
            if sm.get('peaks'):
                sg_txt += ", peaks %d" % sm['peaks']
            if sm.get('warmup'):
                sg_txt += ", run 1 dropped"
            parts.append(sg_txt + ")")
            pwm = (r.get('thermal') or {}).get('pwm_avg')
            if pwm is not None:
                lf = self._live_melt_expect(results + [r])
                if lf is None:
                    h = "heater %.1f%% (fitting line…)" % (pwm * 100)
                elif lf['b'] <= 0:
                    h = ("heater %.1f%% (low-flow line is flat — no "
                         "expected value)" % (pwm * 100))
                else:
                    exp = lf['a'] + lf['b'] * flow
                    h = ("heater %.1f%% vs %.1f%% exp (%.0f%% of line)"
                         % (pwm * 100, exp * 100,
                            (pwm - lf['a']) / (lf['b'] * flow) * 100))
                if pwm >= self.profile.PWM_SATURATED:
                    h += " SATURATED"
                parts.append(h)
            if bdw is not None:
                bdw_after = self._bdwidth_total(bdw)
                cmd_mm = (flow / self.filament_area) * step_duration \
                    * reps_done
                r['bdw'] = self._bdwidth_step(
                    bdw, bdw_before, bdw_after, cmd_mm)
                b = r['bdw']
                if b.get('meas_mm') is not None:
                    parts.append("BDWidth %.0f / %.0f mm (%.0f%%)"
                                 % (b['meas_mm'], b['cmd_mm'],
                                    b['ratio'] * 100))
            if enc is not None:
                r['enc'] = self._encoder_step(enc)
                gap = r['enc'].get('gap_mm')
                if gap is not None:
                    if gap > enc['gap_limit']:
                        parts.append("encoder STOPPED (no pulse for "
                                     "%.0f mm)" % gap)
                    else:
                        parts.append("encoder OK (last pulse %.1f mm "
                                     "ago, limit %.0f mm)"
                                     % (gap, enc['gap_limit']))
            if es:
                parts.append("stopped after rep %d/%d"
                             % (es['reps'], es['of']))
            gcmd.respond_info("  %.1f mm³/s | %s" % (flow,
                                                     " | ".join(parts)))
            return r

        def measure_and_save(flow, phase):
            r = measure_annotated(flow, repeat)
            r['phase'] = phase
            results.append(r)
            self._save_report(
                results, meta, timestamp, None, no_html,
                gcmd=None, announce=False)
            return r

        def check(results):
            """SG-based slip detection. A HARD failure (filament stopped)
            ends the run here instead of returning to the bisection."""
            reason = self._check_triggers_sg(results)
            if reason and getattr(self, '_hard_failure', None):
                results[-1]['hard_fail'] = True
                hard_stop(reason)
            es = (results[-1].get('early_stop') or {}) if results else {}
            if es.get('reason'):
                # A rep already showed the filament stopped; the
                # whole-step checks normally agree, but stop regardless.
                results[-1]['hard_fail'] = True
                hard_stop(es['reason'])
            return reason

        def hard_stop(reason):
            cur = results[-1]['flow'] if results else 0.0
            # Steps that passed: lower flow AND not themselves a hard
            # failure (a jam can sit through several reference steps).
            good = [r['flow'] for r in results[:-1]
                    if r['flow'] < cur - 0.001
                    and not r.get('hard_fail')]
            last_good = max(good) if good else None
            trigger_events_ref.append({
                'phase': results[-1].get('phase', 'coarse') if results
                else 'coarse',
                'flow': cur, 'kind': 'trigger', 'reason': reason,
                'cv': None, 'iqr': None, 'sg_median': None})
            gcmd.respond_info(
                "\n  >>> FILAMENT FAILURE at %.1f mm³/s — %s\n"
                "      Stopping the test: the filament is stripped, "
                "stalled or clogged, so further steps would only grind "
                "it. Reload or trim the filament before the next run."
                % (cur, reason))
            if last_good is None:
                final = None
            else:
                final = {
                    'max_safe': last_good,
                    'verify_cv': 0.0,
                    'quality': 'stopped on filament failure',
                    'stop_reason': reason,
                    'trigger_events': list(trigger_events_ref),
                    'baseline_stats': self._coarse_baseline_stats(results),
                }
            # Save and restore everything BEFORE the end-of-run idle
            # baseline, so cancelling during that wait loses nothing.
            self._save_report(results, meta, timestamp, reason,
                              no_html, gcmd=None, announce=False,
                              final_result=final)
            self.gcode.run_script_from_command("G92 E0")
            self._restore_sgt_if_needed(gcmd, original_sgt, final_sgt,
                                         keep_sgt, save_sgt=save_sgt)
            try:
                self._set_part_cooling_fan(original_fan_pct, gcmd,
                                           quiet=True)
            except Exception:
                pass
            if idle_s > 0:
                gcmd.respond_info(
                    "  Report saved. Measuring the end-of-run idle "
                    "baseline (~%.0f s, no extrusion) for drift "
                    "correction and the melt limit — safe to cancel, "
                    "you only lose those two extras."
                    % (settle_s + idle_s))
            finish_idle()
            if last_good is None:
                gcmd.respond_info(
                    "The filament was not moving from the first step on "
                    "(%.1f mm³/s) — this is a jam or clog, not a flow "
                    "limit. Clear it (cold pull / trim the chewed "
                    "filament) before re-running." % cur)
            else:
                gcmd.respond_info(
                    "\n========== FINAL RESULT ==========\n"
                    "Last step before the filament failure: %.1f mm³/s\n"
                    "(mechanical limit — the extruder could not push "
                    "%.1f)\n"
                    "----------------------------------\n"
                    "Slicer recommendation:\n"
                    "  Conservative (80%%): %.1f mm³/s\n"
                    "  Aggressive (90%%):   %.1f mm³/s\n"
                    "  See the melt-limit line above for print quality.\n"
                    "=================================="
                    % (last_good, cur, last_good * 0.8, last_good * 0.9))
            self._save_report(results, meta, timestamp, reason,
                              no_html, gcmd=gcmd, final_result=final)
            raise _HardStop(reason)

        def finish_idle():
            """Idle heater baseline at the end of the run, so the CSV
            can correct for drift (chamber warm-up) between start and
            end. Recorded only; never affects the result. Also puts the
            BDWidth sensor back in its original mode."""
            if bdw is not None:
                self._bdwidth_restore(gcmd, bdw)
            if enc is not None:
                self._encoder_restore(gcmd, enc)
            if mpc_state is not None:
                self._mpc_restore(gcmd, mpc_state)
            report_slip_onset()
            try:
                if idle_s <= 0:
                    return
                try:
                    ib_end = self._measure_idle_baseline(
                        gcmd, settle_s, idle_s, 'end')
                except Exception as e:
                    gcmd.respond_info("  Idle baseline (end) skipped: %s"
                                      % (e,))
                    return
                if ib_end:
                    meta['idle_pwm_end'] = round(ib_end['pwm'], 4)
                    meta['idle_temp_end'] = round(ib_end['temp'], 2)
                    meta['idle_t_end_s'] = round(
                        ib_end['t_mid'] - t0, 1)
                report_melt_knee()
            finally:
                # Chamber fans ran through the end idle as well.
                self._chamber_fans_restore(gcmd, self._fan_state)

        def report_slip_onset():
            s = self._slip_onset(results)
            if s is None:
                return
            if s['onset'] is None:
                meta['slip_onset'] = ('no slip seen (BDWidth feed stayed '
                                      'within %.0f %% of the low-flow '
                                      'reading up to %.1f mm³/s)'
                                      % (s['frac'] * 100, s['max_tested']))
            else:
                meta['slip_onset'] = (
                    '%.1f mm³/s (feed %.0f %% of the low-flow reading; '
                    'last clean step %s)'
                    % (s['onset'], s['onset_rel'] * 100,
                       ('%.1f mm³/s' % s['last_clean'])
                       if s['last_clean'] is not None else 'none'))
            gcmd.respond_info("Slip start (BDWidth, report only): %s"
                              % meta['slip_onset'])

        def report_melt_knee():
            knee = self._melt_knee(results, meta)
            if knee is None:
                return
            meta['melt_fit_intercept'] = round(knee['intercept'], 4)
            meta['melt_fit_slope'] = round(knee['slope'], 6)
            if knee.get('warn'):
                meta['melt_fit_warning'] = knee['warn']
                gcmd.respond_info("  Melt-limit note: %s" % knee['warn'])
            if knee['knee'] is None:
                meta['melt_limit'] = ('heat per mm³ stayed above %.0f %% '
                                      'of the low-flow line up to '
                                      '%.1f mm³/s'
                                      % (knee['frac'] * 100,
                                         knee['max_tested']))
            else:
                meta['melt_limit'] = ('%s (heat per mm³ fell to '
                                      '%.0f %% of the low-flow line at '
                                      '%.1f mm³/s)'
                                      % ('%.1f mm³/s' % knee['last_ok']
                                         if knee['last_ok'] else 'n/a',
                                         knee['knee_ratio'] * 100,
                                         knee['knee']))
            gcmd.respond_info("\nMelt limit (print quality, report "
                              "only): %s" % meta['melt_limit'])

        # ─── PHASE 0: Reference (baseline) ───
        # Two fixed low-flow steps far below any realistic limit.
        # They feed the trigger baselines (CV, dips, peaks, IQR, gap)
        # with guaranteed slip-free data, so the sweep's very first
        # step is already fully trigger-armed even with a high START.
        if ref_flows:
            gcmd.respond_info(
                "\n>>> Phase 0: Reference baseline <<<\n"
                "  Measuring %s mm³/s — slip-free reference for the "
                "trigger baselines (START=%.0f is far above; disable "
                "with REF_FLOW=0)."
                % (" + ".join("%.0f" % f for f in ref_flows),
                   start_flow))
            for rf in ref_flows:
                measure_and_save(rf, 'ref')
                # Only the hard checks (filament stopped) act here; the
                # SG checks need the sweep history and are ignored.
                check(results)
                self.gcode.run_script_from_command("G4 P500")

        # ─── PHASE 1: Coarse ───
        gcmd.respond_info(
            "\n>>> Phase 1: Coarse Upward Sweep <<<\n"
            "  Stepping up by %.0f mm³/s. Each step: %d reps × %.1f s."
            % (coarse_step, repeat, step_duration))

        flow = start_flow
        low = None
        high = None
        first_trigger_reason = None
        # Track every trigger / borderline / verify-fail event for the
        # HTML "decision trail" panel. Each entry: dict with keys
        # 'phase', 'flow', 'kind' (one of: 'trigger', 'borderline-retest',
        # 'borderline-confirmed', 'verify-fail', 'verify-borderline'),
        # 'reason', 'cv', 'iqr', 'sg_median'.
        trigger_events = trigger_events_ref

        def _record_event(phase, flow, kind, reason, last_result=None):
            ev = {'phase': phase, 'flow': flow,
                  'kind': kind, 'reason': reason,
                  'cv': None, 'iqr': None, 'sg_median': None}
            if last_result is not None:
                rc = last_result.get('run_consistency') or {}
                ev['cv'] = rc.get('sg_cv')
                sg = last_result.get('sg') or {}
                if 'p25' in sg and 'p75' in sg:
                    ev['iqr'] = sg['p75'] - sg['p25']
                if 'median' in sg:
                    ev['sg_median'] = sg['median']
            trigger_events.append(ev)

        while flow <= max_flow + 0.001:
            r = measure_and_save(flow, 'coarse')

            reason = check(results)
            if reason:
                high = flow
                low = flow - coarse_step
                if low < start_flow and ref_flows:
                    # Trigger on the very first sweep step: fall back
                    # to the highest reference step as the safe lower
                    # bound instead of aborting — the refs are proven
                    # slip-free measurements.
                    low = ref_flows[-1]
                first_trigger_reason = reason
                _record_event('coarse', flow, 'trigger', reason,
                              results[-1] if results else None)
                gcmd.respond_info(
                    "  >>> TRIGGER at %.1f mm³/s — %s\n"
                    "      → Safe range narrowed to [%.1f, %.1f] mm³/s"
                    % (flow, reason, low, high))
                break

            low = flow
            self.gcode.run_script_from_command("G4 P500")
            flow += coarse_step

        if high is None:
            gcmd.respond_info(
                "Reached MAX %.1f mm³/s without trigger. Try MAX higher."
                % max_flow)
            finish_idle()
            self._save_report(results, meta, timestamp, None,
                              no_html, gcmd=gcmd)
            self._restore_sgt_if_needed(gcmd, original_sgt, final_sgt,
                                         keep_sgt, save_sgt=save_sgt)
            return

        sweep_floor = ref_flows[-1] if ref_flows else start_flow
        if low < sweep_floor:
            gcmd.respond_info(
                "Trigger fired on first step (%.1f) — lower START." % high)
            finish_idle()
            self._save_report(results, meta, timestamp, first_trigger_reason,
                              no_html, gcmd=gcmd)
            self._restore_sgt_if_needed(gcmd, original_sgt, final_sgt,
                                         keep_sgt, save_sgt=save_sgt)
            return

        # ─── PHASES 2 + 3: Bisection + Verify (with retry) ───
        # Phase 2 narrows the bracket; Phase 3 confirms the result. If
        # Phase 3 fails (verify itself trips a trigger or stays
        # borderline on re-test), we drop back into Phase 2 with a
        # tightened high bound and try again. Up to MAX_VERIFY_FAILURES
        # rounds before giving up — without this cap a hardware /
        # filament issue could loop forever.
        MAX_VERIFY_FAILURES = 3

        # Track which flows we've already re-tested in bisection (set
        # is reset on each bisection re-entry to allow fresh re-tests
        # at the new bracket).
        last_trigger_reason = first_trigger_reason
        verify_failures = 0
        verify_result = None
        verify_cv = 0.0

        while True:
            # ─── PHASE 2: Bisection ───
            if cooldown > 0:
                gcmd.respond_info(
                    "  ... Cool-down: %.0f s ..." % cooldown)
                self.gcode.run_script_from_command(
                    "G4 P%d" % int(cooldown * 1000))

            phase2_label = (("\n>>> Phase 2: Bisection (retry %d) <<<"
                             % verify_failures)
                            if verify_failures > 0
                            else "\n>>> Phase 2: Bisection <<<")
            gcmd.respond_info(
                "%s\n"
                "  Narrowing [%.0f, %.0f] by halving until interval ≤ %.0f. "
                "Up to %d steps.\n"
                "  Borderline measurements (CV 4-7%% or IQR 15-24) are "
                "re-tested once for confirmation."
                % (phase2_label, low, high, min_step, max_bisect))

            bisect_iter = 0
            retested_flows = set()
            while ((high - low) > min_step + 0.001
                   and bisect_iter < max_bisect):
                bisect_iter += 1
                raw_mid = (low + high) / 2.0
                mid = round(raw_mid / min_step) * min_step
                if mid <= low + 0.001 or mid >= high - 0.001:
                    break

                r = measure_and_save(mid, 'bisect')
                reason = check(results)

                # Borderline check: if no clear trigger fired but the
                # data sits in the gray zone, re-measure once before
                # classifying.
                if (reason is None
                        and mid not in retested_flows):
                    borderline_why = self._is_borderline(results)
                    if borderline_why:
                        retested_flows.add(mid)
                        _record_event('bisect', mid, 'borderline-retest',
                                      borderline_why,
                                      results[-1] if results else None)
                        gcmd.respond_info(
                            "  >>> %.1f mm³/s BORDERLINE — %s\n"
                            "      Re-measuring once for confirmation..."
                            % (mid, borderline_why))
                        if cooldown > 0:
                            self.gcode.run_script_from_command(
                                "G4 P%d" % int(cooldown * 1000))
                        # Drop the borderline measurement and replace
                        # with the fresh one.
                        results.pop()
                        r2 = measure_and_save(mid, 'bisect')
                        reason = check(results)
                        if reason is None:
                            re_borderline = self._is_borderline(results)
                            if re_borderline:
                                reason = ("borderline confirmed on re-"
                                          "test: " + re_borderline)
                                _record_event(
                                    'bisect', mid, 'borderline-confirmed',
                                    re_borderline,
                                    results[-1] if results else None)

                if reason:
                    high = mid
                    last_trigger_reason = reason
                    # Avoid double-recording if borderline-confirmed
                    # already recorded above.
                    if not (trigger_events
                            and trigger_events[-1].get('flow') == mid
                            and trigger_events[-1].get('kind')
                            == 'borderline-confirmed'):
                        _record_event('bisect', mid, 'trigger', reason,
                                      results[-1] if results else None)
                    gcmd.respond_info(
                        "  >>> TRIGGER at %.1f — %s\n"
                        "      → [%.1f, %.1f] (%d/%d)"
                        % (mid, reason, low, high,
                           bisect_iter, max_bisect))
                else:
                    low = mid
                    gcmd.respond_info(
                        "  >>> %.1f mm³/s SAFE → [%.1f, %.1f] (%d/%d)"
                        % (mid, low, high, bisect_iter, max_bisect))
                self.gcode.run_script_from_command("G4 P500")

            # ─── PHASE 3: Verify ───
            if cooldown > 0:
                gcmd.respond_info(
                    "  ... Cool-down before verification: %.0f s ..."
                    % cooldown)
                self.gcode.run_script_from_command(
                    "G4 P%d" % int(cooldown * 1000))

            gcmd.respond_info(
                "\n>>> Phase 3: Verification at %.1f mm³/s <<<\n"
                "  Confirming with %d repetitions."
                % (low, verify_repeats))
            verify_result = measure_annotated(low, verify_repeats)
            verify_result['phase'] = 'verify'
            results.append(verify_result)

            # Check verify result against triggers and borderline.
            verify_trigger = check(results)
            if verify_trigger is None:
                verify_borderline = self._is_borderline(results)
                if verify_borderline:
                    # Re-measure verify once at the same flow before
                    # accepting; if it stays borderline, treat as fail.
                    _record_event('verify', low, 'verify-borderline',
                                  verify_borderline,
                                  results[-1] if results else None)
                    gcmd.respond_info(
                        "  >>> Verify at %.1f BORDERLINE — %s\n"
                        "      Re-measuring once for confirmation..."
                        % (low, verify_borderline))
                    if cooldown > 0:
                        self.gcode.run_script_from_command(
                            "G4 P%d" % int(cooldown * 1000))
                    results.pop()
                    verify_result = measure_annotated(low, verify_repeats)
                    verify_result['phase'] = 'verify'
                    results.append(verify_result)
                    verify_trigger = check(results)
                    if verify_trigger is None:
                        re_borderline = self._is_borderline(results)
                        if re_borderline:
                            verify_trigger = ("borderline confirmed on "
                                              "re-test: "
                                              + re_borderline)

            if verify_trigger is None:
                # Verify clean — accept the result.
                break

            # Verify FAILED. Drop back into bisection with a tighter
            # high bound (the just-failed flow becomes the new high).
            verify_failures += 1
            _record_event('verify', low, 'verify-fail', verify_trigger,
                          results[-1] if results else None)
            verify_cv_now = (verify_result.get('run_consistency', {})
                             .get('sg_cv', 0.0))
            gcmd.respond_info(
                "\n  >>> ⚠ VERIFY FAILED at %.1f mm³/s "
                "(CV %.1f%%, attempt %d/%d)\n"
                "      Reason: %s"
                % (low, verify_cv_now, verify_failures,
                   MAX_VERIFY_FAILURES, verify_trigger))

            if verify_failures >= MAX_VERIFY_FAILURES:
                # Give up gracefully — report the next-lower safe step
                # we have evidence for. Search results history for the
                # last bisect/coarse step at flow < low that did NOT
                # trigger.
                fallback = None
                for r in reversed(results[:-1]):  # skip the failed verify
                    if r.get('phase') in ('coarse', 'bisect', 'ref'):
                        if r['flow'] < low - 0.001:
                            fallback = r['flow']
                            break
                if fallback is None:
                    fallback = max(sweep_floor, low - min_step)
                gcmd.respond_info(
                    "  Maximum verify retries reached — falling back "
                    "to last known-good flow %.1f mm³/s." % fallback)
                low = fallback
                last_trigger_reason = verify_trigger
                # Do one final verify at the fallback to populate the
                # final verify_result with this lower value.
                if cooldown > 0:
                    self.gcode.run_script_from_command(
                        "G4 P%d" % int(cooldown * 1000))
                gcmd.respond_info(
                    "\n>>> Phase 3: Final verification at %.1f mm³/s "
                    "<<<\n  Confirming with %d repetitions."
                    % (low, verify_repeats))
                verify_result = measure_annotated(low, verify_repeats)
                verify_result['phase'] = 'verify'
                results.append(verify_result)
                break

            # Tighten bracket: failed flow becomes new high.
            high = low
            # Drop low to the previous safe step we have data for, or
            # one min_step below if no earlier safe step exists.
            new_low = None
            for r in reversed(results[:-1]):  # skip failed verify
                if r.get('phase') in ('coarse', 'bisect', 'ref'):
                    if r['flow'] < high - 0.001:
                        # Was this flow a safe one?  Check it didn't
                        # carry a trigger.
                        # (Heuristic: if the flow appears as a
                        # 'safe' entry — i.e. CV < 4% AND not
                        # borderline by current standards — accept
                        # it as the new low.)
                        rc = r.get('run_consistency') or {}
                        cv = rc.get('sg_cv', 99.0)
                        sgp = r.get('sg') or {}
                        iqr = (sgp.get('p75', 0) - sgp.get('p25', 0)
                               if 'p25' in sgp and 'p75' in sgp else 99)
                        if cv < 4.0 and iqr < 15:
                            new_low = r['flow']
                            break
            if new_low is None:
                new_low = max(sweep_floor, high - 2 * min_step)
            low = new_low
            last_trigger_reason = verify_trigger
            gcmd.respond_info(
                "      → Re-entering bisection with bracket "
                "[%.1f, %.1f]" % (low, high))

            # Sanity: bracket too tight to bisect → just accept low.
            if (high - low) <= min_step + 0.001:
                gcmd.respond_info(
                    "      Bracket already tight, "
                    "accepting %.1f as final." % low)
                # Do one more verify at this low to capture the final
                # measurement metadata.
                if cooldown > 0:
                    self.gcode.run_script_from_command(
                        "G4 P%d" % int(cooldown * 1000))
                gcmd.respond_info(
                    "\n>>> Phase 3: Final verification at %.1f mm³/s "
                    "<<<\n  Confirming with %d repetitions."
                    % (low, verify_repeats))
                verify_result = measure_annotated(low, verify_repeats)
                verify_result['phase'] = 'verify'
                results.append(verify_result)
                break

        verify_cv = (verify_result.get('run_consistency', {}).get('sg_cv', 0)
                     if verify_result.get('run_consistency') else 0)

        max_safe = low
        if verify_cv < 5:
            quality = "excellent (very stable)"
        elif verify_cv < 10:
            quality = "good (stable)"
        elif verify_cv < 20:
            quality = "acceptable (some variation)"
        else:
            quality = "poor (high variation — re-run advised)"

        gcmd.respond_info(
            "\n========== FINAL RESULT ==========\n"
            "Maximum safe volumetric flow: %.1f mm³/s\n"
            "Verification quality: %s (CV = %.1f%%)\n"
            "----------------------------------\n"
            "Slicer recommendation:\n"
            "  Conservative (80%%): %.1f mm³/s   ← recommended\n"
            "  Aggressive (90%%):   %.1f mm³/s   ← only with margin\n"
            "----------------------------------\n"
            "Detailed data: see CSV/HTML report\n"
            "=================================="
            % (max_safe, quality, verify_cv,
               max_safe * 0.8, max_safe * 0.9))

        finish_idle()
        self._save_report(results, meta, timestamp, last_trigger_reason,
                          no_html, gcmd=gcmd,
                          final_result={
                              'max_safe': max_safe,
                              'verify_cv': verify_cv,
                              'quality': quality,
                              'stop_reason': last_trigger_reason,
                              'trigger_events': trigger_events,
                              'baseline_stats':
                                  self._coarse_baseline_stats(results),
                          })
        self.gcode.run_script_from_command("G92 E0")

        # Restore SGT (unless KEEP_SGT=1)
        self._restore_sgt_if_needed(gcmd, original_sgt, final_sgt,
                                     keep_sgt, save_sgt=save_sgt)

        # Restore part-cooling fan to whatever it was before the test
        try:
            self._set_part_cooling_fan(original_fan_pct, gcmd,
                                       quiet=True)
        except Exception:
            pass


def load_config(config):
    return TMCFlowTest(config)
