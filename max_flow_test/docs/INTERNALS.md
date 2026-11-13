# How slip detection works

The plugin samples StallGuard at 50 Hz during each measurement step (5 repeats × 5 s by default). The first 10 % and last 5 % of every run are discarded (spin-up/spin-down transients — SG2 readings are invalid while the motor accelerates or decelerates, which would poison `sg_min`/`sg_max` and the dip/peak counters). From the remaining samples it tracks per-step **median**, **IQR (P25–P75)**, **run-to-run CV**, **collapse-dip count**, **unload-peak count**, and **intra-run trend** (slope of SG over time within a single run).

Slip detection uses **multiple independent triggers** that look for different signatures:

- **Thermal runaway** (all drivers) — the only trigger that does not read StallGuard, and the only one that catches *continuous* grinding: when the hobbed gear chews into filament that no longer advances, the motor load stays high and steady, so medians, dips, peaks and CV all look normal. But melting a given flow absorbs a matching amount of power — if the material stops arriving, that heat sink disappears and the hotend overshoots its target even though the heater was maxed out one step earlier. Fires when the step runs ≥5 °C above target **and** ≥4 °C above the run's own earlier offset (so a hotend that simply runs warm doesn't trip it). Validated: eight clean historical runs peak at 0.7–4.9 °C, a real grinding run reached +17.9 °C.
- **Heater collapse** (all drivers) — the second trigger that does not read StallGuard. Melting more plastic takes more heat, so heater duty must rise with flow. When the commanded flow goes up but the average heater duty drops by more than 10 % (and at least 3 percentage points) versus the step at the highest lower flow, the filament has stopped arriving. This catches what the thermal trigger misses when a strong heater PID holds the overshoot under 5 °C, including an extruder that loses sync at move start and sits still in a spinning field (SG keeps tracking the commanded step rate, so the SG triggers stay quiet). Needs only one lower-flow step, so it is armed from the second step. Observed: a TMC5160 start-of-move stall dropped duty 38 % → 27 % while SG stayed on a smooth curve.
- **BDWidth feed** (any driver, needs a bdwidth sensor) — reads the bdwidth motion counter before and after every step and compares the measured filament feed with the commanded length, scale-corrected by the median of the earlier steps. Fires below `BDW_MIN_RATIO` (default 60 %). Needs no history, so it also covers a stall on the very first step. Deliberately coarse: per-step noise is about ±5 % with occasional dropouts (one −25 % reading seen), so it is a "filament stopped" check, not a slip-onset check.
- **Hard stop** — heater collapse and BDWidth feed mean the filament physically stopped (stall, strip or clog). Either one ends the run on the spot instead of bisecting: every further step would only grind the damaged filament. The result is the last step below the failure (the mechanical limit). Reload or trim the filament before the next run.
- **Melt limit (report only)** — fits a straight line, heater duty = a + b·flow, through the healthy low-flow steps (steps 2–5 in run order, drift-corrected with the start/end idle baselines; failed steps excluded). While the plastic fully melts every step sits on that line; once the melt zone cannot keep up, duty above the intercept falls below the line. The first flow below 85 % of the line is reported as the print-quality limit (console, CSV/HTML header `melt_limit`, per-step `melt_fit_pct` column). Fitting the line instead of dividing by the measured idle absorbs a constant per-rep offset (0.7–3 duty points in real runs) that otherwise makes heat-per-mm³ slide from the first step on. If the intercept is more than 0.05 away from the measured idle, or duty does not rise with flow, the report says so. Never stops the test. On a filament failure the report is saved and SGT/fan restored before the end-of-run idle baseline, so cancelling that wait loses only the drift correction and the melt limit.
- **Collapse dips** (SG2 only) — samples below 25 % of the run median are brief stall events (stick-slip clicks). Median/IQR/CV absorb a handful of dips among hundreds of samples completely, so this is the primary trigger for audible click-stalls. Fires when the step's dip rate exceeds 1 % of samples (≥3 dips) and 3× the prior-step baseline. Logged as `sg_dips` / `sg_dip_rate_pct` in the CSV.
- **Unload peaks** (SG2 only) — the extruder-grinding mirror image of the dip: when the drive gear grinds through the filament the motor is briefly *unloaded* and SG spikes above 145 % of the run median while the median stays smooth. Fires when the peak rate exceeds 0.6 % of samples (≥8 peaks) and 3× the prior-step baseline. Logged as `sg_peaks` / `sg_peak_rate_pct`.
- **SG signal patterns** — snap-back, over-jump, single-step plateau, 2-step cumulative plateau (with saturation-skip and median-baseline)
- **Run-to-run variance** — CV spike, CV jump, rising trend, vs coarse-baseline
- **Sample distribution** — IQR widening (single-step), IQR cumulative growth (vs early-test baseline), IQR vs coarse-baseline, IQR absolute floor
- **Per-run analysis** — single-run outlier detection (warmup-aware), SG max spike for decoupling

The per-step console line also flags **`heater at 100%`** when the heater duty averaged ≥95 % for a whole step: from that flow on, the *hotend* — not the motor — is the limit, and higher numbers are not sustainable in a real print even if no slip is detected.

Trigger maturity: the thermal, dip and CV triggers are active from the 3rd step of the sweep; the distribution triggers (IQR, max-spike, outlier) need 5 steps of history because early steps often sit in the SG saturation region with spiky max/IQR.

Each trigger fires under tighter conditions in **bisection / verify** than in coarse, so the coarse phase stays noise-resistant while the final result is accurate to ±1 mm³/s.

The HTML report's **decision-trail panel** lists every trigger event with the metrics that caused it, so you can see exactly why the plugin chose the value it did.

---

## Plateau and IQR-growth triggers

Two of the most important trigger types deserve more explanation:

### Single-step plateau

Fires when the most recent step's SG-delta is essentially flat compared to the recent baseline (e.g. typical −25 deltas then suddenly +0.5). Catches abrupt plateaus that a 2-step cumulative trigger would smooth over.

The threshold is half the per-driver `PLATEAU_RATIO` — the actual delta must be a much clearer single-step plateau to fire than a borderline 2-step accumulation.

### IQR cumulative growth

Fires when within-step spread (P75 − P25) has gradually doubled relative to the early-test baseline AND now exceeds an absolute floor of 12 raw units. Catches gradual stick-slip onset that single-step ratio triggers don't see.

Requirements:
- ≥5 coarse steps of history (so the baseline is meaningful)
- Last IQR ≥ 12 absolute units
- Last IQR ≥ 2.0 × early baseline (median of 2nd+3rd smallest of first 4 IQRs)

The early baseline uses the median of mid-ranked IQRs (not the smallest, which can be a noise floor). This makes the trigger robust against initial noise.

---

## Warmup-skip

The first repetition of every measurement step shows different behaviour than the others (motor transitions from cold-stop, filament path settles). The plugin detects this drift and excludes run 1 from the median/IQR/CV stats when it deviates significantly from the rest.

The threshold is per-driver:
- TMC5160: 10 % cutoff
- TMC2240: 4 % cutoff (more aggressive — TMC2240 has systematic 3–6 % first-run drift)
- TMC2209: 15 % cutoff (SG4 has higher first-run drift)

---

## Thermal monitoring

During every measurement the plugin captures heater PWM (avg/max), hotend temperature (target/actual/min/drop), and TMC driver thermal flags (otpw / ot). These appear in the CSV and HTML report. **They do NOT trigger anything** — they're context for interpreting the result.

A heuristic 0–100 **thermal stress score** is computed per step combining five signals:

1. **Heater PWM level** (0–30 points) — how hard the heater is working
2. **Temperature drop from target** (0–25 points) — how much actual temp falls below target
3. **PWM rising trend** (0–15 points) — heater taking on more work over the sweep
4. **PWM peak saturation hits** (0–10 points) — saturation events even briefly
5. **Intra-run SG drift** (0–20 points) — load growing within a single run = cold extrusion fingerprint

The chart backgrounds get tinted: **green (0–30 stable)** / **yellow (30–60 moderate)** / **red (60+ likely cold extrusion)**. The first flow where the score crosses 30 is marked as "cold extrusion onset" — purely visual, doesn't affect detection.

> **Why intra-run SG drift matters:** Within a single 5-second run, if SG drifts steadily toward higher load, the filament was getting harder to push as it heated. This is the physical signature of cold extrusion (vs. abrupt motor slip, which jumps suddenly without a gradient). Hover any thermal-chart point to see the per-step drift value.

---

## Per-driver tuning

Each driver family has its own `TriggerProfile` in the source code with all detection thresholds. To tune sensitivity for one driver without affecting the others, edit only its profile class.

- **TMC5160Profile** — validated production baseline (SGT=15, SFILT=1)
- **TMC2240Profile** — inherits TMC5160 base, with TMC2240-specific overrides:
  - `WARMUP_DRIFT_THRESHOLD = 0.04` — catches the systematic first-run drift
  - `PLATEAU_RATIO = 0.2` — the SG2 saturation curve is steeper on TMC2240
- **TMC2209Profile** *(experimental)* — fundamentally different detection strategy from SG2:
  - SG-magnitude triggers (snap-back, plateau, max-spike) **DISABLED** — SG4 doesn't follow the SG2 "smooth load curve" model
  - CV-spike triggers as primary slip indicator (`CV_HIGH_VARIANCE = 12.0` vs. 5.0 for SG2)
  - IQR triggers calibrated for the 0-510 SG4 scale (`IQR_ABSOLUTE_TRIGGER = 50` vs. 25 for SG2)
  - `WARMUP_DRIFT_THRESHOLD = 0.15` — SG4 has higher first-run drift
