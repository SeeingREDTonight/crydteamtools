# Modifications from stock crydteamtools

This branch is based on the `dev` branch of
[Fragmon/crydteamtools](https://github.com/Fragmon/crydteamtools) at
commit `062d1dc` (speed_test v1.6.1, TMC Flow Test v1.2.2). Everything
else is unchanged from upstream. To see every changed line:

```bash
git diff 062d1dc -- speed_test max_flow_test install.sh
```

Developed and tested on an AWD CoreXY (4× TMC5160 on X/Y), Kalico,
TMC5160 extruder at 48 V, MPC hotend control and a BDWidth sensor, with
ABS and PLA.

Install exactly like upstream (`./install.sh`). After replacing a
plugin's `.py` file, restart the Klipper **service**
(`sudo systemctl restart klipper`) — `FIRMWARE_RESTART` does not reload
Python code.

---

## Installer (`install.sh`)

Stock writes the upstream URL and branch `main` into Moonraker's update
manager. Here it uses the clone's own `origin` remote and current
branch, so Moonraker updates come from this fork and branch.

---

## Speed test (`speed_test/`)

| Change | Stock behaviour | Modified behaviour |
|---|---|---|
| `park_x` / `park_y` options | During a single-axis test the other axis stays wherever homing left it (usually an axis end). | The other axis is moved to `park_x` / `park_y` first (while the tested axis is still at home), then the test runs. Keeps the toolhead clear of obstacles near the ends of travel, e.g. front motor towers on AWD printers. |
| `skip_resync_gcode` option | After lost steps the test only re-homes. Re-homing does not restore the phase between the two motors of an AWD belt pair, so later steps can run on desynced motors and read low. | After any step that loses steps, the configured G-code (e.g. `SYNC_MOTORS AXES=x,y` from MRX8024/motors-sync, or `MOTOR_SYNC` from this suite) runs at travel speed/accel, then the next step re-homes. Empty = off. |

Both options are listed (commented out) in `speed_test_settings.cfg`.

---

## Flow test (`max_flow_test/`)

### Why

On a geared extruder with plenty of torque, StallGuard's SG_RESULT
mostly tracks extrusion **speed** (SG × flow stayed constant within a
few % across every run), not load. The stock SG triggers missed real
stalls and strips, and the plateau trigger fired on the normal shape of
the SG curve. The changes below add detection that does not depend on
StallGuard. The stock SG triggers are all still active.

### New stop checks (end the run instead of bisecting)

| Check | Fires when | Notes |
|---|---|---|
| Heater collapse | Average heater duty drops > 10 % and ≥ 3 points vs the step at the highest lower flow, while flow went up | Filament stopped arriving (stall, strip, clog). Caught every stall/strip in testing; never fired on a clean step. |
| BDWidth feed | The bdwidth motion counter saw < `BDW_MIN_RATIO` (default 0.6) of the commanded filament, scale-corrected by earlier steps | Needs a bdwidth sensor; auto-detected. Runs from the first reference step, so a jam is caught immediately. Optical tracking depends on the filament surface, so the limit is deliberately coarse. |

When either fires, the test stops on the spot (bisecting further only
grinds damaged filament), reports the last step that fed cleanly as the
mechanical limit, saves the report, restores SGT / fan / bdwidth, then
measures the end-of-run idle baseline (safe to cancel).

### Melt limit (report only, never stops the test)

Fits a straight line heater-duty = a + b·flow through steps 2–5 (drift
corrected with the idle baselines, failed steps excluded). While the
plastic fully melts every step sits on that line; the first flow below
85 % of it is reported as the **print-quality limit** — console, CSV/HTML
header `melt_limit`, per-step `melt_fit_pct`. This is usually the number
to use in the slicer; the mechanical limit is where the extruder gives
up.

### New parameters for `TMC_FLOW_FIND_MAX`

| Parameter | Default | Purpose |
|---|---|---|
| `TEMP` | 0 | Heat the hotend and wait to within 1 °C (stock commands, not `M109`, which macro packs often override with a loose window). 0 = use the current target. |
| `SETTLE` | 20 | Seconds without extrusion before each idle baseline. |
| `IDLE` | 60 | Seconds of idle heater-duty measurement at the start and end. 0 = off (also disables the melt limit). |
| `BDWIDTH` | auto | bdwidth sensor name, `auto` or `none`. |
| `BDW_MIN_RATIO` | 0.6 | Feed-check limit; 0 = record only. |

The `CRYD_FLOW_TEST` macro does not forward these — call
`TMC_FLOW_FIND_MAX` directly, e.g.

```
TMC_FLOW_FIND_MAX START=25 COARSE_STEP=5 MAX=100 TEMP=240 KEEP_SGT=1
```

Use `COARSE_STEP=5` and `START` ≥ 25 (enables the two reference steps).
With 10 mm³/s steps the stock SG plateau trigger tends to fire on the
normal SG curve.

### Data capture

- New CSV columns: `t_s, idle_pwm, melt_duty, melt_duty_per_flow,
  melt_fit_pct, bdw_cmd_mm, bdw_meas_mm, bdw_ratio`.
- New header lines: idle baselines, fit intercept/slope, `melt_limit`,
  bdwidth mode.
- HTML: new **Melt & feed** tab, full per-step data table with every CSV
  column, melt limit under the result, "Stopped on filament failure"
  status instead of "Verified" for stopped runs.

### Other fixes

- bdwidth is switched to motion-only (and `M221 S100`) for the test and
  its runout logic is pushed out of reach, so it neither changes the
  flow, spams the console nor pauses the printer mid-test; restored at
  the end.
- HTML charts on hidden tabs drew their annotation lines in the wrong
  place (also the stock CV-trigger line); charts are redrawn when their
  tab is opened.

### Validation runs (same printer)

| Material | Temp | Mechanical result | Melt limit | Stopped by |
|---|---|---|---|---|
| PLA (matte black) | 220 °C | 60 | 50 | BDWidth feed |
| PLA (matte black) | 240 °C | 65 | 60 | Strip: heater + feed |
| PLA (matte black) | 260 °C | 75 | 60 | BDWidth feed |
| ABS | 260 °C | 69 (SG plateau) | 65 | — |
| ABS | 280 °C | none ≤ 120 | 65 | — |

Thresholds come from this handful of runs on one printer — treat them as
a starting point, and please share CSVs if a check fires on a step that
was clearly feeding.
