# Modifications from stock crydteamtools

This branch is based on the `dev` branch of
[Fragmon/crydteamtools](https://github.com/Fragmon/crydteamtools) at
commit `062d1dc` (speed_test v1.6.1, TMC Flow Test v1.2.2). Everything
else is unchanged from upstream. To see every changed line:

```bash
git diff 062d1dc -- speed_test max_flow_test pa_test install.sh
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
Re-running it also corrects `origin` / `primary_branch` in an existing
entry.

All macro and settings `.cfg` files go into a `crydteamtools/` folder in
the config directory (stock: the config root), included as
`[include crydteamtools/<file>]`. Existing includes are recognised
wherever they are, also a folder-wide `[include crydteamtools/*.cfg]`.
Files from an older install are moved into the folder and their
includes repointed; a root copy that clashes with a file already in the
folder is archived. Uninstall removes the folder once it is empty.

---

## PA test (`pa_test/`)

| Change | Stock behaviour | Modified behaviour |
|---|---|---|
| Module and section name | Installs as `klippy/extras/pa_test.py` with section `[pa_test]`. That is the name of Klipper/Kalico's own PA tower module, so the installer's symlink replaces the firmware file: `PRINT_PA_TOWER` disappears, the firmware repo shows a modified tracked file, and an existing `[pa_test]` tower section loads the StallGuard probe instead. | Installs as `sg_pa_test.py` with section `[sg_pa_test]`, next to the firmware's module. `./install.sh pa_test` migrates an old install: removes the `pa_test.py` symlink, restores the firmware file with `git checkout`, and renames the section header in the copied `pa_test_settings.cfg`. Commands (`PA_TEST_PROBE`, `PA_TEST_GUI`) are unchanged. |

---

## Speed test (`speed_test/`)

| Change | Stock behaviour | Modified behaviour |
|---|---|---|
| `park_x` / `park_y` options | During a single-axis test the other axis stays wherever homing left it (usually an axis end). | The other axis is moved to `park_x` / `park_y` first (while the tested axis is still at home), then the test runs. Keeps the toolhead clear of obstacles near the ends of travel, e.g. front motor towers on AWD printers. |
| `skip_resync_gcode` option (default `auto`) | After lost steps the test only re-homes. Re-homing does not restore the phase between the two motors of an AWD belt pair, so later steps can run on desynced motors and read low. | After any step that loses steps on a dual-motor axis: accelerometer sync (`SYNC_MOTORS`, MRX8024/motors-sync) if loaded, else StallGuard sync (`MOTOR_SYNC AXIS=…`, this suite's `[motor_sync]`) if enabled, else a console note that the re-sync was skipped and the motors may now be out of sync. Runs at travel speed/accel, then the next step re-homes. Single-motor axes: nothing happens. `none` = off; any other value is run as literal G-code. |

Both options are listed (commented out) in `speed_test_settings.cfg`.

---

## Chopper test (`chopper_test/`)

Stock writes each candidate register value to `stepper_<axis>` only and
judges it by lost steps after a re-home. On an AWD CoreXY four motors
drive every X or Y move, so three kept their normal settings and held
position while the tested one could lose sync audibly — bad candidates
passed.

- Candidates are written to every motor that drives the move (CoreXY:
  all `stepper_x*` / `stepper_y*` with the same driver type; cartesian:
  all motors of the axis). Each driver's own starting values are
  restored at the end, on errors too, and before a re-sync.
- `skip_resync_gcode` (default `auto`, as in the speed test): after a
  measurement that lost steps the dual-motor belts are re-synced before
  the next one, so later candidates are not tested on an out-of-phase
  pair.
- A failed candidate still gets one retry, but the console shows
  `OK on retry (lost N µsteps the first time)` and the CSV records the
  first result (`first_try_lost_usteps`). Stock printed a plain `OK`.
- After the last round of a pass the acceleration is not raised any
  more. Stock announced a raise, never tested it, and started the next
  pass at that untested acceleration. If nothing holds in a later round,
  the next pass starts at the last acceleration where something held.
- `park_x` / `park_y` park the other axis during the test; unset, the
  `[speed_test]` values are used.
- X/Y are always homed before the first measurement, so its reference
  is a homed position.

## Flow test (`max_flow_test/`)

### Why

> **These flow-test changes were made because StallGuard would not
> measure or trigger appropriately on my setup:** a TMC5160 driving the
> extruder with a **0.05 Ω sense resistor**. That resistor is sized for
> high-current motors, so at extruder currents (~0.8 A) the driver works
> in a small part of its current-sense range and StallGuard has little
> resolution. With a sense resistor sized for extruder currents, the
> stock StallGuard detection would most likely work and these changes
> would probably not be needed. They are still useful as an independent
> cross-check, but treat them as a workaround for this kind of driver
> setup, not a replacement for a correctly sized driver.

On my geared extruder with plenty of torque, StallGuard's SG_RESULT
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
| Encoder feed | Used when no bdwidth is found: a Klipper `[filament_motion_sensor]` (BTT SFS etc.) saw no pulse while the extruder pushed more than `ENCODER_GAP` (default 2 × `detection_length`, min 10 mm) | Auto-detected. Encoders keep no counter, so this is a "filament stopped" check, not a ratio. The sensor's runout action is disabled during the test and restored after. |

With no motion sensor at all, the feed checks are skipped (console note)
and heater collapse alone catches stalls and strips.

BDWidth is optical, and some filaments track badly: a glossy black PLA
read 45–71 % of the real feed, jumping ±10 points between steps, while
it fed normally. So the sensor is judged on the first four steps. If the
median reading is below 80 %, or two of the four are more than 10 % off
the median, the feed check and slip start become report-only for that
run (console note, CSV header `bdwidth_untrusted`) and heater collapse
alone stops the test. Until those four steps exist, only a jam-level
reading (under 30 %) stops it, so a jam from the first step is still
caught at once. The scale correction allows 0.4–1.5 (was 0.8–1.2).

The feed checks also run after every rep inside a step (5 reps by
default): when one rep already shows the filament stopped, the step ends
there instead of grinding through the remaining reps. The CSV
`n_repeats` column shows how many reps ran. The heater checks still use
the whole step (one rep is too noisy for them).

When any check fires, the test stops on the spot (bisecting further only
grinds damaged filament), reports the last step before the filament
stopped as the mechanical limit, saves the report, restores SGT / fan / bdwidth, then
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
| `MATERIAL` | (blank) | Filament type (`PLA`, `PETG`, `ABS`, …). Recorded in CSV/HTML and checked by the Compare tab. On a Kalico MPC hotend it also sets MPC's filament density / heat capacity for the test (via `_SET_MPC_MATERIAL` if present, else the table from Kalico's MPC docs, i.e. `MPC_SET HEATER=extruder FILAMENT_DENSITY=… FILAMENT_HEAT_CAPACITY=…`) and restores the previous values afterwards, even on abort. Blank = stock behaviour. |
| `TEMP` | 0 | Heat the hotend and wait to within 1 °C (stock commands, not `M109`, which macro packs often override with a loose window). 0 = use the current target. |
| `SETTLE` | 20 | Seconds without extrusion before each idle baseline. |
| `IDLE` | 60 | Seconds of idle heater-duty measurement at the start and end. 0 = off (also disables the melt limit). |
| `BDWIDTH` | auto | bdwidth sensor name, `auto` or `none`. |
| `BDW_MIN_RATIO` | 0.6 | Feed-check limit; 0 = record only. |
| `ENCODER` | auto | `[filament_motion_sensor]` name, `auto` or `none` (only used without bdwidth). |
| `ENCODER_GAP` | 0 (auto) | mm of extrusion without an encoder pulse that counts as "stopped". |

The `CRYD_FLOW_TEST` macro does not forward these — call
`TMC_FLOW_FIND_MAX` directly, e.g.

```
TMC_FLOW_FIND_MAX START=25 COARSE_STEP=5 MAX=100 TEMP=240 MATERIAL=PLA KEEP_SGT=1
```

Use `COARSE_STEP=5` and `START` ≥ 25 (enables the two reference steps).
With 10 mm³/s steps the stock SG plateau trigger tends to fire on the
normal SG curve.

`BDWIDTH`, `BDW_MIN_RATIO`, `ENCODER`, `ENCODER_GAP`, `SETTLE` and `IDLE`
take their defaults from the config options `bdwidth`, `bdw_min_ratio`,
`encoder`, `encoder_gap`, `idle_settle` and `idle_measure` (listed in the
settings template), so `CRYD_FLOW_TEST` and the GUI — which pass none of
them — can use them too. `TEMP` and `MATERIAL` are command-only.

### Slip start (report only, needs bdwidth)

The 60 % feed check only catches a filament that stopped; partial slip
(e.g. 84 % feed) passes it. The run therefore also reports where slip
starts: the lowest flow from which every higher step reads below 90 % of
the low-flow feed reading (a single noisy dropout lower down is
ignored), plus the last clean step — console, HTML and CSV header
`slip_onset`.

### Data capture

- New CSV columns: `t_s, idle_pwm, melt_duty, melt_duty_per_flow,
  melt_fit_pct, bdw_cmd_mm, bdw_meas_mm, bdw_ratio, enc_gap_mm`.
- New header lines: idle baselines, fit intercept/slope, `melt_limit`,
  bdwidth mode.
- HTML: new **Melt & feed** tab, full per-step data table with every CSV
  column, melt limit under the result, "Stopped on filament failure"
  status instead of "Verified" for stopped runs.
- HTML **Compare runs** tab: load CSVs from other runs (e.g. the same
  filament at other temperatures) to overlay heat-vs-line, BDWidth feed
  and StallGuard, with a table of target temp, melt limit and result per
  run. Computed in the browser from the CSVs (older CSVs work too). The
  CSV that produced the report and files already loaded are rejected;
  a different or unrecorded material is flagged.
- CSV header now also records `result_max_flow` and `result_quality`.
- Console: one line per step with StallGuard (median, run-to-run CV,
  dips), heater duty vs what the low-flow line expects (from step 6 on)
  and BDWidth measured / commanded mm or the encoder state, instead of
  separate lines.

### Chamber fans (config `chamber_fans`)

Optional list of fans to run for the whole test (air filtering), e.g.
`chamber_fans: Nevermore, chamber`. `[fan_generic]` fans run at 100 %;
`[temperature_fan]` fans get their target set to 1 °C (or `min_temp`, if higher) so they run.
Previous speeds / targets are restored at the end, also on abort.

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
