# TMC Flow Test

*Part of [Crydteam Tools](../README.md) — Klipper plugins by Steven (Fragmon).*

**Find your printer's real max flow rate in 30 minutes — without test prints.**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Klipper](https://img.shields.io/badge/Klipper-compatible-green.svg)](https://www.klipper3d.org/)
[![Kalico](https://img.shields.io/badge/Kalico-compatible-brightgreen.svg)](https://github.com/KalicoCrew/kalico)

Plugin by **Steven (Fragmon) — Crydteam**
[![YouTube](https://img.shields.io/badge/YouTube-@crydteamprinting-red?logo=youtube)](https://www.youtube.com/@crydteamprinting)

<p align="center">
  <img src="images/results.png?v=1.0" alt="TMC Flow Test HTML report" width="700">
</p>

---

## What it does

The plugin runs your extruder at progressively higher flow rates and watches the TMC driver's stall signal to find the exact point where the motor starts losing grip. Output: an interactive dashboard showing your max safe flow, the slip onset curve, and recommended slicer values (80% conservative / 90% aggressive).

**Result:** A printer-specific flow curve in 30 minutes, with the same level of detail as a CNC Kitchen video — for free.

## Changes in this fork

Full detail and validation runs: [MODIFICATIONS.md](../MODIFICATIONS.md).

> **Why:** these changes were made because StallGuard would not measure
> or trigger properly on one setup — a TMC5160 driving the extruder with
> a **0.05 Ω sense resistor**, sized for far more current than an
> extruder runs. There, SG_RESULT tracked motor speed rather than load
> and missed stalls. With a sense resistor sized for extruder currents,
> stock StallGuard detection likely works and these extras are a safety
> net rather than a necessity.

What it adds on top of the stock StallGuard detection:

- **Stop on filament failure.** The run stops on the spot, instead of
  bisecting on, when:
  - the **heater duty collapses** (more flow commanded, less heat taken
    — the filament stopped arriving), or
  - a **BDWidth** motion sensor sees under 60 % of the commanded feed, or
  - a **`[filament_motion_sensor]`** (BTT SFS etc.) sees no pulse for
    too long.

  The feed checks run after every rep, so a step ends after the first
  rep that shows the filament stopped instead of grinding through all
  five. The sensors are auto-detected; without one, heater collapse
  alone does the job.
- **Glossy filament guard.** BDWidth is optical and some filaments read
  far low. If it reads under 80 % or erratically on the first four
  steps, its checks become report-only for that run (heater collapse
  still stops a real stall).
- **Melt limit** (report only): where heat per mm³ falls below 85 % of
  the low-flow line — usually the number for the slicer.
- **Slip start** (report only, BDWidth): the flow from which feed stays
  under 90 % — partial slip that does not stop the test.
- **Heat-up, idle baselines and MPC:** `TEMP=` heats and waits within
  1 °C; idle heater duty is measured before and after (drift
  correction); `MATERIAL=` records the filament and sets Kalico MPC's
  filament constants for the run.
- **Chamber fans:** `chamber_fans:` in the settings file runs filter /
  exhaust fans for the whole test.
- **Reports:** one console line per step (StallGuard, heater vs
  expected, feed), extra CSV columns, a **Melt & feed** tab, a full data
  table and a **Compare runs** tab for overlaying other CSVs (e.g. the
  same filament at several temperatures).

## Supported drivers

| Driver | Status | Notes |
|---|---|---|
| **TMC5160** / TMC2130 | ✅ Production-ready | full feature set |
| **TMC2240** | ✅ Production-ready | full feature set |
| **TMC2209** | ⚠️ Experimental | works on some hardware, see [docs/ADVANCED.md](docs/ADVANCED.md#tmc2209-experimental) |

---

## Install

```bash
cd ~
git clone -b awd-sync-and-flow-detection https://github.com/SeeingREDTonight/crydteamtools.git
cd crydteamtools
./install.sh max_flow_test
```

The installer symlinks `tmc_flow_test.py` into Klipper's extras folder and the
optional UI macros into the `crydteamtools/` folder of your config directory
(run it without arguments to pick plugins interactively). It also drops a
fully commented `tmc_flow_test_settings.cfg` into that folder and includes it at
the top of `printer.cfg` — to enable the plugin, uncomment the
`[tmc_flow_test]` line in that file, then `sudo systemctl restart klipper`.

Switching an existing upstream install to this fork: see the
[main README](../README.md#switching-an-existing-install-from-upstream).

**Updates**: `cd ~/crydteamtools && git pull`, then
`sudo systemctl restart klipper` (`FIRMWARE_RESTART` does not reload the
plugin's Python code).

## Configure

In `tmc_flow_test_settings.cfg` (or directly in `printer.cfg`):

```ini
[tmc_flow_test]
extruder_stepper: extruder
filament_diameter: 1.75
melt_zone_length: 42         # see hotend table below
```

Your existing `[tmcXXXX extruder]` section also needs `driver_SFILT: 1` (TMC5160/2240 only). See [docs/ADVANCED.md](docs/ADVANCED.md#configuration) for full driver-specific settings.

Optional, for filter / exhaust fans during the test:

```ini
chamber_fans: Nevermore, chamber
```

`[fan_generic]` fans run at 100 %; `[temperature_fan]` fans get their
target set to 1 °C (or `min_temp`, if higher) so they run. Both are restored afterwards.
Existing installs keep their settings file — add the line by hand.

**Melt zone length**:
- V6 / Revo Six: ~13 mm
- Volcano: ~21 mm
- Mosquito: ~18 mm
- Chube Air: ~42 mm
- Goliath / CHC-XL: ~50 mm

## Control panel (GUI)

New to the plugin? Run **`TMC_FLOW_GUI`** — it writes a beginner-friendly web
page into the Flowtest folder that builds the commands for you: a live
checklist of your driver (and every StallGuard problem it finds), hotend
presets, plain-language help on every field, a runtime estimate, and a
**Copy** / **Send to printer** button for the finished command.

**[How to use it →](docs/GUI.md)**

---

## Run a test

Move the toolhead **above the bed with at least 50 mm of clearance** — the test extrudes 2-5 meters of filament. Then in your printer console:

```
TMC_FLOW_FIND_MAX MAX=120 START=25 COARSE_STEP=5
```

Test takes ~10–15 minutes. CSV and HTML report land in `~/printer_data/config/Flowtest/`. Open the HTML in your browser.

With this fork's options (heat to 240 °C, record PLA, set MPC for it):

```
TMC_FLOW_FIND_MAX START=25 COARSE_STEP=5 MAX=100 TEMP=240 MATERIAL=PLA
```

| Parameter | Default | Purpose |
|---|---|---|
| `TEMP` | 0 | Heat the hotend and wait to within 1 °C. 0 = use the current target. |
| `MATERIAL` | (blank) | Filament type (`PLA`, `PETG`, `ABS`, …). Recorded in the report; on a Kalico MPC hotend also sets MPC's filament constants for the run (via `_SET_MPC_MATERIAL` if you have it) and restores them afterwards. |
| `SETTLE` / `IDLE` | 20 / 60 | Seconds of settle and idle heater measurement before and after the run. `IDLE=0` turns it off (and the melt limit with it). |
| `BDWIDTH` | auto | BDWidth sensor name, `auto` or `none`. |
| `BDW_MIN_RATIO` | 0.6 | Feed-check limit; 0 = record only. |
| `ENCODER` / `ENCODER_GAP` | auto / 0 | `[filament_motion_sensor]` name and the mm without a pulse that counts as stopped (0 = auto). Used only without BDWidth. |

The `CRYD_FLOW_TEST` macro and the GUI do not pass these on — type
`TMC_FLOW_FIND_MAX` directly to use them. Use `COARSE_STEP=5` and
`START` ≥ 25: with 10 mm³/s steps the stock StallGuard plateau trigger
tends to fire on the normal curve.

### Picking the right START / MAX

Set MAX to roughly 1.5× what you think your hotend can flow:

| Hotend type | Command |
|---|---|
| V6 stock 0.4 (~15 mm³/s) | `TMC_FLOW_FIND_MAX MAX=25 START=5 COARSE_STEP=5` |
| Volcano / CHT 0.4 (~25 mm³/s) | `TMC_FLOW_FIND_MAX MAX=40 START=10 COARSE_STEP=5` |
| Rapido HF 0.4 (~30 mm³/s) | `TMC_FLOW_FIND_MAX MAX=60 START=10 COARSE_STEP=5` |
| Rapido HF 0.6 (~50 mm³/s) | `TMC_FLOW_FIND_MAX MAX=80 START=15 COARSE_STEP=5` |
| Goliath / Chc XL 0.8 (~80 mm³/s) | `TMC_FLOW_FIND_MAX MAX=120 START=25 COARSE_STEP=5` |
| Unknown | `TMC_FLOW_FIND_MAX MAX=150 START=10 COARSE_STEP=5` |

If the test reaches MAX without detecting slip, **double MAX** and re-run.

## Reading the report

The HTML report shows:

- **Hero panel** — your max safe flow rate, plus 80%/90% slicer-ready values. When the run stopped on a filament failure it shows the last step before the stop (mechanical limit), with the **melt limit** and **slip start** underneath — for print quality, the melt limit is usually the number to use
- **Insight cards** — quick-glance status: result quality, where slip first triggered, thermal headroom, driver config
- **Three charts** (click through tabs):
  - **StallGuard signal** — the load curve. Sudden flat spots = slip onset. Click any line in the legend to hide/show it.
  - **Thermal profile** — heater PWM, temperature drop, residence time. Tells you if heater could keep up.
  - **Run-to-run variance** — bar chart showing measurement consistency. Coloured bars: green = clean, yellow = borderline, red = slip likely.
- **Melt & feed** tab — heat per mm³ against the low-flow line, and BDWidth feed per step
- **Compare runs** tab — load CSVs from other runs to overlay them, with a table of temperature, material, melt limit and result per run
- **Decision timeline** — phase-by-phase explanation of how the result was chosen
- **Test details** — full data table, configuration snapshot, decision trail

The chart background is tinted in three zones: green (safe), yellow (cold-extrusion suspected), red (slip detected).

## Use the result

Copy the **80% value** from the hero panel into your slicer's max volumetric speed setting. That gives you ~20% safety margin for filament variation, temperature drift, and longer prints. The 90% value works after you've validated it with a few real-world prints.

---

## Troubleshooting

**Test reaches MAX without trigger** → your hotend can flow more than MAX. Double MAX and re-run.

**Trigger fires too early (e.g. flow=50 on a high-flow setup)** → check `driver_SFILT: 1` is in your TMC config, restart Klipper, re-test.

**`Section 'tmc_flow_test' is not a valid config section`** → run `sudo systemctl restart klipper` (NOT `FIRMWARE_RESTART`) and verify the symlink exists in `~/klipper/klippy/extras/tmc_flow_test.py`.

**TMC2240 results below ~70 mm³/s on a fast extruder** → you're in StealthChop mode. Remove `stealthchop_threshold` from your `[tmc2240]` section.

**Test stopped with "BDWidth saw xx % of the expected feed" but the filament was fine** → check the console for the glossy-filament note; on a filament BDWidth reads badly but steadily low the check should switch to report-only after four steps. If it still stops early, run with `BDW_MIN_RATIO=0` (record only) or `BDWIDTH=none`.

**Test stopped on "heater duty fell"** → the filament stopped feeding (strip, stall or clog). Trim past the chewed section and check for a partial clog before re-running.

**More issues** → see [docs/ADVANCED.md#troubleshooting](docs/ADVANCED.md#troubleshooting).

---

## Need more?

- **Control panel (GUI)** → [docs/GUI.md](docs/GUI.md)
- **Driver-specific tuning, chopper mode internals, TMC2209 setup** → [docs/ADVANCED.md](docs/ADVANCED.md)
- **Command reference** → [docs/COMMANDS.md](docs/COMMANDS.md)
- **How slip detection works** → [docs/INTERNALS.md](docs/INTERNALS.md)

---

## Credits & License

Released under **GPL-3.0**.

Inspired by Klipper's StallGuard implementation and [klipper_tmc_autotune](https://github.com/andrewmcgr/klipper_tmc_autotune).

Plugin author: **Steven (Fragmon) — Crydteam** · [YouTube: @crydteamprinting](https://www.youtube.com/@crydteamprinting)

Contributions and feedback welcome — open an issue or PR.
