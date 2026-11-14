> **Modified fork.** This branch adds changes to the speed test and the
> flow test on top of upstream `dev` (commit `062d1dc`). See
> [MODIFICATIONS.md](MODIFICATIONS.md) for what changed and why.

# Crydteam Tools for Klipper

A collection of Klipper diagnostic and tuning plugins by **Steven (Fragmon) — Crydteam**.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Klipper](https://img.shields.io/badge/Klipper-compatible-green.svg)](https://www.klipper3d.org/)
[![Kalico](https://img.shields.io/badge/Kalico-compatible-brightgreen.svg)](https://github.com/KalicoCrew/kalico)
[![YouTube](https://img.shields.io/badge/YouTube-@crydteamprinting-red?logo=youtube)](https://www.youtube.com/@crydteamprinting)

---

## Plugins

| Plugin | What it does | Docs |
| ------ | ------------ | ---- |
| **[Speed Test](speed_test/README.md)** | Finds a motor's real limits: max velocity, max acceleration per speed (limit map), square-corner velocity and the lowest run_current that still holds them — validated with a simulated print. CSV + interactive HTML reports, beginner-friendly control-panel GUI. | [speed_test/](speed_test/README.md) |
| **[TMC Flow Test](max_flow_test/README.md)** | Finds your printer's real max flow rate in ~30 minutes without test prints: runs the extruder at rising flow rates and watches the TMC StallGuard signal for slip onset. Interactive dashboard with recommended slicer values. | [max_flow_test/](max_flow_test/README.md) |
| **[Motor Sync](motor_sync/README.md)** | Synchronizes the two motors of a dual-motor axis (AWD, dual gantry) — **without an accelerometer**: out-of-phase motors fight each other, the extra load shows in StallGuard, and the plugin shifts one motor in microsteps until the fight is gone. Sync approach adapted from [MRX8024/motors-sync](https://github.com/MRX8024/motors-sync). | [motor_sync/](motor_sync/README.md) |
| **[PA Test](pa_test/README.md)** *(prototype)* | Pressure-advance calibration via StallGuard: melt pressure loads the extruder motor, PA errors show up as step-response over/undershoot. Currently ships the feasibility probe; goal is an auto-generated adaptive-PA table (PA vs. flow) for OrcaSlicer. | [pa_test/](pa_test/README.md) |
| **[Chopper Test](chopper_test/README.md)** | Tunes the TMC chopper registers (TOFF, TBL, HSTRT/HEND, TPFD) for **more torque and less heat** — **without an accelerometer**: instead of scoring vibration, it scores the outcome, measuring lost steps under reversal stress and the motor's temperature rise. Sweep order adapted from [MRX8024/chopper-resonance-tuner](https://github.com/MRX8024/chopper-resonance-tuner). | [chopper_test/](chopper_test/README.md) |

More plugins will join this collection — each lives in its own folder with its
own README and docs.

---

## Installation

This fork lives on the `awd-sync-and-flow-detection` branch:

```bash
cd ~
git clone -b awd-sync-and-flow-detection https://github.com/SeeingREDTonight/crydteamtools.git
cd crydteamtools
bash ./install.sh
```

### Switching an existing install from upstream

If `~/crydteamtools` is already a clone of Fragmon's repo:

```bash
cd ~/crydteamtools
git remote set-url origin https://github.com/SeeingREDTonight/crydteamtools.git
git fetch origin
git checkout -B awd-sync-and-flow-detection origin/awd-sync-and-flow-detection
bash ./install.sh            # choose "update"
sudo systemctl restart moonraker
sudo systemctl restart klipper
```

Re-running the installer points the existing `[update_manager crydteamtools]`
entry at this fork and branch, so Mainsail/Fluidd updates come from here.
Your `*_settings.cfg` files are kept; new options (e.g. `chamber_fans`,
`park_x`) are documented in the templates in each plugin folder and can be
copied in by hand.

The installer lets you **pick which plugins to install** (or pass them
directly: `./install.sh all`, `./install.sh speed_test`,
`./install.sh max_flow_test`). It does everything automatically:

- symlinks the plugin files into Klipper's `extras` folder
- symlinks the UI macros into a `crydteamtools/` folder in your config
  directory
- copies a fully commented `…_settings.cfg` template per plugin into the
  same folder (your edits survive updates)
- adds `[include crydteamtools/…]` lines for macros and settings to the
  top of `printer.cfg` — unless they are already included somewhere,
  including a folder-wide `[include crydteamtools/*.cfg]`
- moves files from older installs (config root) into the folder and
  repoints their `[include]` lines; if the folder already has a file of
  the same name, the root copy is archived and the folder's is kept
- registers the repo with Moonraker's update manager

The only manual step left: open the plugin's settings file in
`crydteamtools/` (e.g. `speed_test_settings.cfg` / `tmc_flow_test_settings.cfg`) in the
Mainsail/Fluidd config editor, uncomment the section header
(`[speed_test]` / `[tmc_flow_test]`), then run `sudo systemctl restart klipper`. Every
option is documented in the file with its default value.

Re-running the installer on an existing install asks whether you want to
**update** (re-link the plugin files and add any missing `[include]` lines —
your `*_settings.cfg` files are never overwritten), skip the installed ones,
or abort.

## Uninstalling

```bash
cd ~/crydteamtools
./install.sh uninstall            # interactive selection
./install.sh uninstall all        # or name plugins directly
```

This removes the plugin's symlinks from Klipper's `extras`, deletes its
`[include]` lines from every config file, and — once no Crydteam plugin is
left — removes the `[update_manager crydteamtools]` entry from
`moonraker.conf`. Your settings files are **not deleted**: they are moved to
`~/printer_data/config/archive/`. Afterwards remove the plugin's config
section (e.g. `[speed_test]`) if you had added it to `printer.cfg`, then run
`FIRMWARE_RESTART`.

## Updating

The installer registers the repo with **Moonraker's update manager**, so new
releases appear directly in Mainsail/Fluidd (Machine → Update Manager) —
one-click update, Klipper is restarted automatically.

Manual alternative:

```bash
cd ~/crydteamtools && git pull
sudo systemctl restart klipper
```

Use `sudo systemctl restart klipper`, not `FIRMWARE_RESTART`: a firmware
restart does not reload the plugins' Python code. Re-run `./install.sh` only
if new plugins were added.

### Getting upstream changes

This fork does not follow Fragmon's repo automatically. Upstream changes
are merged into the fork's branch by its maintainer; you then get them with
the normal update above.

---

## Credits & License

Plugins by Steven (Fragmon) — Crydteam · YouTube: [@crydteamprinting](https://www.youtube.com/@crydteamprinting)

Released under the GNU General Public License v3.0. Per-plugin credits are
listed in each plugin's README.

## Contributing

Issues and pull requests welcome. If you've tested on hardware not listed in
the docs, let me know — I'd love to add confirmed-working markers.
