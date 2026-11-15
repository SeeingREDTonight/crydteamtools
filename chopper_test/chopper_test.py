# Chopper Test — TMC chopper tuning by torque and temperature,
# without an accelerometer.
# credits:
#   Steven (Fragmon) — Crydteam
#   YouTube: https://www.youtube.com/@crydteamprinting
#
#   The idea of sweeping the TMC chopper registers (TBL, TOFF, HSTRT,
#   HEND, TPFD) field by field, and the staged order in which they are
#   tuned, come from "chopper-resonance-tuner" by Maksim Bolgov
#   (MRX8024), GPLv3:
#   https://github.com/MRX8024/chopper-resonance-tuner
#
#   The MEASUREMENT is different and needs no extra hardware. That tool
#   scores each register combination by the vibration an ADXL
#   accelerometer records. StallGuard cannot substitute for that,
#   because TBL and TOFF define the very chopper window SG is sampled
#   in — SG readings are not comparable across those settings.
#
#   So this plugin scores the OUTCOME instead of the symptom:
#     • torque      — does the motor still hold a hard reversal-stress
#                     move without losing steps (MCU step counter
#                     compared across a re-home, as in speed_test)?
#     • temperature — how warm does the motor get during a fixed run?
#   Both are what chopper tuning is actually for.
#
# License: GPLv3

import logging
import math
import os
import time

MODULE_NAME = "Chopper Test"
MODULE_VERSION = "0.1.0"

TRINAMIC_DRIVERS = ("tmc2130", "tmc2209", "tmc2240", "tmc5160", "tmc2660")
# Drivers whose CHOPCONF carries a TPFD field (passive fast decay).
# NOT present on TMC2209 (reserved bits) and NOT on TMC2130 (those bits
# are the `sync` field there — writing them would do something else).
TPFD_DRIVERS = ("tmc5160", "tmc2240")

# Internal clock used to convert TPWMTHRS back into a velocity.
TMC_FREQUENCY = {
    "tmc2130": 13200000.,
    "tmc2209": 12000000.,
    "tmc2240": 12500000.,
    "tmc5160": 12500000.,
    "tmc2660": 15000000.,
}

# ─── Safe sweep space (SpreadCycle / chm=0 only) ──────────────────────
# toff=0 disables all bridges ("Driver disable, all bridges off"), and
# toff=1 is only legal with tbl>=2 — starting at 2 removes that ordering
# hazard entirely.
TOFF_VALUES = (2, 3, 4, 5, 6, 7, 8)
TBL_VALUES = (0, 1, 2, 3)
HSTRT_MAX = 7
HEND_MAX = 15
# Datasheet rule (TMC2130/5160/2209, chm=0): effective HEND+HSTRT <= 16.
# Effective HSTRT = reg+1, effective HEND = reg-3, so in the register
# units Klipper writes: hstrt + hend <= 18.
HSTRT_HEND_SUM_MAX = 18
TPFD_VALUES = tuple(range(0, 16, 2))

# Fields this plugin will never touch, and why:
#   vsense    changes the actual motor current by ~2x behind Klipper's
#             back — overcurrent / overheating hazard
#   chm       remaps hstrt/hend/fd3 to a different meaning; SpreadCycle
#             and constant-off-time are two incomparable sweep spaces
#   vhighchm / vhighfs / fd3 / disfdcc
#             inert while chm=0 and thigh=0, or dcStep-only

SETTLE_DELAY_MS = 200
# Rounds per pass, and how much the acceleration rises between them.
ESCALATIONS = 4
ESCALATION_FACTOR = 1.15


class ChopperTest:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.gcode = self.printer.lookup_object('gcode')
        self.reactor = self.printer.get_reactor()

        self.structure = config.get('structure', 'cartesian').lower()
        if self.structure not in ('cartesian', 'corexy'):
            raise config.error(
                "chopper_test: structure must be cartesian or corexy")
        self.default_axis = config.get('default_axis', 'X').upper()
        self.margin = config.getfloat('margin', 20.0, above=0.)
        self.z_pos = config.getfloat('z_pos', 20.0, minval=0.)
        self.max_missed = config.getfloat('max_missed', 1.5, above=0.)
        self.travel_speed = config.getfloat('travel_speed', 100.0, above=0.)
        # Optional Klipper temperature sensor mounted on the motor, e.g.
        # "temperature_sensor motor_x". Without it the thermal stage
        # falls back to the driver's own sensor (TMC2240) or is skipped.
        self.motor_sensor = config.get('motor_sensor', '').strip()
        # Where to park the OTHER axis during a test (e.g. clear of
        # motor towers). Unset = [speed_test]'s park_x / park_y if that
        # plugin is loaded, else the axis stays where homing left it.
        # Distance from the endstop where the stress move starts.
        # 0 = automatic (20 % of the usable travel).
        self.start_offset = config.getfloat('start_offset', 0.0, minval=0.)
        self.park_x = config.getfloat('park_x', None)
        self.park_y = config.getfloat('park_y', None)
        # Re-sync of dual-motor (AWD) belt pairs after lost steps:
        # auto | none | any G-code (same meaning as in [speed_test]).
        self.skip_resync_gcode = config.get(
            'skip_resync_gcode', 'auto').strip() or 'auto'
        self.thermal_seconds = config.getfloat(
            'thermal_seconds', 90.0, minval=10.)

        config_dir = os.path.expanduser('~/printer_data/config')
        if not os.path.isdir(config_dir):
            config_dir = os.path.expanduser('~')
        self.output_dir = os.path.expanduser(config.get(
            'output_dir', os.path.join(config_dir, 'Choppertest')))

        self._tmc_cache = {}
        self._tuned_cache = {}
        self._limits_snapshot = None
        self._runlog = None
        self._runlog_path = None
        self._originals = {}
        self._resync_noted = False

        self.gcode.register_command(
            'CHOPPER_TUNE', self.cmd_CHOPPER_TUNE,
            desc='Tune the TMC chopper registers by measured torque '
                 '(lost steps) and motor temperature — no accelerometer '
                 'needed')
        self.gcode.register_command(
            'CHOPPER_TEST_STATUS', self.cmd_STATUS,
            desc='Show the chopper registers of an axis and whether it '
                 'can be tuned')

    # ─── Run log ────────────────────────────────────────────────────
    # One line per event, written and flushed as it happens, so a run
    # that ends in a homing error or a shutdown still leaves everything
    # up to that point on disk. Mirrored into klippy.log.

    def _log_open(self, timestamp):
        self._runlog = None
        self._runlog_path = None
        self._log_t0 = self.reactor.monotonic()
        try:
            os.makedirs(self.output_dir, exist_ok=True)
            path = os.path.join(self.output_dir,
                                'chopper_%s.log' % timestamp)
            self._runlog = open(path, 'w', encoding='utf-8')
            self._runlog_path = path
        except Exception as e:
            logging.warning("chopper_test: cannot open run log: %s", e)

    def _log(self, event, **kv):
        line = "%8.2f %-9s %s" % (
            self.reactor.monotonic() - getattr(self, '_log_t0', 0.0),
            event,
            " ".join("%s=%s" % (k, kv[k]) for k in kv))
        logging.info("chopper_test: %s", line)
        f = getattr(self, '_runlog', None)
        if f is not None:
            try:
                f.write(line + "\n")
                f.flush()
                os.fsync(f.fileno())
            except Exception:
                pass

    def _log_close(self):
        f = getattr(self, '_runlog', None)
        self._runlog = None
        if f is not None:
            try:
                f.close()
            except Exception:
                pass

    @staticmethod
    def _fmt_combo(combo):
        return ",".join("%s:%s" % kv for kv in sorted(combo.items()))

    # ─── Axis / stepper helpers (same approach as speed_test) ───────

    def _stepper_name(self, axis):
        return 'stepper_' + axis.lower()

    def _get_axis_bounds(self, axis):
        try:
            cfg = self.printer.lookup_object('configfile')
            s = cfg.get_status(self.reactor.monotonic())['settings']
            ax = s[self._stepper_name(axis)]
            lo, hi = float(ax['position_min']), float(ax['position_max'])
        except Exception as e:
            raise self.gcode.error(
                "chopper_test: cannot read %s bounds: %s"
                % (self._stepper_name(axis), e))
        m = min(self.margin, 0.1 * (hi - lo))
        lo += m
        hi -= m
        return lo, hi, (lo + hi) / 2.0, hi - lo

    def _read_mcu_pos(self, axis):
        try:
            kin = self.printer.lookup_object('toolhead').get_kinematics()
            want = self._stepper_name(axis)
            for s in kin.get_steppers():
                if s.get_name() == want:
                    return s.get_mcu_position()
        except Exception:
            return None
        return None

    def _microsteps(self, axis):
        """Live value from the driver, falling back to the config."""
        info = self._lookup_tmc(axis)
        if info:
            try:
                return 256 >> info[1].fields.get_field('mres')
            except Exception:
                pass
        try:
            cfg = self.printer.lookup_object('configfile')
            s = cfg.get_status(self.reactor.monotonic())['settings']
            return int(s[self._stepper_name(axis)].get('microsteps', 16))
        except Exception:
            return 16

    def _lookup_tmc(self, axis):
        if axis in self._tmc_cache:
            return self._tmc_cache[axis]
        name = self._stepper_name(axis)
        for drv in TRINAMIC_DRIVERS:
            tmc = self.printer.lookup_object("%s %s" % (drv, name), None)
            if tmc is not None:
                self._tmc_cache[axis] = (drv, tmc)
                return self._tmc_cache[axis]
        self._tmc_cache[axis] = None
        return None

    def _field(self, axis, name):
        info = self._lookup_tmc(axis)
        if not info:
            return None
        try:
            return info[1].fields.get_field(name)
        except Exception:
            return None

    def _tuned_steppers(self, axis):
        """Every stepper that drives a move on `axis` and has the same
        TMC driver type as stepper_<axis>: [(stepper name, tmc), …].

        A candidate has to reach all of them. On CoreXY every X or Y
        move is driven by both belts, and an AWD machine has two motors
        per belt (stepper_x, stepper_x1, stepper_y, stepper_y1): with
        the candidate on one driver only, the other motors keep their
        old settings, carry the load and hide a bad candidate.
        """
        if axis in self._tuned_cache:
            return self._tuned_cache[axis]
        primary = self._stepper_name(axis)
        info = self._lookup_tmc(axis)
        out = []
        if info:
            prefixes = (('stepper_x', 'stepper_y')
                        if self.structure == 'corexy' else (primary,))
            names = [primary]
            try:
                kin = self.printer.lookup_object('toolhead').get_kinematics()
                for st in kin.get_steppers():
                    n = st.get_name()
                    if n not in names and n.startswith(prefixes):
                        names.append(n)
            except Exception:
                pass
            for n in names:
                tmc = self.printer.lookup_object(
                    "%s %s" % (info[0], n), None)
                if tmc is not None:
                    out.append((n, tmc))
        self._tuned_cache[axis] = out
        return out

    def _read_fields(self, tmc, names):
        vals = {}
        for n in names:
            try:
                vals[n] = tmc.fields.get_field(n)
            except Exception:
                pass
        return vals

    def _set_stepper_field(self, stepper, name, value):
        """Write through Klipper's own SET_TMC_FIELD: it validates the
        field, quantises the value and schedules the write at a print
        time, which a raw register write does not."""
        try:
            self.gcode.run_script_from_command(
                "SET_TMC_FIELD STEPPER=%s FIELD=%s VALUE=%d"
                % (stepper, name, int(value)))
            return True
        except Exception as e:
            logging.warning("chopper_test: SET_TMC_FIELD %s %s=%s failed: "
                            "%s", stepper, name, value, e)
            return False

    def _stealthchop_max_velocity(self, axis):
        """Speed (mm/s) below which StealthChop is used, or None.

        StealthChop applies while TSTEP >= TPWMTHRS, and TSTEP is the
        time between steps — large when slow. So TPWMTHRS is an UPPER
        VELOCITY limit for StealthChop: a large TPWMTHRS confines it to
        a crawl, and TPWMTHRS=0 means no limit (StealthChop always).
        Returns inf when StealthChop is unlimited, 0.0 when it is off.
        """
        tp = self._field(axis, 'tpwmthrs')
        if tp is None:
            return None
        if tp == 0:
            return float('inf')
        info = self._lookup_tmc(axis)
        freq = TMC_FREQUENCY.get(info[0] if info else '', 12.5e6)
        mres = self._field(axis, 'mres')
        try:
            kin = self.printer.lookup_object('toolhead').get_kinematics()
            want = self._stepper_name(axis)
            step_dist = None
            for s in kin.get_steppers():
                if s.get_name() == want:
                    step_dist = s.get_step_dist()
                    break
            if step_dist is None or mres is None:
                return None
            # Klipper converts a velocity to TPWMTHRS with the distance
            # of a 1/256 microstep, so invert exactly that.
            step_dist_256 = step_dist * (256 >> mres) / 256.0
            return freq * step_dist_256 / tp
        except Exception:
            return None

    def _spreadcycle_active(self, axis, velocity=None):
        """Do the chopper registers act at the tested speed?

        Not simply 'en_pwm_mode == 0': a config can enable StealthChop
        while confining it to near-standstill via TPWMTHRS, and then
        every test move still runs in SpreadCycle — tuning is perfectly
        valid there. Judging by en_pwm_mode alone would refuse such a
        machine for no reason.
        """
        en_pwm = self._field(axis, 'en_pwm_mode')       # 5160/2130/2240
        if en_pwm is None:
            # TMC2209: Klipper's field table spells this all-lowercase.
            for name in ('en_spreadcycle', 'en_spreadCycle'):
                v = self._field(axis, name)
                if v is not None:
                    return v == 1
            return None
        if en_pwm == 0:
            return True                     # SpreadCycle at every speed
        if velocity is None:
            return False
        limit = self._stealthchop_max_velocity(axis)
        if limit is None:
            return False
        return velocity > limit

    # ─── Motion + lost-step verdict ─────────────────────────────────

    def _move(self, axis, pos, speed):
        self.gcode.run_script_from_command(
            "G90\nG1 %s%.3f F%.0f" % (axis, pos, max(60.0, speed * 60.0)))

    def _snapshot_limits(self):
        """Remember the LIVE motion limits so they can be put back
        exactly — reading printer.cfg instead would silently discard a
        SET_VELOCITY_LIMIT the user (or a macro) had active."""
        snap = {}
        try:
            st = self.printer.lookup_object('toolhead').get_status(
                self.reactor.monotonic())
            for k in ('max_velocity', 'max_accel',
                      'square_corner_velocity', 'minimum_cruise_ratio'):
                if st.get(k) is not None:
                    snap[k] = float(st[k])
        except Exception:
            pass
        return snap

    def _set_limits(self, velocity=None, accel=None):
        """G1 alone does NOT apply an acceleration — the toolhead uses
        its current limit. Without this the whole escalation would only
        change the move DISTANCE and every candidate would see the same
        acceleration."""
        parts = []
        if velocity is not None:
            parts.append("VELOCITY=%.2f" % velocity)
        if accel is not None:
            parts.append("ACCEL=%.2f" % accel)
            # Klipper derives its own decel limit from this; pin it too
            # so the reversal is as hard as the acceleration.
            parts.append("ACCEL_TO_DECEL=%.2f" % accel)
        if parts:
            self.gcode.run_script_from_command(
                "SET_VELOCITY_LIMIT " + " ".join(parts))

    def _restore_limits(self):
        """Put back exactly what was live before the run.

        ACCEL_TO_DECEL is Klipper's legacy spelling of
        minimum_cruise_ratio: setting it pins the ratio to 0, and
        without restoring it the printer keeps that motion behaviour
        after the test ends.
        """
        snap = self._limits_snapshot or {}
        parts = []
        if 'max_velocity' in snap:
            parts.append("VELOCITY=%.2f" % snap['max_velocity'])
        if 'max_accel' in snap:
            parts.append("ACCEL=%.2f" % snap['max_accel'])
        if 'square_corner_velocity' in snap:
            parts.append("SQUARE_CORNER_VELOCITY=%.3f"
                         % snap['square_corner_velocity'])
        if 'minimum_cruise_ratio' in snap:
            parts.append("MINIMUM_CRUISE_RATIO=%.4f"
                         % snap['minimum_cruise_ratio'])
        if parts:
            self.gcode.run_script_from_command(
                "SET_VELOCITY_LIMIT " + " ".join(parts))

    def _sample_axes(self, axis):
        """Which steppers can lose steps for a move on `axis`.

        On CoreXY both motors drive every X or Y move, so both have to
        be watched and re-homed.
        """
        return ('X', 'Y') if self.structure == 'corexy' else (axis,)

    def _prepare(self, axis):
        """Home what is missing, then lift Z out of the way.

        Without the lift the nozzle would be dragged across the bed for
        the whole sweep — homing Z leaves it at the endstop position.
        """
        kin = self.printer.lookup_object('toolhead').get_kinematics()
        homed = kin.get_status(self.reactor.monotonic())['homed_axes']
        if 'z' not in homed:
            self.gcode.run_script_from_command("G28")
        self.gcode.run_script_from_command(
            "G90\nG1 Z%.3f F%.0f\nM400"
            % (self.z_pos, max(60.0, self.travel_speed * 60.0)))
        # The first measurement compares against the position right
        # after a home, so X/Y are always homed here — also when they
        # already were, because the toolhead may be anywhere.
        self._rehome(axis)

    def _rehome(self, axis):
        """Re-home UNCONDITIONALLY.

        This is what makes the lost-step comparison work at all:
        get_mcu_position() is a commanded step counter, so it only
        reveals a physical loss once a homing move re-anchors the
        commanded-to-physical mapping. A conditional 'home only what is
        unhomed' would be a no-op after the first call and every
        candidate would silently score zero lost steps.
        """
        self.gcode.run_script_from_command(
            "G28 %s\nM400" % " ".join(self._sample_axes(axis)))

    def _homes_to_max(self, axis):
        """True when `axis` (X/Y) homes at its high end.

        Read from the stepper config: homing_positive_dir when given,
        otherwise whichever end position_endstop is nearer to (Klipper's
        own default rule)."""
        if axis.upper() not in ('X', 'Y'):
            return False
        try:
            cfg = self.printer.lookup_object('configfile')
            st = cfg.get_status(self.reactor.monotonic())['settings'][
                'stepper_' + axis.lower()]
            pd = st.get('homing_positive_dir')
            if pd is not None:
                return str(pd).strip().lower() in ('true', '1', 'yes')
            lo, hi = float(st['position_min']), float(st['position_max'])
            return float(st['position_endstop']) > (lo + hi) / 2.0
        except Exception:
            return False

    def _park_pos(self, axis):
        """(other axis, position) to park at during a test on `axis`,
        or None. Own park_x / park_y first, then [speed_test]'s."""
        other = 'Y' if axis == 'X' else 'X'
        attr = 'park_' + other.lower()
        pos = getattr(self, attr, None)
        if pos is None:
            st = self.printer.lookup_object('speed_test', None)
            pos = getattr(st, attr, None) if st is not None else None
        if pos is None:
            return None
        try:
            lo, hi, _, _ = self._get_axis_bounds(other)
            pos = min(max(float(pos), lo), hi)
        except Exception:
            return None
        return other, pos

    def _park_other(self, axis):
        park = self._park_pos(axis)
        if park is not None:
            self._move(park[0], park[1], self.travel_speed)

    def _stress_pattern(self, axis, velocity, accel, repeat):
        """Out-and-back that just reaches V, plus a short cruise. Same
        shape speed_test uses for its fast stall probe: it is the
        cheapest move that still loads the motor fully, and it grinds
        the least when the motor does stall."""
        lo, hi, mid, rng = self._get_axis_bounds(axis)
        dist = (velocity ** 2) / accel + 0.05 * velocity
        dist = max(5.0, min(dist, rng))
        # `start` is the end of the move nearest the endstop; the offset
        # is measured from the end the axis HOMES to, so the re-home
        # after every candidate is a short trip.
        offset = (min(self.start_offset, 0.8 * rng)
                  if self.start_offset > 0 else 0.20 * rng)
        if self._homes_to_max(axis):
            start = hi - offset
            end = start - dist
            if end < lo:
                end = lo
                start = min(hi, end + dist)
        else:
            start = lo + offset
            end = start + dist
            if end > hi:
                end = hi
                start = max(lo, end - dist)
        self._park_other(axis)
        self._move(axis, start, self.travel_speed)
        for _ in range(repeat):
            self._move(axis, end, velocity)
            self._move(axis, start, velocity)
        self.gcode.run_script_from_command("M400")
        return start, end

    def _measure(self, gcmd, axis, velocity, accel, repeat):
        """Run the stress move and report lost steps across a re-home.

        Both samples are taken at the same reference: `before` right
        after the previous home, `after` right after this one.
        """
        sample = self._sample_axes(axis)
        before = {a: self._read_mcu_pos(a) for a in sample}
        t_start = self.reactor.monotonic()
        geo = None
        try:
            self._set_limits(velocity=velocity * 1.2, accel=accel)
            geo = self._stress_pattern(axis, velocity, accel, repeat)
        finally:
            self.gcode.run_script_from_command("M400")
            # Home and travel must not run at the test acceleration.
            self._restore_limits()
        t_moved = self.reactor.monotonic()
        try:
            self._rehome(axis)
        except Exception as e:
            self._log('HOME_FAIL', axis=axis, accel="%.0f" % accel,
                      error=repr(str(e)))
            raise
        worst = 0
        got_any = False
        detail = []
        for a in sample:
            after = self._read_mcu_pos(a)
            if after is None or before.get(a) is None:
                detail.append("%s:n/a" % a)
                continue
            got_any = True
            diff = abs(after - before[a])
            thresh = self.max_missed * self._microsteps(a)
            detail.append("%s:%+d" % (a, after - before[a]))
            if diff > thresh:
                worst = max(worst, diff)
        self._log(
            'MEASURE', axis=axis, v="%.0f" % velocity, a="%.0f" % accel,
            reps=repeat,
            move=("%.1f>%.1f" % geo) if geo else "n/a",
            usteps_after_rehome=",".join(detail),
            limit="%.0f" % (self.max_missed * self._microsteps(axis)),
            lost=worst, move_s="%.1f" % (t_moved - t_start),
            home_s="%.1f" % (self.reactor.monotonic() - t_moved))
        if not got_any:
            raise self.gcode.error(
                "chopper_test: cannot read the stepper MCU position for "
                "%s." % ", ".join(sample))
        return {'lost': worst, 'failed': worst > 0}

    # ─── Sweep ──────────────────────────────────────────────────────

    def _apply(self, axis, combo):
        """Write a register combination, DECREASES FIRST.

        The fields are written one at a time, so the driver briefly
        sees mixed old/new values. Both the old and the new combination
        satisfy hstrt+hend <= 18, but a naive order does not: going
        from (hstrt 5, hend 2) to (hstrt 4, hend 15) passes through
        (5, 15) = 20 and violates the datasheet rule mid-flight.
        Writing every decreasing field before any increasing one keeps
        each intermediate sum at or below max(old_sum, new_sum), so no
        illegal state can exist even for a moment. The same ordering
        protects the toff/tbl pair.
        """
        self._log('APPLY', combo=self._fmt_combo(combo))
        for stepper, tmc in self._tuned_steppers(axis):
            current = self._read_fields(tmc, combo)
            items = sorted(
                combo.items(),
                key=lambda kv: 0 if (current.get(kv[0]) is not None
                                     and kv[1] <= current[kv[0]]) else 1)
            for name, value in items:
                if current.get(name) == value:
                    continue      # already there — skip the round trip
                self._set_stepper_field(stepper, name, value)
        self.gcode.run_script_from_command("G4 P%d" % SETTLE_DELAY_MS)

    # ─── Re-sync of dual-motor belts after lost steps ───────────────

    def _resync_plan(self):
        """(gcode or None, description) — same choice as speed_test."""
        cfg = self.skip_resync_gcode
        if cfg.lower() in ('none', 'off', 'false', '0'):
            return None, 'disabled'
        if cfg.lower() != 'auto':
            return cfg, 'configured G-code'
        axes = []
        try:
            kin = self.printer.lookup_object('toolhead').get_kinematics()
            names = [st.get_name() for st in kin.get_steppers()]
        except Exception:
            names = []
        for ax in ('x', 'y'):
            base = 'stepper_' + ax
            if sum(1 for n in names if n.startswith(base)) > 1:
                axes.append(ax.upper())
        if not axes:
            return None, 'single-motor axes'
        if self.printer.lookup_object('motors_sync', None) is not None:
            return 'SYNC_MOTORS', 'accelerometer sync (motors_sync)'
        if self.printer.lookup_object('motor_sync', None) is not None:
            return ('\n'.join('MOTOR_SYNC AXIS=%s' % a for a in axes),
                    'StallGuard sync (motor_sync)')
        return None, 'unavailable'

    def _resync_after_loss(self, gcmd, axis):
        """Lost steps can leave the two motors of a belt out of phase,
        and homing does not fix that — every later candidate would then
        be tested on a weaker pair. Runs with the candidate still on the
        drivers put aside: the sync is done on the settings the run
        started with. Re-homes afterwards (the sync moves the toolhead).
        """
        script, how = self._resync_plan()
        if script is None:
            if how == 'unavailable' and not self._resync_noted:
                self._resync_noted = True
                gcmd.respond_info(
                    "    NOTE: lost steps on a dual-motor axis, but no "
                    "motor sync is available (no [motors_sync] or "
                    "[motor_sync]). The motors of a belt pair may now "
                    "be out of sync, so later candidates can read low.")
            return
        gcmd.respond_info("    lost steps → re-syncing motors (%s)" % how)
        self._log('RESYNC', how=repr(how))
        live = {}
        for stepper, tmc in self._tuned_steppers(axis):
            live[stepper] = self._read_fields(tmc, self._originals.get(
                stepper, {}))
        self._restore_originals(axis)
        try:
            self.gcode.run_script_from_command(script)
            self.gcode.run_script_from_command("M400")
        except self.gcode.error as e:
            gcmd.respond_info(
                "    WARNING: re-sync failed (%s) — continuing; later "
                "candidates may read low." % (e,))
            self._log('RESYNC_FAIL', error=repr(str(e)))
        self._rehome(axis)
        # Back to the candidate that was being tested.
        for stepper, tmc in self._tuned_steppers(axis):
            want = live.get(stepper) or {}
            cur = self._read_fields(tmc, want)
            for name, value in sorted(
                    want.items(),
                    key=lambda kv: 0 if kv[1] <= cur.get(kv[0], kv[1])
                    else 1):
                if cur.get(name) != value:
                    self._set_stepper_field(stepper, name, value)
        self.gcode.run_script_from_command("G4 P%d" % SETTLE_DELAY_MS)

    def _restore_originals(self, axis):
        """Every tuned driver back to what it had when the run began,
        decreases first (see _apply)."""
        for stepper, tmc in self._tuned_steppers(axis):
            orig = self._originals.get(stepper) or {}
            cur = self._read_fields(tmc, orig)
            for name, value in sorted(
                    orig.items(),
                    key=lambda kv: 0 if kv[1] <= cur.get(kv[0], kv[1])
                    else 1):
                if cur.get(name) != value:
                    self._set_stepper_field(stepper, name, value)
        self.gcode.run_script_from_command("G4 P%d" % SETTLE_DELAY_MS)

    def _measure_scored(self, gcmd, axis, velocity, accel, repeat):
        """_measure plus the one retry a failed candidate gets ('one
        noisy fail must not condemn a setting'), with a motor re-sync
        after every measurement that lost steps. The first result is
        kept in 'first_lost' so a retried pass is never shown as a
        clean one."""
        r = self._measure(gcmd, axis, velocity, accel, repeat)
        if not r['failed']:
            r['first_lost'] = 0
            self._log('VERDICT', result='OK')
            return r
        first = r['lost']
        self._resync_after_loss(gcmd, axis)
        self._log('RETRY', first_lost=first)
        r = self._measure(gcmd, axis, velocity, accel, repeat)
        if r['failed']:
            self._resync_after_loss(gcmd, axis)
        r['first_lost'] = first
        self._log('VERDICT',
                  result='FAIL' if r['failed'] else 'OK_ON_RETRY',
                  lost=r['lost'], first_lost=first)
        return r

    @staticmethod
    def _verdict(r):
        if r['failed']:
            return "lost %d µsteps (and %d on the first try)" % (
                r['lost'], r['first_lost'])
        if r.get('first_lost'):
            return "OK on retry (lost %d µsteps the first time)" % (
                r['first_lost'])
        return "OK"

    def _score_candidates(self, gcmd, axis, field, values, base,
                          velocity, accel, repeat, results):
        """Measure every candidate value of one field at `accel`.

        If more than one survives, the acceleration is raised and only
        the survivors are re-measured — that turns a pass/fail probe
        into a ranking without a full bisection per candidate, and the
        acceleration auto-calibrates to where the motor actually
        discriminates.

        Returns (best_value, accel_used).
        """
        candidates = list(values)
        cur_accel = accel
        for escalation in range(ESCALATIONS):
            survivors = []
            for v in candidates:
                combo = dict(base)
                combo[field] = v
                self._apply(axis, combo)
                r = self._measure_scored(gcmd, axis, velocity, cur_accel,
                                         repeat)
                results.append({
                    'field': field, 'value': v, 'accel': cur_accel,
                    'lost': r['lost'], 'failed': r['failed'],
                    'first_lost': r['first_lost'], 'base': dict(base),
                })
                gcmd.respond_info(
                    "    %s=%-3s @ %.0f mm/s²  →  %s"
                    % (field, v, cur_accel, self._verdict(r)))
                if not r['failed']:
                    survivors.append(v)
            if len(survivors) == 1:
                return survivors[0], cur_accel
            if not survivors:
                # Nothing held: take the value that lost the least, and
                # hand the next pass the last acceleration at which
                # something still held (if there was one).
                rows = [x for x in results
                        if x['field'] == field and x['accel'] == cur_accel]
                best = min(rows, key=lambda x: x['lost'])
                return best['value'], (cur_accel / ESCALATION_FACTOR
                                       if escalation else cur_accel)
            candidates = survivors
            if escalation == ESCALATIONS - 1:
                break
            # Several held — raise the bar and let them fight it out.
            cur_accel *= ESCALATION_FACTOR
            gcmd.respond_info(
                "    %d values held — raising to %.0f mm/s² to separate "
                "them" % (len(survivors), cur_accel))
        # Still tied after the last round: this pass could not separate
        # the values. Never invent a winner from list order — keep what
        # the printer already ran unless it is not among the survivors.
        # cur_accel is the last acceleration that was actually tested;
        # the next pass starts there, not above it.
        incumbent = base.get(field)
        keep = incumbent if incumbent in candidates else candidates[0]
        gcmd.respond_info(
            "    inconclusive — %d values still hold at %.0f mm/s² (the "
            "highest tested); keeping %s=%s"
            % (len(candidates), cur_accel, field, keep))
        return keep, cur_accel

    # ─── Thermal comparison ─────────────────────────────────────────

    def _read_motor_temp(self, axis):
        """Motor temperature, best source available.

        1. a Klipper temperature sensor the user mounted on the motor
        2. the driver's own sensor (TMC2240 reports one via get_status)
        Returns (value, source) or (None, None).
        """
        if self.motor_sensor:
            try:
                obj = self.printer.lookup_object(self.motor_sensor, None)
                if obj is not None:
                    now = self.reactor.monotonic()
                    if hasattr(obj, 'get_temp'):
                        res = obj.get_temp(now)
                        val = res[0] if isinstance(res, tuple) else res
                        return float(val), self.motor_sensor
                    if hasattr(obj, 'last_temp'):
                        return float(obj.last_temp), self.motor_sensor
            except Exception:
                pass
        info = self._lookup_tmc(axis)
        if info:
            try:
                st = info[1].get_status(self.reactor.monotonic())
                t = st.get('temperature')
                if t is not None:
                    return float(t), "%s driver sensor" % info[0]
            except Exception:
                pass
        return None, None

    def _thermal_run(self, gcmd, axis, combo, velocity, accel, seconds):
        """Move continuously for `seconds` and report the temperature
        rise. Returns dict or None when no sensor is available."""
        t0, source = self._read_motor_temp(axis)
        if t0 is None:
            return None
        if combo is None:
            self._restore_originals(axis)   # each driver's own values
        else:
            self._apply(axis, combo)
        lo, hi, mid, rng = self._get_axis_bounds(axis)
        dist = min(rng, max(20.0, (velocity ** 2) / accel + 0.2 * velocity))
        start = mid - dist / 2.0
        end = mid + dist / 2.0
        self._park_other(axis)
        self._move(axis, start, self.travel_speed)
        self.gcode.run_script_from_command("M400")
        deadline = self.reactor.monotonic() + seconds
        peak = t0
        try:
            # Both thermal runs must move identically, so pin the limits
            # here too — otherwise the comparison measures the printer's
            # current accel setting, not the chopper settings.
            self._set_limits(velocity=velocity * 1.2, accel=accel)
            while self.reactor.monotonic() < deadline:
                self._move(axis, end, velocity)
                self._move(axis, start, velocity)
                self.gcode.run_script_from_command("M400")
                t, _ = self._read_motor_temp(axis)
                if t is not None and t > peak:
                    peak = t
        finally:
            self._restore_limits()
        t1, _ = self._read_motor_temp(axis)
        return {'start': t0, 'end': t1 if t1 is not None else peak,
                'peak': peak, 'rise': (t1 if t1 is not None else peak) - t0,
                'source': source}

    # ─── Commands ───────────────────────────────────────────────────

    def cmd_STATUS(self, gcmd):
        axis = gcmd.get('AXIS', self.default_axis).upper()
        info = self._lookup_tmc(axis)
        lines = ["%s v%s — axis %s" % (MODULE_NAME, MODULE_VERSION, axis)]
        if not info:
            lines.append("  no TMC driver found for %s"
                         % self._stepper_name(axis))
            gcmd.respond_info("\n".join(lines))
            return
        drv = info[0]
        lines.append("  driver: %s on %s" % (drv, self._stepper_name(axis)))
        lines.append("  CHOPPER_TUNE writes candidates to: %s"
                     % ", ".join(n for n, _ in self._tuned_steppers(axis)))
        park = self._park_pos(axis)
        if park is not None:
            lines.append("  %s parked at %.1f during the test"
                         % (park[0], park[1]))
        lines.append("  re-sync after lost steps: %s"
                     % self._resync_plan()[1])
        limit = self._stealthchop_max_velocity(axis)
        en_pwm = self._field(axis, 'en_pwm_mode')
        if en_pwm == 0 or (en_pwm is None
                           and self._spreadcycle_active(axis) is True):
            mode = "SpreadCycle at every speed ✓ tunable"
        elif limit in (None, float('inf')):
            mode = ("StealthChop at every speed ✗ — the chopper "
                    "registers do nothing here")
        else:
            mode = ("StealthChop below %.2f mm/s, SpreadCycle above ✓ "
                    "tunable (test moves are far above that)" % limit)
        lines.append("  chopper mode: %s" % mode)
        vals = []
        for f in ('toff', 'tbl', 'hstrt', 'hend', 'tpfd'):
            v = self._field(axis, f)
            if v is not None:
                vals.append("%s=%s" % (f, v))
        lines.append("  registers: %s" % (" ".join(vals) or "n/a"))
        h, e = self._field(axis, 'hstrt'), self._field(axis, 'hend')
        if h is not None and e is not None:
            lines.append("  hstrt+hend = %d (datasheet limit %d)"
                         % (h + e, HSTRT_HEND_SUM_MAX))
        t, src = self._read_motor_temp(axis)
        lines.append("  motor temperature: %s"
                     % ("%.1f °C (%s)" % (t, src) if t is not None
                        else "no sensor — set motor_sensor: in "
                             "[chopper_test] for the thermal stage"))
        lines.append("  tpfd sweep: %s"
                     % ("available" if drv in TPFD_DRIVERS
                        else "not on this driver"))
        gcmd.respond_info("\n".join(lines))

    def cmd_CHOPPER_TUNE(self, gcmd):
        axis = gcmd.get('AXIS', self.default_axis).upper()
        if axis not in ('X', 'Y'):
            raise self.gcode.error("chopper_test: AXIS must be X or Y")
        velocity = gcmd.get_float('VELOCITY', 200.0, above=0.)
        accel = gcmd.get_float('ACCEL', 0.0, minval=0.)
        repeat = gcmd.get_int('REPEAT', 8, minval=2, maxval=40)
        do_tpfd = gcmd.get_int('TPFD', 0, minval=0, maxval=1)
        do_thermal = gcmd.get_int('THERMAL', 1, minval=0, maxval=1)
        thermal_s = gcmd.get_float('THERMAL_SECONDS', self.thermal_seconds,
                                   minval=10.)

        info = self._lookup_tmc(axis)
        if not info:
            raise self.gcode.error(
                "chopper_test: no TMC driver found for %s."
                % self._stepper_name(axis))
        drv = info[0]
        if self._spreadcycle_active(axis, velocity) is False:
            limit = self._stealthchop_max_velocity(axis)
            raise self.gcode.error(
                "chopper_test: %s runs in StealthChop up to %s — the "
                "chopper registers have no effect there, so the test "
                "speed of %.0f mm/s would measure nothing. Either "
                "raise VELOCITY above that, or set "
                "'stealthchop_threshold: 0' in [%s %s] (which keeps "
                "StealthChop for standstill only) and "
                "FIRMWARE_RESTART."
                % (axis,
                   "every speed" if limit in (None, float('inf'))
                   else "%.1f mm/s" % limit,
                   velocity, drv, self._stepper_name(axis)))

        # Original register state — restored no matter how we leave.
        fields = ['toff', 'tbl', 'hstrt', 'hend']
        if drv in TPFD_DRIVERS:
            fields.append('tpfd')
        original = {}
        for f in fields:
            v = self._field(axis, f)
            if v is not None:
                original[f] = v
        if 'toff' not in original:
            raise self.gcode.error(
                "chopper_test: cannot read the chopper registers of %s."
                % drv)

        # Every motor that drives this move gets the candidates, and
        # each one is handed back exactly what it had.
        steppers = self._tuned_steppers(axis)
        self._originals = {n: self._read_fields(tmc, fields)
                           for n, tmc in steppers}
        self._resync_noted = False
        differing = [n for n, _ in steppers
                     if any(self._originals[n].get(f) != v
                            for f, v in original.items())]

        if accel <= 0:
            _, _, _, rng = self._get_axis_bounds(axis)
            accel = max(1000.0, (velocity ** 2) / max(5.0, rng * 0.25))
            gcmd.respond_info(
                "  no ACCEL given — starting the screening at %.0f mm/s² "
                "(derived from the travel); it rises automatically until "
                "the settings separate." % accel)

        gcmd.respond_info(
            "──── CHOPPER_TUNE %s ────\n"
            "  driver %s | start: %s\n"
            "  candidates are written to: %s\n"
            "  stress: %.0f mm/s, %d out-and-back moves per candidate\n"
            "  scoring: lost steps (torque)%s\n"
            "  every register is restored when the run ends."
            % (axis, drv,
               " ".join("%s=%s" % kv for kv in sorted(original.items())),
               ", ".join(n for n, _ in steppers),
               velocity, repeat,
               " + motor temperature" if do_thermal else ""))
        if differing:
            gcmd.respond_info(
                "  NOTE: %s start(s) with different chopper values than "
                "%s. They get the same candidates during the test and "
                "their own values back afterwards."
                % (", ".join(differing), self._stepper_name(axis)))
        park = self._park_pos(axis)
        if park is not None:
            gcmd.respond_info("  %s is parked at %.1f during the test."
                              % (park[0], park[1]))

        results = []
        best = dict(original)
        timestamp = time.strftime('%Y-%m-%d_%H-%M-%S')
        thermal = {}
        self._limits_snapshot = self._snapshot_limits()
        self._log_open(timestamp)
        if self._runlog_path:
            gcmd.respond_info(
                "  run log (written live): %s" % self._runlog_path)
        self._log('START', version=MODULE_VERSION, axis=axis, driver=drv,
                  velocity="%.0f" % velocity, accel="%.0f" % accel,
                  repeat=repeat, tpfd=do_tpfd, thermal=do_thermal,
                  structure=self.structure)
        for n, _tmc in steppers:
            self._log('STEPPER', name=n,
                      start=self._fmt_combo(self._originals.get(n, {})))
        lo_, hi_, _mid, _rng = self._get_axis_bounds(axis)
        self._log('SETUP', bounds="%.1f..%.1f" % (lo_, hi_),
                  homes_to_max=self._homes_to_max(axis),
                  start_offset="%.1f" % self.start_offset,
                  park=("%s%.1f" % park) if park else "none",
                  resync=repr(self._resync_plan()[1]),
                  max_missed=self.max_missed,
                  limits=self._fmt_combo(self._limits_snapshot or {}))
        failed_with = None
        try:
            self._prepare(axis)

            # Pass 1 — TOFF (dominates the chopper frequency)
            self._log('PASS', name='Pass 1: TOFF', accel="%.0f" % accel)
            gcmd.respond_info("\n>>> Pass 1: TOFF (chopper off time) <<<")
            v, accel = self._score_candidates(
                gcmd, axis, 'toff', TOFF_VALUES, best, velocity, accel,
                repeat, results)
            best['toff'] = v
            self._log('PICK', field='toff', value=v, accel="%.0f" % accel)
            gcmd.respond_info("  → toff = %s" % v)

            # Pass 2 — TBL (comparator blank time)
            self._log('PASS', name='Pass 2: TBL', accel="%.0f" % accel)
            gcmd.respond_info("\n>>> Pass 2: TBL (blank time) <<<")
            v, accel = self._score_candidates(
                gcmd, axis, 'tbl', TBL_VALUES, best, velocity, accel,
                repeat, results)
            best['tbl'] = v
            self._log('PICK', field='tbl', value=v, accel="%.0f" % accel)
            gcmd.respond_info("  → tbl = %s" % v)

            # Pass 3 — hysteresis. hstrt and hend are swept together
            # because only their SUM is constrained and they shape the
            # same regulation window.
            self._log('PASS', name='Pass 3: HSTRT / HEND',
                      accel="%.0f" % accel)
            gcmd.respond_info(
                "\n>>> Pass 3: HSTRT / HEND (hysteresis) <<<")
            pairs = []
            for hs in range(0, HSTRT_MAX + 1, 2):
                for he in range(0, HEND_MAX + 1, 3):
                    if hs + he <= HSTRT_HEND_SUM_MAX:
                        pairs.append((hs, he))
            hyst_results = []
            cur_accel = accel
            survivors = pairs
            for rnd in range(ESCALATIONS):
                held = []
                for hs, he in survivors:
                    combo = dict(best)
                    combo['hstrt'] = hs
                    combo['hend'] = he
                    self._apply(axis, combo)
                    r = self._measure_scored(gcmd, axis, velocity,
                                             cur_accel, repeat)
                    hyst_results.append({
                        'field': 'hstrt/hend', 'value': '%d/%d' % (hs, he),
                        'accel': cur_accel, 'lost': r['lost'],
                        'failed': r['failed'],
                        'first_lost': r['first_lost'],
                        'base': dict(best)})
                    gcmd.respond_info(
                        "    hstrt=%d hend=%-2d @ %.0f mm/s²  →  %s"
                        % (hs, he, cur_accel, self._verdict(r)))
                    if not r['failed']:
                        held.append((hs, he))
                if len(held) == 1:
                    survivors = held
                    break
                if not held:
                    rows = [x for x in hyst_results
                            if x['accel'] == cur_accel]
                    bestrow = min(rows, key=lambda x: x['lost'])
                    hs, he = [int(x) for x in bestrow['value'].split('/')]
                    survivors = [(hs, he)]
                    if rnd:
                        cur_accel /= ESCALATION_FACTOR
                    break
                survivors = held
                if rnd == ESCALATIONS - 1:
                    break
                cur_accel *= ESCALATION_FACTOR
                gcmd.respond_info(
                    "    %d pairs held — raising to %.0f mm/s²"
                    % (len(held), cur_accel))
            results.extend(hyst_results)
            accel = cur_accel
            if len(survivors) > 1:
                # Inconclusive — keep the incumbent pair rather than
                # picking whatever happens to be first in the grid.
                incumbent = (original.get('hstrt'), original.get('hend'))
                if incumbent in survivors:
                    survivors = [incumbent]
                gcmd.respond_info(
                    "    inconclusive — %d pairs still hold; keeping "
                    "hstrt=%d hend=%d"
                    % (len(survivors), survivors[0][0], survivors[0][1]))
            best['hstrt'], best['hend'] = survivors[0]
            self._log('PICK', field='hstrt/hend',
                      value="%d/%d" % survivors[0], accel="%.0f" % accel)
            gcmd.respond_info("  → hstrt = %d, hend = %d"
                              % (best['hstrt'], best['hend']))

            # Pass 4 — TPFD (5160/2240 only, opt-in)
            if do_tpfd and drv in TPFD_DRIVERS and 'tpfd' in original:
                self._log('PASS', name='Pass 4: TPFD', accel="%.0f" % accel)
                gcmd.respond_info(
                    "\n>>> Pass 4: TPFD (passive fast decay) <<<")
                v, accel = self._score_candidates(
                    gcmd, axis, 'tpfd', TPFD_VALUES, best, velocity,
                    accel, repeat, results)
                best['tpfd'] = v
                self._log('PICK', field='tpfd', value=v,
                          accel="%.0f" % accel)
                gcmd.respond_info("  → tpfd = %s" % v)

            # Thermal shootout: old vs new, same move, same duration.
            if do_thermal:
                self._log('PASS', name='Thermal',
                          accel="%.0f" % (accel / 1.3))
                gcmd.respond_info(
                    "\n>>> Thermal comparison (%.0f s each) <<<"
                    % thermal_s)
                probe = self._thermal_run(gcmd, axis, None, velocity,
                                          accel / 1.3, thermal_s)
                if probe is None:
                    gcmd.respond_info(
                        "  skipped — no motor temperature source. Add a "
                        "sensor and set 'motor_sensor:' in "
                        "[chopper_test], or use a TMC2240.")
                else:
                    thermal['original'] = probe
                    gcmd.respond_info(
                        "  original: %.1f → %.1f °C (+%.1f) [%s]"
                        % (probe['start'], probe['end'], probe['rise'],
                           probe['source']))
                    gcmd.respond_info("  ... letting the motor cool 60 s")
                    self.gcode.run_script_from_command("G4 P60000")
                    probe2 = self._thermal_run(gcmd, axis, best, velocity,
                                               accel / 1.3, thermal_s)
                    if probe2:
                        thermal['tuned'] = probe2
                        gcmd.respond_info(
                            "  tuned:    %.1f → %.1f °C (+%.1f)"
                            % (probe2['start'], probe2['end'],
                               probe2['rise']))
        except Exception as e:
            # Keep what was measured: the run log is already on disk,
            # and the CSV gets everything up to here.
            failed_with = e
            self._log('ERROR', error=repr(str(e)), measured=len(results))
            self._write_csv(gcmd, timestamp, axis, drv, original, best,
                            results, thermal, velocity, repeat,
                            aborted=str(e))
            raise
        finally:
            # Always hand every driver back exactly as we found it.
            try:
                self._restore_originals(axis)
                self._log('RESTORED', steppers=",".join(
                    n for n, _tmc in steppers))
            finally:
                if failed_with is not None:
                    self._log_close()

        changed = {k: v for k, v in best.items() if original.get(k) != v}
        gcmd.respond_info(
            "\n──── result ────\n"
            "  before: %s\n"
            "  after:  %s\n"
            "  changed: %s"
            % (" ".join("%s=%s" % kv for kv in sorted(original.items())),
               " ".join("%s=%s" % kv for kv in sorted(best.items())),
               " ".join("%s=%s" % kv for kv in sorted(changed.items()))
               or "nothing — your current settings already won"))
        if thermal.get('original') and thermal.get('tuned'):
            d = thermal['tuned']['rise'] - thermal['original']['rise']
            gcmd.respond_info(
                "  temperature rise: %+.1f °C vs the original settings"
                % d)
        if changed:
            gcmd.respond_info(
                "  The registers are back on their original values. To "
                "keep the result, add to %s:\n%s"
                % (" and ".join("[%s %s]" % (drv, n)
                                for n, _ in steppers),
                   "\n".join("    driver_%s: %s" % (k.upper(), v)
                             for k, v in sorted(changed.items()))))
        self._write_csv(gcmd, timestamp, axis, drv, original, best,
                        results, thermal, velocity, repeat)
        self._log('END', result=self._fmt_combo(best),
                  changed=self._fmt_combo(changed) or 'nothing')
        self._log_close()

    def _write_csv(self, gcmd, timestamp, axis, drv, original, best,
                   results, thermal, velocity, repeat, aborted=None):
        try:
            os.makedirs(self.output_dir, exist_ok=True)
            path = os.path.join(self.output_dir,
                                'chopper_%s.csv' % timestamp)
            with open(path, 'w', encoding='utf-8') as f:
                f.write("# Chopper Test v%s results\n" % MODULE_VERSION)
                f.write("# Plugin by Steven (Fragmon) — Crydteam\n")
                f.write("# sweep order adapted from MRX8024/"
                        "chopper-resonance-tuner\n")
                if aborted:
                    f.write("# ABORTED — partial results. Error: %s\n"
                            % aborted.replace("\n", " "))
                f.write("# axis: %s  driver: %s\n" % (axis, drv))
                f.write("# steppers tuned: %s\n" % ", ".join(
                    n for n, _ in self._tuned_steppers(axis)))
                f.write("# stress: %.0f mm/s, %d moves per candidate\n"
                        % (velocity, repeat))
                f.write("# before: %s\n" % " ".join(
                    "%s=%s" % kv for kv in sorted(original.items())))
                f.write("# after: %s\n" % " ".join(
                    "%s=%s" % kv for kv in sorted(best.items())))
                for key, t in sorted(thermal.items()):
                    f.write("# thermal_%s: %.1f -> %.1f C (+%.1f) via %s\n"
                            % (key, t['start'], t['end'], t['rise'],
                               t['source']))
                f.write("field,value,accel_mm_s2,lost_usteps,failed,"
                        "first_try_lost_usteps,base\n")
                for r in results:
                    f.write("%s,%s,%.0f,%d,%d,%d,%s\n" % (
                        r['field'], r['value'], r['accel'], r['lost'],
                        1 if r['failed'] else 0, r.get('first_lost', 0),
                        " ".join("%s=%s" % kv
                                 for kv in sorted(r['base'].items()))))
            gcmd.respond_info("  CSV saved: %s" % path)
        except Exception as e:
            gcmd.respond_info("  CSV write failed: %s" % e)


def load_config(config):
    return ChopperTest(config)
