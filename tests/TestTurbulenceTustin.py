# TestTurbulenceTustin.py
#
# Check that the Tustin turbulence model keeps the intensity requested by
# MIL-F-8785C in free flight, whatever the time step.
#
# Copyright (c) 2026 Alexander Kalmykov
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <http://www.gnu.org/licenses/>
#

import math
from JSBSim_utils import JSBSimTestCase, RunTest


class TestTurbulenceTustin(JSBSimTestCase):
    # sigma_w of severity 4 (the 1e-3 curve of MIL-F-8785C Fig. 7) at the initial
    # altitude, interpolated as in FGWinds. The aircraft only descends from there,
    # so 2*SIGMA_W bounds the whole run.
    SIGMA_W = 10.125  # ft/s
    FREE_FLIGHT_H_SL_FT = 2800  # chosen: the altitude SIGMA_W is read at
    FREE_FLIGHT_SEED = 1  # chosen
    # enum tType in FGWinds.h
    TT_NONE, TT_MILSPEC, TT_TUSTIN = 0, 3, 4

    def testFreeFlightIntensity(self):
        for turb_type in (self.TT_MILSPEC, self.TT_TUSTIN):
            for dt in (1/120, 1/1000):
                with self.subTest(turb_type=turb_type, dt=dt):
                    w_rms = self.runFreeFlight(turb_type, dt)
                    self.assertLess(w_rms, 2.0*self.SIGMA_W)

    # Measured over seeds 1 to 64 without the fix: the last Milspec draw sets the kick,
    # 0.6 to 26 times the cold peak, so one seed can miss it. A plain range, not picked
    # seeds, because each standard library draws another stream.
    MODEL_CHANGE_SEEDS = (1, 2, 3, 4)
    MODEL_CHANGE_H_SL_FT = 5000  # chosen: where the run-away was measured
    MODEL_CHANGE_DT = 1/240  # chosen: the time step of that measurement
    BEFORE_CHANGE_S = 60  # chosen: longer than the time constant L_w/V of about 12 s
    AFTER_CHANGE_S = 30  # measured: the run-away peaks about 7 s after the change
    # Measured over those 64 seeds: at most 2.2 with the fix, median 6.7 without it.
    PEAK_OVER_COLD = 3

    def testModelChange(self):
        for seed in self.MODEL_CHANGE_SEEDS:
            for before, after in (((self.TT_MILSPEC,), self.TT_TUSTIN),
                                  ((self.TT_MILSPEC, self.TT_NONE), self.TT_TUSTIN),
                                  ((self.TT_TUSTIN,), self.TT_MILSPEC)):
                with self.subTest(seed=seed, before=before, after=after):
                    cold_peak = self.peakTurbDown(seed, (self.TT_NONE,), after)
                    self.assertLess(self.peakTurbDown(seed, before, after), self.PEAK_OVER_COLD*cold_peak)

    def peakTurbDown(self, seed, before, after):
        fdm = self.startFreeFlight(seed, self.MODEL_CHANGE_H_SL_FT, self.MODEL_CHANGE_DT)
        for turb_type in before:
            fdm['atmosphere/turb-type'] = turb_type
            t_end = fdm.get_sim_time() + self.BEFORE_CHANGE_S
            while fdm.get_sim_time() < t_end:
                fdm.run()

        fdm['atmosphere/turb-type'] = after
        peak = 0.0
        t_end = fdm.get_sim_time() + self.AFTER_CHANGE_S
        while fdm.get_sim_time() < t_end:
            fdm.run()
            peak = max(peak, abs(fdm['atmosphere/turb-down-fps']))
        return peak

    def startFreeFlight(self, seed, h_sl_ft, dt):
        fdm = self.create_fdm()
        fdm['simulation/randomseed'] = seed
        fdm.load_model('c172p')
        fdm.set_dt(dt)

        # No trim and no autopilot: the aircraft descends by about 1000 ft/min.
        fdm['ic/h-sl-ft'] = h_sl_ft
        fdm['ic/vc-kts'] = 74
        fdm['ic/gamma-deg'] = 0
        fdm['propulsion/engine[0]/set-running'] = 1
        fdm.run_ic()
        fdm['fcs/throttle-cmd-norm'] = 0.7
        fdm['fcs/mixture-cmd-norm'] = 1.0

        fdm['atmosphere/turbulence/milspec/severity'] = 4
        return fdm

    def runFreeFlight(self, turb_type, dt):
        fdm = self.startFreeFlight(self.FREE_FLIGHT_SEED, self.FREE_FLIGHT_H_SL_FT, dt)
        fdm['atmosphere/turb-type'] = turb_type

        sum_w2 = 0.0
        n = 0
        while fdm.get_sim_time() < 60:
            fdm.run()
            sum_w2 += fdm['atmosphere/turb-down-fps']**2
            n += 1

        return math.sqrt(sum_w2/n)

    def testRunICWithAirspeed(self):
        TUSTIN = 4  # enum tType in FGWinds.h
        H_SL_FT = 2800  # chosen: as runFreeFlight
        VC_KTS = 74  # chosen: as runFreeFlight
        SEVERITY = 4  # "moderate" in the table of FGWinds.h
        FLIGHT_S = 1  # chosen: any flight draws turbulence at once

        fdm = self.create_fdm()
        fdm.load_model('c172p')
        fdm['ic/h-sl-ft'] = H_SL_FT
        fdm['ic/vc-kts'] = VC_KTS
        fdm['atmosphere/turbulence/milspec/severity'] = SEVERITY
        fdm['atmosphere/turb-type'] = TUSTIN

        fdm.run_ic()
        self.checkTurbulenceAndVelocityFinite(fdm)
        while fdm.get_sim_time() < FLIGHT_S:
            fdm.run()
        self.checkTurbulenceAndVelocityFinite(fdm)
        for axis in ('north', 'east', 'down'):
            self.assertNotEqual(fdm[f'atmosphere/turb-{axis}-fps'], 0.0)

    def checkTurbulenceAndVelocityFinite(self, fdm):
        for name in ('atmosphere/turb-north-fps', 'atmosphere/turb-east-fps',
                     'atmosphere/turb-down-fps', 'velocities/vt-fps'):
            self.assertTrue(math.isfinite(fdm[name]), name)

RunTest(TestTurbulenceTustin)
