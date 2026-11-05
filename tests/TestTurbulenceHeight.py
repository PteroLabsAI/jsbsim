# TestTurbulenceHeight.py
#
# Check that the height above the ground picks the low altitude model of the Milspec
# turbulence, where MIL-F-8785C defines it, while the table of its high altitude
# model stays read above sea level, and that a jump of that height does not ring
# in the Tustin filters.
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

from JSBSim_utils import JSBSimTestCase, RunTest


class TestTurbulenceHeight(JSBSimTestCase):
    MILSPEC, TUSTIN = 3, 4  # enum tType in FGWinds.h
    SEVERITY = 6  # "severe" in the table of FGWinds.h
    CALM_FPS = 0  # no wind at 20 ft: the low altitude model has no intensity
    WIND_20FT_FPS = 33  # chosen: the wind of the run that found the bug
    LOW_AGL_FT = 30  # chosen: in the low altitude model
    TERRAIN_FT = 6000  # chosen: LOW_AGL_FT above it is in the table's range above sea level
    HIGH_MODEL_TERRAIN_FT = 3000  # chosen: leaves more than 2000 ft above it, in the high model
    # A host with no ground known reports 50 km: FDMGroundCallback::NoGround in PteroSim.
    NO_GROUND_AGL_FT = 164042
    SPEED_KTS = 74  # chosen: as TestTurbulenceTustin
    H_TOLERANCE_FT = 1  # chosen: RunIC does not move the aircraft, this only catches a lost terrain
    TURBULENCE_TOLERANCE_FPS = 1e-6  # measured: 0 between terrains at one altitude
    LOWER_H_SL_FT = 3030  # chosen: the same clipped height as TERRAIN_FT + LOW_AGL_FT, the table reads 21.1 against 23.4
    # Chosen: between the measured 1.11 of the two altitudes and the 1.000005 if the table reads the height.
    MIN_ALTITUDE_RATIO = 1.05

    # A plain range, not picked seeds: each standard library draws another stream.
    JUMP_SEEDS = (1, 2, 3, 4)
    JUMP_DT = 1/240  # chosen: the time step of the measurement
    JUMP_FROM_AGL_FT = 30  # chosen: the low height the jump is measured from
    JUMP_TO_AGL_FT = 130  # measured: the worst of the ordinary jumps without the limit
    BEFORE_JUMP_S = 20  # chosen: over 10 times L_u/V at the lower height, about 1.5 s
    AFTER_JUMP_S = 20  # chosen: over 4 times L_u/V at the upper height, about 4 s
    # Measured over seeds 1 to 128: with the limit at most 1.7 (no ground 1.3), without it 81 seeds above this (no ground 111).
    PEAK_OVER_NO_JUMP = 3

    def testLowOverHighTerrain(self):
        # MIL-F-8785C p. 60: low altitude is "ground level to approximately
        # 2,000 feet AGL". It follows the wind at 20 ft and ignores the severity.
        self.assertEqual(self.turbDown(self.TERRAIN_FT, self.LOW_AGL_FT, self.CALM_FPS), 0.0)
        self.assertGreater(self.turbDown(self.TERRAIN_FT, self.LOW_AGL_FT, self.WIND_20FT_FPS), 0.0)

    def testHighModelTableReadAboveSeaLevel(self):
        # Fig. 7 says only "ALTITUDE": the ground under the aircraft must not move the table.
        h_sl_ft = self.TERRAIN_FT + self.LOW_AGL_FT
        over_sea_level = self.turbDown(0, h_sl_ft, self.CALM_FPS)
        self.assertGreater(over_sea_level, 0.0)
        for terrain_ft in (self.HIGH_MODEL_TERRAIN_FT, h_sl_ft - self.NO_GROUND_AGL_FT):
            with self.subTest(terrain_ft=terrain_ft):
                self.assertAlmostEqual(self.turbDown(terrain_ft, h_sl_ft - terrain_ft, self.CALM_FPS),
                                       over_sea_level, delta=self.TURBULENCE_TOLERANCE_FPS)
        # A ratio, because the first sample goes with the draw of each standard library.
        self.assertGreater(over_sea_level/self.turbDown(0, self.LOWER_H_SL_FT, self.CALM_FPS), self.MIN_ALTITUDE_RATIO)

    def testSecondRunICTakesItsHeight(self):
        # A host runs RunIC at load and again at the real start: in calm the low model has no intensity.
        fdm = self.startFlight(0, self.TERRAIN_FT + self.LOW_AGL_FT, self.CALM_FPS, self.MILSPEC, 'ic/vt-kts')
        fdm['ic/terrain-elevation-ft'] = self.TERRAIN_FT
        fdm['ic/h-agl-ft'] = self.LOW_AGL_FT
        fdm.run_ic()
        fdm.run()
        self.assertEqual(fdm['atmosphere/turb-down-fps'], 0.0)

    def testHeightJumpSmoothed(self):
        # A no-ground answer is measured against the height it left: the limit keeps it out of the high model.
        for jump_to_agl_ft, no_jump_agl_ft in ((self.JUMP_TO_AGL_FT, self.JUMP_TO_AGL_FT),
                                               (self.NO_GROUND_AGL_FT, self.JUMP_FROM_AGL_FT)):
            for seed in self.JUMP_SEEDS:
                with self.subTest(jump_to_agl_ft=jump_to_agl_ft, seed=seed):
                    no_jump_peak = self.peakTurbEast(seed, no_jump_agl_ft, no_jump_agl_ft)
                    self.assertLess(self.peakTurbEast(seed, self.JUMP_FROM_AGL_FT, jump_to_agl_ft),
                                    self.PEAK_OVER_NO_JUMP*no_jump_peak)

    def peakTurbEast(self, seed, before_agl_ft, after_agl_ft):
        fdm = self.startFlight(self.TERRAIN_FT, before_agl_ft, self.WIND_20FT_FPS, self.TUSTIN, 'ic/vc-kts')
        fdm['atmosphere/randomseed'] = seed
        fdm.set_dt(self.JUMP_DT)
        self.flyAt(fdm, before_agl_ft, self.BEFORE_JUMP_S)
        return self.flyAt(fdm, after_agl_ft, self.AFTER_JUMP_S)

    def flyAt(self, fdm, h_agl_ft, seconds):
        peak = 0.0
        t_end = fdm.get_sim_time() + seconds
        while fdm.get_sim_time() < t_end:
            # No autopilot holds the height: the terrain follows the aircraft.
            fdm['ic/terrain-elevation-ft'] = fdm['position/h-sl-ft'] - h_agl_ft
            fdm.run()
            peak = max(peak, abs(fdm['atmosphere/turb-east-fps']))
        # The terrain was set before this step: the aircraft has moved one step from it.
        one_step_ft = abs(fdm['velocities/h-dot-fps'])*fdm.get_delta_t()
        self.assertAlmostEqual(fdm['position/h-agl-ft'], h_agl_ft, delta=self.H_TOLERANCE_FT + one_step_ft)
        return peak

    def turbDown(self, terrain_ft, h_agl_ft, wind_20ft_fps):
        # One true airspeed at every altitude: the first sample goes with the root of V.
        fdm = self.startFlight(terrain_ft, h_agl_ft, wind_20ft_fps, self.MILSPEC, 'ic/vt-kts')
        fdm.run()
        return abs(fdm['atmosphere/turb-down-fps'])

    def startFlight(self, terrain_ft, h_agl_ft, wind_20ft_fps, turb_type, speed_ic):
        fdm = self.create_fdm()
        fdm.load_model('c172p')
        fdm['ic/terrain-elevation-ft'] = terrain_ft
        fdm['ic/h-agl-ft'] = h_agl_ft
        fdm[speed_ic] = self.SPEED_KTS
        fdm.run_ic()

        self.assertAlmostEqual(fdm['position/h-sl-ft'], terrain_ft + h_agl_ft, delta=self.H_TOLERANCE_FT)

        fdm['atmosphere/turbulence/milspec/severity'] = self.SEVERITY
        fdm['atmosphere/turbulence/milspec/windspeed_at_20ft_AGL-fps'] = wind_20ft_fps
        fdm['atmosphere/turb-type'] = turb_type
        return fdm

RunTest(TestTurbulenceHeight)
