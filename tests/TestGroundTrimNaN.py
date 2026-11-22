# TestGroundTrimNaN.py
#
# Check that a trim fails, rather than crashes, when its moment or a part of its
# initial conditions is not a number.
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
from jsbsim import TrimFailureError, TrimMode


class TestGroundTrimNaN(JSBSimTestCase):
    def trim_fails(self, name):
        fdm = self.create_fdm()
        fdm.load_model('c172x')  # its default IC stands on the ground
        fdm[name] = math.nan
        with self.assertRaises(TrimFailureError):
            fdm.do_trim(TrimMode.GROUND)

    # Unless the trim refuses the IC before any model runs, the IC cases but
    # the rate abort a debug build in a propeller table; the rate one trims.
    def test_nan_moment_fails_the_trim(self):
        # A NaN command makes the elevator, hence the pitching moment, NaN.
        self.trim_fails('fcs/elevator-cmd-norm')

    def test_nan_attitude_fails_the_trim(self):
        self.trim_fails('ic/theta-deg')

    def test_nan_speed_fails_the_trim(self):
        self.trim_fails('ic/u-fps')

    def test_nan_wind_fails_the_trim(self):
        self.trim_fails('ic/vw-mag-fps')  # moves only the airspeed

    def test_nan_terrain_fails_the_trim(self):
        self.trim_fails('ic/terrain-elevation-ft')  # moves only the AGL

    def test_nan_rate_fails_the_trim(self):
        # Zeroed in the trim's copy, the IC still holds it for the next run.
        self.trim_fails('ic/p-rad_sec')


RunTest(TestGroundTrimNaN)
