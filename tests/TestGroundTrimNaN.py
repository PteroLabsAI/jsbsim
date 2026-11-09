# TestGroundTrimNaN.py
#
# Check that a ground trim fails, rather than crashes, when its moment or its
# attitude is not a number.
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
    def ground_trim_fails(self, name):
        fdm = self.create_fdm()
        fdm.load_model('c172x')  # its default IC stands on the ground
        fdm[name] = math.nan
        with self.assertRaises(TrimFailureError):
            fdm.do_trim(TrimMode.GROUND)

    def test_nan_moment_fails_the_trim(self):
        # A NaN command makes the elevator, hence the pitching moment, NaN.
        self.ground_trim_fails('fcs/elevator-cmd-norm')

    def test_nan_attitude_fails_the_trim(self):
        # No contact point has a height to pick the closest to the ground by.
        self.ground_trim_fails('ic/theta-deg')


RunTest(TestGroundTrimNaN)
