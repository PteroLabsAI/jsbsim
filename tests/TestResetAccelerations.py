# TestResetAccelerations.py
#
# Check that a reset clears the body angular accelerations, so a model that
# feeds them back into its forces does not start the next run with the last
# run's values.
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

DONT_EXECUTE_RUN_IC = 2  # FGFDMExec::DONT_EXECUTE_RUN_IC
ANGULAR = ['accelerations/pdot-rad_sec2', 'accelerations/qdot-rad_sec2',
           'accelerations/rdot-rad_sec2']


class TestResetAccelerations(JSBSimTestCase):
    def test_reset_clears_angular_accelerations(self):
        fdm = self.create_fdm()
        fdm.load_model('c172x')
        fdm.load_ic('reset01', True)
        fdm.run_ic()
        for name in ANGULAR:
            self.assertNotEqual(fdm[name], 0.0, msg=name)

        fdm.reset_to_initial_conditions(DONT_EXECUTE_RUN_IC)
        for name in ANGULAR:
            self.assertEqual(fdm[name], 0.0, msg=name)


RunTest(TestResetAccelerations)
