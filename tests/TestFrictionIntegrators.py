# TestFrictionIntegrators.py
#
# Check that an aircraft standing on its gear comes to rest whatever scheme
# integrates its body rates: the friction forces at rest have to cancel the
# relative motion in one step of the scheme in use, not of the Euler scheme.
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

FT_TO_M = 0.3048
DT = 1/240
# The c172x's CG stands 4.1 ft over the ground: dropped from just above it.
START_AGL_FT = 5.0
SETTLE_S = 50.0
WATCH_S = 10.0
# simulation/integrator/rate/rotational: FGPropagate::eIntegrateType
SCHEMES = {'Euler': 1, 'AB2': 3, 'AB3': 4, 'AB4': 5}
# With friction aimed at an Euler step, AB2 kept the body rates at 4E-4 rad/s
# and AB3 at 2E-2 rad/s, sliding 1.8 ft; aimed at their own step, every scheme
# here settles under 3E-8 rad/s and 1E-4 ft.
REST_RATE_RAD_S = 1E-6
REST_SLIDE_FT = 1E-3


class TestFrictionIntegrators(JSBSimTestCase):
    def stand(self, scheme):
        fdm = self.create_fdm()
        fdm.set_dt(DT)
        fdm.load_model('c172x')
        fdm['ic/lat-geod-deg'] = 47.0
        fdm['ic/long-gc-deg'] = 8.0
        fdm['ic/terrain-elevation-ft'] = 0.0
        fdm['ic/h-agl-ft'] = START_AGL_FT
        fdm['simulation/integrator/rate/rotational'] = scheme
        fdm.run_ic()
        return fdm

    def rest(self, scheme):
        """Largest body rate (rad/s) and slide (ft) once settled. The model is
        closed before anything is asserted, so a failure leaves no output file
        open in the sandbox."""
        fdm = self.stand(scheme)
        for _ in range(round(SETTLE_S/DT)):
            fdm.run()
        lat_m = fdm['position/distance-from-start-lat-mt']
        lon_m = fdm['position/distance-from-start-lon-mt']
        rate = 0.0
        for _ in range(round(WATCH_S/DT)):
            fdm.run()
            rate = max([rate] + [abs(fdm[f'velocities/{a}-rad_sec'])
                                 for a in 'pqr'])
        slide_ft = math.hypot(
            fdm['position/distance-from-start-lat-mt'] - lat_m,
            fdm['position/distance-from-start-lon-mt'] - lon_m) / FT_TO_M
        self.delete_fdm()
        return rate, slide_ft

    def test_comes_to_rest(self):
        for name, scheme in SCHEMES.items():
            with self.subTest(scheme=name):
                rate, slide_ft = self.rest(scheme)
                self.assertLess(rate, REST_RATE_RAD_S)
                self.assertLess(slide_ft, REST_SLIDE_FT)


RunTest(TestFrictionIntegrators)
