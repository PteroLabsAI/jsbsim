# TestPropellerInPlaneForce.py
#
# Check the in-plane force of a propeller (the H-force of a rotor): a propeller
# moving edgewise through the air is pushed back at its hub by the force that
# its optional C_INPLANE_FORCE table gives, and by no force without the table.
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

PROP = 'propulsion/engine[0]/'
FORCES = [f'forces/fb{axis}-prop-lbs' for axis in 'xyz']
MOMENTS = [f'moments/{axis}-prop-lbsft' for axis in 'lmn']
DIAMETER_FT = 1.0
# Off the CG (the structural origin) on all three axes, so that the in-plane
# force has a moment about each.
HUB_IN = (12.0, -6.0, 4.0)
IN_PER_FT = 12.0
HUB_FT = (-HUB_IN[0]/IN_PER_FT, HUB_IN[1]/IN_PER_FT, -HUB_IN[2]/IN_PER_FT)
C_PER_MU = 0.05  # a linear rotor drag, C_H = c*mu
MU_MAX = 1.0  # end of the C_INPLANE_FORCE table
ALTITUDE_FT = 1000.0
# u, v and a yaw rate that moves the hub by another omega x r
EDGEWISE = (20.0, -10.0, 2.0)
AT_REST = (0.0, 0.0, 0.0)
THROTTLE = 1.0
SPIN_UP_S = 2.0
RTOL = 1E-9


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


class TestPropellerInPlaneForce(JSBSimTestCase):
    def edgewise(self, sense, c_per_mu, throttle=THROTTLE, motion=EDGEWISE):
        inplane = '' if c_per_mu is None else f'''
          <table name="C_INPLANE_FORCE" type="internal">
            <tableData>
              0.0 0.0
              {MU_MAX} {c_per_mu*MU_MAX}
            </tableData>
          </table>'''
        x, y, z = HUB_IN
        with open('edgewise.xml', 'w') as f:
            f.write(f'''<?xml version="1.0"?>
<fdm_config name="edgewise" version="2.0" release="ALPHA">
  <fileheader/>
  <metrics/>
  <mass_balance>
    <ixx unit="SLUG*FT2">1.0</ixx>
    <iyy unit="SLUG*FT2">1.0</iyy>
    <izz unit="SLUG*FT2">1.0</izz>
    <emptywt unit="LBS">10.0</emptywt>
    <location name="CG" unit="IN"><x>0</x><y>0</y><z>0</z></location>
  </mass_balance>
  <propulsion>
    <engine>
      <electric_engine><power unit="WATTS">200.0</power></electric_engine>
      <thruster>
        <propeller version="1.1">
          <ixx unit="SLUG*FT2">0.0001</ixx>
          <diameter unit="FT">{DIAMETER_FT}</diameter>
          <numblades>2</numblades>
          <table name="C_THRUST" type="internal">
            <tableData>
              0.0 0.1
              1.0 0.1
            </tableData>
          </table>
          <table name="C_POWER" type="internal">
            <tableData>
              0.0 0.05
              1.0 0.05
            </tableData>
          </table>{inplane}
        </propeller>
        <location unit="IN"><x>{x}</x><y>{y}</y><z>{z}</z></location>
        <!-- thrust up: the disk is the body x-y plane -->
        <orient unit="DEG"><roll>0</roll><pitch>90</pitch><yaw>0</yaw></orient>
        <sense>{sense}</sense>
      </thruster>
    </engine>
  </propulsion>
  <ground_reactions/>
  <!-- no aerodynamic force, in body axes so that the airspeed is allowed -->
  <aerodynamics><axis name="X"/></aerodynamics>
</fdm_config>''')

        fdm = self.create_fdm()
        fdm.set_aircraft_path('.')
        self.assertTrue(fdm.load_model('edgewise', False))
        fdm['ic/h-agl-ft'] = ALTITUDE_FT
        fdm['ic/u-fps'], fdm['ic/v-fps'], fdm['ic/r-rad_sec'] = motion
        # Held: the forces move nothing, so runs with and without the table
        # see the same air.
        fdm['forces/hold-down'] = 1.0
        fdm.run_ic()
        fdm['fcs/throttle-cmd-norm'] = throttle
        while fdm.get_sim_time() < SPIN_UP_S:
            fdm.run()
        # One more frame without time: its forces belong to the state read.
        fdm.suspend_integration()
        fdm.run()
        return fdm

    def expected_force(self, fdm):
        n = fdm[PROP + 'propeller-rpm'] / 60.0
        uvw = [fdm[f'velocities/{axis}-aero-fps'] for axis in 'uvw']
        pqr = [fdm[f'velocities/{axis}-aero-rad_sec'] for axis in 'pqr']
        u, v, _ = [a + b for a, b in zip(uvw, cross(pqr, HUB_FT))]
        speed = math.hypot(u, v)  # in the disk plane
        mu = speed / (n * DIAMETER_FT)
        self.assertLess(mu, MU_MAX)
        h = C_PER_MU * mu * fdm['atmosphere/rho-slugs_ft3'] * n**2 \
            * DIAMETER_FT**4
        return h, (-h*u/speed, -h*v/speed, 0.0)

    def test_no_table_no_force(self):
        fdm = self.edgewise(1.0, None)
        self.assertGreater(fdm[PROP + 'thrust-lbs'], 0.0)
        self.assertEqual(fdm[PROP + 'propeller-in-plane-force-lbs'], 0.0)
        # Thrust alone, along the body z axis, has no side component at all.
        self.assertEqual(fdm['forces/fby-prop-lbs'], 0.0)

    def test_force_opposes_edgewise_flight(self):
        for sense in (1.0, -1.0):
            with self.subTest(sense=sense):
                bare = self.edgewise(sense, None)
                fdm = self.edgewise(sense, C_PER_MU)
                h, force = self.expected_force(fdm)
                self.assertAlmostEqual(
                    fdm[PROP + 'propeller-in-plane-force-lbs'], h,
                    delta=RTOL*h)
                # The table adds a force and leaves the rest untouched.
                for name in ('thrust-lbs', 'propeller-rpm',
                             'propeller-torque-ftlb'):
                    self.assertEqual(fdm[PROP + name], bare[PROP + name])
                for name, f in zip(FORCES, force):
                    self.assertAlmostEqual(fdm[name] - bare[name], f,
                                           delta=RTOL*h, msg=name)
                for name, m in zip(MOMENTS, cross(HUB_FT, force)):
                    self.assertAlmostEqual(fdm[name] - bare[name], m,
                                           delta=RTOL*h, msg=name)

    def test_finite_at_rest(self):
        for throttle, motion in ((0.0, EDGEWISE), (THROTTLE, AT_REST),
                                 (0.0, AT_REST)):
            with self.subTest(throttle=throttle, motion=motion):
                fdm = self.edgewise(1.0, C_PER_MU, throttle, motion)
                self.assertEqual(fdm[PROP + 'propeller-in-plane-force-lbs'],
                                 0.0)
                for name in FORCES + MOMENTS:
                    self.assertTrue(math.isfinite(fdm[name]), msg=name)


RunTest(TestPropellerInPlaneForce)
