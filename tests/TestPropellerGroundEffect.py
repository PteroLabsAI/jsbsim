# TestPropellerGroundEffect.py
#
# Check the ground effect of a propeller: the thrust of a propeller that blows
# down is multiplied by its optional CT_GROUND_FACTOR table, indexed by the
# distance from the hub to the ground along the slipstream over the radius. A
# propeller that blows sideways or up, or has no table, keeps its thrust.
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
FACTOR = PROP + 'propeller-ground-factor'
THRUST = PROP + 'thrust-lbs'
DIAMETER_FT = 1.0
RADIUS_FT = DIAMETER_FT / 2.0
# Off the CG (the structural origin) on all three axes, so that a tilt moves
# the hub up or down; below the CG, so that the CG stays above the ground.
HUB_IN = (12.0, -6.0, -4.0)
IN_PER_FT = 12.0
HUB_FT = (-HUB_IN[0]/IN_PER_FT, HUB_IN[1]/IN_PER_FT, -HUB_IN[2]/IN_PER_FT)
# Cheeseman-Bennett, T/T_inf = 1/(1 - (R/4z)^2), cut to 1 at z/R = 8
GROUND_TABLE = [(z_r, 1.0/(1.0 - (0.25/z_r)**2)) for z_r in (0.5, 1.0, 2.0)] \
    + [(8.0, 1.0)]
CT = 0.1
LIFT_DEG = 90.0  # orient pitch of a propeller that blows down the body z axis
PUSHER_DEG = 0.0  # blows aft
INVERTED_DEG = -90.0  # blows up the body z axis
TILT_DEG = 45.0
# Every run is at this altitude, the terrain moves: the air is the same.
ALTITUDE_FT = 1000.0
FAR_FT = 900.0  # hub height far above the table
SPIN_UP_S = 1.0
RTOL = 1E-9
# The IC sets the height to some 1E-6 ft, which moves the factor by less.
IC_TOL = 1E-5


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def down_in_body(theta, phi):
    return (-math.sin(theta), math.sin(phi)*math.cos(theta),
            math.cos(phi)*math.cos(theta))


def interp(table, x):
    if x <= table[0][0]:
        return table[0][1]
    for (x0, y0), (x1, y1) in zip(table, table[1:]):
        if x <= x1:
            return y0 + (x - x0)/(x1 - x0)*(y1 - y0)
    return table[-1][1]


class TestPropellerGroundEffect(JSBSimTestCase):
    def load(self, pitch_deg, hub_agl_ft, theta_deg=0.0, phi_deg=0.0, ct=CT,
             table=GROUND_TABLE):
        rows = '\n'.join(f'{z_r!r} {k!r}' for z_r, k in table or ())
        ground = '' if table is None else f'''
          <table name="CT_GROUND_FACTOR" type="internal">
            <tableData>
              {rows}
            </tableData>
          </table>'''
        x, y, z = HUB_IN
        with open('ground.xml', 'w') as f:
            f.write(f'''<?xml version="1.0"?>
<fdm_config name="ground" version="2.0" release="ALPHA">
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
              0.0 {ct}
              1.0 {ct}
            </tableData>
          </table>
          <table name="C_POWER" type="internal">
            <tableData>
              0.0 0.05
              1.0 0.05
            </tableData>
          </table>{ground}
        </propeller>
        <location unit="IN"><x>{x}</x><y>{y}</y><z>{z}</z></location>
        <orient unit="DEG">
          <roll>0</roll><pitch>{pitch_deg}</pitch><yaw>0</yaw>
        </orient>
      </thruster>
    </engine>
  </propulsion>
  <ground_reactions/>
  <aerodynamics><axis name="X"/></aerodynamics>
</fdm_config>''')

        fdm = self.create_fdm()
        fdm.set_aircraft_path('.')
        self.assertTrue(fdm.load_model('ground', False))
        theta, phi = math.radians(theta_deg), math.radians(phi_deg)
        cg_agl_ft = hub_agl_ft + dot(down_in_body(theta, phi), HUB_FT)
        fdm['ic/h-sl-ft'] = ALTITUDE_FT
        fdm['ic/terrain-elevation-ft'] = ALTITUDE_FT - cg_agl_ft
        fdm['ic/theta-rad'], fdm['ic/phi-rad'] = theta, phi
        # Held: nothing moves, so all runs see the same air and RPM.
        fdm['forces/hold-down'] = 1.0
        fdm.run_ic()
        return fdm

    def spin_up(self, *args, **kwargs):
        fdm = self.load(*args, **kwargs)
        fdm['fcs/throttle-cmd-norm'] = 1.0
        while fdm.get_sim_time() < SPIN_UP_S:
            fdm.run()
        # One more frame without time: its thrust belongs to the state read.
        fdm.suspend_integration()
        fdm.run()
        self.assertTrue(math.isfinite(fdm[THRUST]))
        return fdm

    def expected_factor(self, fdm, pitch_deg):
        down = down_in_body(fdm['attitude/theta-rad'], fdm['attitude/phi-rad'])
        pitch = math.radians(pitch_deg)
        # The slipstream leaves the disk opposite the thrust.
        wake = (-math.cos(pitch), 0.0, math.sin(pitch))
        wake_down = dot(down, wake)*math.copysign(1.0, fdm[THRUST])
        if wake_down <= 0.0:
            return 1.0
        hub_agl_ft = fdm['position/h-agl-ft'] - dot(down, HUB_FT)
        return interp(GROUND_TABLE, hub_agl_ft/(wake_down*RADIUS_FT))

    def assertGroundEffect(self, fdm, far, factor):
        self.assertEqual(fdm['atmosphere/rho-slugs_ft3'],
                         far['atmosphere/rho-slugs_ft3'])
        self.assertEqual(fdm[PROP + 'propeller-rpm'],
                         far[PROP + 'propeller-rpm'])
        self.assertEqual(far[FACTOR], 1.0)
        self.assertAlmostEqual(fdm[FACTOR], factor, delta=RTOL*factor)
        self.assertAlmostEqual(fdm[THRUST]/far[THRUST], factor,
                               delta=RTOL*factor)

    def test_no_table_no_effect(self):
        far = self.spin_up(LIFT_DEG, FAR_FT, table=None)
        fdm = self.spin_up(LIFT_DEG, 0.0, table=None)
        self.assertGreater(far[THRUST], 0.0)
        self.assertEqual(fdm[THRUST], far[THRUST])
        self.assertEqual(fdm[FACTOR], 1.0)
        self.assertEqual(far[FACTOR], 1.0)

    def test_lift_propeller(self):
        far = self.spin_up(LIFT_DEG, FAR_FT)
        for z_r in (0.5, 0.75, 1.0, 1.5, 2.0, 5.0, 8.0, 20.0):
            with self.subTest(z_r=z_r):
                fdm = self.spin_up(LIFT_DEG, z_r*RADIUS_FT)
                factor = self.expected_factor(fdm, LIFT_DEG)
                self.assertGroundEffect(fdm, far, factor)
                self.assertAlmostEqual(factor, interp(GROUND_TABLE, z_r),
                                       delta=IC_TOL)

    def test_no_lag(self):
        # No filter on the height: the factor is right from the first frame.
        self.assertEqual(self.load(LIFT_DEG, FAR_FT)[FACTOR], 1.0)
        fdm = self.load(LIFT_DEG, RADIUS_FT)
        self.assertAlmostEqual(fdm[FACTOR], GROUND_TABLE[1][1], delta=IC_TOL)

    def test_tilted(self):
        far = self.spin_up(PUSHER_DEG, FAR_FT)
        for hub_agl_ft in (0.0, RADIUS_FT, 4.0*RADIUS_FT):
            with self.subTest(pusher=hub_agl_ft):
                fdm = self.spin_up(PUSHER_DEG, hub_agl_ft)
                self.assertEqual(fdm[FACTOR], 1.0)
                self.assertEqual(fdm[THRUST], far[THRUST])

        # Tilted by 45 deg, z/R = 2 at a hub height of 2R*cos(45 deg), where
        # the hub height alone would give z/R = 1.41.
        hub_agl_ft = 2.0*RADIUS_FT*math.cos(math.radians(TILT_DEG))
        for pitch_deg, theta_deg, phi_deg in (
                (TILT_DEG, 0.0, 0.0),  # the propeller on the airframe
                (LIFT_DEG, TILT_DEG, 0.0),  # the airframe, nose up
                (LIFT_DEG, 0.0, -TILT_DEG),  # the airframe, rolled left
                (PUSHER_DEG, TILT_DEG, 0.0)):  # a pusher blowing aft and down
            with self.subTest(pitch=pitch_deg, theta=theta_deg, phi=phi_deg):
                far = self.spin_up(pitch_deg, FAR_FT, theta_deg, phi_deg)
                fdm = self.spin_up(pitch_deg, hub_agl_ft, theta_deg, phi_deg)
                factor = self.expected_factor(fdm, pitch_deg)
                self.assertGroundEffect(fdm, far, factor)
                self.assertAlmostEqual(factor, GROUND_TABLE[2][1],
                                       delta=IC_TOL)

        # A pusher nose down blows aft and up.
        fdm = self.spin_up(PUSHER_DEG, 0.0, -TILT_DEG)
        self.assertEqual(fdm[FACTOR], 1.0)

    def test_blowing_up(self):
        for pitch_deg, ct in ((INVERTED_DEG, CT), (LIFT_DEG, -CT)):
            with self.subTest(pitch=pitch_deg, ct=ct):
                far = self.spin_up(pitch_deg, FAR_FT, ct=ct)
                fdm = self.spin_up(pitch_deg, RADIUS_FT, ct=ct)
                self.assertEqual(fdm[FACTOR], 1.0)
                self.assertEqual(fdm[THRUST], far[THRUST])
        # Inverted and reversed, the slipstream goes down again.
        far = self.spin_up(INVERTED_DEG, FAR_FT, ct=-CT)
        fdm = self.spin_up(INVERTED_DEG, RADIUS_FT, ct=-CT)
        factor = self.expected_factor(fdm, INVERTED_DEG)
        self.assertGroundEffect(fdm, far, factor)
        self.assertAlmostEqual(factor, GROUND_TABLE[1][1], delta=IC_TOL)

    def test_on_the_ground(self):
        far = self.spin_up(LIFT_DEG, FAR_FT)
        # Down to the ground and through it, the table holds its first value.
        for hub_agl_ft in (0.0, -RADIUS_FT):
            with self.subTest(hub_agl_ft=hub_agl_ft):
                fdm = self.spin_up(LIFT_DEG, hub_agl_ft)
                self.assertGroundEffect(fdm, far, GROUND_TABLE[0][1])


RunTest(TestPropellerGroundEffect)
