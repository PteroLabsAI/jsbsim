# TestPropellerInertia.py
#
# Check that a propeller spins up with the inertia it is given, however small:
# the rotor and propeller of a small UAV can be well under 1e-6 slug*ft2.
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
SMALL_IXX_SLUG_FT2 = 4.46E-8  # 6.05e-8 kg*m2, a Crazyflie 2.1 Brushless rotor
ZERO_IXX_SLUG_FT2 = 1E-6  # what FGPropeller takes for a zero or missing ixx
# Small enough that the first step from rest, which FGPropeller takes at
# 1 rad/s, stays under the steady RPM.
POWER_W = 0.001
HP_TO_FT_LBS_PER_S = 550.0
STEPS = 10
RTOL = 1E-9


class TestPropellerInertia(JSBSimTestCase):
    def spin_up(self, ixx):
        ixx = '' if ixx is None else f'<ixx unit="SLUG*FT2">{ixx}</ixx>'
        with open('prop.xml', 'w') as f:
            f.write(f'''<?xml version="1.0"?>
<fdm_config name="prop" version="2.0" release="ALPHA">
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
      <electric_engine><power unit="WATTS">{POWER_W}</power></electric_engine>
      <thruster>
        <propeller>
          {ixx}
          <diameter unit="M">0.055</diameter>
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
          </table>
        </propeller>
        <location unit="IN"><x>0</x><y>0</y><z>0</z></location>
        <orient unit="DEG"><roll>0</roll><pitch>0</pitch><yaw>0</yaw></orient>
      </thruster>
    </engine>
  </propulsion>
  <ground_reactions/>
  <aerodynamics/>
</fdm_config>''')

        fdm = self.create_fdm()
        fdm.set_aircraft_path('.')
        self.assertTrue(fdm.load_model('prop', False))
        fdm['forces/hold-down'] = 1.0
        fdm.run_ic()
        fdm['fcs/throttle-cmd-norm'] = 1.0
        return fdm

    def check_inertia(self, fdm, ixx):
        dt = fdm.get_delta_t()
        rpm = fdm[PROP + 'propeller-rpm']
        self.assertEqual(rpm, 0.0)
        for _ in range(STEPS):
            fdm.run()
            power = fdm[PROP + 'power-hp'] * HP_TO_FT_LBS_PER_S \
                - fdm[PROP + 'propeller-power-ftlbps']
            omega = rpm / 60.0 * 2.0 * math.pi
            torque = power / (omega if omega > 0.01 else 1.0)
            new_rpm = fdm[PROP + 'propeller-rpm']
            self.assertAlmostEqual(new_rpm - rpm,
                                   torque / ixx * dt * 60.0 / (2.0 * math.pi),
                                   delta=RTOL * new_rpm)
            rpm = new_rpm
        self.assertGreater(rpm, 0.0)

    def test_small_inertia_is_used(self):
        self.check_inertia(self.spin_up(SMALL_IXX_SLUG_FT2),
                           SMALL_IXX_SLUG_FT2)

    def test_zero_or_missing_inertia(self):
        for ixx in (0.0, None):
            with self.subTest(ixx=ixx):
                self.check_inertia(self.spin_up(ixx), ZERO_IXX_SLUG_FT2)


RunTest(TestPropellerInertia)
