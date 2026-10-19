# TestLGearSideFriction.py
#
# Check that a BOGEY contact reports a side friction coefficient of zero when it
# is not touching the ground, from the moment it is loaded and after a reset.
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

SIDE_FRICTION = 'gear/unit/side_friction_coeff'
LOADS = 5
STEPS = 10
AIRBORNE_FT = 100.0
COMPRESSED_FT = 0.9  # the wheel is 1 ft below the CG
SIDE_SPEED_FPS = 10.0
# Bytes, chosen to cover the size of an FGLGear: blocks of 0xFF freed just
# before a load make an uninitialised double read as NaN rather than as 0.
JUNK_SIZES = range(256, 8192, 64)
JUNK_COPIES = 4


def fill_freed_heap_with_nan():
    junk = [b'\xff' * n for n in JUNK_SIZES for _ in range(JUNK_COPIES)]
    del junk


class TestLGearSideFriction(JSBSimTestCase):
    def test_side_friction_coeff_off_the_ground(self):
        with open('gear.xml', 'w') as f:
            f.write('''<?xml version="1.0"?>
<fdm_config name="gear" version="2.0" release="ALPHA">
  <fileheader/>
  <metrics/>
  <mass_balance>
    <ixx unit="SLUG*FT2">1.0</ixx>
    <iyy unit="SLUG*FT2">1.0</iyy>
    <izz unit="SLUG*FT2">1.0</izz>
    <emptywt unit="LBS">100.0</emptywt>
    <location name="CG" unit="IN"><x>0</x><y>0</y><z>0</z></location>
  </mass_balance>
  <ground_reactions>
    <contact type="BOGEY" name="WHEEL">
      <location unit="FT"><x>0</x><y>0</y><z>-1</z></location>
      <static_friction>0.8</static_friction>
      <dynamic_friction>0.5</dynamic_friction>
      <rolling_friction>0.02</rolling_friction>
      <spring_coeff unit="LBS/FT">1000</spring_coeff>
      <damping_coeff unit="LBS/FT/SEC">100</damping_coeff>
    </contact>
  </ground_reactions>
  <aerodynamics/>
</fdm_config>''')

        for _ in range(LOADS):
            fill_freed_heap_with_nan()
            fdm = self.create_fdm()
            fdm.set_aircraft_path('.')
            self.assertTrue(fdm.load_model('gear', False))
            self.assertEqual(fdm[SIDE_FRICTION], 0.0)

            fdm['ic/h-agl-ft'] = AIRBORNE_FT
            fdm.run_ic()
            self.assertEqual(fdm[SIDE_FRICTION], 0.0)
            for _ in range(STEPS):
                fdm.run()
            self.assertEqual(fdm[SIDE_FRICTION], 0.0)

            # Sliding sideways on the ground gives a slip-based coefficient...
            fdm['ic/h-agl-ft'] = COMPRESSED_FT
            fdm['ic/v-fps'] = SIDE_SPEED_FPS
            fdm.reset_to_initial_conditions(0)
            self.assertNotEqual(fdm[SIDE_FRICTION], 0.0)

            # ...which a reset in the air does not carry over.
            fdm['ic/h-agl-ft'] = AIRBORNE_FT
            fdm['ic/v-fps'] = 0.0
            fdm.reset_to_initial_conditions(0)
            self.assertEqual(fdm[SIDE_FRICTION], 0.0)


RunTest(TestLGearSideFriction)
