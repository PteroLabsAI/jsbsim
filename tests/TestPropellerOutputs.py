# TestPropellerOutputs.py
#
# Check that the propeller outputs that a system can read before the propeller
# is first computed start from zero: the systems run before the propulsion in
# every frame, so the first frame of RunIC reads them as they were initialised.
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

OUTPUTS = ('thrust-coefficient', 'advance-ratio', 'helical-tip-Mach',
           'propeller-power-ftlbps')
PROP = 'propulsion/engine/'
LOADS = 5
STEPS = 10
AIRSPEED_FPS = 20.0
# Bytes, chosen to cover the size of an FGPropeller: blocks of 0xFF and the
# propellers of the previous loads, freed before each load, keep an
# uninitialised double from reading 0 by chance.
JUNK_SIZES = range(256, 8192, 64)
JUNK_COPIES = 4


def fill_freed_heap_with_nan():
    junk = [b'\xff' * n for n in JUNK_SIZES for _ in range(JUNK_COPIES)]
    del junk


# Each switch keeps the value its output had in the first frame of the system.
FIRST_FRAME = '\n'.join(f'''
      <switch name="first-{name}">
        <default value="test/first/{name}"/>
        <test value="{PROP}{name}">test/latched == 0</test>
        <output>test/first/{name}</output>
      </switch>''' for name in OUTPUTS)


class TestPropellerOutputs(JSBSimTestCase):
    def test_outputs_read_in_the_first_frame(self):
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
      <electric_engine><power unit="WATTS">100</power></electric_engine>
      <thruster>
        <propeller>
          <ixx unit="SLUG*FT2">0.001</ixx>
          <diameter unit="M">0.5</diameter>
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
  <system name="first-frame">
    <channel name="first-frame">{FIRST_FRAME}
      <switch name="latch">
        <default value="1"/>
        <output>test/latched</output>
      </switch>
    </channel>
  </system>
  <aerodynamics/>
</fdm_config>''')

        for _ in range(LOADS):
            fill_freed_heap_with_nan()
            fdm = self.create_fdm()
            fdm.set_aircraft_path('.')
            self.assertTrue(fdm.load_model('prop', False))
            fdm['ic/u-fps'] = AIRSPEED_FPS
            fdm['forces/hold-down'] = 1.0
            fdm.run_ic()
            self.check_first_frame(fdm)

            # Spin the propeller up, then check that a reset starts it over.
            fdm['fcs/throttle-cmd-norm'] = 1.0
            for _ in range(STEPS):
                fdm.run()
            for name in OUTPUTS:
                self.assertNotEqual(fdm[PROP + name], 0.0)
            fdm['test/latched'] = 0.0
            fdm.reset_to_initial_conditions(0)
            self.check_first_frame(fdm)

    def check_first_frame(self, fdm):
        for name in OUTPUTS:
            with self.subTest(output=name):
                self.assertEqual(fdm['test/first/' + name], 0.0)


RunTest(TestPropellerOutputs)
