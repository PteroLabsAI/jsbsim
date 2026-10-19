# TestPropellerDiameter.py
#
# Check that a propeller without a <diameter> is refused at load time: thrust
# and power scale with the diameter, which has no meaningful default.
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
from jsbsim import FGLogger, get_logger, set_logger


class LogBuffer(FGLogger):
    def __init__(self):
        super().__init__()
        self.buffer = ''

    def message(self, message: str) -> None:
        self.buffer += message


class TestPropellerDiameter(JSBSimTestCase):
    def test_missing_diameter_is_refused(self):
        with open('prop.xml', 'w') as f:
            f.write('''<?xml version="1.0"?>
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

        log = LogBuffer()
        self.addCleanup(set_logger, get_logger())
        set_logger(log)
        fdm = self.create_fdm()
        fdm.set_aircraft_path('.')
        self.assertFalse(fdm.load_model('prop', False))
        self.assertIn('<diameter> is a mandatory parameter', log.buffer)


RunTest(TestPropellerDiameter)
