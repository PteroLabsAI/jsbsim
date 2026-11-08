# TestFirstPassVRP.py
#
# Check that the systems, which run before FGMassBalance and FGAuxiliary, read
# the visual reference point where the initial conditions put it on the first
# pass of RunIC, its lever arm from the CG included.
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

import xml.etree.ElementTree as et

from JSBSim_utils import JSBSimTestCase, RunTest, FlightModel

# kd 1: a VRP that moves between RunIC's two passes reads as that move over dt.
CLIMB = """
<pid name="test/vrp-climb-fps">
  <input>position/vrp-radius-ft</input>
  <kp>0</kp><ki>0</ki><kd>1</kd>
</pid>"""
# Below the VRP, at the structural origin: a lever arm taken from a CG not yet
# computed, which reads the origin, misses it.
CG_Z_FT = -5.0
# Chosen: a third of the empty weight, far enough below to move the CG by feet.
PILOT_LBS = 1000.0
PILOT_Z_IN = -120.0
STILL_FPS = 1e-3


class TestFirstPassVRP(JSBSimTestCase):
    def load_tripod(self):
        tripod = FlightModel(self, 'tripod')
        cg = et.SubElement(tripod.root.find('mass_balance'), 'location',
                           name='CG', unit='FT')
        for axis, value in zip('xyz', (0.0, 0.0, CG_Z_FT)):
            et.SubElement(cg, axis).text = str(value)
        channel = et.SubElement(et.SubElement(tripod.root, 'system', name='vrp'),
                                'channel', name='vrp')
        channel.append(et.fromstring(CLIMB))
        tripod.tree.write('tripod.xml')
        fdm = tripod.fdm
        fdm.load_model('tripod', False)
        return fdm

    def test_systems_read_the_initial_vrp(self):
        fdm = self.load_tripod()
        fdm.run_ic()
        self.assertAlmostEqual(fdm['test/vrp-climb-fps'], 0.0, delta=STILL_FPS)

    def test_a_mass_moved_while_held_keeps_the_vrp(self):
        # Held, MassBalance keeps the structure, hence the VRP, where it was.
        fdm = self.load_tripod()
        fdm['forces/hold-down'] = 1
        fdm.run_ic()
        fdm.run()
        fdm['inertia/pointmass-weight-lbs[0]'] = PILOT_LBS
        fdm['inertia/pointmass-location-Z-inches[0]'] = PILOT_Z_IN
        fdm.run_ic()
        self.assertAlmostEqual(fdm['test/vrp-climb-fps'], 0.0, delta=STILL_FPS)


RunTest(TestFirstPassVRP)
