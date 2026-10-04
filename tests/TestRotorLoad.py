# TestRotorLoad.py
#
# Check that a rotor linked with <ExternalRPM load="1"> loads its source
# rotor's transmission, so the engine drives both rotors.
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

import math, os, shutil
import xml.etree.ElementTree as et
from JSBSim_utils import JSBSimTestCase, ExecuteUntil, RunTest


class TestRotorLoad(JSBSimTestCase):
    def setUp(self):
        JSBSimTestCase.setUp(self)
        shutil.copytree(self.sandbox.path_to_jsbsim_file('aircraft', 'ah1s'),
                        os.path.join('aircraft', 'ah1s'))

        # Drive the main rotor at constant power with no governor, so that the
        # rotor speed settles where the engine power equals the rotors' power.
        motor = et.Element('electric_engine', name='test_motor')
        et.SubElement(motor, 'power', unit='HP').text = '500.0'
        et.ElementTree(motor).write(self.ah1s_file('Engines', 'test_motor.xml'))

        tree = et.parse(self.ah1s_file('ah1s.xml'))
        root = tree.getroot()
        root.find('propulsion/engine').set('file', 'test_motor')
        for system in root.findall('system'):
            if system.get('file') == 'rpm_governor':
                root.remove(system)
        tree.write(self.ah1s_file('ah1s.xml'))

        # No gear friction: the engine power is then exactly the rotors' power.
        def no_gearloss(root):
            et.SubElement(root, 'gearloss', unit='HP').text = '0.0'
        self.edit_rotor('ah1s_rotor', no_gearloss)

    def ah1s_file(self, *path):
        return os.path.join('aircraft', 'ah1s', *path)

    def edit_rotor(self, name, change):
        path = self.ah1s_file('Engines', name + '.xml')
        tree = et.parse(path)
        change(tree.getroot())
        tree.write(path)

    def link_tail_rotor(self, source, load):
        def change(root):
            extrpm = root.find('ExternalRPM')
            extrpm.text = str(source)
            extrpm.attrib.pop('load', None)
            if load is not None:
                extrpm.set('load', load)
        self.edit_rotor('ah1s_tail_rotor', change)

    def rotor_powers(self, source, load):
        self.link_tail_rotor(source, load)
        fdm = self.create_fdm()
        fdm.set_aircraft_path('aircraft')
        fdm.load_model('ah1s')
        fdm.load_ic('reset00', True)
        fdm['propulsion/engine[1]/x-rpm-dict'] = 6600.0  # only read with ExternalRPM -1
        fdm.run_ic()
        fdm['forces/hold-down'] = 1.0
        fdm['fcs/throttle-cmd-norm'] = 1.0
        fdm['fcs/collective-cmd-norm'] = 0.3
        fdm['fcs/rudder-cmd-norm'] = 0.6
        ExecuteUntil(fdm, 60.0)

        def rotor_hp(i):
            omega = fdm['propulsion/engine[%d]/rotor-rpm' % i] * math.pi / 30.0
            return fdm['propulsion/engine[%d]/torque-lbsft' % i] * omega / 550.0

        engine, main, tail = fdm['propulsion/engine[0]/power-hp'], rotor_hp(0), rotor_hp(1)
        self.delete_fdm()
        self.assertGreater(tail, 0.05 * main)  # well above the 1% tolerance below
        return engine, main, tail

    def testLoadedTailRotor(self):
        engine, main, tail = self.rotor_powers(0, '1')
        self.assertAlmostEqual(engine / (main + tail), 1.0, delta=0.01)

    def testFreeTailRotor(self):
        for load in (None, '0'):
            engine, main, tail = self.rotor_powers(0, load)
            self.assertAlmostEqual(engine / main, 1.0, delta=0.01)

    def testExternalControlIgnoresLoad(self):
        engine, main, tail = self.rotor_powers(-1, '1')
        self.assertAlmostEqual(engine / main, 1.0, delta=0.01)


RunTest(TestRotorLoad)
