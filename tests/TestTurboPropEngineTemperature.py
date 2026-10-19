# TestTurboPropEngineTemperature.py
#
# Check that the turboprop engine temperature starts from a defined value: zero
# when loaded, the ambient temperature once a stopped engine is simulated, and
# a warm-up from zero towards the ITT when the engine is set running.
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

LOADS = 5
AIRBORNE_FT = 1000.0
# Bytes, chosen to cover the size of an FGTurboProp: blocks of 0xFF and the
# engines of the previous load, freed before each load, keep an uninitialised
# double from reading 0 by chance.
JUNK_SIZES = range(256, 8192, 64)
JUNK_COPIES = 4


def fill_freed_heap_with_nan():
    junk = [b'\xff' * n for n in JUNK_SIZES for _ in range(JUNK_COPIES)]
    del junk


class TestTurboPropEngineTemperature(JSBSimTestCase):
    def load_dhc6(self):
        fill_freed_heap_with_nan()
        fdm = self.create_fdm()
        self.assertTrue(fdm.load_model('DHC6'))
        engines = [f'propulsion/engine[{n}]/'
                   for n in range(fdm.get_propulsion().get_num_engines())]
        for engine in engines:
            self.assertEqual(fdm[engine + 'engtemp-c'], 0.0)
        fdm['ic/h-agl-ft'] = AIRBORNE_FT
        fdm['forces/hold-down'] = 1.0
        return fdm, engines

    def test_stopped_engine(self):
        for _ in range(LOADS):
            fdm, engines = self.load_dhc6()
            fdm.run_ic()
            for engine in engines:
                self.assertEqual(fdm[engine + 'engtemp-c'], 0.0)
            fdm.run()
            for engine in engines:
                self.assertEqual(fdm[engine + 'engtemp-c'],
                                 fdm['propulsion/tat-c'])

    def test_running_engine(self):
        for _ in range(LOADS):
            fdm, engines = self.load_dhc6()
            # Runs the engines to a steady state, which warms them up.
            fdm['propulsion/set-running'] = -1
            fdm.run_ic()
            for engine in engines:
                self.assertGreater(fdm[engine + 'engtemp-c'], 0.0)
                self.assertLess(fdm[engine + 'engtemp-c'],
                                fdm[engine + 'itt-c'])


RunTest(TestTurboPropEngineTemperature)
