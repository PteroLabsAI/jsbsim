# TestRotorMomentum.py
#
# Check that the implicit midpoint rule carries the angular momentum H of a
# spinning propeller implicitly: a free body keeps its kinetic energy and
# |J*w + H|, nutates at the frequency the midpoint rule predicts, and agrees
# with explicit Euler as the step shrinks. With the gyroscopic moment of the
# propeller added explicitly, the nutation gains energy every step and grows
# until it overflows.
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

# FGPropagate::eIntegrateType
RECT_EULER = 1
IMPLICIT_MIDPOINT = 10
RATES_HZ = (120, 240, 1000)
SECONDS = 10.0
# Chosen: spins the propeller to 600 to 620 rad/s from 1000 Hz to 120 Hz, H of
# about 60 slug*ft2/s, a nutation of about 60 rad/s: 0.5 rad a step at 120 Hz.
PROP_IXX_SLUGFT2 = 0.1
POWER_W = 2440.0
SPIN_UP_S = 10.0
# Chosen: a small transverse rate on a slow roll.
RATES_RAD_S = (1.0, 0.1, 0.0)
# Chosen: clear of the ground for the fall of the whole run.
ALTITUDE_FT = 30000.0
# (Ixx, Iyy, Izz) slug*ft2 and the propeller orientation (roll, pitch, yaw) deg.
# Chosen: a body symmetric about the propeller axis, and one that is not, the
# propeller tilted off all of its axes.
MODELS = {'axisymmetric': ((0.5, 1.0, 1.0), (0.0, 0.0, 0.0)),
          'tilted': ((0.6, 1.0, 1.3), (10.0, 30.0, 20.0))}
# Requested bound, far above the rounding the midpoint rule leaves.
DRIFT = 1E-9
# Chosen: steps short enough for Euler's error, h*W^2*t/2 at a nutation W, to
# shrink with h. The run covers about 6 rad of nutation.
FINE_RATES_HZ = (4000, 16000)
AGREEMENT_S = 0.1
# Chosen: room for the higher order terms of Euler's error.
RATIO_TOLERANCE = 0.2
INERTIAL_RATES = ('velocities/pi-rad_sec', 'velocities/qi-rad_sec',
                  'velocities/ri-rad_sec')
ANGULAR_MOMENTUM = ('propulsion/angular-momentum-x-slugs_ft2_sec',
                    'propulsion/angular-momentum-y-slugs_ft2_sec',
                    'propulsion/angular-momentum-z-slugs_ft2_sec')


def energy(w, Jw):
    return 0.5*sum(a*b for a, b in zip(w, Jw))


def magnitude(a, b):
    """|a + b|"""
    return math.hypot(*(x + y for x, y in zip(a, b)))


class TestRotorMomentum(JSBSimTestCase):
    def setUp(self):
        JSBSimTestCase.setUp(self)
        for name, ((ixx, iyy, izz), (roll, pitch, yaw)) in MODELS.items():
            with open(name + '.xml', 'w') as f:
                f.write(f'''<?xml version="1.0"?>
<fdm_config name="{name}" version="2.0" release="ALPHA">
  <fileheader/>
  <metrics/>
  <mass_balance>
    <ixx unit="SLUG*FT2">{ixx}</ixx>
    <iyy unit="SLUG*FT2">{iyy}</iyy>
    <izz unit="SLUG*FT2">{izz}</izz>
    <emptywt unit="LBS">10.0</emptywt>
    <location name="CG" unit="IN"><x>0</x><y>0</y><z>0</z></location>
  </mass_balance>
  <propulsion>
    <engine>
      <electric_engine><power unit="WATTS">{POWER_W}</power></electric_engine>
      <thruster>
        <propeller version="1.1">
          <ixx unit="SLUG*FT2">{PROP_IXX_SLUGFT2}</ixx>
          <diameter unit="M">0.5</diameter>
          <table name="C_THRUST" type="internal">
            <tableData>
              0.0 0.0
              1.0 0.0
            </tableData>
          </table>
          <table name="C_POWER" type="internal">
            <tableData>
              0.0 0.0
              1.0 0.0
            </tableData>
          </table>
        </propeller>
        <location unit="IN"><x>0</x><y>0</y><z>0</z></location>
        <orient unit="DEG"><roll>{roll}</roll><pitch>{pitch}</pitch><yaw>{yaw}</yaw></orient>
      </thruster>
    </engine>
  </propulsion>
  <ground_reactions/>
  <aerodynamics/>
</fdm_config>''')

    def start(self, model, hz, scheme):
        """The model released at RATES_RAD_S, its propeller spun up held down
        and then left without power or drag, so H stays constant."""
        fdm = self.create_fdm()
        fdm.set_aircraft_path('.')
        fdm.set_dt(1.0/hz)
        self.assertTrue(fdm.load_model(model, False))
        fdm['simulation/integrator/rate/rotational'] = scheme
        fdm['ic/h-sl-ft'] = ALTITUDE_FT
        fdm['forces/hold-down'] = 1.0
        fdm.run_ic()
        fdm['fcs/throttle-cmd-norm'] = 1.0
        for _ in range(round(SPIN_UP_S*hz)):
            fdm.run()
        fdm['fcs/throttle-cmd-norm'] = 0.0
        fdm['forces/hold-down'] = 0.0
        for prop, rate in zip(('ic/p-rad_sec', 'ic/q-rad_sec', 'ic/r-rad_sec'),
                              RATES_RAD_S):
            fdm[prop] = rate
        fdm.run_ic()
        return fdm

    def rates(self, fdm):
        return [fdm[p] for p in INERTIAL_RATES]

    def momentum(self, fdm, w):
        """J*w of the diagonal inertia."""
        J = [fdm['inertia/i'+a+'-slugs_ft2'] for a in ('xx', 'yy', 'zz')]
        return [j*v for j, v in zip(J, w)]

    def test_keeps_energy_and_momentum(self):
        for model in MODELS:
            for hz in RATES_HZ:
                with self.subTest(model=model, hz=hz):
                    fdm = self.start(model, hz, IMPLICIT_MIDPOINT)
                    w0 = self.rates(fdm)
                    for _ in range(round(SECONDS*hz)):
                        fdm.run()
                    w1 = self.rates(fdm)
                    Jw0, Jw1 = self.momentum(fdm, w0), self.momentum(fdm, w1)
                    self.assertAlmostEqual(energy(w1, Jw1)/energy(w0, Jw0), 1.0,
                                           delta=DRIFT)
                    # H is read at the end: a change would show in |J*w + H|.
                    H = [fdm[p] for p in ANGULAR_MOMENTUM]
                    self.assertAlmostEqual(magnitude(Jw1, H)/magnitude(Jw0, H),
                                           1.0, delta=DRIFT)
                    self.delete_fdm()

    def test_nutation_frequency(self):
        (ia, it, _), _ = MODELS['axisymmetric']
        for hz in RATES_HZ:
            with self.subTest(hz=hz):
                h = 1.0/hz
                fdm = self.start('axisymmetric', hz, IMPLICIT_MIDPOINT)
                H = fdm[ANGULAR_MOMENTUM[0]]
                p, q, r = self.rates(fdm)
                # The transverse rates turn at W in body axes.
                W = (H + (ia - it)*p)/it
                angle = 0.0
                steps = round(SECONDS*hz)
                for _ in range(steps):
                    fdm.run()
                    _, q1, r1 = self.rates(fdm)
                    angle += math.atan2(q*r1 - r*q1, q*q1 + r*r1)
                    q, r = q1, r1
                # The midpoint rule turns this motion by 2*atan(W*h/2)
                # a step, W*h as h goes to 0.
                self.assertAlmostEqual(angle/steps, 2.0*math.atan(0.5*W*h),
                                       delta=DRIFT*W*h)
                self.delete_fdm()

    def test_agrees_with_euler_as_the_step_shrinks(self):
        distances = []
        for hz in FINE_RATES_HZ:
            final = {}
            for scheme in (RECT_EULER, IMPLICIT_MIDPOINT):
                fdm = self.start('tilted', hz, scheme)
                for _ in range(round(AGREEMENT_S*hz)):
                    fdm.run()
                final[scheme] = self.rates(fdm)
                self.delete_fdm()
            distances.append(math.dist(final[RECT_EULER], final[IMPLICIT_MIDPOINT]))
        # Euler is first order: its distance shrinks as the step.
        step_ratio = FINE_RATES_HZ[1]/FINE_RATES_HZ[0]
        self.assertAlmostEqual(distances[0]/distances[1], step_ratio,
                               delta=RATIO_TOLERANCE*step_ratio)


RunTest(TestRotorMomentum)
