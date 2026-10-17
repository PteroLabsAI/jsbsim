# TestStrutDamping.py
#
# Check that the linear damper of a contact is stable whatever its damping: a
# light quadcopter standing on four contacts and kicked in roll comes back to
# rest at a time step where the damper used to flip its roll rate every step.
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

FT_TO_M = 0.3048
N_M_PER_LBS_FT = 14.5939  # JSBSim's own factor (FGXMLElement.cpp)
STANDARD_GRAVITY = 9.80665  # m/s2, only to start close to the rest height

# A 1.4 kg quadcopter (the mass and inertia of the stock F450) on 4 contacts
# at (+/-ARM, +/-ARM, -DROP) from its CG.
MASS_KG = 1.4
IXX_KG_M2 = 0.019
IZZ_KG_M2 = 0.0252
ARM_M = 0.158
DROP_M = 0.060
SPRING_N_M = 3000.0
DT = 1/240
# Roll damping rate 4*c*ARM**2/IXX = 526/s, i.e. 2.2/DT: an explicit damper
# needs less than 2/DT.
STIFF_DAMPING_N_M_S = 100.0
# 0.013/FINE_DT, well under that limit
SOFT_DAMPING_N_M_S = 30.0
FINE_DT = 1/12000
KICK_P_RAD_S = 0.5
SETTLE_S = 1.0
SOFT_RUN_S = 0.25
# The explicit damper flips p by +/-0.24 rad/s every step, forever
REST_P_RAD_S = 1E-3
# Both dampers are first order in time: here they differ by 0.5% (STRUCTURE)
# and 0.7% (BOGEY) of the kick, while a BOGEY damping 13% off, as from a
# wrong strut-to-normal projection, differs by 3.7%.
SOFT_TOLERANCE = 0.015
SPLAY_DEG = 30.0  # leg tilt of the BOGEY variant, so that strut != normal


class TestStrutDamping(JSBSimTestCase):
    def quad(self, damping, dt, contact_type='STRUCTURE', explicit=False):
        k = SPRING_N_M / N_M_PER_LBS_FT
        c = damping / N_M_PER_LBS_FT
        contacts = ''
        for i, (x, y) in enumerate(((ARM_M, ARM_M), (ARM_M, -ARM_M),
                                    (-ARM_M, ARM_M), (-ARM_M, -ARM_M))):
            if explicit:
                # JSBSim leaves <strut_force> functions explicit
                unit = f'contact/unit[{i}]' if contact_type == 'STRUCTURE' \
                    else f'gear/unit[{i}]'
                strut = f'''
      <strut_force><function><sum>
        <product><value>{-k}</value>
          <property>{unit}/compression-ft</property></product>
        <product><value>{-c}</value>
          <property>{unit}/compression-velocity-fps</property></product>
      </sum></function></strut_force>'''
            else:
                strut = f'''
      <spring_coeff unit="LBS/FT">{k}</spring_coeff>
      <damping_coeff unit="LBS/FT/SEC">{c}</damping_coeff>'''
            splay = SPLAY_DEG if y > 0 else -SPLAY_DEG
            contacts += f'''
    <contact type="{contact_type}" name="C{i}">
      <location unit="M"><x>{x}</x><y>{y}</y><z>{-DROP_M}</z></location>
      <orientation unit="DEG"><roll>{splay}</roll><pitch>0</pitch><yaw>0</yaw></orientation>
      <static_friction>0.8</static_friction>
      <dynamic_friction>0.5</dynamic_friction>
      <rolling_friction>0.02</rolling_friction>{strut}
    </contact>'''

        with open('quad.xml', 'w') as f:
            f.write(f'''<?xml version="1.0"?>
<fdm_config name="quad" version="2.0" release="ALPHA">
  <fileheader/>
  <metrics/>
  <mass_balance>
    <ixx unit="KG*M2">{IXX_KG_M2}</ixx>
    <iyy unit="KG*M2">{IXX_KG_M2}</iyy>
    <izz unit="KG*M2">{IZZ_KG_M2}</izz>
    <emptywt unit="KG">{MASS_KG}</emptywt>
    <location name="CG" unit="M"><x>0</x><y>0</y><z>0</z></location>
  </mass_balance>
  <ground_reactions>{contacts}
  </ground_reactions>
  <aerodynamics/>
</fdm_config>''')

        fdm = self.create_fdm()
        fdm.set_aircraft_path('.')
        fdm.set_dt(dt)
        fdm.load_model('quad', False)
        # North pole: no centrifugal nor Coriolis acceleration
        fdm['ic/lat-gc-deg'] = 90.0
        sag_m = MASS_KG * STANDARD_GRAVITY / (4 * SPRING_N_M)
        fdm['ic/h-agl-ft'] = (DROP_M - sag_m) / FT_TO_M
        fdm['ic/p-rad_sec'] = KICK_P_RAD_S
        fdm.run_ic()
        return fdm

    def roll_rates(self, fdm, seconds):
        p = []
        for _ in range(round(seconds / fdm.get_delta_t())):
            fdm.run()
            p.append(fdm['velocities/p-rad_sec'])
        return p

    def test_stiff_damper_comes_to_rest(self):
        for contact_type in ('STRUCTURE', 'BOGEY'):
            with self.subTest(contact_type=contact_type):
                fdm = self.quad(STIFF_DAMPING_N_M_S, DT, contact_type)
                p = self.roll_rates(fdm, SETTLE_S)
                self.assertLess(max(abs(x) for x in p[len(p)//2:]),
                                REST_P_RAD_S)

    def test_spring_is_untouched(self):
        fdm = self.quad(STIFF_DAMPING_N_M_S, DT)
        self.roll_rates(fdm, SETTLE_S)
        # At rest each contact carries a quarter of the weight
        weight = fdm['forces/fbz-weight-lbs']
        k = SPRING_N_M / N_M_PER_LBS_FT
        for i in range(4):
            compression = fdm[f'contact/unit[{i}]/compression-ft']
            self.assertAlmostEqual(4*k*compression/weight, 1.0, delta=1E-5)

    def test_soft_damper_matches_explicit(self):
        for contact_type in ('STRUCTURE', 'BOGEY'):
            with self.subTest(contact_type=contact_type):
                implicit = self.roll_rates(
                    self.quad(SOFT_DAMPING_N_M_S, FINE_DT, contact_type),
                    SOFT_RUN_S)
                explicit = self.roll_rates(
                    self.quad(SOFT_DAMPING_N_M_S, FINE_DT, contact_type, True),
                    SOFT_RUN_S)
                diff = max(abs(a - b) for a, b in zip(implicit, explicit))
                self.assertLess(diff, SOFT_TOLERANCE*KICK_P_RAD_S)


RunTest(TestStrutDamping)
