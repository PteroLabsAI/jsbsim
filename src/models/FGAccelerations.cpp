/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Module:       FGAccelerations.cpp
 Author:       Jon S. Berndt
 Date started: 07/12/11
 Purpose:      Calculates derivatives of rotational and translational rates, and
               of the attitude quaternion.
 Called by:    FGFDMExec

 ------------- Copyright (C) 2011  Jon S. Berndt (jon@jsbsim.org) -------------

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the Free
 Software Foundation; either version 2 of the License, or (at your option) any
 later version.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 details.

 You should have received a copy of the GNU Lesser General Public License along
 with this program; if not, write to the Free Software Foundation, Inc., 59
 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

 Further information about the GNU Lesser General Public License can also be
 found on the world wide web at http://www.gnu.org.

FUNCTIONAL DESCRIPTION
--------------------------------------------------------------------------------
This class encapsulates the calculation of the derivatives of the state vectors
UVW and PQR - the translational and rotational rates relative to the planet
fixed frame. The derivatives relative to the inertial frame are also calculated
as a side effect. Also, the derivative of the attitude quaterion is also
calculated.

HISTORY
--------------------------------------------------------------------------------
07/12/11   JSB   Created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMMENTS, REFERENCES,  and NOTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[1] Stevens and Lewis, "Aircraft Control and Simulation", Second edition (2004)
    Wiley
[2] Richard E. McFarland, "A Standard Kinematic Model for Flight Simulation at
    NASA-Ames", NASA CR-2497, January 1975
[3] Erin Catto, "Iterative Dynamics with Temporal Coherence", February 22, 2005
[4] Mark Harris and Robert Lyle, "Spacecraft Gravitational Torques",
    NASA SP-8024, May 1969

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INCLUDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <algorithm>
#include <limits>
#include "FGAccelerations.h"
#include "FGFDMExec.h"
#include "models/FGPropagate.h"
#include "input_output/FGLog.h"

using namespace std;

namespace JSBSim {

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLASS IMPLEMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

FGAccelerations::FGAccelerations(FGFDMExec* fdmex)
  : FGModel(fdmex, "FGAccelerations")
{
  Debug(0);
  gravTorque = false;

  vPQRidot.InitMatrix();
  vUVWidot.InitMatrix();
  vUVWdot.InitMatrix();
  vBodyAccel.InitMatrix();

  bind();
  Debug(0);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FGAccelerations::~FGAccelerations(void)
{
  Debug(1);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool FGAccelerations::InitModel(void)
{
  if (!FGModel::InitModel()) return false;

  vPQRidot.InitMatrix();
  vUVWidot.InitMatrix();
  vUVWdot.InitMatrix();
  vBodyAccel.InitMatrix();

  return true;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
Purpose: Called on a schedule to calculate derivatives.
*/

bool FGAccelerations::Run(bool Holding)
{
  if (FGModel::Run(Holding)) return true;  // Fast return if we have nothing to do ...
  if (Holding) return false;

  CalculatePQRdot();   // Angular rate derivative
  CalculateUVWdot();   // Translational rate derivative

  if (!FDMExec->GetHoldDown())
    CalculateFrictionForces(in.DeltaT * rate);  // Update rate derivatives with friction forces

  Debug(2);
  return false;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Compute body frame rotational accelerations based on the current body moments
//
// vPQRdot is the derivative of the absolute angular velocity of the vehicle
// (body rate with respect to the ECEF frame), expressed in the body frame,
// where the derivative is taken in the body frame.
// J is the inertia matrix
// Jinv is the inverse inertia matrix
// vMoments is the moment vector in the body frame
// in.vPQRi is the total inertial angular velocity of the vehicle
// expressed in the body frame.
// Reference: See Stevens and Lewis, "Aircraft Control and Simulation",
//            Second edition (2004), eqn 1.5-16e (page 50)

void FGAccelerations::CalculatePQRdot(void)
{
  if (gravTorque) {
    // Compute the gravitational torque
    // Reference: See Harris and Lyle "Spacecraft Gravitational Torques",
    //            NASA SP-8024 (1969) eqn (2) (page 7)
    FGColumnVector3 R = in.Ti2b * in.vInertialPosition;
    double invRadius = 1.0 / R.Magnitude();
    R *= invRadius;
    in.Moment += (3.0 * in.vGravAccel.Magnitude() * invRadius) * (R * (in.J * R));
  }

  // Compute body frame rotational accelerations based on the current body
  // moments and the total inertial angular velocity expressed in the body
  // frame.
//  if (HoldDown && !FDMExec->GetTrimStatus()) {
  if (FDMExec->GetHoldDown()) {
    // The rotational acceleration in ECI is calculated so that the rotational
    // acceleration is zero in the body frame.
    vPQRdot.InitMatrix();
    vPQRidot = vPQRdot - in.vPQRi * (in.Ti2b * in.vOmegaPlanet);
  }
  else {
    vPQRidot = in.Jinv * (in.Moment - in.vPQRi * (in.J * in.vPQRi));
    vPQRdot = vPQRidot + in.vPQRi * (in.Ti2b * in.vOmegaPlanet);
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This set of calculations results in the body and inertial frame accelerations
// being computed.
// Compute body and inertial frames accelerations based on the current body
// forces including centripetal and Coriolis accelerations for the former.
// in.vOmegaPlanet is the Earth angular rate - expressed in the inertial frame -
//   so it has to be transformed to the body frame. More completely,
//   in.vOmegaPlanet is the rate of the ECEF frame relative to the Inertial
//   frame (ECI), expressed in the Inertial frame.
// in.Force is the total force on the vehicle in the body frame.
// in.vPQR is the vehicle body rate relative to the ECEF frame, expressed
//   in the body frame.
// in.vUVW is the vehicle velocity relative to the ECEF frame, expressed
//   in the body frame.
// Reference: See Stevens and Lewis, "Aircraft Control and Simulation",
//            Second edition (2004), eqns 1.5-13 (pg 48) and 1.5-16d (page 50)

void FGAccelerations::CalculateUVWdot(void)
{
  if (FDMExec->GetHoldDown() && !FDMExec->GetTrimStatus())
    vBodyAccel.InitMatrix();
  else
    vBodyAccel = in.Force / in.Mass;

  vUVWdot = vBodyAccel - (in.vPQR + 2.0 * (in.Ti2b * in.vOmegaPlanet)) * in.vUVW;

  // Include Centripetal acceleration.
  vUVWdot -= in.Ti2b * (in.vOmegaPlanet * (in.vOmegaPlanet * in.vInertialPosition));

  if (FDMExec->GetHoldDown()) {
    // The acceleration in ECI is calculated so that the acceleration is zero
    // in the body frame.
    vUVWidot = in.vOmegaPlanet * (in.vOmegaPlanet * in.vInertialPosition);
    vUVWdot.InitMatrix();
  }
  else {
    vUVWdot += in.Tec2b * in.vGravAccel;
    vUVWidot = in.Tb2i * vBodyAccel + in.Tec2i * in.vGravAccel;
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void FGAccelerations::SetHoldDown(bool hd)
{
  if (hd) {
    vUVWidot = in.vOmegaPlanet * (in.vOmegaPlanet * in.vInertialPosition);
    vUVWdot.InitMatrix();
    vPQRdot.InitMatrix();
    vPQRidot = vPQRdot - in.vPQRi * (in.Ti2b * in.vOmegaPlanet);
    // No row is solved held down, and none asks the position a move: not the
    // last free step's either, which the next step's Propagate would take.
    vPseudoUVWi.InitMatrix();
    vPseudoPQRi.InitMatrix();
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Computes the contact forces just before integrating the EOM.
// This routine is using Lagrange multipliers and the projected Gauss-Seidel
// (PGS) method.
// Reference: See Erin Catto, "Iterative Dynamics with Temporal Coherence",
//            February 22, 2005
// In JSBSim there is only one rigid body (the aircraft) and there can be
// multiple points of contact between the aircraft and the ground. As a
// consequence our matrix Jac*M^-1*Jac^T is not sparse and the algorithm
// described in Catto's paper has been adapted accordingly.
// The friction forces are resolved in the body frame relative to the origin
// (Earth center).
// The linear strut dampers are solved here too, as soft constraints along the
// ground normal, so that the friction bounds follow the solved normal force.
// A row may push another body besides the aircraft (LagrangeMultiplier::Other),
// which gives to it too, and hold the contact at that body's face speed; or two
// bodies (Other2), as a body resting on another does.
// The forces aim at one step of the scheme that integrates each rate
// (FGPropagate::Integrate): a multistep scheme takes this step's derivative at
// its weight and adds the past ones', so a force aimed at an Euler step came
// back the next step against itself. The other bodies are moved by Euler's.

namespace {
// The projected Gauss-Seidel's caps, for both solves: as many sweeps as it
// takes, up to what a chain of a heavy aircraft on a light box on a plank
// needs (it passes a push down the chain only as fast as its lightest body
// takes it, and 50 left the box holding a third of the aircraft's weight); a
// solve that converges stops once no row's speed moves by more in a sweep,
// ft/s^2. Rows on one face are redundant and their values drift along what
// moves nothing, which a sum of the values' changes took for motion through all
// 500 sweeps of a parked aircraft. 1E-6 balances the ground reactions to 1E-7
// of the weight as TestGndReactions asks; 1E-5 does not.
constexpr int PGSMaxIterations = 500;
constexpr double PGSTolerance = 1E-6;

// How the bodies rows i and j both push move them together: the dot products
// of their Give on each such body.
double SharedGive(const LagrangeMultiplier& i, const LagrangeMultiplier& j)
{
  const std::uint64_t iOthers[2] = {i.Other, i.Other2};
  const std::uint64_t jOthers[2] = {j.Other, j.Other2};
  const double* iGives[2] = {i.Give, i.Give2};
  const double* jGives[2] = {j.Give, j.Give2};
  double give = 0.;
  for (int p = 0; p < 2; p++)
    for (int q = 0; q < 2; q++)
      if (iOthers[p] && iOthers[p] == jOthers[q])
        for (int k = 0; k < 6; k++)
          give += iGives[p][k]*jGives[q][k];
  return give;
}

// A row as the projected Gauss-Seidel takes it, the rows' matrix never formed:
// how a unit of its value moves its own speed through the aircraft and the
// bodies it names.
struct FrictionRow {
  LagrangeMultiplier M; // its value is where the solve starts, and its answer
  FGColumnVector3 RxU;  // LeverArm x ForceJacobian, ft
  FGColumnVector3 V1;   // the aircraft's acceleration per lbs of value, 1/slug
  FGColumnVector3 V2;   // its angular acceleration per lbs of value, rad/(slug ft)
  double Diag = 0.;     // the row's own speed per lbs of its value, 1/slug
  double Rhs = 0.;      // how much its speed has to change, ft/s^2
  int Normal = -1;      // the row of these whose push bounds it, none if -1
};

// Row m, the aircraft's response to it weighted as its schemes take a step's
// derivatives.
FrictionRow RowOf(const LagrangeMultiplier& m, double mass, const FGMatrix33& Jinv,
                  double translationWeight, double rotationWeight)
{
  FrictionRow row{m};
  row.RxU = m.LeverArm * m.ForceJacobian;
  row.V1 = translationWeight * m.ForceJacobian / mass;
  row.V2 = rotationWeight * (Jinv * row.RxU); // Should be J^-T but J is symmetric and so is J^-1
  row.Diag = DotProduct(row.V1, m.ForceJacobian) + DotProduct(row.V2, row.RxU) + SharedGive(m, m);
  return row;
}

// The bodies a row names besides the aircraft (LagrangeMultiplier::Other and
// Other2), and the coordinates of a body's give.
constexpr size_t RowBodies = 2;
constexpr size_t Dofs = std::extent_v<decltype(LagrangeMultiplier::Give)>;

// The projected Gauss-Seidel over the pushes the rows' values add up to, the
// aircraft's and each body's, in place of their matrix: a sweep costs as the
// rows do, not as their square.
void SolveFrictionRows(vector<FrictionRow>& rows, double dt)
{
  const size_t n = rows.size();
  vector<std::uint64_t> bodies;
  vector<int> bodyAt(RowBodies*n, -1);
  for (size_t i=0; i < n; i++) {
    const std::uint64_t others[RowBodies] = {rows[i].M.Other, rows[i].M.Other2};
    for (size_t p = 0; p < RowBodies; p++) {
      if (!others[p]) continue;
      const auto known = find(bodies.begin(), bodies.end(), others[p]);
      bodyAt[RowBodies*i+p] = static_cast<int>(known - bodies.begin());
      if (known == bodies.end()) bodies.push_back(others[p]);
    }
  }
  // The aircraft's push, lbs and lbs ft, and each body's in its mass-weighted
  // coordinates (LagrangeMultiplier::Give).
  FGColumnVector3 force, moment;
  vector<double> pushed(Dofs*bodies.size(), 0.);
  const auto push = [&](size_t i, double by) {
    const FrictionRow& row = rows[i];
    force += by * row.M.ForceJacobian;
    moment += by * row.RxU;
    const double* gives[RowBodies] = {row.M.Give, row.M.Give2};
    for (size_t p = 0; p < RowBodies; p++)
      if (bodyAt[RowBodies*i+p] >= 0)
        for (size_t k = 0; k < Dofs; k++) pushed[Dofs*bodyAt[RowBodies*i+p]+k] += by * gives[p][k];
  };
  // Row i's speed: its line of the matrix times the values.
  const auto speed = [&](size_t i) {
    const FrictionRow& row = rows[i];
    double s = DotProduct(row.V1, force) + DotProduct(row.V2, moment);
    const double* gives[RowBodies] = {row.M.Give, row.M.Give2};
    for (size_t p = 0; p < RowBodies; p++)
      if (bodyAt[RowBodies*i+p] >= 0)
        for (size_t k = 0; k < Dofs; k++) s += gives[p][k] * pushed[Dofs*bodyAt[RowBodies*i+p]+k];
    return s;
  };
  for (size_t i=0; i < n; i++) push(i, rows[i].M.value);

  vector<double> speeds(n);
  for (size_t i=0; i < n; i++) speeds[i] = speed(i);
  for (int iter=0; iter < PGSMaxIterations; iter++) {
    for (size_t i=0; i < n; i++) {
      FrictionRow& row = rows[i];
      const double lambda0 = row.M.value;
      double lambda = lambda0 + (row.Rhs - speed(i)) / row.Diag;
      // A soft row adds Compliance/dt to its diagonal (Catto, "Soft
      // Constraints", GDC 2011); the branch follows the sign of the solution.
      const double compliance = row.M.Compliance[lambda < 0.0];
      if (compliance > 0.0)
        lambda *= dt*row.Diag / (dt*row.Diag + compliance);
      const double N = row.Normal >= 0 ? rows[row.Normal].M.NormalBase + rows[row.Normal].M.value : 1.0;
      row.M.value = FGJSBBase::Constrain(row.M.Min*N, lambda, row.M.Max*N);
      push(i, row.M.value - lambda0);
    }
    double moved = 0.;
    for (size_t i=0; i < n; i++) {
      const double s = speed(i);
      moved = max(moved, fabs(s - speeds[i]));
      speeds[i] = s;
    }
    if (moved < PGSTolerance) break;
  }
}
}

void FGAccelerations::CalculateFrictionForces(double dt)
{
  vector<LagrangeMultiplier*>& multipliers = *in.MultipliersList;
  size_t n = multipliers.size();

  vFrictionForces.InitMatrix();
  vFrictionMoments.InitMatrix();
  vPseudoUVWi.InitMatrix();
  vPseudoPQRi.InitMatrix();

  // If no gears are in contact with the ground then return
  if (!n) return;

  // With no step, as in a trim, the forces only balance the derivatives.
  const FGPropagate::StepWeights translation = FGPropagate::WeightsOf(dt > 0. ? in.TranslationalRateIntegrator : 0);
  const FGPropagate::StepWeights rotation = FGPropagate::WeightsOf(dt > 0. ? in.RotationalRateIntegrator : 0);

  FGColumnVector3 vdot = vUVWdot;
  FGColumnVector3 wdot = vPQRdot;
  if (dt > 0.) {
    // The multistep part of each step on the inertial derivatives alone: the
    // frame's turning is in the state, which no scheme integrates.
    FGColumnVector3 uvwidot = (translation.Weight - 1.) * vUVWidot;
    FGColumnVector3 pqridot = (rotation.Weight - 1.) * vPQRidot;
    static_assert(std::extent_v<decltype(Inputs::vUVWidotPast)> == std::extent_v<decltype(FGPropagate::StepWeights::Past)>,
                  "A past derivative for each the schemes take");
    for (size_t k = 0; k < std::size(translation.Past); k++) {
      uvwidot += translation.Past[k] * in.vUVWidotPast[k];
      pqridot += rotation.Past[k] * in.vPQRidotPast[k];
    }
    // Zeroes out the relative movement between the aircraft and the ground
    vdot += in.Ti2b * uvwidot + (in.vUVW - in.Tec2b * in.TerrainVelocity) / dt;
    wdot += pqridot + (in.vPQR - in.Tec2b * in.TerrainAngularVel) / dt;
  }

  vector<FrictionRow> rows;
  rows.reserve(n);
  for (size_t i=0; i < n; i++) {
    const LagrangeMultiplier* m = multipliers[i];
    FrictionRow& row = rows.emplace_back(RowOf(*m, in.Mass, in.Jinv, translation.Weight, rotation.Weight));
    // Not a number too: a body's give or the aircraft's mass gone bad.
    if (!(row.Diag > 0.)) {
      LogException err;
      err << "A friction row moves neither the aircraft nor a body it names, or not by a number!";
      throw err;
    }
    row.Rhs = -DotProduct(m->ForceJacobian, vdot + wdot*m->LeverArm);
    if (dt > 0.) // Holds the contact at the speed of the face it pushes
      row.Rhs += m->Bias/dt;
    if (m->Normal) {
      // Looked for back from the row first: a gear and a host list a normal
      // just before its rows.
      const auto back = find(make_reverse_iterator(multipliers.begin() + i), multipliers.rend(), m->Normal);
      const auto normal = back != multipliers.rend() ? back.base() - 1 : find(multipliers.begin() + i, multipliers.end(), m->Normal);
      if (normal != multipliers.end())
        row.Normal = static_cast<int>(normal - multipliers.begin());
      else {
        // A normal solved elsewhere, as a strut's: its push bounds the row as it stands.
        row.M.Min *= m->Normal->NormalBase + m->Normal->value;
        row.M.Max *= m->Normal->NormalBase + m->Normal->value;
      }
    }
  }

  SolveFrictionRows(rows, dt);

  // Calculate the total friction forces and moments

  for (size_t i=0; i < n; i++) {
    multipliers[i]->value = rows[i].M.value;
    vFrictionForces += rows[i].M.value * rows[i].M.ForceJacobian;
    vFrictionMoments += rows[i].M.value * rows[i].RxU;
  }

  FGColumnVector3 accel = vFrictionForces / in.Mass;
  FGColumnVector3 omegadot = in.Jinv * vFrictionMoments;

  vBodyAccel += accel;
  vUVWdot += accel;
  vUVWidot += in.Tb2i * accel;
  vPQRdot += omegadot;
  vPQRidot += omegadot;

  if (dt > 0.) SolvePositionBias(dt);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// The rows' PositionBias solved alone, from rest, on the rows that only push:
// those asked come out, the others are not pushed into. A depth worked off by
// the rows' own speed was kept once they parted: a rod left a wall at 0.29 of
// the speed it met it.
void FGAccelerations::SolvePositionBias(double dt)
{
  vector<const LagrangeMultiplier*> pushers;
  bool asked = false;
  for (const LagrangeMultiplier* m : *in.MultipliersList) {
    const bool pushes = !m->Normal && m->Min == 0. && m->Max == numeric_limits<double>::max() && m->ForceJacobian.Magnitude() > 0.;
    // A move asked of a row this solve does not take, or by no number.
    if (!(m->PositionBias >= 0.) || (m->PositionBias > 0. && !pushes)) {
      LogException err;
      err << "A friction row asks a move of the position that only a push out of a face can take, or not by a number!";
      throw err;
    }
    if (!pushes) continue;
    pushers.push_back(m);
    asked = asked || m->PositionBias > 0.;
  }
  if (!asked) return;

  vector<FrictionRow> rows;
  rows.reserve(pushers.size());
  for (const LagrangeMultiplier* m : pushers) {
    // From rest, hard, the position taking no scheme's weight.
    FrictionRow& row = rows.emplace_back(RowOf(*m, in.Mass, in.Jinv, 1., 1.));
    row.M.value = 0.;
    row.M.Compliance[0] = row.M.Compliance[1] = 0.;
    row.Rhs = m->PositionBias/dt;
  }

  SolveFrictionRows(rows, dt);

  FGColumnVector3 force, moment;
  for (const FrictionRow& row : rows) {
    force += row.M.value * row.M.ForceJacobian;
    moment += row.M.value * row.RxU;
  }
  vPseudoUVWi = in.Tb2i * (dt * force / in.Mass);
  vPseudoPQRi = dt * (in.Jinv * moment);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void FGAccelerations::InitializeDerivatives(void)
{
  // Make an initial run and set past values
  CalculatePQRdot();           // Angular rate derivative
  CalculateUVWdot();           // Translational rate derivative
  CalculateFrictionForces(0.);   // Update rate derivatives with friction forces
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void FGAccelerations::bind(void)
{
  PropertyManager->Tie("accelerations/pdot-rad_sec2", this, eP, &FGAccelerations::GetPQRdot);
  PropertyManager->Tie("accelerations/qdot-rad_sec2", this, eQ, &FGAccelerations::GetPQRdot);
  PropertyManager->Tie("accelerations/rdot-rad_sec2", this, eR, &FGAccelerations::GetPQRdot);

  PropertyManager->Tie("accelerations/pidot-rad_sec2", this, eP, &FGAccelerations::GetPQRidot);
  PropertyManager->Tie("accelerations/qidot-rad_sec2", this, eQ, &FGAccelerations::GetPQRidot);
  PropertyManager->Tie("accelerations/ridot-rad_sec2", this, eR, &FGAccelerations::GetPQRidot);

  PropertyManager->Tie("accelerations/udot-ft_sec2", this, eU, &FGAccelerations::GetUVWdot);
  PropertyManager->Tie("accelerations/vdot-ft_sec2", this, eV, &FGAccelerations::GetUVWdot);
  PropertyManager->Tie("accelerations/wdot-ft_sec2", this, eW, &FGAccelerations::GetUVWdot);

  PropertyManager->Tie("accelerations/uidot-ft_sec2", this, eU, &FGAccelerations::GetUVWidot);
  PropertyManager->Tie("accelerations/vidot-ft_sec2", this, eV, &FGAccelerations::GetUVWidot);
  PropertyManager->Tie("accelerations/widot-ft_sec2", this, eW, &FGAccelerations::GetUVWidot);

  PropertyManager->Tie("accelerations/gravity-ft_sec2", this, &FGAccelerations::GetGravAccelMagnitude);
  PropertyManager->Tie("simulation/gravitational-torque", &gravTorque);
  PropertyManager->Tie("forces/fbx-weight-lbs", this, eX, &FGAccelerations::GetWeight);
  PropertyManager->Tie("forces/fby-weight-lbs", this, eY, &FGAccelerations::GetWeight);
  PropertyManager->Tie("forces/fbz-weight-lbs", this, eZ, &FGAccelerations::GetWeight);

  PropertyManager->Tie("forces/fbx-total-lbs", this, eX, &FGAccelerations::GetForces);
  PropertyManager->Tie("forces/fby-total-lbs", this, eY, &FGAccelerations::GetForces);
  PropertyManager->Tie("forces/fbz-total-lbs", this, eZ, &FGAccelerations::GetForces);
  PropertyManager->Tie("moments/l-total-lbsft", this, eL, &FGAccelerations::GetMoments);
  PropertyManager->Tie("moments/m-total-lbsft", this, eM, &FGAccelerations::GetMoments);
  PropertyManager->Tie("moments/n-total-lbsft", this, eN, &FGAccelerations::GetMoments);

  PropertyManager->Tie("moments/l-gear-lbsft", this, eL, &FGAccelerations::GetGroundMoments);
  PropertyManager->Tie("moments/m-gear-lbsft", this, eM, &FGAccelerations::GetGroundMoments);
  PropertyManager->Tie("moments/n-gear-lbsft", this, eN, &FGAccelerations::GetGroundMoments);
  PropertyManager->Tie("forces/fbx-gear-lbs", this, eX, &FGAccelerations::GetGroundForces);
  PropertyManager->Tie("forces/fby-gear-lbs", this, eY, &FGAccelerations::GetGroundForces);
  PropertyManager->Tie("forces/fbz-gear-lbs", this, eZ, &FGAccelerations::GetGroundForces);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    The bitmasked value choices are as follows:
//    unset: In this case (the default) JSBSim would only print
//       out the normally expected messages, essentially echoing
//       the config files as they are read. If the environment
//       variable is not set, debug_lvl is set to 1 internally
//    0: This requests JSBSim not to output any messages
//       whatsoever.
//    1: This value explicity requests the normal JSBSim
//       startup messages
//    2: This value asks for a message to be printed out when
//       a class is instantiated
//    4: When this value is set, a message is displayed when a
//       FGModel object executes its Run() method
//    8: When this value is set, various runtime state variables
//       are printed out periodically
//    16: When set various parameters are sanity checked and
//       a message is printed out when they go out of bounds

void FGAccelerations::Debug(int from)
{
  if (debug_lvl <= 0) return;

  if (debug_lvl & 1) { // Standard console startup message output
    if (from == 0) { // Constructor

    }
  }
  if (debug_lvl & 2 ) { // Instantiation/Destruction notification
    FGLogging log(LogLevel::DEBUG);
    if (from == 0) log << "Instantiated: FGAccelerations\n";
    if (from == 1) log << "Destroyed:    FGAccelerations\n";
  }
  if (debug_lvl & 4 ) { // Run() method entry print for FGModel-derived objects
  }
  if (debug_lvl & 8 && from == 2) { // Runtime state variables
  }
  if (debug_lvl & 16) { // Sanity checking
  }
  if (debug_lvl & 64) {
    if (from == 0) { // Constructor
    }
  }
}
}
