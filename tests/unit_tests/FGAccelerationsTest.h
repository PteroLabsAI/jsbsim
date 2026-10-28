#include <cxxtest/TestSuite.h>

#include <FGFDMExec.h>
#include <models/FGAccelerations.h>
#include <models/FGPropagate.h>

using namespace JSBSim;

// The aircraft's mass, the speed it slides at and the time step.
constexpr double FGAccelerationsTest_MassSlug = 1.0;
constexpr double FGAccelerationsTest_SlideFtps = 1.0;
constexpr double FGAccelerationsTest_Dt = 1.0 / 120.0;

class FGAccelerationsTest : public CxxTest::TestSuite
{
public:
  ~FGAccelerationsTest() {
    // Avoid constructing `FGLogging` instances in a static instance of
    // FGAccelerationsTest. This is due to the destruction of thread local
    // storage *before* the destruction of static globals.
    JSBSim::FGJSBBase::debug_lvl = 0;
  }

  // A friction row along X at the CG stops an aircraft sliding over the
  // ground; on a body as heavy as it, which gives as much, the two share the
  // stop and the row pushes half as hard; on a face moving along X it holds
  // the aircraft at the face's speed.
  void testRowOnAGivingBody() {
    TS_ASSERT_DELTA(Solve(0, 0.0, 0.0),
                    -FGAccelerationsTest_MassSlug*FGAccelerationsTest_SlideFtps/FGAccelerationsTest_Dt, 1E-9);
    const double bodyGive = 1.0/sqrt(FGAccelerationsTest_MassSlug);
    TS_ASSERT_DELTA(Solve(1, bodyGive, 0.0),
                    -0.5*FGAccelerationsTest_MassSlug*FGAccelerationsTest_SlideFtps/FGAccelerationsTest_Dt, 1E-9);
    const double faceFtps = 0.5*FGAccelerationsTest_SlideFtps;
    TS_ASSERT_DELTA(Solve(1, bodyGive, faceFtps),
                    -0.5*FGAccelerationsTest_MassSlug*(FGAccelerationsTest_SlideFtps - faceFtps)/FGAccelerationsTest_Dt,
                    1E-9);
  }

  // A row that holds two bodies of the aircraft's mass apart at a speed moves
  // both: free, each takes half the push it would alone; the second held still
  // by a row of its own, the first takes all of it, which the old coupling,
  // blind to the second body, gave either way.
  void testRowBetweenTwoBodies() {
    const double give = 1.0/sqrt(FGAccelerationsTest_MassSlug);
    const double alone = FGAccelerationsTest_MassSlug*FGAccelerationsTest_SlideFtps/FGAccelerationsTest_Dt;
    LagrangeMultiplier between{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, -1E9, 1E9, 0.0};
    between.Other = 1;
    between.Give[0] = give;
    between.Other2 = 2;
    between.Give2[0] = -give;
    between.Bias = FGAccelerationsTest_SlideFtps;
    Slide({&between}, FGPropagate::eRectEuler, 1, 0.0);
    TS_ASSERT_DELTA(between.value, 0.5*alone, 1E-9);

    LagrangeMultiplier held{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, -1E9, 1E9, 0.0};
    held.Other = 2;
    held.Give[0] = give;
    between.value = 0.0;
    Slide({&between, &held}, FGPropagate::eRectEuler, 1, 0.0);
    TS_ASSERT_DELTA(between.value, alone, 1E-4);
    TS_ASSERT_DELTA(held.value, alone, 1E-4);
  }

  // Under Adams-Bashforth 2 translational rates a row stops a slide in the
  // step that scheme takes, 3/2 of this step's derivative less half the last:
  // two thirds of Euler's push first, then 8/9 of it with the first's half
  // taken back. Aimed at an Euler step, it pushed alternately too hard and too
  // little.
  void testRowAimsAtTheSchemesStep() {
    const double euler = -FGAccelerationsTest_MassSlug*FGAccelerationsTest_SlideFtps/FGAccelerationsTest_Dt;
    LagrangeMultiplier row{{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, -1E9, 1E9, 0.0};
    std::vector<double> pushes = Slide({&row}, FGPropagate::eAdamsBashforth2, 2, FGAccelerationsTest_SlideFtps);
    TS_ASSERT_DELTA(pushes[0], 2.0/3.0*euler, 1E-6);
    TS_ASSERT_DELTA(pushes[1], 8.0/9.0*euler, 1E-6);
  }

  // A trapezoidal step takes this step's derivative and the last's at half
  // each: aimed at it exactly, a row pushes twice Euler's, then none, then twice
  // again. It is aimed at Euler's step, which settles.
  void testTrapezoidalRatesAimAtEulersStep() {
    const double euler = -FGAccelerationsTest_MassSlug*FGAccelerationsTest_SlideFtps/FGAccelerationsTest_Dt;
    LagrangeMultiplier row{{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, -1E9, 1E9, 0.0};
    for (double push : Slide({&row}, FGPropagate::eTrapezoidal, 3, FGAccelerationsTest_SlideFtps))
      TS_ASSERT_DELTA(push, euler, 1E-6);
  }

private:
  // The row's push, lbs, on the body other giving bodyGive along X at its
  // centre of mass, its face moving at faceFtps.
  double Solve(std::uint64_t other, double bodyGive, double faceFtps) {
    LagrangeMultiplier row{{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, -1E9, 1E9, 0.0};
    row.Other = other;
    row.Give[0] = bodyGive;
    row.Bias = faceFtps;
    return Slide({&row}, FGPropagate::eRectEuler, 1, FGAccelerationsTest_SlideFtps)[0];
  }

  // Each step's push of rows' first on an aircraft sliding along X at
  // slideFtps, its translational rates integrated by scheme, solved steps
  // times as if each step left it as it was.
  std::vector<double> Slide(const std::vector<LagrangeMultiplier*>& rows, int scheme, int steps, double slideFtps) {
    FGFDMExec fdmex;
    auto accelerations = fdmex.GetAccelerations();
    std::vector<LagrangeMultiplier*> list = rows;
    FGAccelerations::Inputs& in = accelerations->in;
    in.J = in.Jinv = in.Ti2b = in.Tb2i = in.Tec2b = in.Tec2i = FGMatrix33(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    in.vUVW = FGColumnVector3(slideFtps, 0.0, 0.0);
    in.Mass = FGAccelerationsTest_MassSlug;
    in.DeltaT = FGAccelerationsTest_Dt;
    in.MultipliersList = &list;
    in.TranslationalRateIntegrator = scheme;
    std::vector<double> pushes;
    for (int step = 0; step < steps; step++) {
      accelerations->Run(false);
      pushes.push_back(rows[0]->value);
    }
    return pushes;
  }
};
