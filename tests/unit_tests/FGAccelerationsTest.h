#include <cxxtest/TestSuite.h>

#include <FGFDMExec.h>
#include <models/FGAccelerations.h>

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

private:
  // The row's push, lbs, on the body other giving bodyGive along X at its
  // centre of mass, its face moving at faceFtps.
  double Solve(std::uint64_t other, double bodyGive, double faceFtps) {
    FGFDMExec fdmex;
    auto accelerations = fdmex.GetAccelerations();
    LagrangeMultiplier row{{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, -1E9, 1E9, 0.0};
    row.Other = other;
    row.Give[0] = bodyGive;
    row.Bias = faceFtps;
    std::vector<LagrangeMultiplier*> rows{&row};

    FGAccelerations::Inputs& in = accelerations->in;
    in.J = in.Jinv = in.Ti2b = in.Tb2i = in.Tec2b = in.Tec2i = FGMatrix33(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    in.vUVW = FGColumnVector3(FGAccelerationsTest_SlideFtps, 0.0, 0.0);
    in.Mass = FGAccelerationsTest_MassSlug;
    in.DeltaT = FGAccelerationsTest_Dt;
    in.MultipliersList = &rows;
    accelerations->Run(false);
    return row.value;
  }
};
