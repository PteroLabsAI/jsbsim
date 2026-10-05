#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <cxxtest/TestSuite.h>

#include <FGFDMExec.h>
#include <initialization/FGInitialCondition.h>
#include <initialization/FGTrim.h>
#include <input_output/FGGroundCallback.h>
#include <models/FGInertial.h>
#include <models/FGPropagate.h>

using namespace JSBSim;

constexpr double degtorad = M_PI / 180.;

// An exactly flat ground: the plane tangent to the ellipsoid at 'origin'. The
// default ellipsoid ground curves away under points that are coplanar in the
// body frame, which hides the rounding a flat ground exposes.
class FlatGroundCallback : public FGGroundCallback
{
public:
  explicit FlatGroundCallback(const FGLocation& origin)
    : point(origin), up(origin.GetTl2ec() * FGColumnVector3(0., 0., -1.)) {}

  double GetAGLevel(double, const FGLocation& location, FGLocation& contact,
                    FGColumnVector3& normal, FGColumnVector3& v,
                    FGColumnVector3& w) const override
  {
    double h = DotProduct(location - point, up);
    contact = location;
    contact -= h * up;
    normal = up;
    v.InitMatrix();
    w.InitMatrix();
    return h;
  }

private:
  FGLocation point;
  FGColumnVector3 up;
};

class FGTrimTest : public CxxTest::TestSuite
{
public:
  // Five feet on a ring under the CG, one of them raised: on a flat ground the
  // aircraft stands level on the four coplanar feet whatever its heading.
  void testGroundTrimOnCoplanarContacts() {
    const int feet = 5;
    const double ringRadius = 0.5; // m
    const double legLength = 1.5;  // m, from the CG down to a foot
    const double raisedBy = 1.0;   // m, how much higher the fifth foot is
    const double latitude = 45.0;  // deg
    const double longitude = 5.0;  // deg
    const double headingStep = 5.0; // deg
    const double levelTolerance = 1.0; // deg
    const double startHeight = 10.0; // ft, the feet start above the ground

    auto dir = std::filesystem::temp_directory_path() / "JSBSimFGTrimTest";
    std::filesystem::create_directories(dir / "coplanar");
    {
      std::ofstream xml(dir / "coplanar" / "coplanar.xml");
      xml << "<fdm_config name=\"coplanar\" version=\"2.0\" release=\"ALPHA\">"
             "<metrics><wingarea unit=\"M2\">1</wingarea>"
             "<wingspan unit=\"M\">1</wingspan><chord unit=\"M\">1</chord>"
             "</metrics>"
             "<mass_balance><emptywt unit=\"KG\">10</emptywt>"
             "<ixx unit=\"KG*M2\">1</ixx><iyy unit=\"KG*M2\">1</iyy>"
             "<izz unit=\"KG*M2\">1</izz>"
             "<location name=\"CG\" unit=\"M\"><x>0</x><y>0</y><z>0</z>"
             "</location></mass_balance><ground_reactions>";
      for (int k = 0; k < feet; ++k) {
        double angle = k * 360.0 / feet * degtorad;
        double z = k == 0 ? raisedBy - legLength : -legLength;
        xml << "<contact type=\"BOGEY\" name=\"FOOT" << k << "\">"
            << "<location unit=\"M\"><x>" << ringRadius * cos(angle) << "</x>"
            << "<y>" << ringRadius * sin(angle) << "</y><z>" << z << "</z>"
            << "</location><static_friction>0.8</static_friction>"
            << "<dynamic_friction>0.5</dynamic_friction>"
            << "<rolling_friction>0.02</rolling_friction>"
            << "<spring_coeff unit=\"N/M\">100000</spring_coeff>"
            << "<damping_coeff unit=\"N/M/SEC\">1000</damping_coeff></contact>";
      }
      xml << "</ground_reactions><aerodynamics/></fdm_config>";
    }

    FGFDMExec fdmex;
    SGPath path(dir.string());
    TS_ASSERT(fdmex.LoadModel(path, path, path, "coplanar"));

    auto planet = fdmex.GetInertial();
    FGLocation origin;
    origin.SetEllipse(planet->GetSemimajor(), planet->GetSemiminor());
    origin.SetPositionGeodetic(longitude * degtorad,
                               latitude * degtorad, 0.0);
    planet->SetGroundCallback(new FlatGroundCallback(origin));

    auto ic = fdmex.GetIC();
    auto propagate = fdmex.GetPropagate();
    for (double psi = 0.0; psi < 360.0; psi += headingStep) {
      // Geodetic, so that the aircraft stands where the plane touches the
      // ellipsoid and the plane is level there.
      ic->SetGeodLatitudeDegIC(latitude);
      ic->SetLongitudeDegIC(longitude);
      ic->SetAltitudeAGLFtIC(startHeight);
      ic->SetPhiDegIC(0.0);
      ic->SetThetaDegIC(0.0);
      ic->SetPsiDegIC(psi);
      TS_ASSERT(fdmex.RunIC());

      std::ostringstream heading;
      heading << "heading " << psi << " deg";
      FGTrim trim(&fdmex, tGround);
      TSM_ASSERT(heading.str(), trim.DoTrim());
      TSM_ASSERT_DELTA(heading.str(), propagate->GetEulerDeg(FGJSBBase::eTht),
                       0.0, levelTolerance);
      TSM_ASSERT_DELTA(heading.str(), propagate->GetEulerDeg(FGJSBBase::ePhi),
                       0.0, levelTolerance);
    }

    std::filesystem::remove_all(dir);
  }
};
