#include <filesystem>
#include <map>
#include <string>
#include <cxxtest/TestSuite.h>

#include <FGFDMExec.h>
#include <input_output/FGXMLFileRead.h>
#include <models/FGPropulsion.h>

using namespace JSBSim;

// The files of an aircraft held in memory by a host, under a directory that
// does not exist on the disk, and every path JSBSim asked the host for.
static std::map<std::string, std::string> hostFiles;
static std::map<std::string, int> hostAsked;

static bool ReadHostFile(const char* path,
                         void (*sink)(void*, const char*, size_t),
                         void* context)
{
  ++hostAsked[path];
  auto file = hostFiles.find(path);
  if (file == hostFiles.end()) return false;
  if (sink) {
    // In two pieces: the reader may hand a file over in several.
    const size_t half = file->second.size() / 2;
    sink(context, file->second.data(), half);
    sink(context, file->second.data() + half, file->second.size() - half);
  }
  return true;
}

class FGXMLFileReadTest : public CxxTest::TestSuite
{
public:
  ~FGXMLFileReadTest() {
    // Avoid constructing `FGLogging` instances in a static instance of
    // FGXMLFileReadTest. This is due to the destruction of thread local
    // storage *before* the destruction of static globals.
    JSBSim::FGJSBBase::debug_lvl = 0;
  }

  void tearDown() {
    FGXMLFileRead::SetFileReader(nullptr);
    hostFiles.clear();
    hostAsked.clear();
  }

  void testAircraftLoadsFromTheHostReader() {
    const auto dir = std::filesystem::temp_directory_path()
                     / "JSBSimFGXMLFileReadTest-not-on-disk";
    TS_ASSERT(!std::filesystem::exists(dir));
    const SGPath root(dir.string());
    const std::string aircraft = (root/"hostcraft"/"hostcraft.xml").utf8Str();
    const std::string motor = (root/"engine"/"motor.xml").utf8Str();
    const std::string nozzle = (root/"engine"/"nozzle.xml").utf8Str();
    const std::string system = (root/"systems"/"sys.xml").utf8Str();

    hostFiles[aircraft] =
      "<fdm_config name=\"hostcraft\" version=\"2.0\" release=\"ALPHA\">"
      "<metrics><wingarea unit=\"M2\">1</wingarea>"
      "<wingspan unit=\"M\">1</wingspan><chord unit=\"M\">1</chord></metrics>"
      "<mass_balance><emptywt unit=\"KG\">10</emptywt>"
      "<ixx unit=\"KG*M2\">1</ixx><iyy unit=\"KG*M2\">1</iyy>"
      "<izz unit=\"KG*M2\">1</izz>"
      "<location name=\"CG\" unit=\"M\"><x>0</x><y>0</y><z>0</z></location>"
      "</mass_balance><ground_reactions/>"
      "<propulsion><engine file=\"motor\">"
      "<location unit=\"M\"><x>0</x><y>0</y><z>0</z></location>"
      "<thruster file=\"nozzle\">"
      "<location unit=\"M\"><x>0</x><y>0</y><z>0</z></location>"
      "</thruster></engine></propulsion>"
      "<system file=\"sys\"/><aerodynamics/></fdm_config>";
    hostFiles[motor] =
      "<electric_engine name=\"motor\"><power unit=\"WATTS\">100</power>"
      "</electric_engine>";
    hostFiles[nozzle] = "<direct name=\"nozzle\"/>";
    hostFiles[system] =
      "<system name=\"sys\"><property value=\"1\">test/from-host</property>"
      "</system>";
    FGXMLFileRead::SetFileReader(ReadHostFile);

    FGFDMExec fdmex;
    TS_ASSERT(fdmex.LoadModel(root, root/"engine", root/"systems",
                              "hostcraft"));

    for (const auto& file: hostFiles)
      TS_ASSERT(hostAsked[file.first] > 0);
    TS_ASSERT_EQUALS(fdmex.GetPropulsion()->GetNumEngines(), 1u);
    auto pm = fdmex.GetPropertyManager();
    TS_ASSERT(pm->HasNode("test/from-host"));
    TS_ASSERT_EQUALS(pm->GetNode("test/from-host")->getDoubleValue(), 1.0);
  }

  void testHostReaderAnswersWhetherAFileExists() {
    const SGPath here(std::filesystem::current_path().string());
    TS_ASSERT(FGXMLFileRead::FileExists(here));

    FGXMLFileRead::SetFileReader(ReadHostFile);
    TS_ASSERT(!FGXMLFileRead::FileExists(here));
    hostFiles[here.utf8Str()] = "";
    TS_ASSERT(FGXMLFileRead::FileExists(here));

    FGXMLFileRead::SetFileReader(nullptr);
    TS_ASSERT(FGXMLFileRead::FileExists(here));
  }

  void testMissingHostFileFailsTheLoad() {
    FGXMLFileRead::SetFileReader(ReadHostFile);
    FGXMLFileRead reader;
    TS_ASSERT(reader.LoadXMLDocument(SGPath("nothing-here.xml"), false)
              == nullptr);
    TS_ASSERT_EQUALS(hostAsked["nothing-here.xml"], 1);
  }
};
