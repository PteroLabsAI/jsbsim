/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Module:       FGXMLFileRead.h
 Author:       Jon S. Berndt
 Date started: 02/04/07
 Purpose:      Shared base class that wraps the XML file reading logic

 ------------- Copyright (C) 2007  Jon S. Berndt (jon@jsbsim.org) -------------

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 details.

 You should have received a copy of the GNU Lesser General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 Place - Suite 330, Boston, MA  02111-1307, USA.

 Further information about the GNU Lesser General Public License can also be found on
 the world wide web at http://www.gnu.org.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INCLUDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <atomic>
#include <memory>
#include <sstream>

#include "FGXMLFileRead.h"
#include "input_output/FGLog.h"
#include "simgear/io/iostreams/sgstream.hxx"

namespace JSBSim {

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLASS IMPLEMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

static std::atomic<FileReader> hostReader{nullptr};

static void AppendTo(void* text, const char* data, size_t size)
{
  static_cast<std::string*>(text)->append(data, size);
}

void FGXMLFileRead::SetFileReader(FileReader reader)
{
  hostReader = reader;
}

bool FGXMLFileRead::FileExists(const SGPath& path)
{
  const FileReader reader = hostReader;
  return reader ? reader(path.utf8Str().c_str(), nullptr, nullptr)
                : path.exists();
}

Element* FGXMLFileRead::LoadXMLDocument(const SGPath& XML_filename,
                                        FGXMLParse& fparse, bool verbose)
{
  SGPath filename(XML_filename);
  if (filename.isNull()) {
    FGLogging log(LogLevel::ERROR);
    log << "No filename given.\n";
    return nullptr;
  }
  if (filename.extension().empty())
    filename.concat(".xml");

  std::unique_ptr<std::istream> input;
  const FileReader reader = hostReader;
  if (reader) {
    std::string text;
    if (reader(filename.utf8Str().c_str(), AppendTo, &text))
      input = std::make_unique<std::istringstream>(std::move(text));
  } else {
    auto file = std::make_unique<sg_ifstream>(filename);
    if (file->is_open())
      input = std::move(file);
  }
  if (!input) {
    if (verbose) {
      FGLogging log(LogLevel::ERROR);
      log << "Could not open file: " << filename << "\n";
    }
    return nullptr;
  }
  readXML(*input, fparse, filename.utf8Str());
  return fparse.GetDocument();
}

} // end namespace JSBSim
