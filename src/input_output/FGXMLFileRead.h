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
SENTRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FGXMLFILEREAD_HEADER_H
#define FGXMLFILEREAD_HEADER_H

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INCLUDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "FGXMLParse.h"
#include "simgear/misc/sg_path.hxx"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FORWARD DECLARATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

namespace JSBSim {

/** A host's own source of the files JSBSim reads, such as an archive in
    memory. It returns false when it has no file at utf8Path. Otherwise it
    hands the content to sink(context, data, size), in one or more pieces, or
    only answers that the file exists when sink is null. */
using FileReader = bool (*)(const char* utf8Path,
                            void (*sink)(void* context, const char* data,
                                         size_t size),
                            void* context);

class JSBSIM_API FGXMLFileRead {
public:
  FGXMLFileRead(void) {}
  ~FGXMLFileRead(void) {}

  /** Every file JSBSim reads afterwards, and every check that one exists,
      goes to reader instead of the disk; nullptr goes back to the disk. Set
      it before any model loads. */
  static void SetFileReader(FileReader reader);

  /** Whether the file exists, asked of the host's reader when one is set. */
  static bool FileExists(const SGPath& path);

  Element* LoadXMLDocument(const SGPath& XML_filename, bool verbose=true)
  {
    return LoadXMLDocument(XML_filename, file_parser, verbose);
  }

  Element* LoadXMLDocument(const SGPath& XML_filename, FGXMLParse& fparse,
                           bool verbose=true);

  void ResetParser(void) {file_parser.reset();}

private:
  FGXMLParse file_parser;
};
}
#endif
