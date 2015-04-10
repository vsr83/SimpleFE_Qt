/* SimpleFE_Qt - A Simple FE Solver with a graphical interface.
   Copyright (C) 2015 Ville Räisänen <vsr at vsr.name>

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MESH_FILE_H
#define MESH_FILE_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <assert.h>

class MeshFile {
 private:
  std::vector<std::string> strlist;
  std::ifstream *inputfile;
 public:
  int num_fields;
  bool error;

  std::vector< std::vector<std::string> > fields;
  std::vector<std::string> field_names;

  std::string version;
  std::string filename;

  MeshFile(const char *filestr);
  ~MeshFile();
};

#endif
