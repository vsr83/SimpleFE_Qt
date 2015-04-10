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

#include "mesh_file.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <assert.h>

MeshFile::MeshFile(const char *filestr) {
  filename = std::string(filestr);
  inputfile = new std::ifstream(filename.c_str(), std::ios_base::in);

  std::string line;
  while(std::getline(*inputfile, line, '\n')) {
    strlist.push_back(line);
  }

  int    field_ind = 0;
  bool   field_on  = false;
  std::string               field_name;
  std::vector <std::string> field_data;
  assert(strlist.size() > 1);

  for (unsigned int ind_line=0; ind_line < strlist.size(); ind_line++) {
    line = strlist[ind_line];

    if (line[0] == '$') {
      if (field_on) {
        // End of a field
        assert(line.length() > 4);
        field_on = false;
        fields.push_back(field_data);
        field_names.push_back(field_name);
        field_ind++;
      } else {
        // Start of a field
        field_on = true;
        assert(line.length() > 1);
        field_name = line.substr(1);
        field_data.clear();
      }
    } else {
      field_data.push_back(line);
      assert(field_on);
    }
  }
  num_fields = field_ind;
}

MeshFile::~MeshFile() {
  delete inputfile;
}

#ifdef MESH_FILE_TEST

int
main(int argc, char **argv) {
  char filename[] = "foo.msh";

  // Test for leakage.
  for (int ind=0; ind<1000; ind++) {
      cout << ind << endl;
      MeshFile meshfile(filename);
  }
}

#endif
