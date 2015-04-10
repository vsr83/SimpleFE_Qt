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

#ifndef MESH_H
#define MESH_H

#include "mesh_file.h"
#include "mesh_element.h"

#include <string>
#include <set>

class Mesh {
 private:
  MeshFile *meshfile;
  std::vector<std::string> * strsplit(const std::string &s, char delim);
 public:
  std::vector <Mesh_Element*> elements;
//    int foobar;

    double * nodes;
    std::set <int> physset;
  std::vector <int> lines;
  std::vector <int> triangles;

  unsigned int num_elements;
  unsigned int num_nodes, num_lines, num_triangles;

  Mesh(MeshFile *meshfile_in);
  Mesh(const Mesh &other);
  ~Mesh();
};

#endif
