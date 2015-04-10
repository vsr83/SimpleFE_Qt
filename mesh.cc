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

#include "mesh.h"

std::vector<std::string> *
Mesh::strsplit(const std::string &s, char delim) {
  std::stringstream ss(s);
  std::string item;

  std::vector<std::string> *elems;
  elems = new std::vector<std::string>;

  while (std::getline(ss, item, delim)) {
    elems->push_back(item);
  }
  return elems;
}

Mesh::Mesh(MeshFile *meshfile_in) {
  meshfile = meshfile_in;

  assert(meshfile->num_fields == 3);
  assert(meshfile->field_names[0] == "MeshFormat");
  assert(meshfile->field_names[1] == "Nodes");
  assert(meshfile->field_names[2] == "Elements");

  std::vector <std::string> *nodefield    = &meshfile->fields[1];

  num_nodes = nodefield->size()-1;
  assert(num_nodes > 0);
  nodes = new double[num_nodes*3];

  std::string nodestr = (*nodefield)[0];
  assert(num_nodes == (unsigned int)atoi(nodestr.c_str()));

  for (unsigned int ind_node=0;ind_node < num_nodes; ind_node++) {
    std::string line;
    line = (*nodefield)[ind_node+1];
    
    std::vector <std::string> *node_strlist = strsplit(line, ' ');
    assert(node_strlist->size() == 4);

    double x = atof((*node_strlist)[1].c_str());
    double y = atof((*node_strlist)[2].c_str());
    double z = atof((*node_strlist)[3].c_str());

    nodes[ind_node*3 + 0] = x;
    nodes[ind_node*3 + 1] = y;
    nodes[ind_node*3 + 2] = z;
    delete node_strlist;
  }

  std::vector <std::string> *elementfield = &meshfile->fields[2];

  num_elements = elementfield->size()-1;
  assert(num_elements > 0);
  std::string elemstr = (*elementfield)[0];
  assert(num_elements == (unsigned int) atoi(elemstr.c_str()));

  num_lines = 0;
  num_triangles = 0;

  for (unsigned int ind_element=0;ind_element < num_elements; ind_element++) {
    std::string elementstr;
    elementstr = (*elementfield)[ind_element+1];

    Mesh_Element *elem = new Mesh_Element(elementstr);
    elements.push_back(elem);

    if (elem->type == ELEMENT_TRIANGLE3 || elem->type == ELEMENT_TRIANGLE6) {
      num_triangles++;
      triangles.push_back(ind_element);
    }
    if (elem->type == ELEMENT_LINE2 || elem->type == ELEMENT_LINE3) {
      num_lines++;
      lines.push_back(ind_element);
    }
  }
}

Mesh::Mesh(const Mesh &other) {
    num_nodes = other.num_nodes;
    num_lines = other.num_lines;
    num_triangles = other.num_triangles;
    num_elements = other.num_elements;

    meshfile = other.meshfile;
    for (unsigned int ind_element=0;ind_element < num_elements; ind_element++) {
        Mesh_Element *elem = new Mesh_Element(*other.elements[ind_element]);
        elements.push_back(elem);
    }

    nodes = new double[num_nodes*3];
    std::copy(other.nodes, other.nodes+num_nodes*3, nodes);
    lines = other.lines;
    triangles = other.triangles;
}

Mesh::~Mesh() {
  delete nodes;
  for (unsigned int ind_element=0;ind_element < num_elements; ind_element++) {
    delete elements[ind_element];
  }

}

#ifdef MESH_TEST

int
main(int argc, char **argv) {
  char filename[] = "foo.msh";

  MeshFile meshfile(filename);

  for (int ind=0; ind < 1000; ind ++){
      cout << ind << endl;
      Mesh mesh(&meshfile);
  }
  for (int ind_field = 0; ind_field < meshfile.num_fields; ind_field++) {
    // cout << "Field " << ind_field << " '" 
    //	 << meshfile.field_names[ind_field] << "'" << endl;

    std::vector<std::string> *fielddata = &meshfile.fields[ind_field];
      
  }
}

#endif
