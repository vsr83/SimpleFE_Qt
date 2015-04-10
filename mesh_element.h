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

#ifndef MESH_ELEMENT_H
#define MESH_ELEMENT_H

#include <iostream>
#include <vector>
#include <sstream>

#include <assert.h>
#include <stdlib.h>

enum {ELEMENT_LINE2         = 1,
      ELEMENT_TRIANGLE3     = 2,
      ELEMENT_QUADRANGLE4   = 3,
      ELEMENT_TETRAHEDRON4  = 4,
      ELEMENT_HEXAHEDRON8   = 5,
      ELEMENT_PRISM6        = 6,
      ELEMENT_PYRAMID5      = 7,
      ELEMENT_LINE3         = 8,
      ELEMENT_TRIANGLE6     = 9,
      ELEMENT_QUADRANGLE9   = 10, 
      ELEMENT_TETRAHEDRON10 = 11,
      ELEMENT_HEXAHEDRON27  = 12,
      ELEMENT_PRISM18       = 13,
      ELEMENT_PYRAMID14     = 14,
      ELEMENT_POINT1        = 15,
      ELEMENT_QUADRANGLE8   = 16,
      ELEMENT_HEXAHEDRON20  = 17,
      ELEMENT_PRISM15       = 18,
      ELEMENT_PYRAMID13     = 19,
      ELEMENT_TRIANGLE9     = 20,
      ELEMENT_TRIANGLE10    = 21,
      ELEMENT_TRIANGLE12    = 22,
      ELEMENT_TRIANGLE15    = 23,
      ELEMENT_TRIANGLE15_5  = 24,
      ELEMENT_TRIANGLE21    = 25,
      ELEMENT_EDGE4         = 26,
      ELEMENT_EDGE5         = 27,
      ELEMENT_EDGE6         = 28,
      ELEMENT_TETRAHEDRON20 = 29,
      ELEMENT_TETRAHEDRON35 = 30,
      ELEMENT_TETRAHEDRON56 = 31,
      //      ELEMENT_TETRAHEDRON64 = 92,
      // ELEMENT_TETRAHEDRON125= 93
};

struct elementdata {
  unsigned int type;
  const char *name;
  unsigned int num_nodes;
  unsigned int order;
  unsigned int hnode, hedges, hfaces, hvolumes;
} static const ElementData[] = {
  0                       , "NULL_ELEMENT",        0,  0,  0, 0, 0, 0,
    ELEMENT_LINE2         , "2-node line",         2,  1,  2, 0, 0, 0,
    ELEMENT_TRIANGLE3     , "3-node triangle",     3,  1,  2, 0, 0, 0,
    ELEMENT_QUADRANGLE4   , "4-node quadrangle",   4,  1,  2, 0, 0, 0,
    ELEMENT_TETRAHEDRON4  , "4-node tetrahedron",  4 , 1,  2, 0, 0, 0,
    ELEMENT_HEXAHEDRON8   , "8-node hexahedron",   8 , 1,  2, 0, 0, 0,
    ELEMENT_PRISM6        , "6-node prism",        6,  1,  2, 0, 0, 0,
    ELEMENT_PYRAMID5      , "5-node pyramid",      5,  1,  2, 0, 0, 0,
    ELEMENT_LINE3         , "3-node line ",        3,  2,  2, 1, 0, 0,
    ELEMENT_TRIANGLE6     , "6-node triangle",     6,  2,  3, 3, 0, 0,
    ELEMENT_QUADRANGLE9   , "9-node quadrangle",   9,  2,  4, 4, 1, 0,
    ELEMENT_TETRAHEDRON10 , "10-node tetrahedron", 10, 2,  4, 6, 0, 0,
    ELEMENT_HEXAHEDRON27  , "27-node hexahedron",  27, 2,  8,12, 6, 1,
    ELEMENT_PRISM18       , "18-node prism",       18, 2,  6, 9, 3, 0,
    ELEMENT_PYRAMID14     , "14-node pyramid",     14, 2,  5, 8, 1, 0,
    ELEMENT_POINT1        , "1-node point",        1,  0,  1, 0, 0, 0,
    ELEMENT_QUADRANGLE8   , "8-node quadrangle",   8,  2,  4, 4, 0, 0,
    ELEMENT_HEXAHEDRON20  , "20-node hexahedron",  20, 2,  8,12, 0, 0,
    ELEMENT_PRISM15       , "15-node prism",       15, 2,  6, 9, 0, 0,
    ELEMENT_PYRAMID13     , "13-node pyramid",     13, 2,  5, 8, 0, 0,
    ELEMENT_TRIANGLE9     , "9-node triangle",     9,  3,  3, 6, 0, 0,
    ELEMENT_TRIANGLE10    , "10-node triangle",    10, 3,  3, 6, 1, 0,
    ELEMENT_TRIANGLE12    , "12-node triangle",    12, 4,  3, 9, 0, 0,
    ELEMENT_TRIANGLE15    , "15-node triangle",    15, 4,  3, 9, 3, 0,
    ELEMENT_TRIANGLE15_5  , "15-node triangle",    15, 5,  3,12, 0, 0,
    ELEMENT_TRIANGLE21    , "21-node triangle",    21, 5,  3,12, 6, 0,
    ELEMENT_EDGE4         , "4-node edge",         4,  3,  2, 2, 0, 0,
    ELEMENT_EDGE5         , "5-node edge",         5,  4,  2, 3, 0, 0,
    ELEMENT_EDGE6         , "6-node edge",         6,  5,  2, 4, 0, 0,
    ELEMENT_TETRAHEDRON20 , "20-node tetrahedron", 20, 3,  4,12, 4, 0,
    ELEMENT_TETRAHEDRON35 , "35-node tetrahedron", 35, 4,  4,18,12, 1,
    ELEMENT_TETRAHEDRON56 , "56-node tetrahedron", 56, 5,  4,24,24, 4
    };

class Mesh_Element {
 private:
  std::vector <int> list_raw;
  std::vector<std::string> *strsplit(const std::string &s, char delim);

 public:
  void disp();

  Mesh_Element(const Mesh_Element &other);
  Mesh_Element(std::string &line);
  ~Mesh_Element();

  unsigned int nnodes;
  unsigned int number;
  unsigned int type;
  unsigned int ntags;
  unsigned int physical;
  unsigned int geometrical;

  unsigned int *tags;
  unsigned int *nodes;
};

#endif
