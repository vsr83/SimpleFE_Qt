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

#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "mesh.h"
#include "element.h"
#include "region.h"
#include <vector>
#include <map>
#include <Eigen/SparseLU>

static double gauss_tri6_u[6] = {0.816847572980459,0.091576213509771,0.091576213509771,
                                 0.108103018168070,0.445948490915965,0.445948490915965};
static double gauss_tri6_v[6] = {0.091576213509771,0.816847572980459,0.091576213509771,
                                 0.445948490915965,0.108103018168070,0.445948490915965};
static double gauss_tri6_w[6] = {0.054975871827661,0.054975871827661,0.054975871827661,
                                 0.111690794839,   0.111690794839,   0.111690794839};

// Assembly objects are used to store the global stiffness and mass matrices
// and all other objects needed for their construction. Since the mesh,
// elements, global_stiff, global_mass objects can be very large, pointers are
// used avoid unnecessary copying.

// The regions map maps physical numbers of the Gmsh mesh to Region objects,
// shich contain material parameters and source densities.
// The elements vector contains the Mesh_Element data structures together with
// the element stiffness and mass matrices.
// The matrices mats_bf[ind_gauss] and mats_bfgrad[ind_gauss] contain the
// values of the basis functions and their gradients at gauss point ind_gauss.

class Assembly {
    Mesh *mesh;
    std::map <int, Region> regions;
    std::vector <Element*> elements;

    const double *gauss_u, *gauss_v, *gauss_w;
    std::vector<Eigen::MatrixXd> mats_bfgrad;
    std::vector<Eigen::MatrixXd> mats_bf;
public:
    Eigen::SparseMatrix<double> *global_stiff, *global_mass;
    Eigen::MatrixXd excitation;
    int num_gauss;
    Assembly(Mesh *_mesh, int _num_gauss, std::map <int, Region> &_regions);
    Assembly(const Assembly &other);
    ~Assembly();
};

#endif // ASSEMBLY_H
