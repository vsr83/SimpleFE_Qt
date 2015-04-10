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

#include "assembly.h"
#include <Eigen/SparseLU>
#include <Eigen/Dense>

Assembly::Assembly(Mesh *_mesh, int _num_gauss, std::map <int, Region> &_regions) {
    mesh = _mesh;
    num_gauss = _num_gauss;
    regions = _regions;

    // Initialize the data structures used in Gaussian integration.

    switch (num_gauss) {
    case 6:
        gauss_u = gauss_tri6_u;
        gauss_v = gauss_tri6_v;
        gauss_w = gauss_tri6_w;
        break;
    default:
        std::cerr << "Number of gauss points not implemented " << num_gauss
                  << std::endl;
        exit(-1);
        break;
    }

    // Construct the Element objects from the Mesh_Element objects, which are
    // triangles.

    for (unsigned int ind_triangle=0; ind_triangle < mesh->num_triangles; ind_triangle++) {
         int ind_elem = mesh->triangles[ind_triangle];
         Mesh_Element *elem = mesh->elements[ind_elem];

         // Find the corner nodes of the triangles.
         double x[elem->nnodes], y[elem->nnodes], z[elem->nnodes];         
         for (unsigned int ind_node=0; ind_node<elem->nnodes; ind_node++) {
             int nodenum = elem->nodes[ind_node];

             // Note that the node, element and triangle numbers used by GMsh differ by
             // one from the corresponding array location:
             x[ind_node] = mesh->nodes[(nodenum-1)*3];
             y[ind_node] = mesh->nodes[(nodenum-1)*3+1];
             z[ind_node] = mesh->nodes[(nodenum-1)*3+2];
         }
        Element *element = new Element(*elem, x, y);
        elements.push_back(element);
    }

    // Tabulate values of the basis functions and their gradients at Gauss
    // points.

    int num_bf = elements[0]->nnodes;
    double val, gx, gy;

    for (int ind_gauss = 0; ind_gauss < num_gauss; ind_gauss++) {
        Eigen::MatrixXd mat_bfgrad(2, num_bf);
        Eigen::MatrixXd mat_bf    (num_bf, 1);

        for (int ind_bf = 0; ind_bf < num_bf; ind_bf++) {
            // Value and the gradient o the basis function no. ind_bf at the
            // Gauss point ind_gauss.
            elements[0]->node_value(ind_bf, gauss_u[ind_gauss], gauss_v[ind_gauss],
                                    val, gx, gy);
            mat_bfgrad(0, ind_bf)     = gx;
            mat_bfgrad(1, ind_bf)     = gy;
            mat_bf(ind_bf, 0)         = val;
        }
        mats_bfgrad.push_back(mat_bfgrad);
        mats_bf.push_back(mat_bf);
    }

    // Compute element stiffness and mass matrices with Gaussian integration.

    Eigen::MatrixXd Se(num_bf, num_bf);
    Eigen::MatrixXd Me(num_bf, num_bf);
    for (unsigned int ind_triangle=0; ind_triangle< mesh->num_triangles; ind_triangle++) {
        Element *elem = elements[ind_triangle];
        Se.setZero();
        Me.setZero();

        Eigen::MatrixXd Jterm(2, 2);
        Eigen::MatrixXd Jinv = elem->Je.inverse();
        Jterm = Jinv*Jinv.transpose();

        double detJ  = elem->Je.determinant();
        double nu    = regions[elem->physical].nu;
        double sigma = regions[elem->physical].sigma;

        for (int ind_gauss=0; ind_gauss <num_gauss; ind_gauss++) {
            Se += fabs(detJ) * gauss_w[ind_gauss] * nu
               * (mats_bfgrad[ind_gauss].transpose() * Jterm * mats_bfgrad[ind_gauss]);
            Me += fabs(detJ) * gauss_w[ind_gauss] * sigma
               * (mats_bf[ind_gauss]*mats_bf[ind_gauss].transpose());
        }
        elem->stiff = Se;
        elem->mass  = Me;
    }

    // A quick way to construct the global sparse stiffness and mass matrices
    // from the element matrices is by a vector of triplets. The triplets
    // contain the global row, column and value corresponding to each element
    // of the element matrices of every triangle.

    std::vector<Eigen::Triplet<double> > triplets_stiff;
    std::vector<Eigen::Triplet<double> > triplets_mass;
    triplets_stiff.reserve(mesh->num_triangles*num_bf*num_bf);
    triplets_mass.reserve (mesh->num_triangles*num_bf*num_bf);

    unsigned int ind_triplet = 0;
    for (unsigned int ind_triangle=0; ind_triangle<mesh->num_triangles; ind_triangle++) {
        Element *elem = elements[ind_triangle];
        Eigen::MatrixXd Se = elem->stiff;
        Eigen::MatrixXd Me = elem->mass;

        int ind_node2, ind_node1;
        for (int ind_bf2=0; ind_bf2<num_bf; ind_bf2++) {
            ind_node2 = elem->nodes[ind_bf2]-1;
            for (int ind_bf1=0; ind_bf1<num_bf; ind_bf1++) {
                ind_node1 = elem->nodes[ind_bf1]-1;
                Eigen::Triplet<double> tri_s(ind_node1, ind_node2, Se(ind_bf1, ind_bf2));
                Eigen::Triplet<double> tri_m(ind_node1, ind_node2, Me(ind_bf1, ind_bf2));

                triplets_stiff.push_back(tri_s);
                triplets_mass.push_back(tri_m);
                ind_triplet++;
            }
        }
    }

    // Construct the global sparse matrices from the triplets.

    global_stiff = new Eigen::SparseMatrix<double>(mesh->num_nodes, mesh->num_nodes);
    global_mass  = new Eigen::SparseMatrix<double>(mesh->num_nodes, mesh->num_nodes);
    global_stiff->setFromTriplets(triplets_stiff.begin(), triplets_stiff.end());
    global_mass->setFromTriplets(triplets_mass.begin(), triplets_mass.end());

    // Compute the excitation vector analytically.

    excitation = Eigen::MatrixXd(mesh->num_nodes, 1);
    excitation.setZero();
    for (unsigned int ind_triangle=0; ind_triangle < mesh->num_triangles; ind_triangle++) {
        Element *elem = elements[ind_triangle];
        double Js = regions[elem->physical].Js;

        for (unsigned int ind_node=0; ind_node<elem->nnodes; ind_node++) {
            if (ind_node < 3) {
                excitation(elem->nodes[ind_node]-1, 0) += elem->area*Js/3;
            }
        }
    }
}

Assembly::Assembly(const Assembly &other) {
    mesh = other.mesh;
    regions = other.regions;

    gauss_u = other.gauss_u;
    gauss_v = other.gauss_v;
    gauss_w = other.gauss_w;

    mats_bfgrad = other.mats_bfgrad;
    mats_bf = other.mats_bf;

    global_stiff = new Eigen::SparseMatrix<double>(*other.global_stiff);
    global_mass  = new Eigen::SparseMatrix<double>(*other.global_mass);

    excitation = other.excitation;
    num_gauss = other.num_gauss;

    for (unsigned int ind_triangle=0; ind_triangle < mesh->num_triangles; ind_triangle++) {
        Element *elem = new Element(*other.elements[ind_triangle]);
        elements.push_back(elem);
    }
}

Assembly::~Assembly() {    
    for (unsigned int ind_triangle=0; ind_triangle < mesh->num_triangles; ind_triangle++) {
        delete elements[ind_triangle];
    }
    delete global_stiff;
    delete global_mass;
}

#ifdef ASSEMBLY_TEST

int
main(int argc, char **argv) {
    double mu0 = 4*M_PI*1e-7;

    Region region_Fe ("Armature", 1/(1000*mu0), 0, 0);
    Region region_Cu ("Coil", 1/(mu0), 0, 1e7);
    Region region_Air("Air", 1/(mu0), 0, 0);
    std::map <int, Region> regions;

    regions[201] = region_Fe;
    regions[202] = region_Fe;
    regions[203] = region_Cu;
    regions[204] = region_Air;
    regions[205] = region_Air;
    regions[206] = region_Air;

    MeshFile *meshfile = new MeshFile("valve.msh");
    Mesh *mesh = new Mesh(meshfile);

    Assembly ass(mesh, 6, regions);

    delete mesh;
    delete meshfile;
}

#endif
