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

#ifndef BASISF_H
#define BASISF_H

#include "mesh_element.h"
#include <vector>
#include <string>
#include <Eigen/Dense>

class Element : public Mesh_Element {
private:
    int num_corners;
    double *corner_x, *corner_y;

    inline double det3(double a, double b, double c,
                       double d, double e, double f,
                       double g, double h, double i) {
        return a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;
    }
    inline double det2(double a, double b,
                       double c, double d) {
        return a*d - b*c;
    }
    void parse_corners(double *_x, double *_y);
public:
    Element(const Mesh_Element &elem, double *_x, double *_y);
    Element(std::string &line, double *_x, double *_y);
    Element(const Element& other);
    ~Element();

    double area; // Surface area of the element in global coordinates.
    int num_dof; // Number of nodes (and DoFs) for the element.

    Eigen::MatrixXd stiff, mass; // Element stiffness and mass matrices.

    // Jacobian matrix for the placement mapping from reference and global
    // coordinates.
    Eigen::MatrixXd Je;

    // Compute value and gradient of the basis function node on the reference
    // element at the reference coordinates (u, v).
    void node_value    (int node, double u, double v,
                        double &val, double &gx, double &gy);

    // Coordinate transformations between the coordinates (u, v) of the
    // reference element and the global coordinates (xp, yp).
    void coord_carttri (double xp, double yp, double &u,  double &v);
    void coord_tricart (double u,  double v,  double &xp, double &yp);
};

#endif // BASISF_H
