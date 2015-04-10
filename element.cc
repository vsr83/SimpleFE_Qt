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

#include <iostream>
#include <vector>
#include "mesh_element.h"
#include "element.h"


void
Element::coord_carttri(double xp, double yp, double &u, double &v) {
    u = det3(1, xp,          yp,
             1, corner_x[0], corner_y[0],
             1, corner_x[2], corner_y[2]) /
        det3(1, corner_x[1], corner_y[1],
             1, corner_x[0], corner_y[0],
             1, corner_x[2], corner_y[2]);
    v = det3(1, xp,          yp,
             1, corner_x[0], corner_y[0],
             1, corner_x[1], corner_y[1]) /
        det3(1, corner_x[2], corner_y[2],
             1, corner_x[0], corner_y[0],
             1, corner_x[1], corner_y[1]);
}

void
Element::coord_tricart(double u, double v, double &xp, double &yp) {
    xp = corner_x[0] + (corner_x[1]-corner_x[0])*u
                     + (corner_x[2]-corner_x[0])*v;
    yp = corner_y[0] + (corner_y[1]-corner_y[0])*u
                     + (corner_y[2]-corner_y[1])*v;
}

void
Element::node_value(int node, double u, double v,
                    double &val, double &gx, double &gy) {
    switch (type) {
    case ELEMENT_TRIANGLE3:
        switch (node) {
        case 0:
            val = 1-u-v;
            gx = -1; gy = -1;
            break;
        case 1:
            val = u;
            gx = 1; gy = 0;
            break;
        case 2:
            val = v;
            gx = 0; gy = 1;
            break;
        default:
            std::cerr << "ELEMENT_TRIANGLE3: Invalid node " << node << std::endl;
            exit(-1);
            break;
        }
        break;
    case ELEMENT_TRIANGLE6:
        switch (node) {
        case 0:
            val = 1-u-v;
            gx = -1; gy = -1;
            break;
        case 1:
            val = u;
            gx = 1; gy = 0;
            break;
        case 2:
            val = v;
            gx = 0; gy = 1;
            break;
        case 3:
            val = (1-u-v)*u;
            gx = 1-v-2*u; gy = -u;
            break;
        case 4:
            val = u*v;
            gx = v; gy = u;
            break;
        case 5:
            val = (1-u-v)*v;
            gx = -v; gy = 1-u-2*v;
            break;
        default:
            std::cerr << "ELEMENT_TRIANGLE6: Invalid node " << node << std::endl;
            break;
        }
        break;
    default:
        std::cerr << "Basis Functions not implemented for element type  "
                  << type << std::endl;
        exit(-1);
        break;
    }
}

Element::Element(const Element &other) : Mesh_Element(other) {
    num_corners = other.num_corners;
    num_dof     = other.num_dof;

    corner_x = new double[num_corners];
    corner_y = new double[num_corners];
    std::copy(other.corner_x, other.corner_x+num_corners, corner_x);
    std::copy(other.corner_y, other.corner_y+num_corners, corner_y);
}

void
Element::parse_corners(double *_x, double *_y) {
    switch (type) {
    case ELEMENT_TRIANGLE3:
        num_corners = 3;
        num_dof = 3;
        break;
    case ELEMENT_TRIANGLE6:
        num_corners = 3;
        num_dof = 6;
        break;
    default:
        std::cerr << "Non-Implemented Element Type: " << type << std::endl;
        exit(-1);
        break;
    }
    corner_x = new double[num_corners];
    corner_y = new double[num_corners];

    for (int ind_corner=0; ind_corner<3; ind_corner++) {
        corner_x[ind_corner] = _x[ind_corner];
        corner_y[ind_corner] = _y[ind_corner];
    }

    switch(num_corners) {
    case 3:
        Je = Eigen::MatrixXd(2, 2);
        Je(0, 0) = corner_x[1]-corner_x[0];
        Je(0, 1) = corner_x[2]-corner_x[0];
        Je(1, 0) = corner_y[1]-corner_y[0];
        Je(1, 1) = corner_y[2]-corner_y[0];

        area =  0.5*fabs((corner_x[2]-corner_x[0])*(corner_y[1]-corner_y[0])
                       - (corner_y[2]-corner_y[0])*(corner_x[1]-corner_x[0]));


        break;
    default:
        std::cerr << "Non-Implemented Number of Corners: " << num_corners
                  << std::endl;
        exit(-1);
        break;
    }
}

Element::Element(const Mesh_Element &elem, double *_x, double *_y) : Mesh_Element(elem){
    parse_corners(_x, _y);
}

Element::Element(std::string &line, double *_x, double *_y) : Mesh_Element(line) {
    parse_corners(_x, _y);
}

Element::~Element() {
    delete [] corner_x;
    delete [] corner_y;
    corner_x = 0;
    corner_y = 0;
}

#ifdef BASISF_TEST

int
main(int argc, char **argv) {
    string testline = string("27568 9 2 2071 2070 584 585 50303 1727 54794 56804");
    double x[] = {0, 1, 0};
    double y[] = {0, 0, 2};

    Element me(testline, x, y);
    me.disp();
    cout << endl;

    double gx, gy, v;

    me.node_value(1, 0.2, 0.2, v, gx, gy);
    cout << v << " " << gx << " " << gy << endl;

    {
    double xx = 1, yy = 0;
    double uu, vv, xp, yp;
    me.coord_carttri(xx, yy, uu, vv);
    me.coord_tricart(uu, vv, xp, yp);
    cout << xx << ", " << yy << "  "
         << uu << ", " << vv << "  "
         << xp << ", " << yp << endl;
    }
    {
    double xx = 0, yy = 1;
    double uu, vv, xp, yp;
    me.coord_carttri(xx, yy, uu, vv);
    me.coord_tricart(uu, vv, xp, yp);
    cout << xx << ", " << yy << "  "
         << uu << ", " << vv << "  "
         << xp << ", " << yp << endl;
    }
    {
    double xx = 0.5, yy = 0.5;
    double uu, vv, xp, yp;
    me.coord_carttri(xx, yy, uu, vv);
    me.coord_tricart(uu, vv, xp, yp);
    cout << xx << ", " << yy << "  "
         << uu << ", " << vv << "  "
         << xp << ", " << yp << endl;
    }

    vector <Element> elements;
    elements.push_back(me);
    elements.push_back(me);
    elements.push_back(me);
    elements.push_back(me);
    elements.push_back(me);

}
#endif
