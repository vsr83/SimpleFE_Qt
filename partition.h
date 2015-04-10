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

#ifndef PARTITION_H
#define PARTITION_H

#include <map>
#include "mesh.h"
#include <Eigen/Sparse>

// The partitions of the stiffness and mass matrices are computed by
// multiplying the global matrices by from both sides with appropriate sparse
// matrices, which select the appropriate rows and columns.
// The constructor for the Partition class is given the mesh object and
// a map, which maps those physical numbers to non-zero partitions, which are
// not put into the default partition.

// - map_nodephys maps the mesh nodes the physical numbers with priority given
//   to physical numbers of the lines over the physical numbers of surfaces.
// - map_physpart maps the physical numbers to partitions.

// - node_parts contains the vector of nodes of each partition.
// - parts_left[col + row*num_partitions] and parts_right[col + row*
//   num_partitions] are the sparse matrices so that the block (row, col) of
//   global matrix G can be obtained as the product:
//   parts_left[col + row*ncols]*G*parts_right[col + row*ncols]
// - parts_size[part] is the number of nodes in the partition part.

class Partition {
private:
    std::map <int, int> map_nodephys, map_physpart;
public:
    std::vector <std::vector<int> > node_parts;
    std::vector <Eigen::SparseMatrix<double> > parts_left, parts_right;
    std::vector <int> parts_size;
    int num_partitions;

    Partition(Mesh *mesh, std::map <int, int> _physmap);
    ~Partition();

    Eigen::SparseMatrix<double> part_left(int row, int col) {
        assert(row >= 0 && row < num_partitions);
        assert(col >= 0 && col < num_partitions);
        return parts_left[col + row*num_partitions];
    }
    Eigen::SparseMatrix<double> part_right(int row, int col) {
        assert(row >= 0 && row < num_partitions);
        assert(col >= 0 && col < num_partitions);
        return parts_right[col + row*num_partitions];
    }
};

#endif // PARTITION_H
