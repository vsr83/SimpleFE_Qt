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

#include "partition.h"
#include <algorithm>

Partition::Partition(Mesh *mesh, std::map <int, int> _physmap) {
    map_physpart   = _physmap;

    num_partitions = 0;
    for (std::map<int, int>::iterator it = _physmap.begin();
         it != _physmap.end(); ++it) {
        int physical  = it->first;
        int part = it->second;

        num_partitions = std::max(num_partitions, part+1);
    }

    // Find out the number of partitions.
/*    num_partitions = (*std::max_element(map_physpart.begin(), map_physpart.end(),
                                        map_physpart.value_comp())).second + 1;
*/
    parts_left.reserve (num_partitions*num_partitions);
    parts_right.reserve(num_partitions*num_partitions);

    // Initialize the node lists for all partitions.
    node_parts = std::vector<std::vector<int> >(num_partitions);

    // Assign all nodes to physical numbers given in elements. Lines are
    // given higher priority.

    for (unsigned int ind_triangle = 0; ind_triangle < mesh->num_triangles; ind_triangle++) {
        Mesh_Element *elem = mesh->elements[mesh->triangles[ind_triangle]];

        for (unsigned int ind_node = 0; ind_node < elem->nnodes; ind_node++) {
            unsigned int node = elem->nodes[ind_node]-1;
            map_nodephys[node] = elem->physical;
        }
    }


    for (unsigned int ind_line = 0; ind_line < mesh->num_lines; ind_line++) {
        Mesh_Element *elem = mesh->elements[mesh->lines[ind_line]];

        for (unsigned int ind_node = 0; ind_node < elem->nnodes; ind_node++) {
            unsigned int node = elem->nodes[ind_node]-1;
            map_nodephys[node] = elem->physical;
        }
    }

    // Use the composite mapping of map_nodephys and map_physpart to map nodes
    // to partitions.

    for (unsigned int ind_node = 0; ind_node < mesh->num_nodes; ind_node++) {
        int phys = map_nodephys[ind_node];
        int part = map_physpart[phys];
        std::cout << phys << " " << part << "/" << num_partitions<<std::endl;
        node_parts[part].push_back(ind_node);
    }

    // Construct the sparse matrices.

    for (int ind_part1 = 0; ind_part1 < num_partitions; ind_part1++) {
        std::vector<int> part1 = node_parts[ind_part1];
        int size_part1 = part1.size();
        parts_size.push_back(size_part1);

        for (int ind_part2 = 0; ind_part2 < num_partitions; ind_part2++) {
            std::vector<int> part2 = node_parts[ind_part2];
            int size_part2 = part2.size();

            Eigen::SparseMatrix<double> Pleft (size_part1, mesh->num_nodes);
            Eigen::SparseMatrix<double> Pright(mesh->num_nodes, size_part2);

            std::vector<Eigen::Triplet<double> > Tleft, Tright;

            for (int ind_node1 = 0; ind_node1 < size_part1; ind_node1++) {
                Eigen::Triplet<double> tleft (ind_node1, part1[ind_node1], 1);
                Tleft.push_back(tleft);
            }
            for (int ind_node2 = 0; ind_node2 < size_part2; ind_node2++) {
                Eigen::Triplet<double> tright (part2[ind_node2], ind_node2, 1);
                Tright.push_back(tright);
            }
            Pleft.setFromTriplets(Tleft.begin(), Tleft.end());
            Pright.setFromTriplets(Tright.begin(), Tright.end());
//          std::cout << "Part (" << ind_part1 << ", " << ind_part2 << ") "
//                    << size_part1 << "x" << size_part2 << std::endl;
            parts_left.push_back(Pleft);
            parts_right.push_back(Pright);
        }
    }
}

Partition::~Partition() {
}

#ifdef TEST_PARTITION
int
main(int argc, char **argv) {
    MeshFile *meshfile = new MeshFile("valve.msh");
    Mesh *mesh = new Mesh(meshfile);
    std::map <int, int> phys_partition;
    phys_partition[206]=1;
    Partition part(mesh, phys_partition);
}
#endif
