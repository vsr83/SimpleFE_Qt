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
#include "region.h"
#include <map>

Region::Region(std::string _name, double _nu, double _sigma, double _Js) {
   name = _name;
   nu = _nu;
   sigma = _sigma;
   Js = _Js;
}

Region::Region() {
   name = "";
   nu = sigma = Js = 0;
}

void
Region::disp() {
    std::cout << name << std::endl;
    std::cout << "nu=" << nu
              << ", sigma=" << sigma
              << ", Js=" << Js
              << std::endl;
}

#ifdef REGION_TEST

int
main(int argc, char ** argv) {
    Region regFe("Core", 1, 2);
    Region regCu("Copper", 2, 3);

    map <int, Region> regions;
    regions[23] = regFe;
    regions[92] = regCu;
    cout << regions.size() << endl;

    cout << "23-" << regions.count(23) << endl;
    cout << "24-" << regions.count(24) << endl;
//    regions[229834].disp();
    cout << regions.size() << endl;

    map <int, Region> tmp;
    tmp = regions;
    tmp[23].disp();
}
#endif
