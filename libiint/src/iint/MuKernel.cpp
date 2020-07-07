/*
 *  src/iint/MuKernel.cpp
 *
 *  Copyright (C) 2020 Mario Prausa
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iint/MuKernel.h>

namespace {
    std::vector<int> prepend(const std::vector<int> &v, int n) {
        std::vector<int> u(n,0);
        u.insert(u.end(),v.begin(),v.end());
        return u;
    }
}

namespace iint {
    MuKernel::MuKernel(int n)
        : _n(n), EllipticKernel(
            prepend({100,45,5},n?n-1:0),
            prepend({-8000,-6400,-1680,0,84,16,1},n?0:1)) {}
}

