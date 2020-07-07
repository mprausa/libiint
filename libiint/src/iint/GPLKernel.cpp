/*
 *  src/iint/GPLKernel.cpp
 *
 *  Copyright (C) 2020 Mario Prausa
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iint/GPLKernel.h>

namespace iint {
    arb::Acb GPLKernel::_calc(const arb::Acb &x, int k) {
        if (x == _a) {
            return (k==-2) ? 1 : 0;
        }

        if (k%2) return 0;
        return -(_a-x).pow(-k/2-1);
    }

    int GPLKernel::_init(const arb::Acb &x) {
        return (x == _a) ? -2 : 0;
    }
}
