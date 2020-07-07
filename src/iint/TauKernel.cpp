/*
 *  src/iint/TauKernel.cpp
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

#include <iint/TauKernel.h>

namespace iint {
    TauKernel::TauKernel() : _reciprocal({0,-20,0,1},{2000,1700,560,85,5}) {
        _singularities = {0.,0.055728090000841214363};
    }

    int TauKernel::_init(const arb::Acb &x) {
        return -_reciprocal.init(x);
    }

    arb::Acb TauKernel::_calc(const arb::Acb &x, int k) {
        auto &data = _points[x];

        assert(k >= data.k0);
        if (k == data.k0) return -arb::Acb::Pi(x.default_prec()).pow(2) / _reciprocal(x,-data.k0);

        arb::Acb res;
        for (int n=1; n<=k-data.k0; ++n) {
            res -= _reciprocal(x,n-data.k0) * data.cache[k-n-data.k0];
        }

        res /= _reciprocal(x,-data.k0);

        //std::cout << "TauKernel -- k = " << k << " res = " << res << std::endl;
        return res;
    }
}
