/*
 *  src/iint/Kernel.cpp
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

#include <iint/Kernel.h>

namespace iint {
    arb::Acb Kernel::operator() (const arb::Acb &x, int k) {
        if (_nocache) return _calc(x,k);
        assert(_points.count(x));
        auto &data = _points[x];

        if (k < data.k0) return arb::Acb(0,x.default_prec());
        if (k-data.k0 < data.cache.size()) return data.cache[k-data.k0];
        while (k-data.k0 > data.cache.size()) (*this)(x,int(data.cache.size())+data.k0);
        assert(k-data.k0 == data.cache.size());

        auto res = _calc(x,k);
        data.cache.push_back(res);
        return res;
    }

    int Kernel::init(const arb::Acb &x) {
        if (_nocache) return _init(x);

        if (_points.count(x)) return _points[x].k0;

        int k0 = _init(x);
        _points[x] = {k0,{}};
        return k0;
    }
}

