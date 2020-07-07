/*
 *  include/iint/Kernel.h
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

#pragma once

#include <arb/Acb.h>
#include <typeinfo>
#include <unordered_map>
#include <complex>

namespace iint {
    // base class for integration kernels
    class Kernel {
        protected:
            bool _nocache = false;

            struct point_data {
                int k0;
                std::vector<arb::Acb> cache;
            };
            std::unordered_map<arb::Acb,point_data> _points;
            std::vector<std::complex<double>> _singularities;
        public:
            // wrapper function to _calc, implements caching
            arb::Acb operator() (const arb::Acb &x, int k);

            // wrapper function to _init, implements caching
            int init(const arb::Acb &x);

            // provide readable string (should be replaced by derived class)
            virtual std::string str() const {
                return typeid(*this).name();
            }
        protected:
            // initialize expansion around x1, return k0
            virtual int _init(const arb::Acb &x1) = 0;

            // calculate expansion coefficient of (x-x1)^(k/2)
            virtual arb::Acb _calc(const arb::Acb &x1, int k) = 0;
    };

    inline std::ostream &operator<<(std::ostream &os, const Kernel &krn) {
        os << krn.str();
        return os;
    }
}
