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
            arb::Acb operator() (const arb::Acb &x, int k);
            int init(const arb::Acb &x);

            virtual std::string str() const {
                return typeid(*this).name();
            }
        protected:
            virtual int _init(const arb::Acb &x) = 0;
            virtual arb::Acb _calc(const arb::Acb &x, int k) = 0;
    };

    inline std::ostream &operator<<(std::ostream &os, const Kernel &krn) {
        os << krn.str();
        return os;
    }
}
