/*
 *  include/iint/ODE.h
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
#include <unordered_map>

namespace iint {
    // solve ordinary differential equation with a power series
    // ansatz (with half-integer powers)
    class ODE {
        protected:
            arb::Acb zero{0};
            std::vector<std::vector<int>> _p;
            int _mhat;

            struct point_data {
                long prec;
                int k0;
                int r;
                std::vector<arb::Acb> a;
                std::vector<std::vector<arb::Acb>> dp;
            };

            std::unordered_map<arb::Acb,point_data> _points;
        public:
            // The ode is defined by 0 = sum_n p_n(x) * Dx^n f(x),
            // where p_n(x) are polynomials. The coefficients of
            // p_n are provided in the std::vector<int> p[n]
            // starting at x^0.
            ODE(const std::vector<std::vector<int>> &p);

            // Initialize power series around x.
            // k0/2 denotes the power of the leading order of the expansion.
            // a are the first few terms of the expansion (initial values).
            // The integer r is used to indicate that C_jk = 0 for j>N-r (N = order of ode).
            void init(const arb::Acb &x, int k0, int r, const std::vector<arb::Acb> &a);

            // return lowest k value of expansion around x
            int start(const arb::Acb &x);

            // calculate expansion term of order k/2
            const arb::Acb &operator() (const arb::Acb &x, int k);
        private:
            arb::Acb C(int j, int k, const point_data &data) const; // C_jk
    };
}
