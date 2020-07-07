/*
 *  include/iint/TRat.h
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

#include <iint/ODE.h>

namespace iint {
    // class for rational functions in t (expanded in x)
    class TRat {
        protected:
            const arb::Acb _zero{0};

            std::vector<int> _numer,_denom;
            ODE _ode;   // Sqrt[1-18*x+x^2]/(2*x)

            struct point_data {
                arb::Acb x;
                int n0;
                bool reciprocal;
                std::vector<arb::Acb> numer,denom;
                int nummin,denmin;
                std::vector<arb::Acb> rat,xrat,cache;
                std::vector<std::vector<arb::Acb>> bell;
            };

            std::unordered_map<arb::Acb,point_data> _points;
        public:
            // Represents numer[t]/denom[t], where t = (1 - 9*x - Sqrt[1 - 18*x + x^2])/(2*x).
            // nummer(denom) is a std::vector<int> of the coefficients of numer[t](denom[t]) starting at t^0
            TRat(const std::vector<int> &numer, const std::vector<int> &denom);

            // initialize expansion around x, return lowest value of n
            int init(const arb::Acb &x);

            // calculate expansion coefficient of (x-x1)^(k/2)
            const arb::Acb &operator() (const arb::Acb &x1, int n);

            void print(std::ostream &os) const;
        private:
            const arb::Acb &xrat_expansion(point_data &data, int n);
            arb::Acb t_expansion(const arb::Acb &x, int n);
            static arb::Acb rat_expansion(point_data &data, int n);
            const arb::Acb &bell(point_data &data, int n, int k);
    };

    inline std::ostream &operator<<(std::ostream &os, const TRat &trat) {
        trat.print(os);
        return os;
    }
}
