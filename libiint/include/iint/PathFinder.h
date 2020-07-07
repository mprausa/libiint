/*
 *  include/iint/PathFinder.h
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

#include <vector>
#include <functional>
#include <complex>
#include <arb/Acb.h>

namespace iint {
    // find matching points along a path
    class PathFinder {
        protected:
            std::vector<arb::Acb> _sings;
        public:
            using mpoint_t = std::array<arb::Acb,3>;    // parameters {x1,x2,x} for IInt::match
            using mpoints_t = std::vector<mpoint_t>;

            using path_fun_t = std::function<arb::Acb(double)>;

            // sings is a list of singularities
            PathFinder(const std::vector<arb::Acb> &sings) : _sings(sings) {}

            // Return matching points.
            // fun(t) parameterizes the path for 0<=t<=1.
            // fun(0) and fun(1) should lie on singularities.
            // For all other values of t, fun(t) must fun(t) must not be a singular point.
            // It is guaranteed that for every returned point |x-x1| <= min(q*r1,a) and |x-x2| <= min(q*r2,a),
            // where r1(r2) is the convergence radius for an expansion around x1(x2).
            mpoints_t operator() (const path_fun_t &fun, double q, const arb::Acb &a = arb::Acb::infty) const;

            // Return matching points for a predefined path from x=1 to x=0 along the euclidean region
            static mpoints_t euclidean(long prec, double q, const arb::Acb &a = arb::Acb::infty);

            // Return matching points for a predefined path from x=1 to x=0 along the physical region
            static mpoints_t physical(long prec, double q, const arb::Acb &a = arb::Acb::infty);
        private:
            arb::Acb radius(const arb::Acb &p) const;
            double scan(const path_fun_t &fun, double t, const arb::Acb &d) const;
            mpoints_t convert(const path_fun_t &fun, const std::vector<double> &path) const;
    };
}
