/*
 *  include/iint/PathFinder.h
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

#pragma once

#include <vector>
#include <functional>
#include <complex>
#include <arb/Acb.h>

namespace iint {
    class PathFinder {
        protected:
            std::vector<arb::Acb> _sings;
        public:
            using mpoint_t = std::array<arb::Acb,3>;
            using mpoints_t = std::vector<mpoint_t>;

            using path_fun_t = std::function<arb::Acb(double)>;
            PathFinder(const std::vector<arb::Acb> &sings) : _sings(sings) {}

            mpoints_t operator() (const path_fun_t &fun, double q, const arb::Acb &a = arb::Acb::infty) const;

            static mpoints_t euclidean(long prec, double q, const arb::Acb &a = arb::Acb::infty);
            static mpoints_t physical(long prec, double q, const arb::Acb &a = arb::Acb::infty);
        private:
            arb::Acb radius(const arb::Acb &p) const;
            double scan(const path_fun_t &fun, double t, const arb::Acb &d) const;
            mpoints_t convert(const path_fun_t &fun, const std::vector<double> &path) const;
    };
}
