/*
 *  include/iint/ODE.h
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

#include <arb/Acb.h>
#include <unordered_map>

namespace iint {
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
            ODE(const std::vector<std::vector<int>> &p);

            void init(const arb::Acb &x, int k0, int r, const std::vector<arb::Acb> &a);
            int start(const arb::Acb &x);
            const arb::Acb &operator() (const arb::Acb &x, int k);
        private:
            arb::Acb C(int j, int k, const point_data &data) const;
    };
}
