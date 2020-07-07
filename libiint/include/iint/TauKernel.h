/*
 *  include/iint/TauKernel.h
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

#include <iint/EllipticKernel.h>

namespace iint {
    // tau[x] = ((5*I)*Pi*(4 + t)*(5 + t)*(4 + t*(6 + t)))/(t*(-20 + t^2)) / psi^2
    class TauKernel : public Kernel {
        protected:
            EllipticKernel _reciprocal; // -Pi^2/tau[x]
        public:
            TauKernel();

            virtual std::string str() const {
                return "tau";
            }
        protected:
            virtual int _init(const arb::Acb &x);
            virtual arb::Acb _calc(const arb::Acb &x, int k);
    };
}

