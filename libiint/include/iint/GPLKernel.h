/*
 *  include/iint/GPLKernel.h
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

#include <iint/Kernel.h>

namespace iint {
    class GPLKernel : public Kernel {
        protected:
            arb::Acb _a;
        public:
            GPLKernel(const arb::Acb &a) : _a(a) {
                _nocache = true;
            }

            virtual std::string str() const {
                return "omega("+_a.str()+")";
            }
        protected:
            virtual arb::Acb _calc(const arb::Acb &x, int k);
            virtual int _init(const arb::Acb &x);
    };
}
