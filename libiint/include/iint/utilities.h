/*
 *  include/iint/utilities.h
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

namespace iint {
    // add small positive imaginary part to x if x is real
    inline arb::Acb xeps(const arb::Acb &x) {
        if (!x.imag().contains_zero()) return x;
        long prec = x.default_prec();
        arb::Acb eps;
        acb_mul_2exp_si(eps.get(),arb::Acb(1,prec).get(),-2*prec);
        return x + arb::Acb::I*eps;
    }
}
