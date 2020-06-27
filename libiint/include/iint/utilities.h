#pragma once

#include <arb/Acb.h>

namespace iint {
    inline arb::Acb xeps(const arb::Acb &x) {
        if (!x.imag().contains_zero()) return x;
        long prec = x.default_prec();
        arb::Acb eps;
        acb_mul_2exp_si(eps.get(),arb::Acb(1,prec).get(),-2*prec);
        return x + arb::Acb::I*eps;
    }
}
