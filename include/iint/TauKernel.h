#pragma once

#include <iint/EllipticKernel.h>

namespace iint {
    // tau[x] = ((5*I)*Pi*(4 + t)*(5 + t)*(4 + t*(6 + t)))/(t*(-20 + t^2)) / psi^2
    class TauKernel : public Kernel {
        protected:
            EllipticKernel _reciprocal; // I*Pi/tau[x]
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

