#pragma once

#include <iint/ODE.h>
#include <iint/TRat.h>
#include <iint/Kernel.h>

namespace iint {
    // t = (1 - 9*x - Sqrt[1 - 18*x + x^2])/(2*x)
    // z = (t*(4 + t)^5)/((4 + 6*t + t^2)^2*(20 + 8*t + t^2))

    //        / Sqrt[2]*Pi*Hypergeometric2F1[1/4,3/4,1,z]                                               ;  Re[x] < 7-4*Sqrt[3]
    // psi = <
    //        \ Sqrt[2]*Pi*Hypergeometric2F1[1/4,3/4,1,z] + 2*I*Pi*Hypergeometric2F1[1/4,3/4,1,1-z]     ;  Re[x] >= 7-4*Sqrt[3]

    // ell[x] = I*Pi * numer[t]/denom[t] * (20 + 8*t + t^2)/(4 + 6*t + t^2) * psi^2
    class EllipticKernel : public Kernel {
        protected:
            TRat _trat;
            ODE _phi;   // phi = (20 + 8*t + t^2)/(4 + 6*t + t^2) * psi^2
        public:
            EllipticKernel(const std::vector<int> &numer, const std::vector<int> &denom);
        protected:
            virtual arb::Acb _calc(const arb::Acb &x, int k);
            virtual int _init(const arb::Acb &x);
    };
}
