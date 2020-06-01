#pragma once

#include <iint/ODE.h>


namespace iint {
    // t = (1 - 9*x - Sqrt[1 - 18*x + x^2])/(2*x)
    // psi^2 = 2*Pi^2*Hypergeometric2F1[1/4,3/4,1,(t*(4 + t)^5)/((4 + 6*t + t^2)^2*(20 + 8*t + t^2))]^2

    class PsiSquared {
        protected:
            ODE _phi;   // phi = (20 + 8*t + t^2)/(4 + 6*t + t^2) * psi^2
        public:
            PsiSquared();

            int init(const arb::Acb &x);
            arb::Acb &operator() (const arb::Acb &x, int k);
    };
}
