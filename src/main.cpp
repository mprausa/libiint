#include <iint/GPLKernel.h>
#include <iint/IInt.h>
#include <iint/ODE.h>
#include <iint/PsiSquared.h>
#include <iint/TRat.h>
#include <memory>
#include <iostream>

#if 0
int main() {
    auto w0 = std::make_shared<iint::GPLKernel>(0);
    auto w1 = std::make_shared<iint::GPLKernel>(1);
    auto wn = std::make_shared<iint::GPLKernel>(-1);
    long prec = 100;
    arb::Acb zero(0,prec);
    arb::Acb one(1,prec);
    arb::Acb x1(.25,prec);
    arb::Acb x2(.5,prec);

    auto hpl = iint::IInt::fetch({w1,w0,w1,wn,w1,w0});

    hpl->match(zero,x2);
    hpl->match(x2,one);

    std::cout << ">> " << (*hpl)(zero,x1) << std::endl;
    std::cout << ">> " << (*hpl)(x2,x1-x2) << std::endl;
    std::cout << ">> " << (*hpl)(one,x1-one) << std::endl;

    return 0;
}
#elif 0
int main() {
    iint::ODE ellipticK({{-1},{4,-8},{0,4,-4}});

    const long prec=100;
    arb::Acb x(.5,prec);

    ellipticK.init(x,0,0,{
        arb::Acb("1.854074677301371918433850347195260046217598823521766905585928045056021776838119978357271861650371897",prec),
        arb::Acb(0,prec),
        arb::Acb("0.8472130847939790866064991234821916364814459103269421850605793726597340048341347597232002939946112299",prec)
    });

    for (int k=0; k<=100; ++k) {
        std::cout << "a[" << k << "] = " << ellipticK(x,k) << std::endl;
    }

    return 0;
}
#elif 0
int main() {
    const long prec=100;

    iint::PsiSquared psi2;
    arb::Acb x(.5,prec);

    psi2.init(x);

    return 0;
}
#else
int main() {
    const long prec=100;
    //auto x = arb::Acb(2,prec)/99;
    arb::Acb x(0,prec);

    iint::TRat trat({-2,-3,-4,-15,-7,8},{0,0,-1,2,-2,-1,7,6});

    std::cout << trat << std::endl;
    int n0 = trat.init(x);

    for (int n=n0; n<=20; ++n) {
        std::cout << "trat[" << n << "] = " << trat(x,n) << std::endl;
    }
}

#endif
