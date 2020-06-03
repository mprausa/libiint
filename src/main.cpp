#include <iint/GPLKernel.h>
#include <iint/IInt.h>
#include <iint/ODE.h>
#include <iint/EllipticKernel.h>
#include <iint/TauKernel.h>
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

    iint::EllipticKernel ell({-2,-3,-4,-15,-7,8},{0,0,-1,2,-2,-1,7,6});
    arb::Acb x(.5,prec);

    int n0 = ell.init(x);
    for (int n=n0; n<10; ++n) {
        std::cout << "ell[" << n << "] = " << ell(x,n) << std::endl;
    }

    return 0;
}
#elif 0
int main() {
    const long prec=300;
    //auto x = arb::Acb(2,prec)/99;
    arb::Acb x(0,prec);
    //auto x = 9 - 4*arb::Acb(5,prec).sqrt();

    iint::TRat trat({0,0,-2,-3,-4,-15,-7,8},{-1,2,-2,-1,7,6});

    std::cout << trat << std::endl;
    int n0 = trat.init(x);

    for (int n=n0; n<=20; ++n) {
        std::cout << "trat[" << n << "] = " << trat(x,n) << std::endl;
    }
}
#elif 0
int main() {
    const long prec=300;
    arb::Acb x(0,prec);

    iint::TauKernel tau;
    int n0 = tau.init(x);

    for (int n=n0; n<=10; ++n) {
        std::cout << "tau[" << n << "] = " << tau(x,n) << std::endl;
    }

    return 0;
}
#elif 1
int main() {
    const long prec = 100;

    auto tau = std::make_shared<iint::TauKernel>();
    arb::Acb zero(0,prec);
    arb::Acb x1(0.03125,prec);
    //arb::Acb x1(0.015625,prec);

    auto itau = iint::IInt::fetch({tau});

    std::cout << *itau << " @ " << x1 << std::endl;
    std::cout << ">> " << (*itau)(zero,x1) << std::endl;

    return 0;
}

#endif
