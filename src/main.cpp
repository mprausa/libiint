#include <iint/GPLKernel.h>
#include <iint/IInt.h>
#include <iint/ODE.h>
#include <iint/EllipticKernel.h>
#include <iint/TauKernel.h>
#include <iint/TRat.h>
#include <iint/Matching.h>
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

    for (int n=1; n<32; ++n) {
        auto x = arb::Acb(0.03125,prec)*n;

        std::cout << "x=" << x << " => " << (*hpl)(x) << std::endl;
    }

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
    arb::Acb x(1,prec);

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
    const long prec = 300;

    auto tau = std::make_shared<iint::TauKernel>();
    arb::Acb zero(0,prec);
    arb::Acb one(1,prec);
    auto sing = 9-4*arb::Acb(5,prec).sqrt();

    auto itau = iint::IInt::fetch({tau});

    auto points = iint::Matching::points3(zero,sing,one,.03125);

    for (auto &x : points) {
        std::cout << "matching " << x[0] << " -> " << x[1] << " @ " << x[2] << std::endl;
        itau->match(x[0],x[1],x[2]);
    }

    arb::Acb x(.9375,prec);

    std::cout << *itau << " @ " << x << ": " << (*itau)(x) << std::endl;
    return 0;

#if 0
    for (int n=0; n<=10; ++n) {
        std::cout << "itau[" << n << "] = " << (*itau)(x1,n,0) << std::endl;
    }

    std::cout << (*itau)(x1,x1/16) << std::endl;

    return 0;

    itau->match(zero,x1/4,x1/8);
    itau->match(x1/4,x1/2,3*x1/8);
    itau->match(x1/2,3*x1/4,5*x1/8);
    itau->match(3*x1/4,x1,7*x1/8);

    std::cout << (*itau)(x1,x1/16) << std::endl;
#endif

    return 0;
}
#elif 0
int main() {
    const long prec = 100;
    arb::Acb x1(0.015625,prec);
    arb::Acb zero(0,prec);

    auto tau = std::make_shared<iint::TauKernel>();
    auto w0 = std::make_shared<iint::GPLKernel>(0);
    auto w1 = std::make_shared<iint::GPLKernel>(1);
    auto wn = std::make_shared<iint::GPLKernel>(-1);

    auto iint1 = iint::IInt::fetch({tau,w1,w0,wn,w1});
    auto iint2 = iint::IInt::fetch({w1,tau,w0,wn,w1});
    auto iint3 = iint::IInt::fetch({w1,w0,tau,wn,w1});
    auto iint4 = iint::IInt::fetch({w1,w0,wn,tau,w1});
    auto iint5 = iint::IInt::fetch({w1,w0,wn,w1,tau});

    auto res = (*iint1)(zero,x1) + (*iint2)(zero,x1) + (*iint3)(zero,x1) + (*iint4)(zero,x1) + (*iint5)(zero,x1);

    std::cout << ">> " << res << std::endl;

    return 0;
}
#endif
