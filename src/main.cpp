#include <iint/GPLKernel.h>
#include <iint/IInt.h>
#include <iint/ODE.h>
#include <iint/EllipticKernel.h>
#include <iint/TauKernel.h>
#include <iint/KappaKernel.h>
#include <iint/TRat.h>
#include <iint/Matching.h>
#include <memory>
#include <iostream>

int main() {
    const long prec = 300;

    arb::Acb zero(0,prec);
    arb::Acb one(1,prec);
    auto sing = 9-4*arb::Acb(5,prec).sqrt();

    #if 0
        auto x0 = zero;
        auto points = iint::Matching::points3(zero,sing,one,.03125);
        arb::Acb x(.9375,prec);
    #else
        auto x0 = one;
        auto points = iint::Matching::points3(one,sing,zero,.03125);
        arb::Acb x(.015625,prec);
    #endif

    auto tau = std::make_shared<iint::TauKernel>();
    auto w0 = std::make_shared<iint::GPLKernel>(0);
    auto w1 = std::make_shared<iint::GPLKernel>(1);
    auto wn = std::make_shared<iint::GPLKernel>(-1);
    auto kappa = std::make_shared<iint::KappaKernel>();

    auto itau = iint::IInt::fetch({tau},x0);
    auto hpl = iint::IInt::fetch({kappa,w0,wn,w1},x0);

    auto ii1 = iint::IInt::fetch({tau,kappa,w0,wn,w1},x0);
    auto ii2 = iint::IInt::fetch({kappa,tau,w0,wn,w1},x0);
    auto ii3 = iint::IInt::fetch({kappa,w0,tau,wn,w1},x0);
    auto ii4 = iint::IInt::fetch({kappa,w0,wn,tau,w1},x0);
    auto ii5 = iint::IInt::fetch({kappa,w0,wn,w1,tau},x0);

    for (auto &x : points) {
        std::cout << "matching " << x[0] << " -> " << x[1] << " @ " << x[2] << std::endl;

        ii1->match(x[0],x[1],x[2]);
        ii2->match(x[0],x[1],x[2]);
        ii3->match(x[0],x[1],x[2]);
        ii4->match(x[0],x[1],x[2]);
        ii5->match(x[0],x[1],x[2]);
    }


    std::cout << "@ " << x << ":" << std::endl;

    auto vtau = (*itau)(x);
    auto vhpl = (*hpl)(x);

    auto v1 = (*ii1)(x);
    auto v2 = (*ii2)(x);
    auto v3 = (*ii3)(x);
    auto v4 = (*ii4)(x);
    auto v5 = (*ii5)(x);

    std::cout << "  " << (*itau) << " = " << vtau << std::endl;
    std::cout << "  " << (*hpl) << " = " << vhpl << std::endl;
    std::cout << "  " << (*ii1) << " = " << v1 << std::endl;
    std::cout << "  " << (*ii2) << " = " << v2 << std::endl;
    std::cout << "  " << (*ii3) << " = " << v3 << std::endl;
    std::cout << "  " << (*ii4) << " = " << v4 << std::endl;
    std::cout << "  " << (*ii5) << " = " << v5 << std::endl;
    std::cout << "  shuffle: " << (v1 + v2 + v3 + v4 + v5 - vtau*vhpl) << std::endl;

    return 0;
}
