#include <iint/GPLKernel.h>
#include <iint/IInt.h>
#include <iint/ODE.h>
#include <iint/EllipticKernel.h>
#include <iint/TauKernel.h>
#include <iint/KappaKernel.h>
#include <iint/MuKernel.h>
#include <iint/TRat.h>
#include <iint/PathFinder.h>
#include <memory>
#include <iostream>
#include <fstream>

int main() {
    iint::verbose = true;

    const long prec = 300;
    double q = 0.015625;

    arb::Acb zero(0,prec);
    arb::Acb one(1,prec);
    auto sing = 9-4*arb::Acb(5,prec).sqrt();
    auto x0 = one;

    iint::PathFinder pf({-one,zero,one,sing});

#if 0
    arb::Acb x(.015625,prec);
    auto points = iint::PathFinder::euclidean(prec,q);
#else
    arb::Acb x(-.625,prec);
    auto points = iint::PathFinder::physical(prec,q);
#endif

    auto tau = std::make_shared<iint::TauKernel>();
    auto omega0 = std::make_shared<iint::GPLKernel>(0);
    auto omega1 = std::make_shared<iint::GPLKernel>(1);
    auto omegaN1 = std::make_shared<iint::GPLKernel>(-1);
    auto kappa = std::make_shared<iint::KappaKernel>();
    auto mu2 = std::make_shared<iint::MuKernel>(2);

    auto itau = iint::IInt::fetch({tau},x0);
    auto hpl = iint::IInt::fetch({kappa,omega0,mu2,omega1},x0);

    auto ii1 = iint::IInt::fetch({tau,kappa,omega0,mu2,omega1},x0);
    auto ii2 = iint::IInt::fetch({kappa,tau,omega0,mu2,omega1},x0);
    auto ii3 = iint::IInt::fetch({kappa,omega0,tau,mu2,omega1},x0);
    auto ii4 = iint::IInt::fetch({kappa,omega0,mu2,tau,omega1},x0);
    auto ii5 = iint::IInt::fetch({kappa,omega0,mu2,omega1,tau},x0);

    {
        size_t cnt=1;
        for (auto &x : points) {
            std::cout << "[" << cnt++ << "/" << points.size() << "] matching " << x[0] << " -> " << x[1] << " @ " << x[2] << std::endl;

            ii1->match(x[0],x[1],x[2]);
            ii2->match(x[0],x[1],x[2]);
            ii3->match(x[0],x[1],x[2]);
            ii4->match(x[0],x[1],x[2]);
            ii5->match(x[0],x[1],x[2]);
        }
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

