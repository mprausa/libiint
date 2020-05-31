#include <iint/GPLKernel.h>
#include <iint/IInt.h>
#include <memory>
#include <iostream>

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
