#include <iint/KernelFactory.h>
#include <iint/GPLKernel.h>
#include <iint/TauKernel.h>
#include <iint/MuKernel.h>
#include <iint/KappaKernel.h>

namespace iint {
    std::unordered_map<std::string,std::shared_ptr<iint::Kernel>> KernelFactory::_kernels;

    std::shared_ptr<Kernel> KernelFactory::get(const std::string &s, long prec) {
        auto &krn = _kernels[s];
        if (krn) return krn;

        if (s == "tau") {
            krn = std::make_shared<iint::TauKernel>();
        } else if (s == "kappa") {
            krn = std::make_shared<iint::KappaKernel>();
        } else if (s.size() > 4 && s.substr(0,3) == "mu(" && s.back() == ')') {
            int n = std::stoi(s.substr(3,s.size()-4));
            krn = std::make_shared<iint::MuKernel>(n);
        } else if (s.size() > 7 && s.substr(0,6) == "omega(" && s.back() == ')') {
            arb::Acb a(s.substr(6,s.size()-7),prec);
            krn = std::make_shared<iint::GPLKernel>(a);
        } else {
            assert(false);
        }

        return krn;
    }
}

