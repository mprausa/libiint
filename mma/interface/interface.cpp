#include <mathlink.h>
#include <iint/IInt.h>
#include <iint/GPLKernel.h>
#include <iint/TauKernel.h>
#include <iint/MuKernel.h>
#include <iint/KappaKernel.h>
#include <iint/Matching.h>
#include <iostream>

static long prec = 100;
static std::unordered_map<std::string,std::shared_ptr<iint::Kernel>> kernels;
static std::vector<std::shared_ptr<iint::IInt>> iints;

static std::shared_ptr<iint::Kernel> get_kernel(const std::string &s) {
    auto &krn = kernels[s];
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
        krn = nullptr;
    }

    return krn;
}

void IIntInit(const char *sprec) {
    prec = std::atol(sprec);
    std::cout << "setting default precision to " << prec << std::endl;
}

void IIntCreate(const char *ckernels, const char *cx0) {
    std::string s;
    std::vector<std::shared_ptr<iint::Kernel>> kernels;

    for (size_t n=0; ckernels[n]; ++n) {
        auto &c = ckernels[n];

        if (c == ',') {
            kernels.push_back(get_kernel(s));
            s.clear();
        } else {
            s += c;
        }
    }
    kernels.push_back(get_kernel(s));

    arb::Acb x0(cx0,prec);

    auto iint = iint::IInt::fetch(kernels,x0);

    size_t id = iints.size();
    iints.push_back(iint);

    std::cout << "identifier for " << *iint << " is II" << id << std::endl;

    if (!MLPutString(stdlink,("II"+std::to_string(id)).c_str())) {
        std::cerr << "MLPutString - error." << std::endl;
    }
}

void IIntMatch(const char *cx1, const char *cx2, const char *cx3, const char *cq) {
    arb::Acb x1(cx1,prec);
    arb::Acb x2(cx2,prec);
    arb::Acb x3(cx3,prec);
    double q = std::stod(cq);

    auto points = iint::Matching::points3(x1,x2,x3,q);

    std::cout << "matiching " << x1 << " -> " << x2 << " -> " << x3 << " (q=" << q << ")" << std::endl;

    for (auto &x : points) {
        for (auto &i : iints) {
            i->match(x[0],x[1],x[2]);
        }
    }
    std::cout << "done." << std::endl;
}

void IIntEvaluate(const char *cid, const char *cx) {
    std::string sid = cid;

    std::cout << "evaluating " << sid << std::endl;

    if (sid.substr(0,2) != "II") {
        if (!MLPutString(stdlink,"$Failed")) {
            std::cerr << "MLPutString - error." << std::endl;
        }
        return;
    }

    size_t id = std::stol(sid.substr(2));
    if (id >= iints.size()) {
        if (!MLPutString(stdlink,"$Failed")) {
            std::cerr << "MLPutString - error." << std::endl;
        }
        return;
    }

    arb::Acb x(cx,prec);

    auto res = (*iints[id])(x);

    if (!MLPutString(stdlink,res.mma().c_str())) {
        std::cerr << "MLPutString - error." << std::endl;
    }
}

int main(int argc, char **argv) {
    return MLMain(argc,argv);
}

