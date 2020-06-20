#include <mathlink.h>
#include <iint/IInt.h>
#include <iint/GPLKernel.h>
#include <iint/TauKernel.h>
#include <iint/MuKernel.h>
#include <iint/KappaKernel.h>
#include <iint/PathFinder.h>
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

void IIntInit(int iprec, int iverbose) {
    prec = iprec;
    iint::verbose = (bool)iverbose;

    if (iint::verbose) std::cout << "setting default precision to " << prec << std::endl;
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

    if (iint::verbose) std::cout << "identifier for " << *iint << " is II" << id << std::endl;

    if (!MLPutString(stdlink,("II"+std::to_string(id)).c_str())) {
        std::cerr << "MLPutString - error." << std::endl;
    }
}

void IIntMatchEuclidean(const char *cq, const char *ca) {
    double q = std::stod(cq);
    auto a = (std::string(ca) == "default") ? arb::Acb::infty : arb::Acb(ca,prec);

    auto points = iint::PathFinder::euclidean(prec,q,a);

    size_t cnt=1;
    for (auto &x : points) {
        if (iint::verbose) std::cout << "[" << cnt << "/" << points.size() << "] matching " << x[0] << " -> " << x[1] << " @ " << x[2] << std::endl;

        for (auto &i : iints) {
            if (MLAbort) {
                if (iint::verbose) std::cout << "aborted.";
                return;
            }
            i->match(x[0],x[1],x[2]);
        }
        ++cnt;
    }
    if (iint::verbose) std::cout << "done." << std::endl;
}

void IIntMatchPhysical(const char *cq, const char *ca) {
    double q = std::stod(cq);
    auto a = (std::string(ca) == "default") ? arb::Acb::infty : arb::Acb(ca,prec);

    auto points = iint::PathFinder::physical(prec,q,a);

    size_t cnt=1;
    for (auto &x : points) {
        if (iint::verbose) std::cout << "[" << cnt << "/" << points.size() << "] matching " << x[0] << " -> " << x[1] << " @ " << x[2] << std::endl;

        for (auto &i : iints) {
            if (MLAbort) {
                if (iint::verbose) std::cout << "aborted.";
                return;
            }
            i->match(x[0],x[1],x[2]);
        }
        ++cnt;
    }
    if (iint::verbose) std::cout << "done." << std::endl;
}

void IIntEvaluate(const char *cid, const char *cx) {
    std::string sid = cid;

    if (iint::verbose) std::cout << "evaluating " << sid << std::flush;

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

    if (iint::verbose) std::cout << " => " << res << std::endl;

    if (!MLPutString(stdlink,res.mma().c_str())) {
        std::cerr << "MLPutString - error." << std::endl;
    }
}

int main(int argc, char **argv) {
    return MLMain(argc,argv);
}

