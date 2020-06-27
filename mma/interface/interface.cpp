#include <mathlink.h>
#include <iint/IInt.h>
#include <iint/KernelFactory.h>
#include <iint/PathFinder.h>
#include <iostream>
#include <fstream>

static long prec = 100;
static std::unordered_map<std::string,std::shared_ptr<iint::Kernel>> kernels;
static std::vector<std::shared_ptr<iint::IInt>> iints;

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
            kernels.push_back(iint::KernelFactory::get(s,prec));
            s.clear();
        } else {
            s += c;
        }
    }
    kernels.push_back(iint::KernelFactory::get(s,prec));

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

void IIntSave(const char *fn) {
    if (iint::verbose) std::cout << "saving constants to " << fn << std::endl;
    std::ofstream file(fn);
    file << iint::IInt::store() << std::endl;
    file.close();
}

void IIntLoad(const char *fn) {
    if (iint::verbose) std::cout << "loading constants from " << fn << std::endl;
    iint::IInt::restore(YAML::LoadFile(fn),prec);
}

void IIntEvaluate(const char *cid, const char *cx) {
    std::string sid = cid;

    if (iint::verbose) std::cout << "evaluating " << sid << " @ " << cx << std::flush;

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

