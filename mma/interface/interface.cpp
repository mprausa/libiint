/*
 *  mma/interface/interface.cpp
 *
 *  Copyright (C) 2020 Mario Prausa
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
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

void IIntSeries(const char *cid, const char *cvar, const char *cx, int order) {
    std::string sid = cid;
    std::string sx = cx;
    std::string svar = cvar;

    if (iint::verbose) std::cout << "series expansion of " << sid << " around " << cx << " upto order " << order << std::endl;

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

    auto iint = iints[id];
    arb::Acb x(sx,prec);

    int n0 = iint->start(x);
    int mhat = iint->maxlog(x);

    std::string res = "O["+svar+","+sx+"]^"+std::to_string(order+1);

    for (int n=n0; n<=2*order+1; ++n) {
        for (int m=0; m<=mhat; ++m) {
            auto v = (*iint)(x,n,m);
            if (v.contains_zero()) continue;
            res += " + (" + v.mma() + ")*(" + svar + "-(" + sx + "))^("+std::to_string(n)+"/2)*Log["+svar+"-("+sx+")]^"+std::to_string(m);
        }
    }

    if (!MLPutString(stdlink,res.c_str())) {
        std::cerr << "MLPutString - error." << std::endl;
    }
}

int main(int argc, char **argv) {
    return MLMain(argc,argv);
}

