/*
 *  src/iint/PathFinder.cpp
 *
 *  Copyright (C) 2020 Mario Prausa
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iint/PathFinder.h>
#include <iostream>
#include <limits>

namespace iint {
    PathFinder::mpoints_t PathFinder::operator() (const path_fun_t &fun, double q, const arb::Acb &a) const {
        double t=0.;
        std::vector<double> path;
        double t1;

        {
            auto d = radius(fun(1.))*q;
            if (d > a) d = a;
            t1 = 1.-scan([&fun](double t) {return fun(1.-t);},0.,d);
        }

        while (t<t1 || (path.size() % 2)) {
            path.push_back(t);
            auto d = radius(fun(t))*q;
            if (d > a) d = a;
            t = scan(fun,t,d);
        }

        path.push_back(1.);

        return convert(fun,path);
    }

    arb::Acb PathFinder::radius(const arb::Acb &p) const {
        auto r = arb::Acb::infty;

        for (auto s : _sings) {
            auto r0 = (p-s).abs();
            if (r0.contains_zero()) continue;
            if (r0 < r) r = r0;
        }

        return r;
    }

    double PathFinder::scan(const path_fun_t &fun, double t, const arb::Acb &d) const {
        constexpr double eps = 1e-15;
        auto p0 = fun(t);

        double a,b;

        for (b=t; b<1.; b+=.01) {
            auto p = fun(b);
            if ((p-p0).abs() < d) {
                a = b;
            } else {
                break;
            }
        }

        if (b >= 1.) b = 1.;

        while (std::abs(a-b) > eps) {
            double c = (a+b)/2;
            auto p = fun(c);
            if ((p-p0).abs() < d) {
                a = c;
            } else {
                b = c;
            }
        }


        return (a+b)/2;
    }

    PathFinder::mpoints_t PathFinder::convert(const path_fun_t &fun, const std::vector<double> &path) const {
        assert(path.size() % 2);

        mpoints_t points;

        for (size_t n=0,sz=path.size(); n<sz-2; n+=2) {
            auto &t1 = path[n];
            auto &tm = path[n+1];
            auto &t2 = path[n+2];

            points.push_back({fun(t1),fun(t2),fun(tm)});
        }

        return points;
    }

    PathFinder::mpoints_t PathFinder::euclidean(long prec, double q, const arb::Acb &a) {
        arb::Acb zero(0,prec);
        arb::Acb one(1,prec);
        auto sing = 9-4*arb::Acb(5,prec).sqrt();

        iint::PathFinder pf({-one,zero,one,sing});

        auto points = pf([&sing,&one](double t) {return (sing-one)*t + one;},q,arb::Acb::infty);
        auto points1 = pf([&sing](double t) {return sing*(1.-t);},q,arb::Acb::infty);
        points.insert(points.end(),points1.begin(),points1.end());

        return points;
    }

    PathFinder::mpoints_t PathFinder::physical(long prec, double q, const arb::Acb &a) {
        arb::Acb zero(0,prec);
        arb::Acb one(1,prec);
        auto sing = 9-4*arb::Acb(5,prec).sqrt();

        iint::PathFinder pf({-one,zero,one,sing});

        auto points = pf([&prec](double t) {return arb::Acb(t,prec).exp_pi_i();},q,arb::Acb::infty);
        auto points1 = pf([&prec](double t) {return arb::Acb(-1.+t,prec);},q,arb::Acb::infty);
        points.insert(points.end(),points1.begin(),points1.end());

        return points;
    }
}
