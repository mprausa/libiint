/*
 *  include/iint/IInt.h
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

#pragma once

#include <iint/Kernel.h>
#include <memory>
#include <unordered_map>
#include <iostream>

namespace iint {
    extern bool verbose;

    // class for iterated integrals
    class IInt {
        protected:
            struct cdbl_less {
                bool operator() (std::complex<double> a, std::complex<double> b) {
                    constexpr double eps = 1e-16;

                    if (b.real() - a.real() > eps) return true;
                    if (a.real() - b.real() > eps) return false;
                    if (b.imag() - a.imag() > eps) return true;
                    return false;
                }
            };

        public:
            using kernels_t = std::vector<std::shared_ptr<Kernel>>;
        protected:
            kernels_t _kernels;     // kernels f_1,...,f_N
            arb::Acb _x0;

            // _constants[x] contains the constant term of the expansion around x
            std::unordered_map<arb::Acb,arb::Acb> _constants;
            std::shared_ptr<IInt> _subiint = nullptr;

            struct cache_key_s {
                arb::Acb x;
                int n,m;

                bool operator==(const cache_key_s &other) const {
                    return x == other.x && n == other.n && m == other.m;
                }
            };

            struct cache_hasher {
                static std::size_t combine(std::size_t seed, std::size_t hash) {
                        return seed ^ (hash + 0x9e3779b9 + (seed << 6) + (seed >> 2));
                }

                std::size_t operator()(const cache_key_s &key) const {
                    std::size_t hash = std::hash<arb::Acb>()(key.x);
                    hash = combine(hash,std::hash<int>()(key.n));
                    hash = combine(hash,std::hash<int>()(key.m));
                    return hash;
                }
            };

            std::unordered_map<cache_key_s,arb::Acb,cache_hasher> _cache;

            struct args_hasher {
                static std::size_t combine(std::size_t seed, std::size_t hash) {
                        return seed ^ (hash + 0x9e3779b9 + (seed << 6) + (seed >> 2));
                }

                std::size_t operator()(const std::pair<kernels_t,arb::Acb> &args) const {
                    std::size_t hash = std::hash<size_t>()(args.first.size());
                    for (auto &k : args.first) {
                        hash = combine(hash,std::hash<std::shared_ptr<Kernel>>()(k));
                    }

                    hash = combine(hash,std::hash<arb::Acb>()(args.second));

                    return hash;
                }

            };

            static std::unordered_map<std::pair<kernels_t,arb::Acb>,std::shared_ptr<IInt>,args_hasher> _iints;
            static arb::Acb _zero,_one;
        public:
            IInt(const kernels_t &kernels, const arb::Acb &x0);

            // match the constant term of the expansion around x2
            // by evaluating the expansions around x1 and x2 at x
            void match(const arb::Acb &x1, const arb::Acb &x2, const arb::Acb &x);

            int start(const arb::Acb &x);   // lowest n-value of expansion around x
            int maxlog(const arb::Acb &x);  // largest m-value of expansion around x

            // return series coefficient of (x-x1)^(n/2)*Log[x-x1]^m
            const arb::Acb &operator() (const arb::Acb &x1, int n, int m);

            // evaluate series around x at x+delta
            arb::Acb operator() (const arb::Acb &x, const arb::Acb &delta);

            // Evaluate iterative integral at x.
            // This function looks for the closest point x1 where a constant term
            // is available and calls (*this)(x1,x-x1).
            arb::Acb operator() (const arb::Acb &x);

            // New IInt objects should always be created using the fetch function.
            // Provides caching
            static std::shared_ptr<IInt> fetch(const kernels_t &kernels, const arb::Acb &x0);

            // string representation of IInt
            std::string str() const {
                if (_kernels.empty()) return "II(;"+_x0.mma()+")";

                std::string s="II(";

                for (auto &k : _kernels) {
                    s += k->str()+",";
                }
                s.pop_back();

                s += ";"+_x0.mma()+")";
                return s;
            }

            // store all constants in a YAML node
            YAML::Node store_constants();

            // store all constants of all IInts in a YAML node
            static YAML::Node store();

            // restore constants from YAML node
            void restore_constants(const YAML::Node &node);

            // restore constants of all IInts in YAML node
            static void restore(const YAML::Node &node, long prec);
    };

    inline std::ostream &operator<<(std::ostream &os, const IInt &iint) {
        os << iint.str();
        return os;
    }
}
