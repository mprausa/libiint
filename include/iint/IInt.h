#pragma once

#include <iint/Kernel.h>
#include <memory>
#include <unordered_map>

namespace iint {
    class IInt {
        public:
            using kernels_t = std::vector<std::shared_ptr<Kernel>>;
        protected:
            kernels_t _kernels;
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

            struct kernels_hasher {
                static std::size_t combine(std::size_t seed, std::size_t hash) {
                        return seed ^ (hash + 0x9e3779b9 + (seed << 6) + (seed >> 2));
                }

                std::size_t operator()(const kernels_t &kernels) const {
                    std::size_t hash = std::hash<size_t>()(kernels.size());
                    for (auto &k : kernels) {
                        hash = combine(hash,std::hash<std::shared_ptr<Kernel>>()(k));
                    }

                    return hash;
                }

            };

            static std::unordered_map<kernels_t,std::shared_ptr<IInt>,kernels_hasher> _iints;
            static arb::Acb _zero,_one;
        public:
            IInt(const kernels_t &kernels);

            void match(const arb::Acb &x1, const arb::Acb &x2, const arb::Acb &x);

            int start(const arb::Acb &x);
            int maxlog(const arb::Acb &x);
            const arb::Acb &operator() (const arb::Acb &x, int n, int m);
            arb::Acb operator() (const arb::Acb &x, const arb::Acb &delta);
            arb::Acb operator() (const arb::Acb &x);

            static std::shared_ptr<IInt> fetch(const kernels_t &kernels);

            std::string str() const {
                if (_kernels.empty()) return "IInt[]";

                std::string s="IInt[";

                for (auto &k : _kernels) {
                    s += k->str()+",";
                }
                s.pop_back();
                s += "]";
                return s;
            }
    };

    inline std::ostream &operator<<(std::ostream &os, const IInt &iint) {
        os << iint.str();
        return os;
    }
}
