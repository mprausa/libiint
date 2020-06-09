#include <iint/MuKernel.h>

namespace {
    std::vector<int> prepend(const std::vector<int> &v, int n) {
        std::vector<int> u(n,0);
        u.insert(u.end(),v.begin(),v.end());
        return u;
    }
}

namespace iint {
    MuKernel::MuKernel(int n)
        : _n(n), EllipticKernel(
            prepend({100,45,5},n?n-1:0),
            prepend({-8000,-6400,-1680,0,84,16,1},n?0:1)) {}
}

