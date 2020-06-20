#pragma once

#include <iint/Kernel.h>

namespace iint {
    class KernelFactory {
        protected:
            static std::unordered_map<std::string,std::shared_ptr<iint::Kernel>> _kernels;
        public:
            static std::shared_ptr<Kernel> get(const std::string &s, long prec);
    };
}

