/*
 *  include/iint/KernelFactory.h
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

namespace iint {
    class KernelFactory {
        protected:
            static std::unordered_map<std::string,std::shared_ptr<iint::Kernel>> _kernels;
        public:
            // retrieve integration kernel by string representation
            static std::shared_ptr<Kernel> get(const std::string &s, long prec);
    };
}

