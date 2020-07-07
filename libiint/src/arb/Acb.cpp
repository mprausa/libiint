/*
 *  src/arb/Acb.cpp
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

#include <arb/Acb.h>
#include <cassert>

namespace YAML {
    long convert<arb::Acb>::default_prec = 100;
}

namespace arb {
    Acb Acb::infty;
    Acb Acb::nan;
    Acb Acb::I;

    static class __Initializer {
        public:
            __Initializer() {
                arb_pos_inf(acb_realref(Acb::infty.get()));
                arb_indeterminate(acb_realref(Acb::nan.get()));
                acb_onei(Acb::I.get());
            }
    } __initializer;

    Acb Acb::Pi(long prec) {
        Acb pi;
        arb_const_pi(acb_realref(pi.get()),prec);
        pi.update_default_prec(prec);

        return pi;
    }

    Acb Acb::fac(unsigned long n, long prec) {
        Acb res;
        arb_fac_ui(acb_realref(res.get()),n,prec);
        res.update_default_prec(prec);
        return res;
    }

    Acb Acb::binom(unsigned long n, unsigned long k, long prec) {
        Acb res;
        arb_bin_uiui(acb_realref(res.get()),n,k,prec);
        res.update_default_prec(prec);
        return res;
    }

    Acb Acb::Zeta(const Acb &s, long prec) {
        Acb res;
        acb_zeta(res.get(),s.get(),prec);
        res.update_default_prec(prec);
        return res;
    }

    Acb Acb::Gamma(const Acb &s, long prec) {
        Acb res;
        acb_gamma(res.get(),s.get(),prec);
        res.update_default_prec(prec);
        return res;
    }

    Acb Acb::PolyLog(long s, const Acb &z, long prec) {
        Acb res;
        acb_polylog_si(res.get(),s,z.get(),prec);
        res.update_default_prec(prec);
        return res;
    }

    int Acb::xalloc1 = std::ios_base::xalloc();
    int Acb::xalloc2 = std::ios_base::xalloc();

    Acb::Acb(const std::string &str, long prec) : defPrec(prec) {
        acb_init(acb);
        size_t sz = str.size();

        if (sz < 3) {
            arb_set_str(acb_realref(acb),str.c_str(),prec);
            return;
        }

        if (str[sz-2] == '*' && str[sz-1] == 'I') {
            size_t plus = str.find(" + ");
            if (plus == std::string::npos) {
                arb_set_str(acb_imagref(acb),std::string(str,0,sz-2).c_str(),prec);
                return;
            }

            arb_set_str(acb_realref(acb),std::string(str,0,plus).c_str(),prec);
            arb_set_str(acb_imagref(acb),std::string(str,plus+3,sz-plus-5).c_str(),prec);
        } else {
            arb_set_str(acb_realref(acb),str.c_str(),prec);
        }
    }

    std::string Acb::mma(long digits) const {
        auto arb2mma = [](const arb_t arb, long digits0) {
            if (arb_is_zero(arb)) return std::string();

            long bits = arb_rel_accuracy_bits(arb);
            if (bits < 0 && arb_contains_zero(arb)) {
                bits = mag_is_zero(arb_radref(arb)) ? ARF_PREC_EXACT : -MAG_EXP(arb_radref(arb));
            }
            if (bits == ARF_PREC_EXACT) {
                fmpz_t n,d;

                fmpz_init(n);
                fmpz_init(d);

                arf_get_fmpz_2exp(n,d,arb_midref(arb));

                long e = fmpz_get_si(d);

                bool isfrac = e < 0;
                if (isfrac) e = -e;

                fmpz_set_ui(d,2);

                fmpz_pow_ui(d,d,e);

                if (!isfrac) fmpz_mul(n,n,d);

                char *cstr;

                cstr = fmpz_get_str(NULL,10,n);
                std::string s = cstr;
                flint_free(cstr);

                if (isfrac) {
                    cstr = fmpz_get_str(NULL,10,d);
                    s += "/" + std::string(cstr);
                    flint_free(cstr);
                }

                fmpz_clear(n);
                fmpz_clear(d);

                return s;
            }

            long digits = 0.30103f*bits;

            if (digits0 < digits) digits = digits0;

            char *cstr = arb_get_str(arb,(long)digits,ARB_STR_NO_RADIUS);
            std::string s = (cstr[0] == '[') ? "0" : cstr;
            flint_free(cstr);

            auto epos = s.find('e');
            std::string expo;

            if (epos != std::string::npos) {
                expo = s.substr(epos+1);
                s = s.substr(0,epos);
            }

            s += "`" + std::to_string(digits);

            if (!expo.empty()) s += "*^" + expo;

            return s;
        };


        std::string real = arb2mma(acb_realref(acb),digits);
        std::string imag = arb_contains_zero(acb_imagref(acb)) ? "" : arb2mma(acb_imagref(acb),digits);

        if (real.empty() && imag.empty()) return "0";

        if (!imag.empty()) {
            if (imag == "1") {
                imag = "I";
            } else if (imag == "-1") {
                imag = "-I";
            } else {
                imag = imag+"*I";
            }
        }

        if (imag.empty()) return real;
        if (real.empty()) return imag;

        if (imag[0] != '-') imag = "+" + imag;

        return real + imag;
    }
}


