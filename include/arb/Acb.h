#pragma once

#include <acb.h>
#include <fmpqxx.h>
#include <yaml-cpp/yaml.h>
#include <sstream>
#include <complex>
#include <cassert>

namespace arb {
    class Acb;
    inline std::ostream &operator<<(std::ostream &os, const Acb &acb);

    class Acb {
        friend std::ostream &operator<<(std::ostream &os, const Acb &acb);
        protected:
            long defPrec;
            acb_t acb;
            static int xalloc1,xalloc2,xalloc3;
        public:
            static Acb infty;
            static Acb nan;
            static Acb I;
            static Acb Pi(long prec);
            static Acb Zeta(const Acb &s, long prec);
            static Acb PolyLog(long s, const Acb &z, long prec);
            static Acb fac(unsigned long n, long prec);
            static Acb binom(unsigned long n, unsigned long k, long prec);

            Acb() : defPrec(ARF_PREC_EXACT) {
                acb_init(acb);
            }

            Acb(long n, long prec=ARF_PREC_EXACT) : defPrec(prec) {
                acb_init(acb);
                acb_set_si(acb,n);
            }

            Acb(unsigned long n, long prec=ARF_PREC_EXACT) : defPrec(prec) {
                acb_init(acb);
                acb_set_ui(acb,n);
            }

            Acb(int n, long prec=ARF_PREC_EXACT) : defPrec(prec) {
                acb_init(acb);
                acb_set_si(acb,n);
            }

            Acb(unsigned int n, long prec=ARF_PREC_EXACT) : defPrec(prec) {
                acb_init(acb);
                acb_set_ui(acb,n);
            }

            Acb(double d, long prec=ARF_PREC_EXACT) : defPrec(prec) {
                acb_init(acb);
                acb_set_d(acb,d);
            }

            Acb(std::complex<double> cmplx, long prec=ARF_PREC_EXACT) : defPrec(prec) {
                acb_init(acb);
                acb_set_d_d(acb,cmplx.real(),cmplx.imag());
            }

            Acb(const flint::fmpqxx &q, long prec) : defPrec(prec) {
                acb_init(acb);
                acb_set_fmpq(acb, q._data().inner, prec);
            }

            Acb(const std::string &str, long prec);

            Acb(const Acb &other) : defPrec(other.defPrec) {
                acb_init(acb);
                acb_set(acb,other.acb);
            }

            Acb(const Acb &a, const Acb &b) {
                acb_init(acb);

                acb_union(acb,a.acb,b.acb,ARF_PREC_EXACT);
                update_default_prec(std::min(a.defPrec,b.defPrec));
            }

            ~Acb() {
                acb_clear(acb);
            }

            acb_t &get() {
                return acb;
            }

            const acb_t &get() const {
                return acb;
            }

            Acb &operator=(const Acb &other) {
                defPrec = other.defPrec;
                acb_set(acb,other.acb);
                return *this;
            }

            explicit operator std::complex<double>() const {
                return std::complex<double>(
                    arf_get_d(arb_midref(acb_realref(acb)),ARF_RND_NEAR),
                    arf_get_d(arb_midref(acb_imagref(acb)),ARF_RND_NEAR));
            }

            explicit operator long() const {
                return arf_get_si(arb_midref(acb_realref(acb)),ARF_RND_NEAR);
            }

            std::string str(long digits=LONG_MAX) const {
                std::string res;

                if (digits == LONG_MAX) {
                    digits = prec();
                    if (digits == LONG_MAX || digits < 0) {
                        digits = (defPrec<LONG_MAX) ? defPrec : 100;
                    }

                    digits *= 0.30103f;
                }

                if (!arb_is_nonzero(acb_imagref(acb)) || arb_is_nonzero(acb_realref(acb))) {
                    char *cstr;
                    if (arb_is_int(acb_realref(acb))) {
                        fmpz_t v;
                        fmpz_init(v);

                        arf_get_fmpz(v,arb_midref(acb_realref(acb)),ARF_RND_NEAR);
                        cstr = fmpz_get_str(NULL,10,v);

                        fmpz_clear(v);
                    } else {
                        cstr = arb_get_str(acb_realref(acb),digits,0);
                    }

                    res = cstr;
                    flint_free(cstr);
                }

                if (arb_is_nonzero(acb_imagref(acb))) {
                    char *cstr;
                    if (arb_is_int(acb_realref(acb))) {
                        fmpz_t v;
                        fmpz_init(v);

                        arf_get_fmpz(v,arb_midref(acb_realref(acb)),ARF_RND_NEAR);
                        cstr = fmpz_get_str(NULL,10,v);

                        fmpz_clear(v);
                    } else {
                        cstr = arb_get_str(acb_realref(acb),digits,0);
                    }

                    if (!res.empty()) res += " + ";
                    res += std::string(cstr)+"*I";
                    flint_free(cstr);
                }

                return res;
            }

            std::string mma(long digits=LONG_MAX) const;

            long prec() const {
                long bits = acb_rel_accuracy_bits(acb);
                if (bits < 0 && contains_zero()) {
                    long rexp = mag_is_zero(arb_radref(acb_realref(acb))) ? -ARF_PREC_EXACT : MAG_EXP(arb_radref(acb_realref(acb)));
                    long iexp = mag_is_zero(arb_radref(acb_imagref(acb))) ? -ARF_PREC_EXACT : MAG_EXP(arb_radref(acb_imagref(acb)));
                    return -std::max(rexp,iexp);
                }
                return bits;
            }

            long default_prec() const {
                return defPrec;
            }

            void update_default_prec(long prec0) {
                defPrec = std::max(prec(),prec0);
                if (defPrec == ARF_PREC_EXACT) defPrec = prec0;
            }

            void trim() {
                acb_trim(acb,acb);
            }

            Acb real() const {
                Acb res = *this;
                arb_zero(acb_imagref(res.acb));
                return res;
            }

            Acb imag() const {
                Acb res = *this;
                arb_zero(acb_realref(res.acb));
                arb_swap(acb_realref(res.acb),acb_imagref(res.acb));
                return res;
            }

            Acb abs() const {
                long prec0 = prec();
                if (prec0 == ARF_PREC_EXACT) prec0 = defPrec;

                Acb res;

                if (arb_contains_zero(acb_realref(acb))) {
                    arb_abs(acb_realref(res.acb),acb_imagref(acb));
                } else if (arb_contains_zero(acb_imagref(acb))) {
                    arb_abs(acb_realref(res.acb),acb_realref(acb));
                } else {
                    acb_abs(acb_realref(res.acb),acb,prec0);
                }

                res.update_default_prec(defPrec);
                return res;
            }

            Acb round() const {
                fmpz_t r,i;

                fmpz_init(r);
                fmpz_init(i);

                arf_get_fmpz(r,arb_midref(acb_realref(acb)),ARF_RND_NEAR);
                arf_get_fmpz(i,arb_midref(acb_imagref(acb)),ARF_RND_NEAR);

                Acb res;
                res.defPrec = defPrec;

                arb_set_fmpz(acb_realref(res.acb),r);
                arb_set_fmpz(acb_imagref(res.acb),i);

                fmpz_clear(r);
                fmpz_clear(i);

                return res;
            }

            Acb ceil() const {
                if (contains_int()) return round();
                Acb res;
                arb_ceil(acb_realref(res.acb),acb_realref(acb),ARF_PREC_EXACT);
                arb_ceil(acb_imagref(res.acb),acb_imagref(acb),ARF_PREC_EXACT);
                res.update_default_prec(defPrec);
                return res;
            }

            Acb floor() const {
                if (contains_int()) return round();
                Acb res;
                arb_floor(acb_realref(res.acb),acb_realref(acb),ARF_PREC_EXACT);
                arb_floor(acb_imagref(res.acb),acb_imagref(acb),ARF_PREC_EXACT);
                res.update_default_prec(defPrec);
                return res;
            }

            bool contains_int() const {
                return acb_contains_int(acb);
            }

            bool is_positive() const {
                return arb_is_positive(acb_realref(acb));
            }

            bool overlaps(const Acb &other) const {
                return acb_overlaps(acb,other.acb);
            }

            bool is_nan() const {
                return arf_is_nan(arb_midref(acb_realref(acb))) || arf_is_nan(arb_midref(acb_imagref(acb)));
            }

            bool is_infty() const {
                return arf_is_inf(arb_midref(acb_realref(acb))) || arf_is_inf(arb_midref(acb_imagref(acb)));
            }

            Acb operator-() const {
                Acb res;

                acb_neg(res.acb,acb);
                res.defPrec = defPrec;

                return res;
            }

            Acb operator+(const Acb &other) const {
                Acb res;

                acb_add(res.acb,acb,other.acb,ARF_PREC_EXACT);
                acb_trim(res.acb,res.acb);
                res.update_default_prec(std::min(defPrec,other.defPrec));

                return res;
            }

            Acb operator-(const Acb &other) const {
                Acb res;

                acb_sub(res.acb,acb,other.acb,ARF_PREC_EXACT);
                acb_trim(res.acb,res.acb);
                res.update_default_prec(std::min(defPrec,other.defPrec));

                return res;
            }

            Acb operator*(const Acb &other) const {
                Acb res;

                acb_mul(res.acb,acb,other.acb,ARF_PREC_EXACT);
                acb_trim(res.acb,res.acb);
                res.update_default_prec(std::min(defPrec,other.defPrec));

                return res;
            }

            Acb operator*(const flint::fmpqxx &q) const {
                long prec0 = prec();
                if (prec0 == ARF_PREC_EXACT) prec0 = defPrec;
                if (prec0 < ARF_PREC_EXACT) prec0 += 16;

                return (*this) * Acb(q,prec0);
            }

            Acb operator/(const Acb &other) const {
                Acb res;

                long prec0 = other.prec();
                if (prec0 == ARF_PREC_EXACT) prec0 = std::min(defPrec,other.defPrec);
                if (prec0 < ARF_PREC_EXACT) prec0 += 16;

                acb_div(res.acb,acb,other.acb,prec0);
                acb_trim(res.acb,res.acb);
                res.update_default_prec(std::min(defPrec,other.defPrec));

                return res;
            }

            Acb &operator+=(const Acb &other) {
                acb_add(acb,acb,other.acb,ARF_PREC_EXACT);
                acb_trim(acb,acb);
                update_default_prec(std::min(defPrec,other.defPrec));

                return *this;
            }

            Acb &operator-=(const Acb &other) {
                acb_sub(acb,acb,other.acb,ARF_PREC_EXACT);
                acb_trim(acb,acb);
                update_default_prec(std::min(defPrec,other.defPrec));

                return *this;
            }

            Acb &operator*=(const Acb &other) {
                acb_mul(acb,acb,other.acb,ARF_PREC_EXACT);
                acb_trim(acb,acb);
                update_default_prec(std::min(defPrec,other.defPrec));

                return *this;
            }

            Acb &operator/=(const Acb &other) {
                long prec0 = other.prec();
                if (prec0 == ARF_PREC_EXACT) prec0 = std::min(defPrec,other.defPrec);
                if (prec0 < ARF_PREC_EXACT) prec0 += 16;

                acb_div(acb,acb,other.acb,prec0);
                acb_trim(acb,acb);
                update_default_prec(std::min(defPrec,other.defPrec));

                return *this;
            }

            bool contains_zero() const {
                if (is_nan()) return false;
                return acb_contains_zero(acb);
            }

            bool operator==(const Acb &other) const {
                return !acb_ne(acb,other.acb);
            }

            bool operator!=(const Acb &other) const {
                return !(*this == other);
            }

            bool operator<(const Acb &other) const {
                if (arb_lt(acb_realref(acb),acb_realref(other.acb))) return true;
                if (arb_gt(acb_realref(acb),acb_realref(other.acb))) return false;
                return arb_lt(acb_imagref(acb),acb_imagref(other.acb));
            }

            bool operator>(const Acb &other) const {
                return other < *this;
            }

            bool operator<=(const Acb &other) const {
                return (*this < other) || (*this == other);
            }

            bool operator>=(const Acb &other) const {
                return other <= *this;
            }

            Acb pow(long n) const {
                Acb res;

                long prec0 = n<0 ? defPrec : ARF_PREC_EXACT;
                if (prec0 < ARF_PREC_EXACT) prec0 += 4;

                acb_pow_si(res.acb,acb,n,prec0);
                res.update_default_prec(defPrec);

                return res;
            }

            Acb pow(const Acb &e) const {
                long prec0;

                if (prec() == ARF_PREC_EXACT && e.prec() == ARF_PREC_EXACT && e.contains_int() && (e.is_positive() || e.contains_zero())) {
                    prec0 = ARF_PREC_EXACT;
                } else {
                    prec0 = std::min(prec(),e.prec());
                    if (prec0 == ARF_PREC_EXACT) {
                        prec0 = std::min(defPrec,e.defPrec);
                    }
                    assert(prec0 < ARF_PREC_EXACT);
                }

                Acb res;

                acb_pow(res.acb,acb,e.acb,prec0);
                res.update_default_prec(std::min(defPrec,e.defPrec));

                return res;
            }

            Acb sqrt() const {
                Acb res;

                long prec0 = prec();
                if (prec0 == ARF_PREC_EXACT) prec0 = defPrec;

                acb_sqrt(res.acb,acb,prec0);
                res.update_default_prec(defPrec);

                return res;
            }

            Acb log() const {
                assert(defPrec != ARF_PREC_EXACT);

                Acb res;
                long prec0 = defPrec+4;

                acb_log(res.acb,acb,prec0);
                res.update_default_prec(defPrec);

                return res;
            }

            Acb exp_pi_i() const {
                assert(defPrec != ARF_PREC_EXACT);

                Acb res;
                long prec0 = defPrec+4;

                acb_exp_pi_i(res.acb,acb,prec0);
                res.update_default_prec(defPrec);

                return res;
            }

            //output manipulators
            static std::ios_base &mma(std::ios_base &os) {
                os.iword(xalloc1) = (os.iword(xalloc1) << 1) | 1;
                return os;
            }

            static std::ios_base &standard(std::ios_base &os) {
                os.iword(xalloc1) = os.iword(xalloc1) << 1;
                return os;
            }

            static std::ios_base &rev(std::ios_base &os) {
                os.iword(xalloc1) = os.iword(xalloc1) >> 1;
                return os;
            }

            struct digits {
                int digs;
                digits(long digs) : digs(digs) {}
                void manip(std::ostream &os) const {
                    os.iword(Acb::xalloc2) = digs;
                }
            };

            static std::ios_base &all(std::ios_base &os) {
                os.iword(xalloc2) = 0;
                return os;
            }
    };

    inline std::ostream &operator<<(std::ostream &os, const Acb &acb) {
        long digits = os.iword(Acb::xalloc2);
        if (!digits) digits = LONG_MAX;
        os << ((os.iword(Acb::xalloc1) & 1) ? acb.mma(digits) : acb.str(digits));
        return os;
    }

    inline std::ostream &operator<<(std::ostream &os, const Acb::digits &m) {
        m.manip(os);
        return os;
    }

    inline Acb operator+(long n, const Acb &acb) {
        return Acb(n)+acb;
    }

    inline Acb operator-(long n, const Acb &acb) {
        return Acb(n)-acb;
    }

    inline Acb operator*(long n, const Acb &acb) {
        return Acb(n)*acb;
    }

    inline Acb operator/(long n, const Acb &acb) {
        return Acb(n)/acb;
    }

}

namespace std {
    template<> struct hash<arb::Acb> {
        size_t operator()(const arb::Acb &acb) const noexcept {
            auto cmplx = std::complex<double>(acb);

            size_t rhash = hash<double>()(cmplx.real());
            size_t ihash = hash<double>()(cmplx.imag());

            return rhash ^ (ihash + 0x9e3779b9 + (rhash << 6) + (rhash >> 2));
        }
    };
}

namespace YAML {
    template<> struct convert<arb::Acb> {
        static long default_prec;

        static Node encode(const arb::Acb &acb) {
            Node node;
            node = acb.str();
            return node;
        }

        static bool decode(const Node &node, arb::Acb &res) {
            if (!node.IsScalar()) return false;
            res = arb::Acb(node.as<std::string>(),default_prec);
            return true;
        }
    };
}

