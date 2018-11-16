#ifndef BASIS_H
#define BASIS_H

#include <array>
#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cassert>

using Vec3d   = std::array<double, 3>;
using cdouble = std::complex<double>;


class Shell {
   public:
    static const std::array<int, 11> crt_siz;
    static const std::array<int, 11> labels;
    static const std::map<char, int> charmap;

    enum Shells { S = 0,
                  P = 1,
                  D = 2,
                  F = 3,
                  G = 4,
                  H = 5,
                  I = 6,
                  K = 7,
                  L = 8
    };

    friend Shell char2shell(const char &c);
    friend char shell2char(const Shell &shell);
    friend int shell2int(const Shell &shell);

    Shell(const Shell &shell) = default;
    Shell(const int &shell) {
        shl = static_cast<Shells>(shell);
    }

   private:
    Shells shl;
};

class GTOPW {
   public:
    GTOPW(const std::vector<double> &exponents,
          const std::vector<cdouble> &coefficients,
          const Shell &shell,
          const Vec3d kvec)
        : shl(shell) {
        assert(exponents.size() == coefficients.size());

        size  = exponents.size();
        exps  = exponents;
        coefs = coefficients;
        k     = kvec;
    };

    GTOPW() : exps(), coefs(), shl(0), size(0) {
        k = {0, 0, 0};
    };

    friend std::ostream &operator<<(std::ostream &os, const GTOPW &rhs);
    bool read(std::istream &is);
    int functions_number() const;
    inline void set_kvec(const Vec3d &kvec) { k = kvec; };

   private:
    std::vector<double> exps;
    std::vector<cdouble> coefs;
    Shell shl;
    Vec3d k;
    int size;
};

class Basis {
   public:
    Basis() : gtopws(), label(), charge(0.) {
        position = {0, 0, 0};
    };

    bool read(std::istream &is);
    friend std::ostream &operator<<(std::ostream &os, const Basis &rhs);
    int functions_number() const;
    inline void set_position(const Vec3d &pos) { position = pos; };
    inline void set_label(const std::string &lb) { label = lb; };

   private:
    std::vector<GTOPW> gtopws;
    std::string label;
    double charge;
    Vec3d position;
};

#endif