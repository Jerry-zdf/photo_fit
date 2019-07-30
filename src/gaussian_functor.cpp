#include "gaussian_functor.h"

#include <cassert>
#include <cmath>

using std::exp;
using std::pow;

Gaussian_fit::Gaussian_fit(const int &gauss_cout, const int &l, const Eigen::ArrayXd &grid, const Eigen::ArrayXcd &function_values)
    : Eigen::DenseFunctor<double>(4, 2 * grid.size()), _yv(function_values), _rv(grid), _l(l), _gauss_cout(gauss_cout) {
    assert(_rv.size() == _yv.size());
    assert(_gauss_cout <= _rv.size());
}

int Gaussian_fit::operator()(const InputType &x, ValueType &fvec) const {
    assert(x.size() == this->inputs());

    fvec.resize(this->values());
    fvec.head(_rv.size()) = _yv.real();
    fvec.tail(_rv.size()) = _yv.imag();

    auto ls = generate_least_squares(x);

    for (int k = 0; k < _gauss_cout; ++k) {
        const double ed = x(0) * pow(x(1), k) * (1. + x(2) * pow((double)k / _gauss_cout, x(3)));
        fvec.head(_rv.size()).array() -= ls(k) * exp(-ed * pow(_rv, 2)) * pow(_rv, _l);
        fvec.tail(_rv.size()).array() -= ls(_gauss_cout + k) * exp(-ed * pow(_rv, 2)) * pow(_rv, _l);
    }
    return 0;
}

Eigen::VectorXd Gaussian_fit::generate_least_squares(const Eigen::Vector4d &x) const {
    Eigen::MatrixXd mat(_rv.size(), _gauss_cout);

    for (int k = 0; k < _gauss_cout; ++k) {
        const double ed    = x(0) * pow(x(1), k) * (1. + x(2) * pow((double)k / _gauss_cout, x(3)));
        mat.col(k).array() = exp(-ed * pow(_rv, 2)) * pow(_rv, _l);
    }

    Eigen::VectorXd sol;
    sol.resize(2 * _gauss_cout);
    sol.head(_gauss_cout) = (mat.transpose() * mat).ldlt().solve(mat.transpose() * _yv.real().matrix());
    sol.tail(_gauss_cout) = (mat.transpose() * mat).ldlt().solve(mat.transpose() * _yv.imag().matrix());
    //    sol.head(_gauss_cout) = mat.colPivHouseholderQr().solve(yv.real().matrix());
    //    sol.tail(_gauss_cout) = mat.colPivHouseholderQr().solve(yv.imag().matrix());
    return sol;
}
