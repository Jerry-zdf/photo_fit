#include "GaussianFitFunctor.h"

#include <cmath>
#include <cassert>


using std::pow;
using std::exp;


GaussianFit::GaussianFit(int _GaussNumber, DoubleData _rv, ComplexData _yv)
    : Eigen::DenseFunctor<double>(4, 2 * _rv.size()) {
    assert(_rv.size() == _yv.size());
    assert(_GaussNumber <= _rv.size());

    GaussNumber = _GaussNumber;
    rv = _rv;
    yv = _yv;
}


void GaussianFit::setL(int _l) {
    assert(_l >= 0);
    l = _l;
}


int GaussianFit::operator()(const InputType &x, ValueType &fvec) const {
    assert(x.size() == this->inputs());
    assert(fvec.size() == this->values());

    Eigen::VectorXd ls;
    this->generateLeastSqares(x, ls);

    double ed;
    fvec.head(rv.size()) = yv.real();
    fvec.tail(rv.size()) = yv.imag();

    for(int k = 0; k < GaussNumber; ++k) {
        ed = x(0) * pow(x(1), k) * (1. + x(2) * pow((double) k / GaussNumber, x(3)));

        fvec.head(rv.size()).array() -= ls(k)               * exp(-ed * pow(rv, 2)) * pow(rv, l);
        fvec.tail(rv.size()).array() -= ls(GaussNumber + k) * exp(-ed * pow(rv, 2)) * pow(rv, l);
    }
    return 0;
}


void GaussianFit::generateLeastSqares(const Parameters& x, Eigen::VectorXd &sol) const {
    Eigen::MatrixXd mat(rv.size(), GaussNumber);
    double ed;

    for(int k = 0; k < GaussNumber; ++k)
    {
        ed = x(0) * pow(x(1), k) * (1. + x(2) * pow((double) k / GaussNumber, x(3)));
        mat.col(k).array() = exp(-ed * pow(rv, 2)) * pow(rv, l);
    }

    sol.resize(2 * GaussNumber);
    sol.head(GaussNumber) = (mat.transpose() * mat).ldlt().solve(mat.transpose() * yv.real().matrix());
    sol.tail(GaussNumber) = (mat.transpose() * mat).ldlt().solve(mat.transpose() * yv.imag().matrix());
    //    sol.head(GaussNumber) = mat.colPivHouseholderQr().solve(yv.real().matrix());
    //    sol.tail(GaussNumber) = mat.colPivHouseholderQr().solve(yv.imag().matrix());
}


