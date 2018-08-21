#ifndef GAUSSIAN_FIT_FUNCTION_H
#define GAUSSIAN_FIT_FUNCTION_H

#include "unsupported/Eigen/LevenbergMarquardt"
#include <Eigen/Dense>


struct GaussianFit : Eigen::DenseFunctor<double> {

public:

    typedef Eigen::ArrayXcd ComplexData;
    typedef Eigen::ArrayXd DoubleData;
    typedef Eigen::Vector4d Parameters;

    GaussianFit(int _GaussNumber, DoubleData _rv, ComplexData _yv);
    GaussianFit(const GaussianFit&) = default;

    int operator()(const InputType &x, ValueType &fvec) const;

    void generateLeastSqares(const Parameters &x, Eigen::VectorXd& sol) const;
    void setL(int _l);

private:

    ComplexData yv;
    DoubleData rv;
    int l = 0;
    int GaussNumber;
};

#endif // GAUSSIAN_FIT_FUNCTION_H
