#pragma once

#include <Eigen/Dense>
#include <eigen3/unsupported/Eigen/LevenbergMarquardt>

struct Gaussian_fit : Eigen::DenseFunctor<double> {
   public:
    Gaussian_fit(const int &gauss_cout, const int &l, const Eigen::ArrayXd &grid, const Eigen::ArrayXcd &function_values);
    Gaussian_fit(const Gaussian_fit &) = default;

    int operator()(const InputType &x, ValueType &fvec) const;

    Eigen::VectorXd generate_least_squares(const Eigen::Vector4d &x) const;

   private:
    const Eigen::ArrayXcd _yv;
    const Eigen::ArrayXd _rv;
    const int _l;
    const int _gauss_cout;
};
