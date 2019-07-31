#include "procedures.h"

#include <iomanip>
#include <iostream>

#include <eigen3/unsupported/Eigen/LevenbergMarquardt>

#include <ecf/ECF.h>

#include "ecf_extension/EvalOp.h"
#include "ecf_extension/MyFloatingPoint.h"

using namespace Eigen;
using namespace std;

VectorXd solve_levenberg_marquardt(const Gaussian_fit& fit, VectorXd& vec_to_minimzie) {
    cout << " Starting solution:\n " << setprecision(15) << vec_to_minimzie.transpose() << "\n\n";

    NumericalDiff<Gaussian_fit> num_diff(fit);
    LevenbergMarquardt<NumericalDiff<Gaussian_fit>> lm(num_diff);

    lm.setMaxfev(10000);
    lm.setXtol(1.0e-20);
    VectorXd fvec(fit.values());
    fit(vec_to_minimzie, fvec);

    cout << scientific;

    cout << " Max fun eval:              " << setw(20) << setprecision(12) << lm.maxfev() << "\n"
         << " x tolerance:               " << setw(20) << setprecision(12) << lm.xtol() << "\n"
         << " Init funct val:            " << setw(20) << setprecision(12) << fvec.dot(fvec) << "\n\n"
         << " Minimizing...\n\n";

    const int ret = lm.minimize(vec_to_minimzie);
    fit(vec_to_minimzie, fvec);
    cout << " Function value at minimun: " << setw(20) << setprecision(12) << fvec.dot(fvec) << "\n\n";

    cout << " Iterations count: " << lm.iterations() << "\n"
         << " Routine status:   " << ret << "\n\n";
    cout << " Solution:\n " << setprecision(15) << vec_to_minimzie.transpose() << "\n\n";

    return fit.generate_least_squares(vec_to_minimzie);
}

VectorXd get_evolutionary_minimizer(const Gaussian_fit& fit, char* ecf_command[]) {
    StateP state(new State);
    MyFloatingPointP gen = static_cast<MyFloatingPointP>(new MyFloatingPoint);
    state->setEvalOp(new EvalOp(fit));
    state->addGenotype(gen);
    state->initialize(3, ecf_command);
    state->run();

    MyFloatingPoint* best_gene =
        static_cast<MyFloatingPoint*>(state->getHoF()->getBest()[0]->getGenotype().get());

    return Map<VectorXd>(best_gene->realValue.data(), best_gene->realValue.size());
}

Gaussian_fit get_gaussian_functor_from_file(const string& file_name, const Control_data& control, const int& l) {
    ifstream input(control.in_path + "/" + file_name);

    vector<double> rvec;
    vector<complex<double>> yvec;

    if (input.is_open()) {
        string buff;
        while (getline(input, buff)) {
            stringstream line(buff);
            double rn, real, imag;
            line >> rn >> real >> imag;
            rvec.push_back(rn);
            yvec.push_back(real + 1i * imag);
        }
    } else
        throw runtime_error("Unable to open input_file " + file_name + " !");

    input.close();

    auto rv = Map<ArrayXd, Eigen::Unaligned>(rvec.data(), rvec.size());
    auto yv = Map<ArrayXcd, Eigen::Unaligned>(yvec.data(), yvec.size());

    return Gaussian_fit(control.contraction_size, l, rv, yv);
}

GTOPW_contraction make_gtopw_contraction(const Control_data& control, const VectorXd& ls_sol,
                                         const VectorXd& lm_sol, const double& kval, const int& l) {
    assert(ls_sol.size() == 2 * control.contraction_size);
    assert(lm_sol.size() == 4);

    VectorXd expv(control.contraction_size);
    expv(0) = lm_sol(0);

    for (int k = 1; k < expv.size(); ++k) {
        expv(k) = lm_sol(0) * pow(lm_sol(1), k) * (1. + lm_sol(2) * pow((double)k / expv.size(), lm_sol(3)));
    }
    GTOPW_contraction contr;
    contr.gtopws.reserve(control.contraction_size);
    contr.shl = int_to_shell(l);

    const auto kvec = control.k_dir * kval;

    for (int i = 0; i < control.contraction_size; ++i) {
        GTOPW_primitive prim;
        prim.exp  = expv(i);
        prim.coef = ls_sol.head(expv.size())(i) + 1i * ls_sol.tail(expv.size())(i);
        prim.coef /= sqrt(pow(expv(i), 1.5 + l));  //normalize the coefficients
        prim.k[0] = kvec(0);
        prim.k[1] = kvec(1);
        prim.k[2] = kvec(2);
        contr.gtopws.push_back(prim);
    }
    return contr;
}