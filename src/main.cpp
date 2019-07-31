#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>

#include <ecf/ECF.h>

#include "ecf_extension/EvalOp.h"
#include "ecf_extension/MyFloatingPoint.h"

#include "basis.h"
#include "control_data.h"
#include "gaussian_functor.h"
#include "utils.h"

using namespace Eigen;
using namespace std;

double Norm(double a, int l) {
    double norm_sqrt = 1.0 / pow(a, 1.5 + l);
    return sqrt(norm_sqrt);
}

int main(int argc, char *argv[]) {
    Clock clk;

    if (argc != 4) {
        cout << " Usage: ./photo_fit <config_file> <evol_parameters_file> <k_val>\n";
        return EXIT_SUCCESS;
    }

    const string kstr = argv[3];

    const auto control = Control_data::parse_input_file(argv[1]);

    const string input_with_k  = regex_replace(control.input_file_pattern, regex("<k>"), kstr);
    const string output_with_k = regex_replace(control.output_file_pattern, regex("<k>"), kstr);

    cout << control
         << "\n\n";

    cout << " Reading following files:\n";

    vector<string> input_files;
    regex l_reg("<l>");
    for (int l = 0; l <= control.max_l; ++l) {
        input_files.emplace_back(regex_replace(input_with_k, l_reg, std::to_string(l)));
        cout << " " << input_files.back() << '\n';
    }
    cout << "\n\n";

    const double kval = std::stod(kstr);

    Atom atom;
    atom.label = control.basis_name;

    for (int l = 0; l <= control.max_l; ++l) {
        ifstream input(control.in_path + "/" + input_files[l]);

        vector<double> rvec;
        vector<complex<double>> yvec;
        string buff;

        if (input.is_open()) {
            while (getline(input, buff)) {
                stringstream line(buff);
                double rn, real, imag;
                line >> rn >> real >> imag;
                rvec.push_back(rn);
                yvec.push_back(complex<double>(real, imag));
            }
        } else
            throw runtime_error("Unable to open input_file !");

        input.close();

        ArrayXd rv  = Map<ArrayXd, Eigen::Unaligned>(rvec.data(), rvec.size());
        ArrayXcd yv = Map<ArrayXcd, Eigen::Unaligned>(yvec.data(), yvec.size());

        Gaussian_fit fit(control.contraction_size, l, rv, yv);

        cout << "\n\n"
             << "======================================"
             << "    STARTING EVOLUTIONARY ROUTINE     "
             << "======================================"
             << "\n\n";

        StateP state(new State);
        MyFloatingPointP gen = static_cast<MyFloatingPointP>(new MyFloatingPoint);
        state->setEvalOp(new EvalOp(fit));
        state->addGenotype(gen);
        char *ecf_command[] = {"./photo_fit", argv[2]};
        state->initialize(3, ecf_command);
        state->run();

        MyFloatingPoint *best_gene =
            static_cast<MyFloatingPoint *>(state->getHoF()->getBest()[0]->getGenotype().get());

        VectorXd lm_sol(best_gene->realValue.size());
        for (uint i = 0; i < best_gene->realValue.size(); ++i)
            lm_sol[i] = best_gene->realValue[i];

        cout << "\n\n"
             << "======================================"
             << " STARTING LEVENBERG-MARQUARDT ROUTINE "
             << "======================================"
             << "\n\n";

        cout << " Starting solution:\n " << lm_sol.transpose() << "\n";

        NumericalDiff<Gaussian_fit> num_diff(fit);
        LevenbergMarquardt<NumericalDiff<Gaussian_fit>> lm(num_diff);

        lm.setMaxfev(10000);
        lm.setXtol(1.0e-20);
        VectorXd fvec(fit.values());
        fit(lm_sol, fvec);

        cout << " Max fun eval:     " << lm.maxfev() << "\n";
        cout << " x tolerance:      " << lm.xtol() << "\n";
        cout << " Init funct val:   " << fvec.dot(fvec) << "\n";

        int ret = lm.minimize(lm_sol);
        cout << " Iterations count: " << lm.iterations() << "\n";
        cout << " Routine status:   " << ret << "\n\n";
        cout << setprecision(15);
        cout << " Solution:\n " << lm_sol.transpose() << "\n";

        const auto ls_sol = fit.generate_least_squares(lm_sol);
        fit(lm_sol, fvec);
        cout << " Function value at minimun: " << fvec.dot(fvec) << "\n";

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
            prim.coef *= Norm(expv(i), l);
            prim.k[0] = kvec(0);
            prim.k[1] = kvec(1);
            prim.k[2] = kvec(2);
            contr.gtopws.push_back(prim);
        }

        atom.contractions.push_back(contr);
    }

    ofstream output(control.out_path + "/" + output_with_k);
    output << atom << '\n';
    output.close();

    cout << " Wall time: " << setprecision(5) << fixed << clk << "\n\n";
    return EXIT_SUCCESS;
}
