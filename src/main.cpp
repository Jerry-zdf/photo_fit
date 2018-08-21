#include <ecf/ECF.h>

#include "MyFloatingPoint.h"
#include "EvalOp.h"
#include "basis.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <ctime>

using namespace Eigen;
using namespace std;

double Norm(double a, int l)
{
    double norm_sqrt = 1.0 / pow(a, 1.5 + l);
    return sqrt(norm_sqrt);
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << " Usage: ./comb_fit <input_file> \n";
        exit(EXIT_SUCCESS);
    }

    clock_t t_start = clock();

    int gauss_n = 10;

    string input_file = argv[1];
    size_t l_pos = input_file.find("l");
    size_t l2_pos = input_file.find(".dat");
    string l_number = input_file.substr(l_pos + 1, l2_pos - l_pos - 1);
    int l = atoi(l_number.c_str());

    size_t k_pos = input_file.find("k");
    size_t k2_pos = input_file.find("_l");
    string kstr = input_file.substr(k_pos + 1, k2_pos - k_pos - 1);
    double kval = atof(kstr.c_str());

    ifstream input;
    input.open(input_file);

    double rn, real, imag;
    vector<double> rvec;
    vector<complex<double>> yvec;
    string buff;

    if (input.is_open())
    {
        while (getline(input, buff))
        {
            stringstream line(buff);
            line >> rn >> real >> imag;
            rvec.push_back(rn);
            yvec.push_back(complex<double>(real, imag));
        }
    }
    else
    {
        cerr << " Unable to open input_file !\n";
        exit(EXIT_FAILURE);
    }
    input.close();

    ArrayXd rv = Map<ArrayXd, Eigen::Unaligned>(rvec.data(), rvec.size());
    ArrayXcd yv = Map<ArrayXcd, Eigen::Unaligned>(yvec.data(), yvec.size());

    GaussianFit fit(gauss_n, rv, yv);
    fit.setL(l);

    cout << "\n\n"
         << "============================="
         << " STARTING EVOLUTIONARY ROUTINE "
         << "============================="
         << "\n\n";

    StateP state(new State);
    MyFloatingPointP gen = static_cast<MyFloatingPointP>(new MyFloatingPoint);
    state->setEvalOp(new EvalOp(fit));
    state->addGenotype(gen);
    char *ecf_command[2] = {"./comb_fit", "parameters.txt"};
    state->initialize(3, ecf_command);
    state->run();

    MyFloatingPoint *best_gene =
        static_cast<MyFloatingPoint *>(state->getHoF()->getBest()[0]->getGenotype().get());

    VectorXd lm_sol(best_gene->realValue.size());
    for (uint i = 0; i < best_gene->realValue.size(); ++i)
        lm_sol[i] = best_gene->realValue[i];

    cout << "\n\n"
         << "============================="
         << " STARTING LEVENBERG-MARQUARDT ROUTINE "
         << "============================="
         << "\n\n";

    cout << " Starting solution:\n " << lm_sol.transpose() << "\n";

    NumericalDiff<GaussianFit> numDiff(fit);
    LevenbergMarquardt<NumericalDiff<GaussianFit>> lm(numDiff);

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

    VectorXd ls_sol;
    fit.generateLeastSqares(lm_sol, ls_sol);
    fit(lm_sol, fvec);
    cout << " Function value at minimun: " << fvec.dot(fvec) << "\n";

    VectorXd expv(gauss_n);
    expv(0) = lm_sol(0);

    for (int k = 1; k < expv.size(); ++k)
    {
        expv[k] = lm_sol(0) * pow(lm_sol(1), k) * (1. + lm_sol(2) * pow((double)k / expv.size(), lm_sol(3)));
    }

    size_t dir_pos = input_file.find("input/");
    input_file.erase(input_file.begin(), input_file.begin() + 6 + dir_pos);

    string output_file = "output2/fit_z1_k" + kstr + ".dat";
    ofstream output(output_file, ios::app);
    if (l == 0)
        output << "CONT 0   0.000000 0.000000 0.000000\n";

    vector<double> exps(gauss_n);
    vector<cdouble> coefs(gauss_n);

    for (int i = 0; i < gauss_n; ++i)
    {
        exps[i] = expv[i];
        coefs[i] = std::complex<double>(ls_sol.head(expv.size())[i], ls_sol.tail(expv.size())[i]);
        coefs[i] *= Norm(exps[i], l);
    }

    GTOPW orbital(exps, coefs, Shell(l), {kval, 0., 0.});
    output << orbital;
    output.close();

    cout << endl
         << setprecision(3);
    clock_t t_end = clock();
    float duration = t_end - t_start;
    cout << " CPU time: " << duration / CLOCKS_PER_SEC << "\n";
}
