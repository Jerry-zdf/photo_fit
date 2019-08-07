#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>

#include "basis.h"
#include "control_data.h"
#include "procedures.h"
#include "utils.h"

using namespace Eigen;
using namespace std;

#ifdef FIT_DEBUG
int main() {
    int argc = 4;
    char* argv[4];
    argv[0] = "./photo_fit";
    argv[1] = "config.txt";
    argv[2] = "parameters.txt";
    argv[3] = "0.05";
#else
int main(int argc, char* argv[]) {
#endif

    Clock clk;

    if (argc != 4) {
        cout << " Usage: ./photo_fit <config_file> <evol_parameters_file> <k_val>\n";
        return EXIT_SUCCESS;
    }

    cout << scientific;
    const auto control = Control_data::parse_input_file(argv[1]);

    stringstream ss;
    ss << fixed << setprecision(control.k_precision) << stod(argv[3]);
    const string kstr = ss.str();
    const double kval = stod(kstr);

    const string input_with_k  = regex_replace(control.input_file_pattern, regex("<k>"), kstr);
    const string output_with_k = regex_replace(control.output_file_pattern, regex("<k>"), kstr);

    cout << control
         << "\n\n";

    Atom atom;
    int min_l = 0;
    {
        ifstream out_file(control.out_path + "/" + output_with_k);
        if (out_file) {
            cout << " Found old output file.  We will make use of that\n";
            atom.read(out_file);
            Basis bs;
            bs.atoms.push_back(atom);
            min_l = shell_to_int(bs.get_max_shell()) + 1;
        }
    }

    cout << " Reading following files:\n";

    vector<string> input_files;
    for (int l = min_l; l <= control.max_l; ++l) {
        input_files.emplace_back(regex_replace(input_with_k, regex("<l>"), std::to_string(l)));
        cout << " " << input_files.back() << '\n';
    }
    cout << "\n\n";

    atom.label = control.basis_name;

    for (int l = min_l; l <= control.max_l; ++l) {
        const auto fit = get_gaussian_functor_from_file(input_files[l - min_l], control, l);

        cout << "\n\n"
             << "    STARTING EVOLUTIONARY ROUTINE     "
             << "\n\n";

        char* ecf_command[] = {argv[0], argv[2]};
        VectorXd lm_sol     = get_evolutionary_minimizer(fit, ecf_command);

        cout
            << "\n\n"
            << " STARTING LEVENBERG-MARQUARDT ROUTINE "
            << "\n\n";

        const auto ls_sol = solve_levenberg_marquardt(fit, lm_sol);

        atom.contractions.emplace_back(make_gtopw_contraction(control, ls_sol, lm_sol, kval, l));
    }

    ofstream output(control.out_path + "/" + output_with_k, output.trunc);
    output << atom << '\n';
    output.close();

    cout << " Wall time: " << setprecision(5) << fixed << clk << "\n\n";
    return EXIT_SUCCESS;
}
