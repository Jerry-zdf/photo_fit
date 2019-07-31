#pragma once

#include <string>

#include <eigen3/Eigen/Dense>

#include "basis.h"
#include "control_data.h"
#include "gaussian_functor.h"

Eigen::VectorXd solve_levenberg_marquardt(const Gaussian_fit& fit, Eigen::VectorXd& vec_to_minimzie);

Eigen::VectorXd get_evolutionary_minimizer(const Gaussian_fit& fit, char* ecf_command[]);

Gaussian_fit get_gaussian_functor_from_file(const std::string& file_name, const Control_data& control, const int& l);

GTOPW_contraction make_gtopw_contraction(const Control_data& control, const Eigen::VectorXd& ls_sol,
                                         const Eigen::VectorXd& lm_sol, const double& kval, const int& l);