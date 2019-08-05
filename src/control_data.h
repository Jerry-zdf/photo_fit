#pragma once

#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <eigen3/Eigen/Dense>

class Control_data {
   public:
    std::string basis_name{"CONT"};
    std::string input_file_pattern{"z1_k<k>_l<l>.dat"};
    std::string output_file_pattern{"fit_z1_k<k>.dat"};
    std::string out_path{"./output"};
    std::string in_path{"./input"};

    Eigen::Vector3d k_dir{0.0, 0.0, 1.0};

    int k_precision{3};
    int contraction_size{10};
    int max_l{6};

    bool use_k{true};

    static Control_data parse_input_file(const std::string &input_file);

   private:
    static std::map<std::string, std::vector<std::string>> read_keys(std::ifstream &input_file);
};

std::ostream &operator<<(std::ostream &os, const Control_data &rhs);
