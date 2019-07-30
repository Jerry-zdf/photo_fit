#include "control_data.h"

#include <regex>

Control_data Control_data::parse_input_file(const std::string &input_file) {
    std::ifstream file(input_file);
    if (!file.is_open())
        throw std::runtime_error("Input file is not open!");

    const auto keys = read_keys(file);
    Control_data cd;

    auto set_unique_string = [&](const std::string &key, std::string &val) {
        const auto search = keys.find(key);
        if (search != keys.end())
            val = search->second.at(0);
    };

    set_unique_string("INPUT_FILE", cd.input_file_pattern);
    set_unique_string("OUTPUT_FILE", cd.output_file_pattern);
    set_unique_string("OUTPUT_PATH", cd.out_path);
    set_unique_string("INPUT_PATH", cd.in_path);
    set_unique_string("BASIS_NAME", cd.basis_name);

    auto set_unique_int = [&](const std::string &key, int &val) {
        const auto search = keys.find(key);
        if (search != keys.end())
            val = std::stoi(search->second.at(0));
    };

    set_unique_int("CONTRACTION_SIZE", cd.contraction_size);
    set_unique_int("K_PRECISION", cd.k_precision);
    set_unique_int("MAX_L", cd.max_l);

    {
        const auto search = keys.find("K_DIRECTION");
        if (search != keys.end()) {
            cd.k_dir(0) = std::stod(search->second.at(0));
            cd.k_dir(1) = std::stod(search->second.at(1));
            cd.k_dir(2) = std::stod(search->second.at(2));
        }
        cd.k_dir /= cd.k_dir.norm();
    }

    file.close();
    return cd;
}

std::map<std::string, std::vector<std::string>> Control_data::read_keys(std::ifstream &file) {
    file.seekg(0, std::ios::beg);

    std::string line;
    std::map<std::string, std::vector<std::string>> keys;

    while (std::getline(file, line)) {
        if (line.empty())
            continue;

        const std::regex reg("\\s+");

        std::sregex_token_iterator beg(line.begin(), line.end(), reg, -1);
        std::sregex_token_iterator end;

        const std::string key = *beg;
        const std::vector<std::string> vec(++beg, end);
        keys.emplace(std::make_pair(key, vec));
    }

    return keys;
}

std::ostream &operator<<(std::ostream &os, const Control_data &rhs) {
    os << "============================================================================\n";
    os << "BASIS_NAME                           " << rhs.basis_name << '\n';
    os << "CONTRACTION_SIZE                     " << rhs.contraction_size << '\n';
    os << "============================================================================\n";
    os << "INPUT_FILE                           " << rhs.input_file_pattern << '\n';
    os << "OUTPUT_FILE                          " << rhs.output_file_pattern << '\n';
    os << "OUTPUT_PATH                          " << rhs.out_path << '\n';
    os << "INPUT_PATH                           " << rhs.in_path << '\n';
    os << "============================================================================\n";
    os << "K_PRECISION                          " << rhs.k_precision << '\n';
    os << "K_DIRECTION                          " << rhs.k_dir.transpose() << '\n';
    os << "MAX_L                                " << rhs.max_l << '\n';
    os << "============================================================================\n";
    return os;
}