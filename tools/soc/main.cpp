#include <CppLibrary/argparse.hpp>
#include <CppLibrary/utility.hpp>

#include <vibron/wfn.hpp>

void perturb(const std::vector<std::vector<std::pair<double, vibron::Wfn *>>> & energies_wfns,
const double & x, const double & y, const double & z);

argparse::ArgumentParser parse_args(const size_t & argc, const char ** & argv) {
    CL::utility::echo_command(argc, argv, std::cout);
    std::cout << '\n';
    argparse::ArgumentParser parser("Spin-orbit coupling for 2 nonrelativistic doublets");

    // vibronic wave function definition
    parser.add_argument("-w","--wfn",       '+', false, "vibronic wave function definition files of each irreducible");
    parser.add_argument("-v","--vibration", '+', false, "vibrational basis definition files");
    parser.add_argument("-l","--level",     '+', false, "vibronic level list files of each irreducible");

    // spin-orbit coupling
    parser.add_argument("-x","--x", 1, false, "x component of Breit-Pauli spin-orbit coupling");
    parser.add_argument("-y","--y", 1, false, "y component of Breit-Pauli spin-orbit coupling");
    parser.add_argument("-z","--z", 1, false, "z component of Breit-Pauli spin-orbit coupling");

    parser.parse_args(argc, argv);
    return parser;
}

int main(size_t argc, const char ** argv) {
    std::cout << "Spin-orbit coupling for 2 nonrelativistic doublets\n"
              << "Yifan Shen 2021\n\n";
    argparse::ArgumentParser args = parse_args(argc, argv);
    CL::utility::show_time(std::cout);
    std::cout << '\n';

    auto wfn_files = args.retrieve<std::vector<std::string>>("wfn");
    auto vib_files = args.retrieve<std::vector<std::string>>("vibration");
    std::vector<std::shared_ptr<vibron::Options>> ops(wfn_files.size());
    for (size_t i = 0; i < ops.size(); i++) {
        ops[i] = std::make_shared<vibron::Options>(wfn_files[i], vib_files);
        std::cout << "Vibronic wave function definition for irreducible " << i + 1 << ":\n";
        ops[i]->pretty_print(std::cout);
        std::cout << '\n';
    }

    auto level_files = args.retrieve<std::vector<std::string>>("level");
    if (level_files.size() != wfn_files.size()) throw std::invalid_argument(
    "inconsistent number between vibronic wave function definitions and vibronic level lists");
    std::vector<std::vector<std::pair<double, vibron::Wfn *>>> energies_wfns(level_files.size());
    for (size_t i = 0; i < level_files.size(); i++) {
        std::ifstream ifs; ifs.open(level_files[i]);
        while (true) {
            std::string line;
            std::getline(ifs, line);
            if (! ifs.good()) break;
            auto strs = CL::utility::split(line);
            std::pair<double, vibron::Wfn *> energy_wfn;
            energy_wfn.first = std::stod(strs[0]);
            energy_wfn.second = new vibron::Wfn(ops[i]);
            energy_wfn.second->read(strs[1]);
            energies_wfns[i].push_back(energy_wfn);
        }
        ifs.close();
    }

    double x = args.retrieve<double>("x"),
           y = args.retrieve<double>("y"),
           z = args.retrieve<double>("z");

    perturb(energies_wfns, x, y, z);

    std::cout << '\n';
    CL::utility::show_time(std::cout);
    std::cout << "Mission success\n";
}