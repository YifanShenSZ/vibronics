#include <CppLibrary/argparse.hpp>
#include <CppLibrary/chemistry.hpp>

#include <tchem/utility.hpp>
#include <tchem/intcoord.hpp>

#include <plot/wfn.hpp>

argparse::ArgumentParser parse_args(const size_t & argc, const char ** & argv) {
    CL::utility::echo_command(argc, argv, std::cout);
    std::cout << '\n';
    argparse::ArgumentParser parser("Plot vibronic wave function along 2 directions");

    // Vibronic wave function definition
    parser.add_argument("-w","--wfn",         1, false, "vibronic wave function definition file");
    parser.add_argument("-v","--vibration", '+', false, "vibrational basis definition files");
    parser.add_argument("-f","--frequency", '+', false, "frequency of each normal mode");
    parser.add_argument("-p","--prefix",      1, false, "prefix of the seed or the check point vector to continue with");

    // Coordiante definition
    parser.add_argument("-F","--format", 1, false, "internal coordinate definition format (Columbus7, default)");
    parser.add_argument("-i","--IC",     1, false, "internal coordinate definition file");
    parser.add_argument("-s","--SAS",    1, false, "symmetry adaptation and scale definition file");
    parser.add_argument("-l","--L^-1", '+', false, "The transformation matrix from internal coordiante to normal coordinate");
    parser.add_argument("-g","--geom",   1, false, "geometry");
    parser.add_argument("-x","--x",      1, false, "x direction");
    parser.add_argument("-y","--y",      1, false, "y direction");

    // Optional argumentes
    parser.add_argument("-1","--stepx", 2, true, "x direction number of steps and step length (default = 10, 0.01)");
    parser.add_argument("-2","--stepy", 2, true, "y direction number of steps and step length (default = 10, 0.01)");

    parser.parse_args(argc, argv);
    return parser;
}

int main(size_t argc, const char ** argv) {
    std::cout << "Plot vibronic wave function along 2 directions\n"
              << "Yifan Shen 2021\n\n";
    argparse::ArgumentParser args = parse_args(argc, argv);
    CL::utility::show_time(std::cout);
    std::cout << '\n';

    auto wfn_file  = args.retrieve<std::string>("wfn");
    auto vib_files = args.retrieve<std::vector<std::string>>("vibration");
    auto op = std::make_shared<vibron::Options>(wfn_file, vib_files);

    auto freq_files = args.retrieve<std::vector<std::string>>("frequency");
    std::vector<std::vector<double>> freqs(freq_files.size());
    for (size_t i = 0; i < freq_files.size(); i++)
    freqs[i] = CL::utility::read_vector(freq_files[i]);

    plot::Wfn wfn(op, freqs);
    auto prefix = args.retrieve<std::string>("prefix");
    wfn.read(prefix);

    std::string format = args.retrieve<std::string>("format"),
                IC     = args.retrieve<std::string>("IC"),
                SAS    = args.retrieve<std::string>("SAS");
    auto sasicset = std::make_shared<tchem::IC::SASICSet>(format, IC, SAS);

    auto Linv_files = args.retrieve<std::vector<std::string>>("L^-1");
    std::vector<at::Tensor> Linvs(Linv_files.size());
    for (size_t i = 0; i < Linv_files.size(); i++) {
        Linvs[i] = tchem::utility::read_vector(Linv_files[i]);
        int64_t size = sqrt(Linvs[i].numel());
        Linvs[i] = Linvs[i].view({size, size});
    }

    std::string geom_file = args.retrieve<std::string>("geom");
    CL::chem::xyz<double> geom(geom_file, true);
    std::vector<double> coords = geom.coords();
    at::Tensor r0 = at::from_blob(coords.data(), coords.size(), at::TensorOptions().dtype(torch::kFloat64));
    
    std::string x_file = args.retrieve<std::string>("x"),
                y_file = args.retrieve<std::string>("y");
    at::Tensor x = tchem::utility::read_vector(x_file),
               y = tchem::utility::read_vector(y_file);
    if (x.size(0) != r0.size(0)) throw std::invalid_argument(
    "x must be a 3 * NAtoms dimensional vector");
    if (y.size(0) != r0.size(0)) throw std::invalid_argument(
    "y must be a 3 * NAtoms dimensional vector");
    x /= x.norm();
    y /= y.norm();

    int64_t NSteps_x = 10, NSteps_y = 10;
    double dx = 0.01, dy = 0.01;
    if (args.gotArgument("stepx")) {
        auto NSteps_d = args.retrieve<std::vector<std::string>>("stepx");
        NSteps_x = std::stoi(NSteps_d[0]);
        dx = std::stod(NSteps_d[1]);
    }
    if (args.gotArgument("stepy")) {
        auto NSteps_d = args.retrieve<std::vector<std::string>>("stepy");
        NSteps_y = std::stoi(NSteps_d[0]);
        dy = std::stod(NSteps_d[1]);
    }

    std::vector<std::ofstream> ofs(op->NStates);
    for (size_t i = 0; i < op->NStates; i++) ofs[i].open("state-" + std::to_string(i + 1) + ".txt");

    for (int64_t i = -NSteps_x; i <= NSteps_x; i++)
    for (int64_t j = -NSteps_y; j <= NSteps_y; j++) {
        at::Tensor r = r0 + i * dx * x + j * dy * y;
        auto qs = (*sasicset)(sasicset->tchem::IC::IntCoordSet::operator()(r));
        std::vector<std::vector<double>> Qs(qs.size());
        for (size_t irred = 0; irred < qs.size(); irred++) {
            at::Tensor Q = Linvs[irred].mv(qs[irred]);
            Qs[irred].resize(Q.numel());
            std::memcpy(Qs[irred].data(), Q.data_ptr<double>(), Q.numel() * sizeof(double));
        }
        auto value = wfn(Qs);
        for (size_t state = 0; state < op->NStates; state++) ofs[state] << value[state] << '\n';
    }

    for (auto & of : ofs) of.close();

    std::cout << '\n';
    CL::utility::show_time(std::cout);
    std::cout << "Mission success\n";
}