#include <CppLibrary/argparse.hpp>
#include <CppLibrary/utility.hpp>

#include <vibron/options.hpp>

#include <Lanczos/mv.hpp>
#include <Lanczos/Lanczos.hpp>

argparse::ArgumentParser parse_args(const size_t & argc, const char ** & argv) {
    CL::utility::echo_command(argc, argv, std::cout);
    std::cout << '\n';
    argparse::ArgumentParser parser("Vibronics version 2");

    // required arguments
    parser.add_argument("-w","--wfn",         1, false, "vibronic wave function definition file");
    parser.add_argument("-v","--vibration", '+', false, "vibrational basis definition files");
    parser.add_argument("-H","--Hd",        '+', false, "anharmonic diabatic Hamiltonian definition files");
    parser.add_argument("-f","--frequency", '+', false, "frequency of each normal mode");
    parser.add_argument("-p","--prefix",      1, false, "prefix of the seed or the check point vector to continue with");

    // optional arguments
    parser.add_argument("-a","--alpha",         1, true, "alpha (diag) to continue with");
    parser.add_argument("-b","--beta",          1, true, "beta (subdiag) to continue with");
    parser.add_argument("--wprefix",            1, true, "prefix of the w vector to continue with");
    parser.add_argument("-m","--max_iteration", 1, true, "max number of iterations to perform (default = 100)");

    parser.parse_args(argc, argv);
    return parser;
}

int main(size_t argc, const char ** argv) {
    std::cout << "Vibronics version 2\n"
              << "This version reorthogonalizes the Krylov vectors on hard disk\n"
              << "Yifan Shen 2021\n\n";
    argparse::ArgumentParser args = parse_args(argc, argv);
    CL::utility::show_time(std::cout);
    std::cout << '\n';

    auto wfn_file  = args.retrieve<std::string>("wfn");
    auto vib_files = args.retrieve<std::vector<std::string>>("vibration");
    auto op = std::make_shared<vibron::Options>(wfn_file, vib_files);
    op->pretty_print(std::cout);
    std::cout << std::endl;

    auto Hd_files = args.retrieve<std::vector<std::string>>("Hd");
    auto Hd = std::make_shared<Lanczos::Hd>(op->NStates, op->NModes, Hd_files);
    Hd->pretty_print(std::cout);
    std::cout << std::endl;

    auto freq_files = args.retrieve<std::vector<std::string>>("frequency");
    auto mvkernel = std::make_shared<Lanczos::MVKernel>(Hd, op, freq_files);
    Lanczos::Kernel kernel(mvkernel);
    kernel.pretty_print(std::cout);
    std::cout << std::endl;

    vibron::Wfn v1(op), v2(op), w(op);

    auto prefix = args.retrieve<std::string>("prefix");
    std::vector<std::vector<std::fstream>> fs(op->NSegs);
    for (size_t i = 0; i < op->NSegs; i++) {
        fs[i].resize(op->NStates);
        for (size_t j = 0; j < op->NStates; j++) {
            std::string file = prefix + "-" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + ".wfn";
            fs[i][j].open(file, std::fstream::in | std::fstream::out | std::fstream::binary);
            if (! fs[i][j].good()) throw CL::utility::file_error(file);
        }
    }

    size_t start;
    std::ofstream alpha_ofs, beta_ofs;
    if (args.gotArgument("alpha") && args.gotArgument("beta") && args.gotArgument("wprefix")) {
        std::cout << "Continue from check point\n";
        auto alpha_file = args.retrieve<std::string>("alpha"),
              beta_file = args.retrieve<std::string>( "beta");
        start = CL::utility::NLines(alpha_file);
        alpha_ofs.open(alpha_file, std::ofstream::app);
         beta_ofs.open( beta_file, std::ofstream::app);
        auto wprefix = args.retrieve<std::string>("wprefix");
        w.read(wprefix);
    }
    else {
        std::cout << "Initial iteration\n";
        start = 1;
        alpha_ofs.open("alpha.txt");
         beta_ofs.open( "beta.txt");
        CL::utility::show_time(std::cout);
        v1.read(fs);
        double alpha = kernel.initialize(v1, w);
        alpha_ofs << alpha << '\n';
    }
    std::cout << std::endl;

    size_t max_iteration = 100;
    if (args.gotArgument("max_iteration")) max_iteration = args.retrieve<size_t>("max_iteration");
    for (size_t i = start; i < start + max_iteration; i += 2) {
        double alpha, beta;
        std::cout << "Iteration " << i + 1 << std::endl;
        CL::utility::show_time(std::cout);
        std::tie(alpha, beta) = kernel.iterate(v1, v2, w);
        alpha_ofs << alpha << '\n';
         beta_ofs << beta  << '\n';
        for (auto & seg : fs) for (auto & state : seg) state.seekg(0);
        for (size_t j = 0; j < i; j++) {
            v1.read(fs);
            double coeff = v1.dot(w);
            w.sub_(coeff, v1);
        }
        v2.write(fs);
        std::cout << "Iteration " << i + 2 << std::endl;
        CL::utility::show_time(std::cout);
        std::tie(alpha, beta) = kernel.iterate(v2, v1, w);
        alpha_ofs << alpha << '\n';
         beta_ofs << beta  << '\n';
        // Reorthogonalize the Krylov vectors
        for (auto & seg : fs) for (auto & state : seg) state.seekg(0);
        for (size_t j = 0; j < i + 1; j++) {
            v2.read(fs);
            double coeff = v2.dot(w);
            w.sub_(coeff, v2);
        }
        v1.write(fs);
    }
    std::cout << '\n';

    alpha_ofs.close();
     beta_ofs.close();
    for (auto & seg : fs) for (auto & state : seg) state.close();

    std::cout << '\n';
    CL::utility::show_time(std::cout);
    std::cout << "Mission success\n";
}