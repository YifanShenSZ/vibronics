#include <CppLibrary/argparse.hpp>
#include <CppLibrary/utility.hpp>

#include <vibron/options.hpp>

#include <Lanczos/mv.hpp>
#include <Lanczos/Lanczos.hpp>

argparse::ArgumentParser parse_args(const size_t & argc, const char ** & argv) {
    CL::utility::echo_command(argc, argv, std::cout);
    std::cout << '\n';
    argparse::ArgumentParser parser("Vibronics version 0");

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
    std::cout << "Vibronics version 0\n"
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
    v1.read(prefix);

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
        CL::utility::show_time(std::cout);
        std::cout << "Iteration 0\n";
        start = 0;
        alpha_ofs.open("alpha.txt");
         beta_ofs.open( "beta.txt");
        double alpha = kernel.initialize(v1, w);
        alpha_ofs << alpha << '\n';
    }
    std::cout << std::endl;

    size_t max_iteration = 100;
    if (args.gotArgument("max_iteration")) max_iteration = args.retrieve<size_t>("max_iteration");
    for (size_t i = start; i < start + max_iteration; i += 2) {
        double alpha, beta;
        CL::utility::show_time(std::cout);
        std::cout << "Iteration " << i + 1 << '\n' << std::endl;
        std::tie(alpha, beta) = kernel.iterate(v1, v2, w);
        alpha_ofs << alpha << '\n';
         beta_ofs << beta  << '\n';
        CL::utility::show_time(std::cout);
        std::cout << "Iteration " << i + 2 << '\n' << std::endl;
        std::tie(alpha, beta) = kernel.iterate(v2, v1, w);
        alpha_ofs << alpha << '\n';
         beta_ofs << beta  << '\n';
    }

    // Write checkpoint
    std::cout << "Restart Lanczos vectors are written to v*.wfn and w*.wfn\n";
    v1.write("v");
    w .write("w");

    alpha_ofs.close();
     beta_ofs.close();

    std::cout << '\n';
    CL::utility::show_time(std::cout);
    std::cout << "Mission success\n";
}