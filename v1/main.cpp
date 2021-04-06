#include <CppLibrary/argparse.hpp>
#include <CppLibrary/utility.hpp>

#include <vibron/options.hpp>

#include <Lanczos/mv.hpp>
#include <Lanczos/Lanczos.hpp>

argparse::ArgumentParser parse_args(const size_t & argc, const char ** & argv) {
    CL::utility::echo_command(argc, argv, std::cout);
    std::cout << '\n';
    argparse::ArgumentParser parser("Vibronics version 1");

    // required arguments
    parser.add_argument("-w","--wfn",         1, false, "vibronic wave function definition file");
    parser.add_argument("-v","--vibration", '+', false, "vibrational basis definition files");
    parser.add_argument("-H","--Hd",        '+', false, "anharmonic diabatic Hamiltonian definition files");
    parser.add_argument("-f","--frequency", '+', false, "frequency of each normal mode");
    parser.add_argument("-p","--prefix",      1, false, "prefix of the seed or the check point vector to continue with");

    // optional arguments
    parser.add_argument("-m","--max_iteration", 1, true, "max number of iterations to perform (default = 100)");

    parser.parse_args(argc, argv);
    return parser;
}

int main(size_t argc, const char ** argv) {
    std::cout << "Vibronics version 1\n"
              << "This version reorthogonalizes the Krylov vectors in RAM\n"
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
    auto freq_files = args.retrieve<std::vector<std::string>>("frequency");
    auto mvkernel = std::make_shared<Lanczos::MVKernel>(Hd, op, freq_files);
    Lanczos::Kernel kernel(mvkernel);
    kernel.pretty_print(std::cout);
    std::cout << std::endl;

    size_t max_iteration = 100;
    if (args.gotArgument("max_iteration")) max_iteration = args.retrieve<size_t>("max_iteration");
    vibron::Wfn w(op);
    std::vector<vibron::Wfn *> wfns(max_iteration);
    for (auto & wfn : wfns) wfn = new vibron::Wfn(op);

    auto prefix = args.retrieve<std::string>("prefix");
    std::vector<std::vector<std::ifstream>> ifs(op->NSegs);
    for (size_t i = 0; i < op->NSegs; i++) {
        ifs[i].resize(op->NStates);
        for (size_t j = 0; j < op->NStates; j++) {
            std::string file = prefix + "-" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + ".wfn";
            ifs[i][j].open(file, std::ifstream::binary);
            if (! ifs[i][j].good()) throw CL::utility::file_error(file);
        }
    }
    wfns[0]->read(ifs);
    for (auto & seg : ifs) for (auto & state : seg) state.close();

    size_t start;
    std::ofstream alpha_ofs, beta_ofs;
    alpha_ofs.open("alpha.txt");
     beta_ofs.open( "beta.txt");
    std::cout << "Iteration 0\n";
    CL::utility::show_time(std::cout);
    double alpha = kernel.initialize(*wfns[0], w);
    alpha_ofs << alpha << '\n';
    std::cout << std::endl;

    for (size_t i = 1; i < max_iteration; i++) {
        double alpha, beta;
        std::cout << "Iteration " << i << std::endl;
        CL::utility::show_time(std::cout);
        std::tie(alpha, beta) = kernel.iterate(*wfns[i - 1], *wfns[i], w);
        alpha_ofs << alpha << '\n';
         beta_ofs << beta  << '\n';
        // Reorthogonalize the Krylov vectors
        for (size_t j = 0; j < i; j++) {
            double coeff = wfns[j]->dot(w);
            w.sub_(coeff, *wfns[j]);
        }
    }

    alpha_ofs.close();
     beta_ofs.close();

    std::cout << '\n';
    CL::utility::show_time(std::cout);
    std::cout << "Mission success\n";
}