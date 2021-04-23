#include <fstream>
#include <cstring>
#include <memory>

#include <CppLibrary/argparse.hpp>
#include <CppLibrary/utility.hpp>
#include <CppLibrary/linalg.hpp>

#include <vibron/options.hpp>

namespace { extern "C" {
    int32_t my_dsteqr_(double * diag, double * subdiag, double * eigvec, const int32_t & N);
} }

void output_spectrum(const std::vector<double> & energy, const CL::utility::matrix<double> & eigvec,
const double & last_beta, const double & threshold);

void output_wfn(const std::vector<size_t> & vec_requests, const CL::utility::matrix<double> & eigvec,
const std::shared_ptr<vibron::Options> & op, const std::string & prefix);

argparse::ArgumentParser parse_args(const size_t & argc, const char ** & argv) {
    CL::utility::echo_command(argc, argv, std::cout);
    std::cout << '\n';
    argparse::ArgumentParser parser("Diagonalize Lanczos tridiagonal matrix");

    // required arguments
    parser.add_argument("-a","--alpha", 1, false, "alpha (diag)");
    parser.add_argument("-b","--beta",  1, false, "beta (subdiag)");

    // optional arguments
    parser.add_argument("-t","--threshold", 1, true, "intensity threshold (default = 1e-6)");
    parser.add_argument("-V","--vector",  '+', true, "eigenvectors to compute");

    // wave function definition, required if eigenvector requested
    parser.add_argument("-w","--wfn",         1, true, "vibronic wave function definition file");
    parser.add_argument("-v","--vibration", '+', true, "vibrational basis definition files");
    parser.add_argument("-p","--prefix",      1, true, "prefix of the stored Krylov vectors");

    parser.parse_args(argc, argv);
    return parser;
}

int main(size_t argc, const char ** argv) {
    std::cout << "Vibronics version 0\n"
              << "Yifan Shen 2021\n\n";
    argparse::ArgumentParser args = parse_args(argc, argv);
    CL::utility::show_time(std::cout);
    std::cout << '\n';

    auto alpha_file = args.retrieve<std::string>("alpha"),
          beta_file = args.retrieve<std::string>( "beta");
    auto alpha = CL::utility::read_vector(alpha_file),
          beta = CL::utility::read_vector( beta_file);
    if (alpha.size() != beta.size() + 1) throw std::invalid_argument(
    "inconsistent number of elements between alpha and beta");

    int32_t N = alpha.size();
    double * eigvec_ptr = new double[N * N];
    int32_t info = my_dsteqr_(alpha.data(), beta.data(), eigvec_ptr, N);
    if (info != 0) std::cerr << "Warning: diagonalization of the tridiagonal matrix failed, infomation = " << info << '\n';

    CL::utility::matrix<double> eigvec(N);
    for (auto & row : eigvec) {
        std::memcpy(row.data(), eigvec_ptr, N);
        eigvec_ptr += N;
        if (row[0] < 0.0) row *= -1.0;
    }

    double threshold = 1e-6;
    if (args.gotArgument("threshold")) threshold = args.retrieve<double>("threshold");
    output_spectrum(alpha, eigvec, beta.back(), threshold);

    if (args.gotArgument("vector")) {
        std::cout << '\n';

        auto vec_requests = args.retrieve<std::vector<size_t>>("vector");

        auto wfn_file  = args.retrieve<std::string>("wfn");
        auto vib_files = args.retrieve<std::vector<std::string>>("vibration");
        auto op = std::make_shared<vibron::Options>(wfn_file, vib_files);

        auto prefix = args.retrieve<std::string>("prefix");
        output_wfn(vec_requests, eigvec, op, prefix);
    }

    std::cout << '\n';
    CL::utility::show_time(std::cout);
    std::cout << "Mission success\n";
}