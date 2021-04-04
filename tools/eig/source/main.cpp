#include <iomanip>

#include <CppLibrary/argparse.hpp>
#include <CppLibrary/utility.hpp>
#include <CppLibrary/linalg.hpp>

namespace { extern "C" {
    int32_t my_dsteqr_(double * diag, double * subdiag, double * eigvec, const int32_t & N);
} }

argparse::ArgumentParser parse_args(const size_t & argc, const char ** & argv) {
    CL::utility::echo_command(argc, argv, std::cout);
    std::cout << '\n';
    argparse::ArgumentParser parser("Diagonalize Lanczos tridiagonal matrix");

    // required arguments
    parser.add_argument("-a","--alpha", 1, false, "alpha (diag)");
    parser.add_argument("-b","--beta",  1, false, "beta (subdiag)");

    // optional arguments
    parser.add_argument("-v","--vector", 0, true, "require eigenvectors");

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

    auto energy = alpha / 4.556335830019422e-6;
    CL::utility::matrix<double> eigvec(N);
    size_t count = 0;
    for (auto & row : eigvec) for (auto & el : row) {
        el = eigvec_ptr[count];
        count++;
    }
    std::vector<double> amplitude(N), strength(N);
    for (size_t i = 0; i < N; i++) {
        amplitude[i] = eigvec[i][0];
        strength [i] = amplitude[i] * amplitude[i];
    }
    auto normalized_strength = strength / *std::max_element(strength.begin(), strength.end());

    std::cout << "Energy / cm^-1    amplitude     strength    normalized strength\n";
    for (size_t i = 0; i < N; i++)
    std::cout << std::fixed << std::setw(14) << std::setprecision(2) << energy   [i] << "    "
              << std::fixed << std::setw( 9) << std::setprecision(6) << amplitude[i] << "    "
              << std::fixed << std::setw( 9) << std::setprecision(6) << strength [i] << "    "
              << std::fixed << std::setw( 9) << std::setprecision(6) << normalized_strength[i] << '\n';
}