#include <iostream>

#include <CppLibrary/utility.hpp>

#include <Lanczos/Lanczos.hpp>

int main() {
    std::vector<std::string> vib_files = {"vibration1.in"};
    auto op = std::make_shared<vibron::Options>("symmetry.in", vib_files);
    op->pretty_print(std::cout);

    vibron::Wfn v1(op), v2(op), w(op);

    std::vector<std::string> Hd_files({"Hd_1-1.in", "Hd_1-2.in", "Hd_2-2.in"});
    auto Hd = std::make_shared<Lanczos::Hd>(2, op->NModes, Hd_files);
    std::vector<std::string> freq_files = {"frequency.in"};
    auto mvkernel = std::make_shared<Lanczos::MVKernel>(Hd, op, freq_files);

    v1 = 0.0;
    v1[{0, 0}][0] = 1.0 / sqrt(2.0);
    v1[{0, 1}][0] = 1.0 / sqrt(2.0);

    std::ofstream alpha_fs, beta_fs;
    alpha_fs.open("alpha.txt");
     beta_fs.open( "beta.txt");

    Lanczos::Kernel kernel(mvkernel);
    std::cout << "Iteration " << 1 << '\n';
    CL::utility::show_time(std::cout);
    double alpha = kernel.initialize(v1, w);
    alpha_fs << alpha << '\n';

    for (size_t i = 1; i < 100; i += 2) {
        double alpha, beta;
        std::cout << "Iteration " << i + 1 << '\n';
        CL::utility::show_time(std::cout);
        std::tie(alpha, beta) = kernel.iterate(v1, v2, w);
        alpha_fs << alpha << '\n';
         beta_fs << beta  << '\n';
         std::cout << "Iteration " << i + 2 << '\n';
        CL::utility::show_time(std::cout);
        std::tie(alpha, beta) = kernel.iterate(v2, v1, w);
        alpha_fs << alpha << '\n';
         beta_fs << beta  << '\n';
    }

    alpha_fs.close();
     beta_fs.close();
}