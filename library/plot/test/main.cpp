#include <iostream>

#include <CppLibrary/utility.hpp>

#include <plot/wfn.hpp>

int main() {
    std::vector<std::string> vib_files = {"vib.in"};
    auto op = std::make_shared<vibron::Options>("wfn.in", vib_files);

    auto freq = CL::utility::read_vector("freq.in");
    std::vector<std::vector<double>> freqs = {freq};
    plot::Wfn v(op, freqs);

    v.vibron::Wfn::operator=(0.0);
    v.select(0, 0, 0) = 1.0 / sqrt(2.0);
    v.select(0, 1, 1) = 1.0 / sqrt(2.0);

    std::ofstream ofs_lower, ofs_upper;
    ofs_lower.open("lower.txt");
    ofs_upper.open("upper.txt");

    double sqrt1 = sqrt(freq[0]),
           sqrt2 = sqrt(freq[1]);
    std::vector<double> Q = {0.0, 0.0};
    std::vector<std::vector<double>> Qs = {Q};
    for (int64_t i = -10; i <= 10; i++)
    for (int64_t j = -10; j <= 10; j++) {
        Qs[0][0] = i * 0.2 / sqrt1;
        Qs[0][1] = j * 0.2 / sqrt2;
        auto wfn = v(Qs);
        ofs_lower << wfn[0] << '\n';
        ofs_upper << wfn[1] << '\n';
    }

    ofs_lower.close();
    ofs_upper.close();
}