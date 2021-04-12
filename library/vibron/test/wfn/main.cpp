#include <iostream>

#include <vibron/wfn.hpp>

int main() {
    std::vector<std::string> vib_files = {"vibration1.in", "vibration2.in"};
    auto op = std::make_shared<vibron::Options>("symmetry.in", vib_files);
    op->pretty_print(std::cout);

    vibron::Wfn v(op);

    std::cout << v[0].size() - 2 << ' '
              << v[{0, 1}][2].item<double>() - v.select(0, 1, 2) << '\n';
}