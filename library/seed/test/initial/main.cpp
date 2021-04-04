#include <seed/initial.hpp>

int main() {
    std::vector<std::string> vib_files = {"vibration1.in"};
    auto op = std::make_shared<vibron::Options>("symmetry.in", vib_files);
    op->pretty_print(std::cout);

    vibron::Wfn wfn(op);

    std::vector<std::vector<size_t>> phonons(1);
    phonons[0] = {0, 0, 0};
    seed::Initial init_bias({1.0, 1.0}, op, phonons);

    double norm = init_bias.generate_seed(wfn);
    std::cout << "Norm of the seed vector = " << norm << '\n';
}