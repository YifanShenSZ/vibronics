#include <iostream>

#include <Lanczos/Hamiltonian.hpp>

int main() {
    std::vector<std::string> Hd_files({"Hd_1-1.in", "Hd_1-2.in", "Hd_2-2.in"});
    Lanczos::Hd Hd(2, {3}, Hd_files);

    std::cerr << Hd.max_order(0, 0) << ' '
              << Hd.max_order(0, 1) << ' '
              << Hd.max_order(0, 2) << '\n';
    for (size_t i = 0; i < Hd.NStates(); i++)
    for (size_t j = i; j < Hd.NStates(); j++)
    std::cerr << "Hd " << i << ' ' << j << '\n'
              << "max excitation " << Hd[{i, j}]->max_excitation() << '\n'
              << Hd[{i, j}]->excitation(2)[0]->second[0].irred << ' '
              << Hd[{i, j}]->excitation(2)[0]->second[0].mode  << ' '
              << Hd[{i, j}]->excitation(2)[0]->second[0].order << '\n'
              << Hd[{i, j}]->excitation(2)[0]->second[1].irred << ' '
              << Hd[{i, j}]->excitation(2)[0]->second[1].mode  << ' '
              << Hd[{i, j}]->excitation(2)[0]->second[1].order << '\n';
}