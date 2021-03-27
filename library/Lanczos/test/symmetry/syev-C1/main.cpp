#include <iostream>

#include <Lanczos/mv.hpp>

int main() {
    std::vector<std::string> vib_files = {"vibration1.in"};
    auto op = std::make_shared<vibron::Options>("symmetry.in", vib_files);
    op->pretty_print(std::cout);

    vibron::Wfn vi(op), vj(op), Hvj(op);

    std::vector<std::string> Hd_files({"Hd_1-1.in", "Hd_1-2.in", "Hd_2-2.in"});
    auto Hd = std::make_shared<Lanczos::Hd>(2, op->NModes, Hd_files);
    std::vector<std::string> freq_files = {"frequency.in"};
    Lanczos::MVKernel mvkernel(Hd, op, freq_files);
    mvkernel.pretty_print(std::cout);

    size_t NVibron = op->NVibron();
    std::cout << "There are " << NVibron << " vibronic basis functions\n";
    at::Tensor Hvibronic = at::eye(NVibron, c10::TensorOptions().dtype(torch::kFloat64));

    size_t i = 0;
    for (size_t iseg = 0; iseg < op->NSegs; iseg++)
    for (size_t istate = 0; istate < op->NStates; istate++)
    for (size_t ivib = 0; ivib < op->stops[iseg][istate] - op->starts[iseg][istate]; ivib++) {
        vi = 0.0;
        vi[{iseg, istate, ivib}] = 1.0;
        size_t j = 0;
        for (size_t jseg = 0; jseg < op->NSegs; jseg++)
        for (size_t jstate = 0; jstate < op->NStates; jstate++)
        for (size_t jvib = 0; jvib < op->stops[jseg][jstate] - op->starts[jseg][jstate]; jvib++) {
            vj = 0.0;
            vj[{jseg, jstate, jvib}] = 1.0;
            mvkernel(vj, Hvj);
            Hvibronic[i][j] = vi.dot(Hvj);
            j++;
        }
        i++;
    }
    std::cout << "Vibronic Hamiltonian:\n"
              << Hvibronic << '\n';

    at::Tensor energy, vibronic_state;
    std::tie(energy, vibronic_state) = Hvibronic.symeig(true);
    energy /= 4.556335830019422e-6;
    std::cout << energy << '\n';
    vibronic_state.transpose_(0, 1);
    at::Tensor population = vibronic_state[0] * vibronic_state[0];
    for (size_t i = 0; i < 3; i++)
    std::cout << population[i].item<double>() << '\n';
}