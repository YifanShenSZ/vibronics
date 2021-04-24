#include <vibron/wfn.hpp>

// M. S. Schuurman, D. E. Weinberg, D. R. Yarkony, J. Chem. Phys. 2007, 127, 104309 (https://doi.org/10.1063/1.2764052)
void perturb(const std::vector<std::vector<std::pair<double, vibron::Wfn *>>> & energies_wfns,
const double & x, const double & y, const double & z) {
    int64_t Neig = 0;
    for (const auto & irred : energies_wfns) Neig += irred.size();
    // Prepare generalized Ham reduction factor matrix O
    at::Tensor O = at::zeros({Neig, Neig}, c10::TensorOptions().dtype(torch::kFloat64));
    size_t i = 0;
    for (const auto & i_irred : energies_wfns)
    for (const auto & i_energy_wfn : i_irred) {
        const auto & i_wfn = i_energy_wfn.second;
        const auto & i_op = i_wfn->options();
        size_t j = 0;
        for (const auto & j_irred : energies_wfns)
        for (const auto & j_energy_wfn : j_irred) {
            const auto & j_wfn = j_energy_wfn.second;
            const auto & j_op = j_wfn->options();
            // Ham reduction between i-th lower state and j-th upper state
            if (i_op->vib_irreds[0] == j_op->vib_irreds[1])
            for (size_t seg = 0; seg < i_op->NSegs; seg++)
            O[i][j] += (*i_wfn)[{seg, 0}].dot((*j_wfn)[{seg, 1}]);
            // Ham reduction between i-th upper state and j-th lower state
            if (i_op->vib_irreds[1] == j_op->vib_irreds[0])
            for (size_t seg = 0; seg < i_op->NSegs; seg++)
            O[i][j] -= (*i_wfn)[{seg, 1}].dot((*j_wfn)[{seg, 0}]);
            j++;
        }
        i++;
    }
    // Prepare diagonal nonrelativistic eigenvalue matrix
    at::Tensor Enr = O.new_zeros({Neig, Neig});
    size_t count = 0;
    for (const auto & irred : energies_wfns)
    for (const auto & energy_wfn : irred) {
        Enr[count][count].fill_(energy_wfn.first);
        count++;
    }
    // Construct real-valued perturbative Hamiltonian
    at::Tensor H = O.new_zeros({4 * Neig, 4 * Neig});
    H.slice(0, 0, Neig).slice(1, 0, Neig).copy_(Enr);
    H.slice(0, 0, Neig).slice(1,     Neig, 2 * Neig).copy_(-x * O);
    H.slice(0, 0, Neig).slice(1, 2 * Neig, 3 * Neig).copy_(-y * O);
    H.slice(0, 0, Neig).slice(1, 3 * Neig, 4 * Neig).copy_(-z * O);
    H.slice(0, Neig, 2 * Neig).slice(1, Neig, 2 * Neig).copy_(Enr);
    H.slice(0, Neig, 2 * Neig).slice(1, 2 * Neig, 3 * Neig).copy_(-z * O);
    H.slice(0, Neig, 2 * Neig).slice(1, 3 * Neig, 4 * Neig).copy_( y * O);
    H.slice(0, 2 * Neig, 3 * Neig).slice(1, 2 * Neig, 3 * Neig).copy_(Enr);
    H.slice(0, 2 * Neig, 3 * Neig).slice(1, 3 * Neig, 4 * Neig).copy_(-x * O);
    H.slice(0, 3 * Neig, 4 * Neig).slice(1, 3 * Neig, 4 * Neig).copy_(Enr);
    // Diagonalize the real-valued perturbative Hamiltonian
    at::Tensor energy, state;
    std::tie(energy, state) = H.symeig(true);

    std::cout << "The relativistic energy levels can be found in soc-level.txt\n"
              << "The complex -> real arithmetic introduces 2-fold degeneracy\n"
              << "Kramers degeneracy theorem introduces 2-fold degeneracy for doublets\n"
              << "So for our 2 nonrelativistic doublets, there will be 4-fold degeneracy\n";
    std::ofstream ofs; ofs.open("soc-level.txt");
    double E0 = energy[0].item<double>() / 4.556335830019422e-6;
    ofs << "Ground state energy = " << std::fixed << std::setprecision(2) << E0 << " cm^-1\n"
        << "ΔE / cm^-1    E / Hartree\n";
    for (size_t i = 0; i < 4 * Neig; i += 4) {
        double E = energy[i].item<double>();
        ofs << std::setw(10) << std::fixed << std::setprecision(2) << E / 4.556335830019422e-6 - E0 
            << std::setw(25) << std::scientific << std::setprecision(15) << E
            << '\n';
    }
}