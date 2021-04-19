#include <seed/final.hpp>

namespace seed {

Final::Final() {}
Final::Final(const std::vector<double> & _dipole, const std::shared_ptr<vibron::Options> & _op,
const at::Tensor & init_freq, const at::Tensor & final_freq,
const at::Tensor & tran_matrix, const at::Tensor & shift_vector)
: dipole_(_dipole), op_(_op) {
    if (_dipole.size() != _op->NStates) throw std::invalid_argument(
    "seed::Final: inconsistent number of electronic states in dipole");
    ifoverlap_ = std::make_shared<IFOverlap>(_op, init_freq, final_freq, tran_matrix, shift_vector);
}
Final::~Final() {}

void Final::pretty_print(std::ostream & stream) const {
    ifoverlap_->pretty_print(stream);
}

// Generate seed vector, return the norm of the seed vector and the corresponding unit vector
double Final::generate_seed(vibron::Wfn & wfn) const {
    if (op_ != wfn.options()) throw std::invalid_argument(
    "seed::Final::generate_seed: inconsistent vibronic wave function definition");
    size_t intdim = op_->intdim();
    #pragma omp parallel for
    for (size_t seg = 0; seg < op_->NSegs; seg++)
    for (size_t state = 0; state < op_->NStates; state++) {
        const auto & vib_set = op_->vib_sets[op_->vib_irreds[state]];
        for (size_t vib = 0; vib < op_->stops[seg][state] - op_->starts[seg][state]; vib++) {
            std::vector<std::vector<size_t>> C1_phonons(1);
            C1_phonons[0].resize(intdim);
            size_t count = 0;
            for (const auto & irred : (*vib_set)[op_->vib_index_abs(seg, state, vib)].phonons())
            for (const size_t & phonon : irred) {
                C1_phonons[0][count] = phonon;
                count++;
            }
            wfn.select(seg, state, vib) = (*ifoverlap_)(C1_phonons);
        }
        wfn[{seg, state}].mul_(dipole_[state]);
    }
    double norm = sqrt(wfn.dot(wfn));
    wfn *= 1.0 / norm;
    return norm;
}

} // namespace seed