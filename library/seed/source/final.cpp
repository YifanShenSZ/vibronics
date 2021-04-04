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

// Generate seed vector, return the norm of the seed vector and the corresponding unit vector
double Final::generate_seed(vibron::Wfn & wfn) const {
    size_t intdim = op_->intdim();
    #pragma omp parallel for
    for (size_t iseg = 0; iseg < op_->NSegs; iseg++)
    for (size_t istate = 0; istate < op_->NStates; istate++)
    for (size_t ivib = 0; ivib < op_->stops[iseg][istate] - op_->starts[iseg][istate]; ivib++) {
        size_t abs_ivib = op_->vib_index_abs(iseg, istate, ivib);
        std::vector<std::vector<size_t>> C1_phonons(1);
        C1_phonons[0].resize(intdim);
        size_t count = 0;
        for (const auto & irred : op_->vib_sets[op_->vib_irreds[istate]][abs_ivib].phonons())
        for (const size_t & phonon : irred) {
            C1_phonons[0][count] = phonon;
            count++;
        }
        wfn.select(iseg, istate, ivib) = (*ifoverlap_)(C1_phonons);
    }
    double norm = sqrt(wfn.dot(wfn));
    wfn *= 1.0 / norm;
    return norm;
}

} // namespace seed