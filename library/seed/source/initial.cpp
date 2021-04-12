#include <seed/initial.hpp>

namespace seed {

Initial::Initial() {}
Initial::Initial(const std::vector<double> & _dipole,
const std::shared_ptr<vibron::Options> & _op, const std::vector<std::vector<size_t>> & _phonons)
: dipole_(_dipole), op_(_op), phonons_(_phonons) {
    if (_dipole.size() != _op->NStates) throw std::invalid_argument(
    "seed::Initial: inconsistent number of electronic states in dipole");
    if (_phonons.size() != _op->NIrreds) throw std::invalid_argument(
    "seed::Initial: inconsistent number of irreducibles in phonons");
    for (size_t i = 0; i < _op->NIrreds; i++)
    if (_phonons[i].size() != _op->NModes[i]) throw std::invalid_argument(
    "seed::Initial: inconsistent number of normal modes in irreducible " + std::to_string(i));
}
Initial::~Initial() {}

// Generate seed vector, return the norm of the seed vector and the corresponding unit vector
double Initial::generate_seed(vibron::Wfn & wfn) const {
    wfn = 0.0;
    size_t irred = op_->vib_irred(phonons_);
    int64_t abs_vib = op_->vib_sets[irred]->index_vibration(phonons_);
    if (abs_vib < 0) throw std::invalid_argument(
    "seed::Initial::generate_seed: the initial vibrational state is not included in the vibrational basis set");
    for (size_t istate = 0; istate < op_->NStates; istate++)
    if (op_->vib_irreds[istate] == irred) {
        size_t iseg, ivib;
        std::tie(iseg, ivib) = op_->vib_index(istate, abs_vib);
        wfn.select(iseg, istate, ivib) = dipole_[istate];
    }
    double norm = sqrt(wfn.dot(wfn));
    wfn *= 1.0 / norm;
    return norm;
}

} // namespace seed