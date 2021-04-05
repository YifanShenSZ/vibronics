#include <Lanczos/mv.hpp>

namespace {
    size_t NDiff(const std::vector<std::vector<size_t>> & a, const std::vector<std::vector<size_t>> & b) {
        size_t NDiff = 0;
        for (size_t i = 0; i < a.size(); i++)
        for (size_t j = 0; j < a[i].size(); j++)
        if (a[i][j] != b[i][j])
        NDiff += 1;
        return NDiff;
    }
}

namespace Lanczos {

SegStateVibValue::SegStateVibValue() {}
SegStateVibValue::SegStateVibValue(const size_t & _seg, const size_t & _state, const size_t & _vib, const double & _value)
: seg(_seg), state(_state), vib(_vib), value(_value) {}
SegStateVibValue::~SegStateVibValue() {}





// Construct `alloweds_` given constructed `Hd_`, `op_`, `integrator_`
void MVKernel::construct_nonzero() {
    // Shape vector
    alloweds_.resize(op_->NSegs);
    for (size_t iseg = 0; iseg < op_->NSegs; iseg++) {
        alloweds_[iseg].resize(op_->NStates);
        for (size_t istate = 0; istate < op_->NStates; istate++)
        alloweds_[iseg][istate].resize(op_->stops[iseg][istate] - op_->starts[iseg][istate]);
    }
    // Loop over rows of the vibronic Hamiltonian matrix
    #pragma omp parallel for
    for (size_t iseg = 0; iseg < op_->NSegs; iseg++)
    for (size_t istate = 0; istate < op_->NStates; istate++)
    for (size_t ivib = 0; ivib < op_->stops[iseg][istate] - op_->starts[iseg][istate]; ivib++) {
        // diagonal
        double value = this->Hdelement(iseg, istate, ivib, iseg, istate, ivib)
                     // same vibration can have Hd constant term
                     + (*Hd_)[{istate, istate}]->constant();
        // diagonal has harmonic oscillator eigenvalue term
        size_t iirred = op_->vib_irreds[istate],
               ivib_abs = ivib + op_->starts[iseg][istate];
        const auto & iphonons = op_->vib_sets[iirred][ivib_abs].phonons();
        for (size_t irred = 0; irred < op_->NIrreds; irred++)
        for (size_t mode = 0; mode < op_->NModes[irred]; mode++)
        value += integrators_[irred][mode].frequency() * (0.5 + iphonons[irred][mode]);
        // append
        if (abs(value) > 1e-15)
        alloweds_[iseg][istate][ivib].push_back(SegStateVibValue(iseg, istate, ivib, value));

        // same vibration but different electronic state
        for (size_t jstate = 0; jstate < op_->NStates; jstate++)
        if (jstate != istate && op_->vib_irreds[istate] == op_->vib_irreds[jstate]) {
            double value = this->Hdelement(iseg, istate, ivib, iseg, jstate, ivib)
                         // same vibration can have Hd constant term
                         + (*Hd_)[{istate, jstate}]->constant();
            if (abs(value) > 1e-15)
            alloweds_[iseg][istate][ivib].push_back(SegStateVibValue(iseg, jstate, ivib, value));
        }

        // remaining off-diagonals
        // Loop over excitations
        for (size_t excitation = 1; excitation <= Hd_->max_excitation(); excitation++) {
            // Select basic excited modes by indices
            std::vector<size_t> excited_indices(excitation);
            std::iota(excited_indices.begin(), excited_indices.end(), 0);
            // Map indices to modes of irreducibles
            std::vector<std::pair<size_t, size_t>> excited_modes(excitation);
            for (size_t i = 0; i < excitation; i++)
            excited_modes[i] = op_->vib_irred_mode(excited_indices[i]);
            generate_all(iseg, istate, ivib, excited_modes);
            // Loop over possible excited modes as an excitation-nary counter, with ascending digits
            while (true) {
                excited_indices.back() += 1;
                // Carry to former digit
                for (size_t i = excitation - 1; i > 0; i--)
                if (excited_indices[i] > op_->intdim() + i - excitation) {
                    excited_indices[i - 1] += 1;
                    excited_indices[i] = 0;
                }
                // Guarantee ascendance
                for (size_t i = 0; i < excitation - 1; i++)
                if (excited_indices[i] >= excited_indices[i + 1])
                excited_indices[i + 1] = excited_indices[i] + 1;
                // Finish when counter overflows
                if (excited_indices.back() >= op_->intdim()) break;
                // Map indices to modes of irreducibles
                std::vector<std::pair<size_t, size_t>> excited_modes(excitation);
                for (size_t i = 0; i < excitation; i++)
                excited_modes[i] = op_->vib_irred_mode(excited_indices[i]);
                generate_all(iseg, istate, ivib, excited_modes);
            }
        }
    }
}

// Supports `construct_nonzero`
// Given excited modes, generate all possible vibrations
void MVKernel::generate_all(const size_t & iseg, const size_t & istate, const size_t & ivib,
const std::vector<std::pair<size_t, size_t>> & excited_modes) {
    const auto & iphonons = op_->vib_sets[op_->vib_irreds[istate]][ivib + op_->starts[iseg][istate]].phonons();
    std::vector<size_t> min_phonons(excited_modes.size()), max_phonons(excited_modes.size());
    for (size_t i = 0; i < excited_modes.size(); i++) {
        const auto & mode = excited_modes[i];
        int64_t iphonon = iphonons[mode.first][mode.second];
        int64_t lower = iphonon - Hd_->max_order(mode),
                upper = iphonon + Hd_->max_order(mode);
        min_phonons[i] = lower > 0 ? lower : 0;
        const size_t & vib_max = op_->max_phonons[mode.first][mode.second];
        max_phonons[i] = upper < vib_max ? upper : vib_max;
    }
    // basic case: the excited modes have lowest possible phonons
    auto jphonons = iphonons;
    for (size_t i = 0; i < excited_modes.size(); i++) {
        const auto & mode = excited_modes[i];
        jphonons[mode.first][mode.second] = min_phonons[i];
    }
    if (NDiff(iphonons, jphonons) == excited_modes.size()) {
        size_t jirred = op_->vib_irred(jphonons);
        int64_t jvib_abs = op_->vib_sets[jirred].index_vibration(jphonons);
        if (jvib_abs >= 0)
        for (size_t jstate = 0; jstate < op_->NStates; jstate++)
        if (op_->vib_irreds[jstate] == jirred) {
            auto jseg_jvib = op_->vib_index(jstate, jvib_abs);
            size_t jseg = jseg_jvib.first, jvib = jseg_jvib.second;
            double value = this->Hdelement(iseg, istate, ivib, jseg, jstate, jvib, excited_modes);
            // append
            if (abs(value) > 1e-15)
            alloweds_[iseg][istate][ivib].push_back(SegStateVibValue(jseg, jstate, jvib, value));
        }
    }
    // Loop over possible phonons as a excited_modes.size()-nary counter
    while (true) {
        const auto & mode = excited_modes[0];
        jphonons[mode.first][mode.second] += 1;
        // Carry to latter digit
        for (size_t i = 0; i < excited_modes.size() - 1; i++) {
            const auto & mode1 = excited_modes[i];
            const auto & mode2 = excited_modes[i + 1];
            if (jphonons[mode1.first][mode1.second] > max_phonons[i]) {
                jphonons[mode1.first][mode1.second] = min_phonons[i];
                jphonons[mode2.first][mode2.second] += 1;
            }
        }
        // Finish when counter overflows
        if (jphonons[excited_modes.back().first][excited_modes.back().second] > max_phonons.back()) break;
        // Compute Hd element
        if (NDiff(iphonons, jphonons) == excited_modes.size()) {
            size_t jirred = op_->vib_irred(jphonons);
            int64_t jvib_abs = op_->vib_sets[jirred].index_vibration(jphonons);
            if (jvib_abs >= 0)
            for (size_t jstate = 0; jstate < op_->NStates; jstate++)
            if (op_->vib_irreds[jstate] == jirred) {
                auto jseg_jvib = op_->vib_index(jstate, jvib_abs);
                size_t jseg = jseg_jvib.first, jvib = jseg_jvib.second;
                double value = this->Hdelement(iseg, istate, ivib, jseg, jstate, jvib, excited_modes);
                // append
                if (abs(value) > 1e-15)
                alloweds_[iseg][istate][ivib].push_back(SegStateVibValue(jseg, jstate, jvib, value));
            }
        }
    }
}

// Supports `construct_nonzero`
// Return <basis[iseg, istate, ivib]|Hd|basis[jseg, jstate, jvib]>
double MVKernel::Hdelement(
const size_t & iseg, const size_t & istate, const size_t & ivib,
const size_t & jseg, const size_t & jstate, const size_t & jvib) const {
    double Hdelement = 0.0;
    const auto & iphonons = op_->vib_sets[op_->vib_irreds[istate]][ivib + op_->starts[iseg][istate]].phonons();
    const auto & jphonons = op_->vib_sets[op_->vib_irreds[jstate]][jvib + op_->starts[jseg][jstate]].phonons();
    for (const auto & coeff_sap : *(*Hd_)[{istate, jstate}]) {
        double Hdterm = 1.0;
        for (const auto & monomial : coeff_sap.second)
        Hdterm *= integrators_[monomial.irred][monomial.mode]({
            iphonons[monomial.irred][monomial.mode],
            monomial.order,
            jphonons[monomial.irred][monomial.mode]
        });
        // Assume other normal modes integrate to 1, i.e. i and j share same phonons in the remaining modes
        Hdelement += coeff_sap.first * Hdterm;
    }
    return Hdelement;
}
double MVKernel::Hdelement(
const size_t & iseg, const size_t & istate, const size_t & ivib,
const size_t & jseg, const size_t & jstate, const size_t & jvib,
const std::vector<std::pair<size_t, size_t>> & excited_modes) const {
    double Hdelement = 0.0;
    const auto & iphonons = op_->vib_sets[op_->vib_irreds[istate]][ivib + op_->starts[iseg][istate]].phonons();
    const auto & jphonons = op_->vib_sets[op_->vib_irreds[jstate]][jvib + op_->starts[jseg][jstate]].phonons();
    // i and j have different phonons is the excited modes, so the Hd term must include these modes
    for (size_t ex = excited_modes.size(); ex <= (*Hd_)[{istate, jstate}]->max_excitation(); ex++)
    for (const auto & coeff_sap : (*Hd_)[{istate, jstate}]->excitation(ex))
    if (coeff_sap->second >= excited_modes) {
        double Hdterm = 1.0;
        for (const auto & monomial : coeff_sap->second)
        Hdterm *= integrators_[monomial.irred][monomial.mode]({
            iphonons[monomial.irred][monomial.mode],
            monomial.order,
            jphonons[monomial.irred][monomial.mode]
        });
        // i and j share same phonons in the remaining modes
        Hdelement += coeff_sap->first * Hdterm;
    }
    return Hdelement;
}

MVKernel::MVKernel() {}
MVKernel::MVKernel(const std::shared_ptr<Hd> & _Hd, const std::shared_ptr<vibron::Options> & _op,
const std::vector<std::string> & frequency_files) : Hd_(_Hd), op_(_op) {
    if (_Hd->NStates() != _op->NStates) throw std::invalid_argument(
    "Lanczos::MVKernel::MVKernel: inconsistent number of electronic states between Hd and wfn");
    // Read frequencies and construct integrators
    if (frequency_files.size() != _op->NIrreds) throw std::invalid_argument(
    "Lanczos::MVKernel::MVKernel: one frequency input file per irreducible");
    integrators_.resize(_op->NIrreds);
    for (size_t i = 0; i < _op->NIrreds; i++) {
        integrators_[i].resize(_op->NModes[i]);
        std::ifstream ifs; ifs.open(frequency_files[i]);
        if (! ifs.good()) throw CL::utility::file_error(frequency_files[i]);
        for (size_t j = 0; j < _op->NModes[i]; j++) {
            double frequency;
            ifs >> frequency;
            if (! ifs.good()) throw CL::utility::file_error(frequency_files[i]);
            integrators_[i][j] = vibron::Integrator(frequency, _op->max_phonons[i][j], _Hd->max_order(i, j));
        }
        ifs.close();
    }
    // Infer the rest
    this->construct_nonzero();
}
MVKernel::~MVKernel() {}

void MVKernel::pretty_print(std::ostream & ostream) const {
    for (size_t iseg = 0; iseg < op_->NSegs; iseg++)
    for (size_t istate = 0; istate < op_->NStates; istate++)
    for (size_t ivib = 0; ivib < op_->stops[iseg][istate] - op_->starts[iseg][istate]; ivib++)
    ostream << "Segment " << iseg << ", state " << istate << ", vibration " << ivib << ": "
            << alloweds_[iseg][istate][ivib].size() << " non-trivial columns\n";
}

void MVKernel::operator()(const vibron::Wfn & wfn, vibron::Wfn & Hwfn) const {
    if (op_ != wfn.options() || op_ != Hwfn.options()) throw std::invalid_argument(
    "Lanczos::MVKernel::operator(): wave functions must share a same definition");
    #pragma omp parallel for
    for (size_t iseg = 0; iseg < op_->NSegs; iseg++)
    for (size_t istate = 0; istate < op_->NStates; istate++)
    for (size_t ivib = 0; ivib < op_->stops[iseg][istate] - op_->starts[iseg][istate]; ivib++) {
        Hwfn.select(iseg, istate, ivib) = 0.0;
        for (const auto & element : alloweds_[iseg][istate][ivib])
        Hwfn.select(iseg, istate, ivib) += element.value * wfn.select(element.seg, element.state, element.vib);
    }
}

} // namespace Lanczos