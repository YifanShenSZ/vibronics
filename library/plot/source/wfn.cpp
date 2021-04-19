#include <plot/wfn.hpp>

namespace plot {

Wfn::Wfn() {}
Wfn::Wfn(const std::shared_ptr<vibron::Options> & _op, const std::vector<std::vector<double>> & _frequencies)
: vibron::Wfn(_op), frequencies_(_frequencies) {}
Wfn::~Wfn() {}

// Given normal coordinate Q, return vibronic wave function value
std::vector<double> Wfn::operator()(const std::vector<std::vector<double>> & Qs) const {
    if (Qs.size() != op_->NIrreds) throw std::invalid_argument(
    "plot::Wfn::operator(): inconsistent number of irreducible representations");
    for (size_t i = 0; i < op_->NIrreds; i++)
    if (Qs[i].size() != op_->NModes[i]) throw std::invalid_argument(
    "plot::Wfn::operator(): inconsistent number of normal modes");
    std::vector<double> wfn(op_->NStates, 0.0);
    std::vector<std::vector<double>> wfns(op_->NSegs);
    #pragma omp parallel for
    for (size_t seg = 0; seg < op_->NSegs; seg++) {
        auto & wfn = wfns[seg];
        wfn.resize(op_->NStates);
        std::fill(wfn.begin(), wfn.end(), 0.0);
        for (size_t state = 0; state < op_->NStates; state++) {
            const auto & vib_set = op_->vib_sets[op_->vib_irreds[state]];
            for (size_t vib = 0; vib < op_->stops[seg][state] - op_->starts[seg][state]; vib++)
            wfn[state] += data_ptrs_[seg][state][vib]
                        * (*vib_set)[vib + op_->starts[seg][state]].value(frequencies_, Qs);
        } 
    }
    for (size_t state = 0; state < op_->NStates; state++)
    for (size_t seg = 0; seg < op_->NSegs; seg++)
    wfn[state] += wfns[seg][state];
    return wfn;
}

} // namespace plot