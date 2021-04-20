#include "HdKernel.hpp"

HdKernel::HdKernel() {}
HdKernel::HdKernel(const std::vector<std::vector<double>> & _freqs, const std::shared_ptr<Lanczos::Hd> & _Hanharmonic)
: freqs_(_freqs), Hanharmonic_(_Hanharmonic) {}
HdKernel::~HdKernel() {}

at::Tensor HdKernel::operator()(const std::vector<std::vector<double>> & Qs) const {
    // harmonic term
    double harmonicity = 0.0;
    for (size_t irred = 0; irred < Qs.size(); irred++)
    for (size_t mode = 0; mode < Qs[irred].size(); mode++) {
        double temp = freqs_[irred][mode] * Qs[irred][mode];
        harmonicity += temp * temp;
    }
    harmonicity /= 2.0;
    // anharmonic term
    int64_t NStates = Hanharmonic_->NStates();
    at::Tensor Hd = at::empty({NStates, NStates}, c10::TensorOptions().dtype(torch::kFloat64));
    auto anharmonicity = (*Hanharmonic_)(Qs);
    // combine them together
    for (size_t i = 0; i < NStates; i++) {
        Hd[i][i] = harmonicity + anharmonicity[i][i];
        for (size_t j = i + 1; j < NStates; j++) Hd[i][j].fill_(anharmonicity[i][j]);
    }
    // output
    return Hd;
}