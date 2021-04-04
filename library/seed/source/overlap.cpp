#include <CppLibrary/math.hpp>

#include <seed/overlap.hpp>

namespace seed {

IFOverlap::IFOverlap() {}
// Let the initial-state normal coordiante q and the final-state normal coordinate Q be related by:
//    q = T . Q + d
// Reference: M. S. Schuurman and D. R. Yarkony, J. Chem. Phys. 2008, 128, 044119 (https://doi.org/10.1063/1.2826380)
// The notations here follow the reference
IFOverlap::IFOverlap(const std::shared_ptr<vibron::Options> & _op,
const at::Tensor & init_freq, const at::Tensor & final_freq,
const at::Tensor & T, const at::Tensor & d)
: op_(_op) {
    size_t intdim = _op->intdim();
    if (init_freq.sizes().size() != 1) throw std::invalid_argument(
    "seed::IFOverlap: initial-state frequency must be a vector");
    if (init_freq.size(0) != intdim) throw std::invalid_argument(
    "seed::IFOverlap: initial-state frequency must share a same dimension with the internal coordinate system");
    if (final_freq.sizes().size() != 1) throw std::invalid_argument(
    "seed::IFOverlap: final-state frequency must be a vector");
    if (final_freq.size(0) != intdim) throw std::invalid_argument(
    "seed::IFOverlap: final-state frequency must share a same dimension with the internal coordinate system");
    if (T.sizes().size() != 2) throw std::invalid_argument(
    "seed::IFOverlap: T must be a matrix");
    if (T.size(0) != T.size(1)) throw std::invalid_argument(
    "seed::IFOverlap: T must be a square matrix");
    if (T.size(0) != intdim) throw std::invalid_argument(
    "seed::IFOverlap: T must share a same dimension with the internal coordinate system");
    if (d.sizes().size() != 1) throw std::invalid_argument(
    "seed::IFOverlap: d must be a vector");
    if (d.size(0) != intdim) throw std::invalid_argument(
    "seed::IFOverlap: d must share a same dimension with the internal coordinate system");
    // Prepare the coefficient matrices and vector
    at::Tensor A = 0.5 * (final_freq.diag() + T.transpose(0, 1).mm(init_freq.diag().mm(T)));
    at::Tensor Acholesky = A.cholesky();
    at::Tensor Ainv = at::cholesky_inverse(Acholesky);
    at::Tensor SqrtDetA = Acholesky.diag().prod();
    at::Tensor SqrtBeta = at::sqrt(final_freq).diag();
    at::Tensor Abar = Ainv.mm(SqrtBeta);
    at::Tensor Atilde = SqrtBeta.mm(Abar);
    at::Tensor b = 0.5 * T.transpose(0, 1).mv(init_freq * d);
    at::Tensor bT_Abar = at::matmul(b, Abar);
    at::Tensor Atilde_1 = Atilde - at::eye(intdim, c10::TensorOptions().dtype(torch::kFloat64));
    // Generate the necessary final-state-biased vibrational basis functions
    std::vector<size_t> max_phonons(op_->intdim());
    size_t count = 0;
    for (const auto & irred : op_->max_phonons)
    for (const size_t & phonon : irred) {
        max_phonons[count] = phonon;
        count++;
    }
    vib_set_ = std::make_shared<vibron::VibrationSet>(max_phonons);
    // Recurse over all the generated vibrational basis functions
    integrals_.resize(vib_set_->size());
    integrals_[0] = 1.0;
    for (size_t i = 1; i < vib_set_->size(); i++) {
        const vibron::Vibration & vib = (*vib_set_)[i];
        std::vector<size_t> phonons = vib.phonons()[0];
        // Use the first non-zero mode as k
        size_t k;
        for (k = 0; k < phonons.size(); k++) if (phonons[k] > 0) break;
        std::vector<size_t> n = phonons;
        n[k]--;
        int64_t index_n = vib_set_->index_vibration(vibron::Vibration({n}));
        if (index_n < 0) throw std::runtime_error(
        "seed::IFOverlap: recurse out of the generated vibration set");
        integrals_[i] = bT_Abar[k].item<double>() * integrals_[index_n];
        for (size_t l = 0; l < intdim; l++) {
            std::vector<size_t> nl = n;
            if (nl[l] == 0) continue; // * 0 has no contribution
            nl[l]--;
            int64_t index_nl = vib_set_->index_vibration(vibron::Vibration({nl}));
            if (index_nl < 0) throw std::runtime_error(
            "seed::IFOverlap: recurse out of the generated vibration set");
            integrals_[i] += Atilde_1[k][l].item<double>() * n[l] * integrals_[index_nl];
        }
        integrals_[i] *= 2.0;
    }
    // C(n) -> <initial vibrational state|final-state-biased vibrational basis function>
    double Gtilde = (at::sqrt(T.det().abs_() * at::sqrt((init_freq * final_freq).prod())) / SqrtDetA
                  * at::exp(-0.5 * d.dot(init_freq * d) + b.dot(Ainv.mv(b)))).item<double>();
    integrals_[0] *= Gtilde;
    for (size_t i = 1; i < vib_set_->size(); i++) {
        const std::vector<size_t> & phonons = (*vib_set_)[i].phonons()[0];
        double coeff = 1.0;
        for (const size_t & phonon : phonons) if (phonon > 0)
        coeff *= pow(2, phonon) * CL::math::dFactorial(phonon);
        integrals_[i] *= Gtilde / sqrt(coeff);
    }
}
IFOverlap::~IFOverlap() {}

const double & IFOverlap::operator()(const size_t & index) const {return integrals_[index];}
const double & IFOverlap::operator()(const vibron::Vibration &vibration) const {
    int64_t index = vib_set_->index_vibration(vibration);
    if (index < 0) throw std::invalid_argument(
    "seed::IFOverlap::operator(): vibration out of range");
    return integrals_[index];
}

}; // namespace seed