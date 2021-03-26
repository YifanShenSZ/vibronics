#include <vibron/wfn.hpp>

namespace vibron {

Wfn::Wfn() {}
Wfn::Wfn(const std::shared_ptr<Options> & _op) : op_(_op) {
    data_.resize(_op->NSegs);
    data_ptrs_.resize(_op->NSegs);
    for (size_t i = 0; i < _op->NSegs; i++) {
        data_[i].resize(_op->NStates);
        data_ptrs_[i].resize(_op->NStates);
    }
    for (size_t i = 0; i < _op->NSegs; i++)
    for (size_t j = 0; j < _op->NStates; j++) {
       data_[i][j] = at::empty(_op->stops[i][j] - _op->starts[i][j],
                     c10::TensorOptions().dtype(torch::kFloat64));
       data_ptrs_[i][j] = data_[i][j].data_ptr<double>();
    }
}
Wfn::~Wfn() {}

const std::shared_ptr<Options> & Wfn::options() const {return op_;}

const std::vector<at::Tensor> & Wfn::operator[](const size_t & seg) const {return data_[seg];}
at::Tensor Wfn::operator[](const CL::utility::triple<size_t, size_t, size_t> & seg_state_vib) const {
    return data_[seg_state_vib.first][seg_state_vib.second][seg_state_vib.third];
}

double & Wfn::select(const size_t & seg, const size_t & state, const size_t & vib) {
    return data_ptrs_[seg][state][vib];
}
const double & Wfn::select(const size_t & seg, const size_t & state, const size_t & vib) const {
    return data_ptrs_[seg][state][vib];
}

double Wfn::dot(const Wfn & other) const {
    std::vector<double> products(op_->NSegs, 0.0);
    #pragma omp parallel for
    for (size_t i = 0; i < op_->NSegs; i++)
    for (size_t j = 0; j < op_->NStates; j++)
    products[i] += data_[i][j].dot(other.data_[i][j]).item<double>();
    return std::accumulate(products.begin(), products.end(), 0.0);
}
void Wfn::operator=(const double & scalar) {
    #pragma omp parallel for
    for (size_t i = 0; i < op_->NSegs; i++)
    for (size_t j = 0; j < op_->NStates; j++)
    data_[i][j].fill_(scalar);
}
void Wfn::mul(const double & scalar, Wfn & result) const {
    #pragma omp parallel for
    for (size_t i = 0; i < op_->NSegs; i++)
    for (size_t j = 0; j < op_->NStates; j++)
    result.data_[i][j].copy_(scalar * data_[i][j]);
}
void Wfn::operator*=(const double & scalar) {
    #pragma omp parallel for
    for (size_t i = 0; i < op_->NSegs; i++)
    for (size_t j = 0; j < op_->NStates; j++)
    data_[i][j].mul_(scalar);
}
void Wfn::sub_(const double & c, const Wfn & sub) {
    #pragma omp parallel for
    for (size_t i = 0; i < op_->NSegs; i++)
    for (size_t j = 0; j < op_->NStates; j++)
    data_[i][j] -= c * sub.data_[i][j];
}
void Wfn::sub2_(const double & c1, const Wfn & sub1, const double & c2, const Wfn & sub2) {
    #pragma omp parallel for
    for (size_t i = 0; i < op_->NSegs; i++)
    for (size_t j = 0; j < op_->NStates; j++)
    data_[i][j] -= c1 * sub1.data_[i][j] + c2 * sub2.data_[i][j];
}

} // namespace vibron