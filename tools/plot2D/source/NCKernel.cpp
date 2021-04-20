#include "NCKernel.hpp"

NCKernel::NCKernel() {}
NCKernel::NCKernel(const std::shared_ptr<tchem::IC::SASICSet> & _sasicset, const std::vector<at::Tensor> & _Linvs)
: sasicset_(_sasicset), Linvs_(_Linvs) {}
NCKernel::~NCKernel() {}

std::vector<std::vector<double>> NCKernel::operator()(const at::Tensor & x) const {
    if (x.sizes().size() != 1) throw std::invalid_argument(
    "NCKernel::operator(): x must be a vector");
    std::vector<at::Tensor> qs(sasicset_->NIrreds());
    // x is internal coordinate vector q
    if (x.size(0) == sasicset_->intdim()) {
        auto NSASICs = sasicset_->NSASICs();
        size_t start = 0;
        for (size_t i = 0; i < qs.size(); i++) {
            size_t stop = start + NSASICs[i];
            qs[i] = x.slice(0, start, stop);
            start = stop;
        }
    }
    // x is Cartesian coordinate vector r
    else qs = (*sasicset_)(sasicset_->tchem::IC::IntCoordSet::operator()(x));
    std::vector<std::vector<double>> Qs(qs.size());
    for (size_t i = 0; i < qs.size(); i++) {
        at::Tensor Q = Linvs_[i].mv(qs[i]);
        Qs[i].resize(Q.numel());
        std::memcpy(Qs[i].data(), Q.data_ptr<double>(), Q.numel() * sizeof(double));
    }
    return Qs;
}