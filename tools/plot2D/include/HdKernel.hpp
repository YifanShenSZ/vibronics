#ifndef HdKernel_hpp
#define HdKernel_hpp

#include <torch/torch.h>

#include <Lanczos/Hamiltonian.hpp>

class HdKernel {
    private:
        std::vector<std::vector<double>> freqs_;
        std::shared_ptr<Lanczos::Hd> Hanharmonic_;
    public:
        HdKernel();
        HdKernel(const std::vector<std::vector<double>> & _freqs, const std::shared_ptr<Lanczos::Hd> & _Hanharmonic);
        ~HdKernel();

        at::Tensor operator()(const std::vector<std::vector<double>> & Qs) const;
};

#endif