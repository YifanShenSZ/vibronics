#ifndef NCKernel_hpp
#define NCKernel_hpp

#include <tchem/intcoord.hpp>

class NCKernel {
    private:
        std::shared_ptr<tchem::IC::SASICSet> sasicset_;
        std::vector<at::Tensor> Linvs_;
    public:
        NCKernel();
        NCKernel(const std::shared_ptr<tchem::IC::SASICSet> & _sasicset, const std::vector<at::Tensor> & _Linvs);
        ~NCKernel();

        std::vector<std::vector<double>> operator()(const at::Tensor & x) const;
};

#endif