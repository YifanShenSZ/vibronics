#ifndef Lanczos_Lanczos_hpp
#define Lanczos_Lanczos_hpp

#include <vibron/wfn.hpp>

#include <Lanczos/mv.hpp>

namespace Lanczos {

class Kernel {
    private:
        std::shared_ptr<MVKernel> mvkernel_;
    public:
        Kernel();
        Kernel(const std::shared_ptr<MVKernel> & _mvkernel);
        ~Kernel();

        void pretty_print(std::ostream & stream) const;

        // w = H . v
        // alpha = w . v
        // w = w - alpha * v
        // return alpha, modify w
        double initialize(const vibron::Wfn & v, vibron::Wfn & w) const;

        // beta = ||w||
        // v = w / beta
        // w = H . v
        // alpha = w . v
        // w = w - alpha * v - beta * vold
        // return (alpha, beta), modify v and w
        std::tuple<double, double> iterate(const vibron::Wfn & vold, vibron::Wfn & v, vibron::Wfn & w) const;
};

} // namespace Lanczos

#endif