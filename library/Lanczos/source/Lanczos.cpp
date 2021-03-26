#include <Lanczos/Lanczos.hpp>

namespace Lanczos {

Kernel::Kernel() {}
Kernel::Kernel(const std::shared_ptr<MVKernel> & _mvkernel) : mvkernel_(_mvkernel) {}
Kernel::~Kernel() {}

// w = H . v
// alpha = w . v
// w = w - alpha * v
// return alpha, modify w
double Kernel::initialize(const vibron::Wfn & v, vibron::Wfn & w) const {
    (*mvkernel_)(v, w);
    double alpha = w.dot(v);
    w.sub_(alpha, v);
    return alpha;
}

// beta = ||w||
// v = w / beta
// w = H . v
// alpha = w . v
// w = w - alpha * v - beta * vold
// return (alpha, beta), modify v and w
std::tuple<double, double> Kernel::iterate(const vibron::Wfn & vold, vibron::Wfn & v, vibron::Wfn & w) const {
    double beta = sqrt(w.dot(w));
    w.mul(1.0 / beta, v);
    (*mvkernel_)(v, w);
    double alpha = w.dot(v);
    w.sub2_(alpha, v, beta, vold);
    return std::make_tuple(alpha, beta);
}

} // namespace Lanczos