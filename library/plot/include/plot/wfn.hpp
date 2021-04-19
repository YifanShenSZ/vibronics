#ifndef plot_wfn_hpp
#define plot_wfn_hpp

#include <vibron/wfn.hpp>

namespace plot {

// To evaluate vibronic wave function value
class Wfn : public vibron::Wfn {
    private:
        std::vector<std::vector<double>> frequencies_;
    public:
        Wfn();
        Wfn(const std::shared_ptr<vibron::Options> & _op, const std::vector<std::vector<double>> & _frequencies);
        ~Wfn();

        // Given normal coordinate Q, return vibronic wave function value
        std::vector<double> operator()(const std::vector<std::vector<double>> & Qs) const;
};

} // namespace plot

#endif