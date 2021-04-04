#ifndef seed_overlap_hpp
#define seed_overlap_hpp

#include <vibron/wfn.hpp>

namespace seed {

// When adopting final-state normal modes to define vibration,
// it is not trivial to expand the initial vibrational state, i.e.
// <initial vibrational state|final-state-biased vibrational basis function>
// has to be evaluated from a complicated recursion
// Currently only support <0|
class IFOverlap {
    private:
        // vibronic wave function definition
        std::shared_ptr<vibron::Options> op_;

        // Given op_->max_phonons, generate a vibration set including all possible vibrational basis functions
        // The normal modes of different irreducibles are concatenated
        std::shared_ptr<vibron::VibrationSet> vib_set_;
        // Store all possible integrals g
        std::vector<double> integrals_;
    public:
        IFOverlap();
        // Let the initial-state normal coordiante q and the final-state normal coordinate Q be related by:
        //    q = T . Q + d
        IFOverlap(const std::shared_ptr<vibron::Options> & _op,
                  const at::Tensor & init_freq, const at::Tensor & final_freq,
                  const at::Tensor & T, const at::Tensor & d);
        ~IFOverlap();

        const double & operator()(const size_t & index) const;
        const double & operator()(const vibron::Vibration &vibration) const;
};

}; // namespace seed

#endif