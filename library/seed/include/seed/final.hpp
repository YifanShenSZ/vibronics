#ifndef seed_final_hpp
#define seed_final_hpp

#include <vibron/wfn.hpp>

#include <seed/overlap.hpp>

namespace seed {

// Use final-state normal modes to define vibration
class Final {
    private:
        // one component or norm of transition dipole
        std::vector<double> dipole_;
        // vibronic wave function definition
        std::shared_ptr<vibron::Options> op_;
        // <initial vibrational state|final-state-biased vibrational basis function>
        std::shared_ptr<IFOverlap> ifoverlap_;
    public:
        Final();
        Final(const std::vector<double> & _dipole, const std::shared_ptr<vibron::Options> & _op,
              const at::Tensor & init_freq, const at::Tensor & final_freq,
              const at::Tensor & tran_matrix, const at::Tensor & shift_vector);
        ~Final();

        // Generate seed vector, return the norm of the seed vector and the corresponding unit vector
        double generate_seed(vibron::Wfn & wfn) const;
};

} // namespace seed

#endif