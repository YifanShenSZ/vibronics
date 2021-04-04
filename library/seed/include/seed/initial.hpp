#ifndef seed_initial_hpp
#define seed_initial_hpp

#include <vibron/wfn.hpp>

namespace seed {

// Use initial-state normal modes to define vibration
class Initial {
    private:
        // One component or norm of transition dipole
        std::vector<double> dipole_;
        // vibronic wave function definition
        std::shared_ptr<vibron::Options> op_;
        // initial vibrational state phonons (under harmonic approximation)
        std::vector<std::vector<size_t>> phonons_;
    public:
        Initial();
        Initial(const std::vector<double> & _dipole,
                const std::shared_ptr<vibron::Options> & _op, const std::vector<std::vector<size_t>> & _phonons);
        ~Initial();

        // Generate seed vector, return the norm of the seed vector and the corresponding unit vector
        double generate_seed(vibron::Wfn & wfn) const;
};

} // namespace seed

#endif