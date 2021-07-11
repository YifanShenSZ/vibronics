#ifndef Lanczos_mv_hpp
#define Lanczos_mv_hpp

#include <vibron/harmonic.hpp>
#include <vibron/wfn.hpp>

#include <Lanczos/Hamiltonian.hpp>

namespace Lanczos {

struct SegStateVibValue {
    uint16_t seg, state;
    uint32_t vib;
    double value;

    SegStateVibValue();
    SegStateVibValue(const uint16_t & _seg, const uint16_t & _state, const uint32_t & _vib, const double & _value);
    ~SegStateVibValue();
};

class MVKernel {
    private:
        // diabatic electronic Hamiltonian definition
        std::shared_ptr<Hd> Hd_;
        // vibronic wave function definition
        std::shared_ptr<vibron::Options> op_;
        // vibronic Hamiltonian matrix element (integral) evaluators
        std::vector<std::vector<vibron::Integrator>> integrators_;

        // the non-trivial vibronic Hamiltonian matrix elements
        // basis[i, j, k] can only have coupling with basis in alloweds_[i][j][k]
        // let element = alloweds_[i][j][k][l]
        // element.value = <basis[i, j, k]|Hd|basis[element.seg, element.state, element.vib]>
        std::vector<std::vector<std::vector<std::vector<SegStateVibValue>>>> alloweds_;

        // Construct `alloweds_` given constructed `Hd_`, `op_`, `integrator_`
        void construct_nonzero();

        // Supports `construct_nonzero`
        // Given excited modes, generate all possible vibrations
        void generate_all(const size_t & iseg, const size_t & istate, const size_t & ivib,
        const std::vector<std::pair<size_t, size_t>> & excited_modes);

        // Supports `construct_nonzero`
        // Return <basis[iseg, istate, ivib]|Hd|basis[jseg, jstate, jvib]>
        double Hdelement(const size_t & iseg, const size_t & istate, const size_t & ivib,
                         const size_t & jseg, const size_t & jstate, const size_t & jvib) const;
        double Hdelement(const size_t & iseg, const size_t & istate, const size_t & ivib,
                         const size_t & jseg, const size_t & jstate, const size_t & jvib,
                         const std::vector<std::pair<size_t, size_t>> & excited_modes) const;
    public:
        MVKernel();
        MVKernel(const std::shared_ptr<Hd> & _Hd, const std::shared_ptr<vibron::Options> & _op,
        const std::vector<std::string> & frequency_files);
        ~MVKernel();

        void pretty_print(std::ostream & stream) const;

        void operator()(const vibron::Wfn & wfn, vibron::Wfn & Hwfn) const;
};

} // namespace Lanczos

#endif