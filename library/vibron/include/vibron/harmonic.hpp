#ifndef vibron_harmonic_hpp
#define vibron_harmonic_hpp

#include <cstddef>
#include <cstdint>
#include <vector>

#include <CppLibrary/utility.hpp>

namespace vibron {

struct BraOrderKet {
    size_t bra, order, ket;

    BraOrderKet();
    BraOrderKet(const size_t & _bra, const size_t & _order, const size_t & _ket);
    BraOrderKet(const std::initializer_list<size_t> & items);
    ~BraOrderKet();
};

// Compute 1-dimensional harmonic oscillator integratals <m|Q^a|n>
class Integrator {
    private:
        // Harmonic frequency
        double frequency_;
        // Highest phonon, highest polynomial order
        size_t phonon_, order_;
        // integrals_[a][m][n] = <m|Q^a|n>
        std::vector<std::vector<std::vector<double>>> integrals_;
    public:
        Integrator();
        Integrator(const double & _frequency, const size_t & _phonon, const size_t & _order);
        ~Integrator();

        const double & frequency() const;
        const size_t & phonon() const;
        const size_t & order() const;

        double operator()(const BraOrderKet & braorderket) const;
};

} // namespace vibron

#endif