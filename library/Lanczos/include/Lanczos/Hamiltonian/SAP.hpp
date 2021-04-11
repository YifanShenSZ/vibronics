#ifndef Lanczos_Hamiltonian_SAP_hpp
#define Lanczos_Hamiltonian_SAP_hpp

#include <vector>
#include <string>

#include <Lanczos/Hamiltonian/monomial.hpp>

namespace Lanczos {

// A symmetry adapted polynomial
class SAP {
    private:
        // The monomials constituting the polynomial
        std::vector<Monomial> coords_;
    public:
        SAP();
        SAP(const std::string & line);
        ~SAP();

        size_t NMonomials() const;
        const Monomial & operator[](const size_t & index) const;
        std::vector<Monomial>::const_iterator begin() const noexcept;
        std::vector<Monomial>::const_iterator end() const noexcept;

        // Check if the monomials match the given normal modes
        bool operator==(const std::vector<std::pair<size_t, size_t>> & irred_modes) const;
        // Check if the monomials include the given normal modes
        bool operator>=(const std::vector<std::pair<size_t, size_t>> & irred_modes) const;
};

} // namespace Lanczos

#endif