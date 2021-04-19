#ifndef Lanczos_Hamiltonian_SAPSet_hpp
#define Lanczos_Hamiltonian_SAPSet_hpp

#include <Lanczos/Hamiltonian/SAP.hpp>

namespace Lanczos {

// An Hd elements is a linear combination of symmetry adapted polynomials
class SAPSet {
    private:
        // 0th order term (constant) is more convenient to be treated separately
        double constant_ = 0.0;
        // 1st and higher order terms
        std::vector<std::pair<double, SAP>> terms_;

        // Highest possible excitation difference can be coupled by this Hd element
        // i.e. max number of monomials among the symmetry adapted polynomials 
        size_t max_excitation_;
        // A view to `terms_` grouped by the number of monomials
        std::vector<std::vector<const std::pair<double, SAP> *>> excitations_;

        // Construct `excitation_` and `excitations_` given constructed `vibrations_`
        void construct_excitation();
    public:
        SAPSet();
        SAPSet(const std::string & Hd_file);
        ~SAPSet();

        const double & constant() const;

        size_t size() const;
        const std::pair<double, SAP> & operator[](const size_t & index) const;
        std::vector<std::pair<double, SAP>>::const_iterator begin() const noexcept;
        std::vector<std::pair<double, SAP>>::const_iterator end() const noexcept;

        const size_t & max_excitation() const;
        const std::vector<const std::pair<double, SAP> *> & excitation(const size_t & ex) const;

        // Return the symmetry adapted polynomial set value given normal coordinate Q
        double operator()(const std::vector<std::vector<double>> & Q) const;
};

} // namespace Lanczos

#endif