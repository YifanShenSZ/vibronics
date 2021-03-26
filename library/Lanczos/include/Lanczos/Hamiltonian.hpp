// The vibronic Hamiltonian consists of 3 parts:
//     Hvib = T + Hharmonic + Hd
// the vibrational basis functions are the eigen functions of T + Hharmonic
// Hd is the diabatic electronic Hamiltonian without Hharmonic

#ifndef Lanczos_Hamiltonian_hpp
#define Lanczos_Hamiltonian_hpp

#include <CppLibrary/utility.hpp>

namespace Lanczos {

struct Monomial {
    size_t irred, mode, order;

    Monomial();
    Monomial(const size_t & _irred, const size_t & _mode, const size_t & _order);
    ~Monomial();

    bool operator<(const Monomial & other) const;
    bool operator==(const std::pair<size_t, size_t> & irred_mode) const;
    bool operator!=(const std::pair<size_t, size_t> & irred_mode) const;
};

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

        bool same_modes(const std::vector<std::pair<size_t, size_t>> & excited_modes) const;
};

// An Hd elements is a linear combination of symmetry adapted polynomials
class SAPSet {
    private:
        // 0th order term (constant) is more convenient to be treated separately
        double constant_;
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
};

// An Hd matrix is a collection of elements
class Hd {
    private:
        size_t NStates_;
        std::vector<size_t> NModes_;
        CL::utility::matrix<SAPSet *> Hd_;

        // Highest possible excitation difference can be coupled by this Hd matrix
        size_t max_excitation_;
        // The highest order of each normal mode
        std::vector<std::vector<size_t>> max_orders_;

        // Construct `max_excitation_`, `max_orders_`
        void construct_max();
    public:
        Hd();
        Hd(const size_t & _NStates, const std::vector<size_t> & _NModes, const std::vector<std::string> & Hd_files);
        ~Hd();

        const size_t & NStates() const;
        const size_t & max_excitation() const;
        const size_t & max_order(const size_t & irred, const size_t & mode) const;
        const size_t & max_order(const std::pair<size_t, size_t> & irred_mode) const;

        const SAPSet * operator[](const std::pair<size_t, size_t> & indices) const;
};

} // namespace Lanczos

#endif