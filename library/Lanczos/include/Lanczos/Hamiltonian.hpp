// The vibronic Hamiltonian consists of 3 parts:
//     Hvib = T + Hharmonic + Hanharmonic
// the vibrational basis functions are the eigen functions of T + Hharmonic
// Hanharmonic = diabatic electronic Hamiltonian - Hharmonic
// for convenience we also call Hanharmonic by "Hd" here, since T + Hharmonic is trivial

#ifndef Lanczos_Hamiltonian_hpp
#define Lanczos_Hamiltonian_hpp

#include <set>

#include <CppLibrary/utility.hpp>
#include <CppLibrary/linalg.hpp>

#include <Lanczos/Hamiltonian/SAPSet.hpp>

namespace Lanczos {

// An Hd matrix is a collection of elements
class Hd {
    private:
        size_t NStates_;
        std::vector<size_t> NModes_;
        CL::utility::matrix<SAPSet *> Hd_;

        // The highest order of each normal mode
        std::vector<std::vector<size_t>> max_orders_;

        // Highest possible excitation difference can be coupled by this anharmonic diabatic electronic Hamiltonian
        size_t max_excitation_;
        // Possbile ways to couple different modes grouped by excitation
        // for triple and higher excitations only, since usualy
        // singles and doubles are all coupled but only a few triples and higher
        std::vector<std::set<std::vector<std::pair<size_t, size_t>>> *> excitations_;

        // Construct `max_orders_`, `max_excitation_`, `excitations_`
        void construct_max();
    public:
        Hd();
        Hd(const size_t & _NStates, const std::vector<size_t> & _NModes, const std::vector<std::string> & Hd_files);
        ~Hd();

        const size_t & NStates() const;
        const size_t & max_order(const size_t & irred, const size_t & mode) const;
        const size_t & max_order(const std::pair<size_t, size_t> & irred_mode) const;
        const size_t & max_excitation() const;
        const std::set<std::vector<std::pair<size_t, size_t>>> * excitation(const size_t & index) const; 

        void pretty_print(std::ostream & stream) const;

        const SAPSet * operator[](const std::pair<size_t, size_t> & indices) const;
};

} // namespace Lanczos

#endif