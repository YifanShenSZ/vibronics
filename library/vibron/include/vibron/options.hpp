#ifndef vibron_options_hpp
#define vibron_options_hpp

#include <CppLibrary/utility.hpp>

#include <vibron/vibration.hpp>

namespace vibron {

struct Options {
    private:
        // Construct `max_phonons`
        void construct_phonons();
        // Construct `NSegs`, `starts`, `stops`
        void construct_segmentation();
    public:
        // basic information
        size_t NIrreds;
        CL::utility::matrix<size_t> product_table;
        std::vector<size_t> NModes;
        size_t NStates;
        std::vector<size_t> vib_irreds;
        // vibrational basis function details
        std::vector<VibrationSet> vib_sets;
    
        // The highest phonon of each normal mode
        std::vector<std::vector<size_t>> max_phonons;
    
        // segmentation
        size_t NSegs;
        // k in wfn[i, j, k] corresponds to VibrationSet[k + starts[i][j]] 
        std::vector<std::vector<size_t>> starts, stops;

        Options();
        Options(const std::string & wfn_file, const std::vector<std::string> & vib_files);
        ~Options();

        void pretty_print(std::ostream & stream) const;

        size_t intdim() const;
        // Given the index of a normal mode, return its irreducible and index within the irreducible
        std::pair<size_t, size_t> irred_mode(const size_t & index) const;
        // Given the irreducible and the index within the irreducible of a normal mode, return its index
        size_t index(const size_t & irred, const size_t & mode) const;

        size_t NVibron() const;

        // Given a vibrational basis function defined by phonons, return its irreducible
        size_t determine_irreducible(const std::vector<std::vector<size_t>> & phonons) const;

        // Given segment, electronic state, vibrational index in this segment,
        // return the vibrational index in the vibrational basis function set
        size_t abs_vib_index(const size_t & seg, const size_t & state, const size_t & vib) const;
        // Given electronic state, vibrational index in the vibrational basis function set
        // return segment and vibrational index in this segment
        std::pair<size_t, size_t> vib_index(const size_t & state, const size_t & abs_vib) const;
};

} // namespace vibron

#endif