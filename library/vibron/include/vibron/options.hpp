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
        // number of irreducible representations
        size_t NIrreds;
        // point group product table
        CL::utility::matrix<size_t> product_table;
        // number of normal modes per irreducible
        std::vector<size_t> NModes;
        // number of electronic states
        size_t NStates;
        // vibrational irreducible of each electronic state
        std::vector<size_t> vib_irreds;
        // vibrational basis function details of each irreducible
        std::vector<VibrationSet *> vib_sets;

        // The highest phonon of each normal mode
        std::vector<std::vector<uint16_t>> max_phonons;

        // number of segmentations
        size_t NSegs;
        // k in wfn[i, j, k] corresponds to VibrationSet[k + starts[i][j]] 
        std::vector<std::vector<size_t>> starts, stops;

        Options();
        Options(const std::string & wfn_file, const std::vector<std::string> & vib_files);
        ~Options();

        void pretty_print(std::ostream & stream) const;

        // internal coordinate dimension = sum(NModes)
        size_t intdim() const;
        // Given the index of a normal mode, return its irreducible and index within the irreducible
        std::pair<size_t, size_t> vib_irred_mode(const size_t & index) const;
        // Given the irreducible and the index within the irreducible of a normal mode, return its index
        size_t vib_C1_index(const size_t & irred, const size_t & mode) const;

        // number of vibronic basis functions
        size_t NVibron() const;

        // Given a vibrational basis function defined by phonons, return its irreducible
        size_t vib_irred(const std::vector<std::vector<size_t>> & phonons) const;
        size_t vib_irred(const std::vector<std::vector<uint16_t>> & phonons) const;

        // Given segment, electronic state, vibrational index in this segment,
        // return the vibrational index in the vibrational basis function set
        size_t vib_index_abs(const size_t & seg, const size_t & state, const size_t & vib) const;
        size_t vib_index_abs(const CL::utility::triple<size_t, size_t, size_t> & seg_state_vib) const;
        // Given electronic state, vibrational index in the vibrational basis function set
        // return segment and vibrational index in this segment
        std::pair<size_t, size_t> vib_index(const size_t & state, const size_t & abs_vib) const;

        // Given the index in the concatenated wave function vector,
        // return segment, electronic state, vibrational index in this segment
        CL::utility::triple<size_t, size_t, size_t> seg_state_vib(const size_t & index) const;
};

} // namespace vibron

#endif