#ifndef vibron_vibration_set_hpp
#define vibron_vibration_set_hpp

#include <vibron/vibration/vibration.hpp>

namespace vibron {

class VibrationSet {
    private:
        // Vibrational basis functions constituting the set, requirements:
        //     1. excitations are sorted ascendingly
        //     2. same excitation terms have the excited modes sorted ascendingly,
        //        where the comparison is made from the first mode to the last
        //     3. same excited modes terms have phonons sorted ascendingly,
        //        where the comparison is made from the last mode to the first
        // e.g. 1-irreducible 3-modes up to 2 phonons per mode:
        //     |0, 0, 0>
        //     |1, 0, 0>, |2, 0, 0>, |0, 1, 0>, |0, 2, 0>, |0, 0, 1>, |0, 0, 2>
        //     |1, 1, 0>, |2, 1, 0>, |1, 2, 0>, |2, 2, 0>
        //     |1, 0, 1>, |2, 0, 1>, |1, 0, 2>, |2, 0, 2>
        //     |0, 1, 1>, |0, 2, 1>, |0, 1, 2>, |0, 2, 2>
        //     |1, 1, 1>, |2, 1, 1>, |1, 2, 1>, |2, 2, 1>, |1, 1, 2>, |2, 1, 2>, |1, 2, 2>, |2, 2, 2>
        std::vector<Vibration> vibrations_;

        // Highest phonons among the vibrational basis functions
        std::vector<std::vector<uint16_t>> max_phonons_;
        // Highest excitation among the vibrational basis functions
        uint16_t max_excitation_;
        // A view to `vibrations_` grouped by excitation
        std::vector<std::vector<const Vibration *>> excitations_;

        // Construct `max_phonons_`, `max_excitation_` and `excitations_` based on constructed `vibrations_`
        void construct_exciations_();

        // Support `VibrationSet(const std::vector<size_t> & max_phonons)`
        void generate_all_(const std::vector<size_t> & excited_modes, const std::vector<size_t> & max_phonons);

        // Support `index_vibration`
        // Given a vibrational basis function, try to locate its index within [lower, upper]
        // index = -1 if not found
        void bisect_(const Vibration & vibration, const size_t & lower, const size_t & upper, int64_t & index) const;
    public:
        VibrationSet();
        // `vib_file` contains one Vibration per NIrreds lines
        VibrationSet(const std::string & vib_file, const size_t & NIrreds);
        // Generate all possible vibrational basis functions given the max phonon of each normal mode, assume C1 symmetry
        VibrationSet(const std::vector<size_t> & max_phonons);
        ~VibrationSet();

        size_t size() const;
        const std::vector<std::vector<uint16_t>> & max_phonons() const;
        const uint16_t & max_phonon(const size_t & irred, const size_t & mode) const;
        const uint16_t & max_phonon(const std::pair<size_t, size_t> & irred_mode) const;
        const uint16_t & max_excitation() const;
        // A read-only accessor to vibrations_[index]
        const Vibration & operator[](const size_t & index) const;

        std::vector<Vibration>::const_iterator begin() const noexcept {return vibrations_.begin();}
        std::vector<Vibration>::const_iterator end() const noexcept {return vibrations_.end();}

        void pretty_print(std::ostream & stream) const;

        // Given a vibrational basis function, return its index in this vibration set
        // Return -1 if not found
        int64_t index_vibration(const Vibration & vibration) const;
};

} // namespace vibron

#endif