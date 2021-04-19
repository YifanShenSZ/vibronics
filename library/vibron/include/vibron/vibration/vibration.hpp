#ifndef vibron_vibration_vibration_hpp
#define vibron_vibration_vibration_hpp

namespace vibron {

// A vibrational basis function made up with
// a direct product of 1D harmonic oscillator wave functions
class Vibration {
    private:
        // The j-th normal mode in i-th irreducible representation
        // contributes | phonons[i][j] > to the direct product
        std::vector<std::vector<size_t>> phonons_;

        // Number of excited normal modes
        size_t excitation_;
        // The irreducible and the index of each excited normal mode
        std::vector<std::pair<size_t, size_t>> excited_modes_;

        // Construct `excitation_` and `excited_modes_` given constructed `phonons_`
        void construct_excitation();
    public:
        Vibration();
        Vibration(const std::vector<std::vector<size_t>> & _phonons);
        // Each line contains the phonons of all normal modes in one irreducible
        Vibration(const std::vector<std::string> & lines);
        ~Vibration();

        const std::vector<std::vector<size_t>> & phonons() const;
        const size_t & excitation() const;
        const std::vector<std::pair<size_t, size_t>> & excited_modes() const;

        void pretty_print(std::ostream & stream) const;

        // Given frequency `w` and normal coordiante `Q`, return harmonic oscillator wave function value
        double value(const std::vector<std::vector<double>> & ws, const std::vector<std::vector<double>> & Qs) const;
};

} // namespace vibron

#endif