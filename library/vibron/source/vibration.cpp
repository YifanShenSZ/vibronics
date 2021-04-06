#include <fstream>
#include <numeric>

#include <CppLibrary/utility.hpp>

#include <vibron/vibration.hpp>

namespace vibron {

void Vibration::construct_excitation() {
    excitation_ = 0;
    excited_modes_.clear();
    for (size_t irred = 0; irred < phonons_.size(); irred++)
    for (size_t mode = 0; mode < phonons_[irred].size(); mode++)
    if (phonons_[irred][mode] > 0) {
        excitation_++;
        excited_modes_.push_back(std::pair<size_t, size_t>(irred, mode));
    }
    excited_modes_.shrink_to_fit();
}

Vibration::Vibration() {}
Vibration::Vibration(const std::vector<std::vector<size_t>> & _phonons)
: phonons_(_phonons) {this->construct_excitation();}
Vibration::Vibration(const std::vector<std::string> & lines) {
    phonons_.resize(lines.size());
    for (size_t i = 0; i < lines.size(); i++) {
        std::vector<std::string> strs = CL::utility::split(lines[i]);
        phonons_[i].resize(strs.size());
        for (size_t j = 0; j < strs.size(); j++)
        phonons_[i][j] = std::stoul(strs[j]);
    }
    this->construct_excitation();
}
Vibration::~Vibration() {}

const std::vector<std::vector<size_t>> & Vibration::phonons() const {return phonons_;}
const size_t & Vibration::excitation() const {return excitation_;}
const std::vector<std::pair<size_t, size_t>> & Vibration::excited_modes() const {return excited_modes_;}

void Vibration::pretty_print(std::ostream & stream) const {
    for (size_t i = 0; i < phonons_.size(); i++) {
        stream << "Irreducible " << i + 1 << ":\n";
        for (const auto & phonon : phonons_[i]) stream << phonon << ' ';
        stream << '\n';
    }
}

Vibration Vibration::create(const size_t & irred, const size_t & mode, const int64_t & phonon) const {
    if (irred >= this->phonons_.size()) throw std::invalid_argument(
        "vibron::vib::Vibration::create: irred out of range"
    );
    if (mode >= this->phonons_[irred].size()) throw std::invalid_argument(
        "vibron::vib::Vibration::create: mode out of range"
    );
    Vibration vibration(this->phonons_);
    vibration.phonons_[irred][mode] += phonon;
    if (vibration.phonons_[irred][mode] < 0) throw std::invalid_argument(
        "vibron::vib::Vibration::create: phonon out of range"
    );
    return vibration;
}





// Construct `max_phonons_`, `max_excitation_` and `excitations_` based on constructed `vibrations_`
void VibrationSet::construct_exciations_() {
    // quick return if possible
    if (vibrations_.empty()) {
        max_excitation_ = 0;
        excitations_.resize(max_excitation_ + 1);
        return;
    }
    // Find out the highest phonons and excitation among the vibrational basis functions
    max_phonons_ = vibrations_[0].phonons();
    max_excitation_ = 0;
    for (const auto & vibration : vibrations_) {
        // Compare phonons
        const auto & phonons = vibration.phonons();
        for (size_t i = 0; i < phonons.size(); i++)
        for (size_t j = 0; j < phonons[i].size(); j++)
        if (phonons[i][j] > max_phonons_[i][j])
        max_phonons_[i][j] = phonons[i][j];
        // Compare excitation
        if (vibration.excitation() > max_excitation_)
        max_excitation_ = vibration.excitation();
    }
    // Construct a view to `vibrations_` grouped by excitation
    excitations_.clear();
    excitations_.resize(max_excitation_ + 1);
    for (const auto & vibration : vibrations_)
    excitations_[vibration.excitation()].push_back(& vibration);
    for (auto & excitation : excitations_) excitation.shrink_to_fit();
}

// Support `VibrationSet(const std::vector<size_t> & max_phonons)`
void VibrationSet::generate_all_(const std::vector<size_t> & excited_modes, const std::vector<size_t> & max_phonons) {
    // quick return if there are no such excitations
    for (const size_t & excited_mode : excited_modes)
    if (max_phonons[excited_mode] == 0) return;
    // basic case: the excited modes has |1>, others |0>
    std::vector<size_t> phonons(max_phonons.size(), 0);
    for (const size_t & excited_mode : excited_modes) phonons[excited_mode] = 1;
    vibrations_.push_back(Vibration({phonons}));
    // Loop as a excited_modes.size()-nary counter
    while (true) {
        phonons[excited_modes[0]]++;
        // Carry to latter digit
        for (size_t i = 1; i < excited_modes.size(); i++)
        if (phonons[excited_modes[i - 1]] > max_phonons[excited_modes[i - 1]]) {
            phonons[excited_modes[i - 1]] = 1;
            phonons[excited_modes[i]]++;
        }
        // Finish when counter overflows
        if (phonons[excited_modes.back()] > max_phonons[excited_modes.back()]) break;
        vibrations_.push_back(Vibration({phonons}));
    }
}

// Support `index_vibration`
// Given a vibrational basis function, try to locate its index within [lower, upper]
// index = -1 if not found
void VibrationSet::bisect_(const Vibration & vibration, const size_t & lower, const size_t & upper, int64_t & index) const {
    // Final round
    if (upper - lower == 1) {
        // Try lower
        bool match = true;
        const auto & phonons = vibration.phonons(),
                   & phonons_lower = vibrations_[lower].phonons();
        for (size_t i = 0; i < phonons.size(); i++)
        for (size_t j = 0; j < phonons[i].size(); j++)
        if (phonons[i][j] != phonons_lower[i][j]) {
            match = false;
            break;
        }
        if (match) {
            index = lower;
            return;
        }
        // Try upper
        match = true;
        const auto & phonons_upper = vibrations_[upper].phonons();
        for (size_t i = 0; i < phonons.size(); i++)
        for (size_t j = 0; j < phonons[i].size(); j++)
        if (phonons[i][j] != phonons_upper[i][j]) {
            match = false;
            break;
        }
        if (match) {
            index = upper;
            return;
        }
        // Neither
        index = -1;
    }
    // Normal bisection
    else {
        // Try bisection
        size_t bisection = (lower + upper) / 2;
        bool match = true;
        // 1st, compare excited modes
        const auto & modes = vibration.excited_modes(),
                   & modes_ref = vibrations_[bisection].excited_modes();
        size_t i;
        for (i = 0; i < modes.size(); i++)
        if (modes[i] != modes_ref[i]) {
            match = false;
            break;
        }
        // 2nd, compare phonons
        if (match) {
            const auto & phonons = vibration.phonons(),
                       & phonons_ref = vibrations_[bisection].phonons();
            int64_t i, j;
            for (i = phonons.size() - 1; i >= 0; i--) {
                for (j = phonons[i].size() - 1; j >= 0; j--)
                if (phonons[i][j] != phonons_ref[i][j]) {
                    match = false;
                    break;
                }
                if (! match) break;
            }
            if (match) {
                index = bisection;
                return;
            }
            // Next range
            else {
                if (phonons[i][j] > phonons_ref[i][j]) bisect_(vibration, bisection, upper, index);
                else                                   bisect_(vibration, lower, bisection, index);
            }
        }
        // Next range
        else {
            if (modes[i] > modes_ref[i]) bisect_(vibration, bisection, upper, index);
            else                         bisect_(vibration, lower, bisection, index);
        }
    }
}

VibrationSet::VibrationSet() {}
VibrationSet::VibrationSet(const std::string & vib_file, const size_t & NIrreds) {
    std::ifstream ifs; ifs.open(vib_file);
    if (! ifs.good()) throw CL::utility::file_error(vib_file);
    while (true) {
        std::vector<std::string> lines(NIrreds);
        for (std::string & line : lines) std::getline(ifs, line);
        if (! ifs.good()) break;
        vibrations_.push_back(Vibration(lines));
    }
    ifs.close();
    vibrations_.shrink_to_fit();
    this->construct_exciations_();
}
// Generate all possible vibrational basis functions given the max phonon of each normal mode, assume C1 symmetry
VibrationSet::VibrationSet(const std::vector<size_t> & max_phonons) {
    size_t intdim = max_phonons.size();
    size_t size = 1, max_excitation = 0;
    std::vector<size_t> possible_modes;
    for (size_t i = 0; i < max_phonons.size(); i++) {
        size *= max_phonons[i] + 1;
        if (max_phonons[i] > 0) {
            max_excitation++;
            possible_modes.push_back(i);
        }
    }
    possible_modes.shrink_to_fit();
    // |0>
    vibrations_.push_back(Vibration({std::vector<size_t>(intdim, 0)}));
    // excited ones
    for (size_t excitation = 1; excitation < max_excitation + 1; excitation++) {
        // basic case: the leading excitation modes are excited
        std::vector<size_t> excited_indices(excitation);
        std::iota(excited_indices.begin(), excited_indices.end(), 0);
        std::vector<size_t> excited_modes(excitation);
        for (size_t i = 0; i < excitation; i++) excited_modes[i] = possible_modes[excited_indices[i]];
        generate_all_(excited_modes, max_phonons);
        // Loop over possible excited modes as an excitation-nary counter, with ascending digits
        while (true) {
            excited_indices.back()++;
            // Carry to former digit
            for (int64_t i = - 1; i > -excitation; i--)
            if (excited_indices[excitation + i] >= max_excitation + 1 + i) {
                excited_indices[excitation + i - 1]++;
                excited_indices[excitation + i] = 0;
            }
            // Guarantee ascendance
            for (size_t i = 0; i < excitation - 1; i++)
            if (excited_indices[i] >= excited_indices[i + 1])
            excited_indices[i + 1] = excited_indices[i] + 1;
            // Finish when counter overflows
            if (excited_indices.back() >= max_excitation) break;
            for (size_t i = 0; i < excitation; i++) excited_modes[i] = possible_modes[excited_indices[i]];
            generate_all_(excited_modes, max_phonons);
        }
    }
    vibrations_.shrink_to_fit();
    this->construct_exciations_();
}
VibrationSet::~VibrationSet() {}

size_t VibrationSet::size() const {return vibrations_.size();}
const std::vector<std::vector<size_t>> & VibrationSet::max_phonons() const {return max_phonons_;}
const size_t & VibrationSet::max_phonon(const size_t & irred, const size_t & mode) const {return max_phonons_[irred][mode];}
const size_t & VibrationSet::max_phonon(const std::pair<size_t, size_t> & irred_mode) const {return max_phonons_[irred_mode.first][irred_mode.second];}
const size_t & VibrationSet::max_excitation() const {return max_excitation_;}
// A read-only accessor to vibrations_[index]
const Vibration & VibrationSet::operator[](const size_t & index) const {return vibrations_[index];}

void VibrationSet::pretty_print(std::ostream & stream) const {
    stream << "Highest phonons among the vibrational basis functions:\n";
    for (const auto & irred : max_phonons_) {
        for (const auto & phonon : irred) stream << phonon << ", ";
        stream << '\n';
    }
}

// Given a vibrational basis function, return its index in this vibration set
// Return -1 if not found
int64_t VibrationSet::index_vibration(const Vibration & vibration) const {
    if (vibrations_.empty()) return -1; // quick return if there is no vibrational basis function
    size_t lower = 0;
    for (size_t i = 0; i < vibration.excitation(); i++) lower += excitations_[i].size();
    size_t upper = lower + excitations_[vibration.excitation()].size() - 1;
    int64_t index;
    this->bisect_(vibration, lower, upper, index);
    return index;
}

} // namespace vibron