#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <cmath>

#include <CppLibrary/utility.hpp>
#include <CppLibrary/math.hpp>

#include <vibron/vibration/vibration.hpp>

namespace vibron {

void Vibration::construct_excitation() {
    excitation_ = 0;
    excited_modes_.clear();
    for (uint16_t irred = 0; irred < phonons_.size(); irred++)
    for (uint16_t mode = 0; mode < phonons_[irred].size(); mode++)
    if (phonons_[irred][mode] > 0) {
        excitation_++;
        excited_modes_.push_back(std::pair<uint16_t, uint16_t>(irred, mode));
    }
    excited_modes_.shrink_to_fit();
}

Vibration::Vibration() {}
Vibration::Vibration(const std::vector<std::vector<uint16_t>> & _phonons)
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

const std::vector<std::vector<uint16_t>> & Vibration::phonons() const {return phonons_;}
const uint16_t & Vibration::excitation() const {return excitation_;}
const std::vector<std::pair<uint16_t, uint16_t>> & Vibration::excited_modes() const {return excited_modes_;}

void Vibration::pretty_print(std::ostream & stream) const {
    for (size_t i = 0; i < phonons_.size(); i++) {
        stream << "Irreducible " << i + 1 << ":\n";
        for (const auto & phonon : phonons_[i]) stream << phonon << ' ';
        stream << '\n';
    }
}

// Given frequency `w` and normal coordiante `Q`, return harmonic oscillator wave function value
double Vibration::value(const std::vector<std::vector<double>> & ws, const std::vector<std::vector<double>> & Qs) const {
    if (ws.size() != Qs.size() || Qs.size() != phonons_.size()) throw std::invalid_argument(
    "vibron::Wfn::operator(): inconsistent number of irreducible representations");
    for (size_t i = 0; i < ws.size(); i++)
    if (ws[i].size() != Qs[i].size() || Qs[i].size() != phonons_[i].size()) throw std::invalid_argument(
    "vibron::Wfn::operator(): inconsistent number of normal modes");
    double value = 1.0;
    size_t intdim = 0, total_phonon = 0;
    for (size_t irred = 0; irred < ws.size(); irred++) {
        intdim += ws[irred].size();
        for (size_t mode = 0; mode < ws[irred].size(); mode++) {
            const double & w = ws[irred][mode];
            const double & Q = Qs[irred][mode];
            const size_t & n = phonons_[irred][mode];
            total_phonon += n;
            value *= pow(w, 0.25) / sqrt(CL::math::dFactorial(n)) * exp(-0.5 * w * Q * Q) * std::hermite(n, sqrt(w) * Q);
        }
    }
    return value / pow(M_PI, (double)intdim / 4.0) / pow(2.0, (double)total_phonon / 2.0);
}

} // namespace vibron