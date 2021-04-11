#include <CppLibrary/utility.hpp>

#include <Lanczos/Hamiltonian.hpp>

namespace Lanczos {

// Construct `max_orders_`, `max_excitation_`, `excitations_`
void Hd::construct_max() {
    // Count the highest order of each normal mode
    max_orders_.resize(NModes_.size());
    for (size_t i = 0; i < NModes_.size(); i++) {
        max_orders_[i].resize(NModes_[i]);
        std::fill(max_orders_[i].begin(), max_orders_[i].end(), 0);
    }
    for (size_t i = 0; i < NStates_; i++)
    for (size_t j = i; j < NStates_; j++)
    for (const auto & coeff_sap : *Hd_[i][j])
    for (const auto & monomial : coeff_sap.second)
    if (monomial.order > max_orders_[monomial.irred][monomial.mode])
    max_orders_[monomial.irred][monomial.mode] = monomial.order;
    // Count the highest possible excitation difference can be coupled by this anharmonic diabatic electronic Hamiltonian
    max_excitation_ = 0;
    for (size_t i = 0; i < Hd_.size(); i++)
    for (size_t j = i; j < Hd_.size(); j++)
    if (Hd_[i][j]->max_excitation() > max_excitation_)
    max_excitation_ = Hd_[i][j]->max_excitation();
    // Count the possbile ways to couple different modes
    excitations_.resize(max_excitation_ + 1);
    for (size_t excitation = 3; excitation <= max_excitation_; excitation++) {
        std::vector<std::vector<std::pair<size_t, size_t>>> excitation_vector;
        for (size_t i = 0; i < Hd_.size(); i++)
        for (size_t j = i; j < Hd_.size(); j++)
        for (const auto & sap : Hd_[i][j]->excitation(excitation)) {
            std::vector<std::pair<size_t, size_t>> excited_modes;
            for (const auto & coord : sap->second) excited_modes.push_back({coord.irred, coord.mode});
            excitation_vector.push_back(excited_modes);
        }
        excitations_[excitation] = new std::set<std::vector<std::pair<size_t, size_t>>>(excitation_vector.begin(), excitation_vector.end());
    }
}

Hd::Hd() {}
Hd::Hd(const size_t & _NStates, const std::vector<size_t> & _NModes, const std::vector<std::string> & Hd_files)
: NStates_(_NStates), NModes_(_NModes) {
    if (Hd_files.size() != ( NStates_ + 1) *  NStates_ / 2) throw std::invalid_argument(
    "Lanczos::Hd: The number of input files must equal to the number of upper triangle elements");
    Hd_.resize( NStates_);
    size_t count = 0;
    for (size_t i = 0; i <  NStates_; i++)
    for (size_t j = i; j <  NStates_; j++) {
        Hd_[i][j] = new SAPSet(Hd_files[count]);
        count++;
    }
    this->construct_max();
}
Hd::~Hd() {}

const size_t & Hd::NStates() const {return NStates_;}
const size_t & Hd::max_order(const size_t & irred, const size_t & mode) const {return max_orders_[irred][mode];}
const size_t & Hd::max_order(const std::pair<size_t, size_t> & irred_mode) const {return max_orders_[irred_mode.first][irred_mode.second];}
const size_t & Hd::max_excitation() const {return max_excitation_;}
const std::set<std::vector<std::pair<size_t, size_t>>> * Hd::excitation(const size_t & index) const {return excitations_[index];}

void Hd::pretty_print(std::ostream & stream) const {
    stream << "Summary of this anharmonic diabatic electronic Hamiltonian\n";
    stream << "Highest possible excitation difference to couple = " << max_excitation_ << '\n';
    stream << "The highest order of each normal mode:\n";
    for (const auto & irred : max_orders_) {
        stream << "    ";
        for (const auto & max_order : irred) stream << max_order << ' ';
        stream << '\n';
    }
    stream << "The number of expansion terms\n";
    for (size_t i = 0; i < Hd_.size(); i++)
    for (size_t j = i; j < Hd_.size(); j++)
    stream << "    Hd" << i << ',' << j << ": " << Hd_[i][j]->size() << '\n';
}

const SAPSet * Hd::operator[](const std::pair<size_t, size_t> & indices) const {
    size_t row = std::min(indices.first, indices.second),
           col = std::max(indices.first, indices.second);
    if (col >= Hd_.size()) throw std::invalid_argument(
    "Lanczos::Hd::operator[]: index out of range");
    return Hd_[row][col];
}

} // namespace Lanczos