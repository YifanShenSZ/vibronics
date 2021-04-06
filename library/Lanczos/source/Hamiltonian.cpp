#include <fstream>

#include <CppLibrary/utility.hpp>

#include <Lanczos/Hamiltonian.hpp>

namespace Lanczos {

Monomial::Monomial() {}
Monomial::Monomial(const size_t & _irred, const size_t & _mode, const size_t & _order)
: irred(_irred), mode(_mode), order(_order) {}
Monomial::~Monomial() {}

bool Monomial::operator<(const Monomial & other) const {
    if (irred < other.irred) return true;
    else {
        if (mode < other.mode) return true;
        else {
            if (order < other.order) return true;
            else                     return false;
        }
    }
}
bool Monomial::operator==(const std::pair<size_t, size_t> & irred_mode) const {
    return irred == irred_mode.first && mode == irred_mode.second;
}
bool Monomial::operator!=(const std::pair<size_t, size_t> & irred_mode) const {
    return irred != irred_mode.first || mode != irred_mode.second;
}





SAP::SAP() {}
SAP::SAP(const std::string & line) {
    // Read the symmetry adapted polynomial in expanded form, e.g. x * x rather than x^2
    std::vector<std::string> strs = CL::utility::split(line);
    std::vector<std::pair<size_t, size_t>> coords_temp(std::stoul(strs[0]));
    for (size_t i = 0; i < coords_temp.size(); i++) {
        std::vector<std::string> irred_mode = CL::utility::split(strs[i + 1], ',');
        coords_temp[i].first  = std::stoul(irred_mode[0]) - 1;
        coords_temp[i].second = std::stoul(irred_mode[1]) - 1;
    }
    std::sort(coords_temp.begin(), coords_temp.end());
    // Count the unique monomials
    std::pair<size_t, size_t> coord_old(-1, -1);
    for (const auto & coord : coords_temp)
    if (coord != coord_old) {
        coords_.push_back(Monomial(coord.first, coord.second, 1));
        coord_old = coord;
    }
    else {
        coords_.back().order += 1;
    }
    coords_.shrink_to_fit();
}
SAP::~SAP() {}

size_t SAP::NMonomials() const {return coords_.size();}
const Monomial & SAP::operator[](const size_t & index) const {return coords_[index];}
std::vector<Monomial>::const_iterator SAP::begin() const noexcept {return coords_.begin();}
std::vector<Monomial>::const_iterator SAP::end() const noexcept {return coords_.end();}

// Check if the monomials match the given normal modes
bool SAP::operator==(const std::vector<std::pair<size_t, size_t>> & irred_modes) const {
    if (coords_.size() != irred_modes.size()) return false;
    auto coords_copy = coords_;
    std::sort(coords_copy.begin(), coords_copy.end());
    auto irred_modes_copy = irred_modes;
    std::sort(irred_modes_copy.begin(), irred_modes_copy.end());
    bool same = true;
    for (size_t i = 0; i < coords_.size(); i++)
    if (coords_copy[i] != irred_modes_copy[i]) {
        same = false;
        break;
    }
    return same;
}
// Check if the monomials include the given normal modes
bool SAP::operator>=(const std::vector<std::pair<size_t, size_t>> & irred_modes) const {
    if (coords_.size() < irred_modes.size()) return false;
    bool include = true;
    for (size_t i = 0; i < irred_modes.size(); i++) {
        size_t j;
        for (j = 0; j < coords_.size(); j++)
        if (coords_[j] == irred_modes[i]) break;
        if (j >= coords_.size()) {
            include = false;
            break;
        }
    }
    return include;
}





void SAPSet::construct_excitation() {
    max_excitation_ = 0;
    for (const auto & term : terms_)
    if (term.second.NMonomials() > max_excitation_)
    max_excitation_ = term.second.NMonomials();
    excitations_.clear();
    excitations_.resize(max_excitation_ + 1);
    for (const auto & term : terms_)
    excitations_[term.second.NMonomials()].push_back(& term);
    for (auto & excitation : excitations_) excitation.shrink_to_fit();
}

SAPSet::SAPSet() {}
SAPSet::SAPSet(const std::string & Hd_file) {
    std::ifstream ifs; ifs.open(Hd_file);
    while (true) {
        std::string line1, line2;
        std::getline(ifs, line1);
        if (! ifs.good()) break;
        std::getline(ifs, line2);
        if (! ifs.good()) throw CL::utility::file_error(Hd_file);
        auto strs = CL::utility::split(line1);
        if (std::stoul(strs[0]) == 0) constant_ = std::stod(line2);
        else terms_.push_back(std::pair<double, SAP>(std::stod(line2), SAP(line1)));
    }
    ifs.close();
    terms_.shrink_to_fit();
    this->construct_excitation();
}
SAPSet::~SAPSet() {}

const double & SAPSet::constant() const {return constant_;}

size_t SAPSet::size() const {return terms_.size();}
const std::pair<double, SAP> & SAPSet::operator[](const size_t & index) const {return terms_[index];}
std::vector<std::pair<double, SAP>>::const_iterator SAPSet::begin() const noexcept {return terms_.begin();}
std::vector<std::pair<double, SAP>>::const_iterator SAPSet::end() const noexcept {return terms_.end();}

const size_t & SAPSet::max_excitation() const {return max_excitation_;}
const std::vector<const std::pair<double, SAP> *> & SAPSet::excitation(const size_t & ex) const {return excitations_[ex];}





// Construct `max_excitation_`, `max_orders_`
void Hd::construct_max() {
    // Count the highest possible excitation difference can be coupled by this Hd matrix
    max_excitation_ = 0;
    for (size_t i = 0; i < Hd_.size(); i++)
    for (size_t j = i; j < Hd_.size(); j++)
    if (Hd_[i][j]->max_excitation() > max_excitation_)
    max_excitation_ = Hd_[i][j]->max_excitation();
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
const size_t & Hd::max_excitation() const {return max_excitation_;}
const size_t & Hd::max_order(const size_t & irred, const size_t & mode) const {return max_orders_[irred][mode];}
const size_t & Hd::max_order(const std::pair<size_t, size_t> & irred_mode) const {return max_orders_[irred_mode.first][irred_mode.second];}

const SAPSet * Hd::operator[](const std::pair<size_t, size_t> & indices) const {
    size_t row = std::min(indices.first, indices.second),
           col = std::max(indices.first, indices.second);
    if (col >= Hd_.size()) throw std::invalid_argument(
    "Lanczos::Hd::operator[]: index out of range");
    return Hd_[row][col];
}

} // namespace Lanczos