#include <CppLibrary/utility.hpp>

#include <Lanczos/Hamiltonian/SAP.hpp>

namespace Lanczos {

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

// Return the symmetry adapted polynomial value given normal coordinate Q
double SAP::operator()(const std::vector<std::vector<double>> & Q) const {
    double result = 1.0;
    for (const Monomial & coord : coords_) result *= coord(Q);
    return result;
}

} // namespace Lanczos