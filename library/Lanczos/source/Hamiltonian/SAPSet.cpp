#include <fstream>

#include <CppLibrary/utility.hpp>

#include <Lanczos/Hamiltonian/SAPSet.hpp>

namespace Lanczos {

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

// Return the symmetry adapted polynomial set value given normal coordinate Q
double SAPSet::operator()(const std::vector<std::vector<double>> & Q) const {
    double result = constant_;
    for (const auto & term : terms_) result += term.first * term.second(Q);
    return result;
}

} // namespace Lanczos