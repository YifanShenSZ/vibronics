#include <Lanczos/Hamiltonian/monomial.hpp>

namespace Lanczos {

Monomial::Monomial() {}
Monomial::Monomial(const size_t & _irred, const size_t & _mode, const size_t & _order)
: irred(_irred), mode(_mode), order(_order) {}
Monomial::~Monomial() {}

bool Monomial::operator<(const Monomial & other) const {
    if (irred < other.irred) return true;
    else if (irred > other.irred) return false;
    else {
        if (mode < other.mode) return true;
        else if (mode > other.mode) return false;
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

} // namespace Lanczos