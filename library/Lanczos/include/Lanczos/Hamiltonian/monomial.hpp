#ifndef Lanczos_Hamiltonian_monomial_hpp
#define Lanczos_Hamiltonian_monomial_hpp

#include <cstddef>
#include <vector>
#include <utility>

namespace Lanczos {

struct Monomial {
    size_t irred, mode, order;

    Monomial();
    Monomial(const size_t & _irred, const size_t & _mode, const size_t & _order);
    ~Monomial();

    bool operator<(const Monomial & other) const;
    bool operator==(const std::pair<size_t, size_t> & irred_mode) const;
    bool operator!=(const std::pair<size_t, size_t> & irred_mode) const;

    // Return the monomial value given normal coordinate Q
    double operator()(const std::vector<std::vector<double>> & Q) const;
};

} // namespace Lanczos

#endif