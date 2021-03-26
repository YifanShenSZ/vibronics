#include <cmath>

#include <vibron/harmonic.hpp>

namespace vibron {

BraOrderKet::BraOrderKet() {}
BraOrderKet::BraOrderKet(const size_t & _bra, const size_t & _order, const size_t & _ket)
: bra(_bra), order(_order), ket(_ket) {}
BraOrderKet::BraOrderKet(const std::initializer_list<size_t> & items) {
    std::vector<size_t> vector(items);
    bra   = vector[0];
    order = vector[1];
    ket   = vector[2];
}
BraOrderKet::~BraOrderKet() {}





Integrator::Integrator() {}
Integrator::Integrator(const double & _frequency, const size_t & _phonon, const size_t & _order)
: frequency_(_frequency), phonon_(_phonon), order_(_order) {
    if (_order < 1) return; // <m|Q^0|n> is trivial
    integrals_.resize(_order + 1);
    for (size_t a = 1; a <= _order; a++) integrals_[a].resize(_phonon + 1);
    // The recursion relation is:
    // <m|Q^a|n> = sqrt(n + 1) / sqrt(2w) * <m|Q^(a - 1)|n + 1>
    //           + sqrt(n    ) / sqrt(2w) * <m|Q^(a - 1)|n - 1>
    // Let coeff = 1 / sqrt(2w), so
    // <m|Q^a|n> = coeff * (
    //                 sqrt(n + 1) * <m|Q^(a - 1)|n + 1>
    //               + sqrt(n    ) * <m|Q^(a - 1)|n - 1>
    //             )
    double coeff = 1.0 / sqrt(2.0 * _frequency);
    // basic case: a = 1
    integrals_[1][0].resize(2);
    integrals_[1][0][0] = 0.0;
    integrals_[1][0][1] = coeff;
    for (size_t m = 1; m <= _phonon; m++) {
        integrals_[1][m].resize(m + 1 + 1);
        std::fill(integrals_[1][m].begin(), integrals_[1][m].end(), 0.0);
        integrals_[1][m][m - 1] = coeff * sqrt(m);
        integrals_[1][m][m + 1] = coeff * sqrt(m + 1);
    }
    // recursion
    for (size_t a = 2; a <= _order; a++)
    for (size_t m = 0; m <= _phonon; m++) {
        integrals_[a][m].resize(m + a + 1);
        std::fill(integrals_[a][m].begin(), integrals_[a][m].end(), 0.0);
        // |m + a> can only be annihilated
        size_t n = m + a;
        integrals_[a][m][n] = coeff * sqrt(n) * integrals_[a - 1][m][n - 1];
        // Non-trivial n starts from m - a
        if (m >= a) {
            // |m - a> can only be created
            size_t n = m - a;
            integrals_[a][m][n] = coeff * sqrt(n + 1) * integrals_[a - 1][m][n + 1];
            // normal ones
            for (size_t n = m - a + 2; n <= m + a - 2; n += 2)
            integrals_[a][m][n] = coeff * (
                  sqrt(n + 1) * integrals_[a - 1][m][n + 1]
                + sqrt(n    ) * integrals_[a - 1][m][n - 1]
            );
        }
        else {
            // Non-trivial n starts from 0
            if (a % 2 == m % 2) {
                // |0> can only be created
                integrals_[a][m][0] = coeff * integrals_[a - 1][m][1];
                // normal ones
                for (size_t n = 2; n <= m + a - 2; n += 2)
                integrals_[a][m][n] = coeff * (
                      sqrt(n + 1) * integrals_[a - 1][m][n + 1]
                    + sqrt(n    ) * integrals_[a - 1][m][n - 1]
                );
            }
            // Non-trivial n starts from 1
            else
            for (size_t n = 1; n <= m + a - 2; n += 2)
            integrals_[a][m][n] = coeff * (
                  sqrt(n + 1) * integrals_[a - 1][m][n + 1]
                + sqrt(n    ) * integrals_[a - 1][m][n - 1]
            );
        }
    }
}
Integrator::~Integrator() {}

const double & Integrator::frequency() const {return frequency_;}
const size_t & Integrator::phonon() const {return phonon_;}
const size_t & Integrator::order() const {return order_;}

double Integrator::operator()(const BraOrderKet & braorderket) const {
    if (abs(braorderket.bra - braorderket.ket) > braorderket.order) return 0.0;
    else {
        if (braorderket.order < 1) return 1.0;
        else return integrals_[braorderket.order][braorderket.bra][braorderket.ket];
    }
}

} // namespace vibron