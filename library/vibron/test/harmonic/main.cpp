#include <iostream>

#include <vibron/harmonic.hpp>

int main() {
    vibron::Integrator integrator(0.5, 5, 4);

    std::cout << "Harmonic oscillator integrals:\n";
    std::cout << integrator({3, 3, 2}) - 15.5885 << ' '
              << integrator({5, 4, 5}) - 183 << '\n';
}