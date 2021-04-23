#include <iomanip>
#include <fstream>

#include <CppLibrary/utility.hpp>
#include <CppLibrary/linalg.hpp>

void output_spectrum(const std::vector<double> & energy, const CL::utility::matrix<double> & eigvec,
const double & last_beta, const double & threshold) {
    size_t N = energy.size();
    
    std::vector<double> convergence(N);
    for (size_t i = 0; i < N; i++) convergence[i] = std::abs(last_beta * eigvec[i].back());

    std::vector<double> amplitude(N), intensity(N);
    for (size_t i = 0; i < N; i++) {
        amplitude[i] = eigvec[i][0];
        intensity [i] = amplitude[i] * amplitude[i];
    }

    auto Ewavenumber = (energy - energy[0]) / 4.556335830019422e-6;

    double max_intensity = *std::max_element(intensity.begin(), intensity.end());
    auto normalized_intensity = intensity / max_intensity;

    std::cout << "The energy levels and spectral intensities can be found in spectrum.txt\n";

    std::ofstream ofs; ofs.open("spectrum.txt");
    ofs << "Ground state energy = " << std::fixed << std::setprecision(2) << energy[0] / 4.556335830019422e-6 << " cm^-1\n"
        << "Max intensity = " << std::fixed << std::setprecision(6) << max_intensity << '\n'
        << "Number    ΔE / cm^-1    convergence    amplitude    intensity    normalized    E / Hartree\n";
    for (size_t i = 0; i < N; i++) if (intensity[i] > threshold)
    ofs << std::setw(6) << i + 1
        << std::setw(14) << std::fixed << std::setprecision(2) << Ewavenumber[i]
        << std::setw(15) << std::scientific << std::setprecision(5) << convergence[i]
        << std::setw(13) << std::fixed << std::setprecision(7) << amplitude[i]
        << std::setw(13) << std::fixed << std::setprecision(7) << intensity[i]
        << std::setw(14) << std::fixed << std::setprecision(8) << normalized_intensity[i]
        << std::setw(25) << std::scientific << std::setprecision(15) << energy[i]
        << '\n';
    ofs.close();
}