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

    auto E = energy - energy[0];

    double max_intensity = *std::max_element(intensity.begin(), intensity.end());
    auto normalized_intensity = intensity / max_intensity;

    std::cout << "The energy levels and spectral intensities can be found in spectrum.txt\n";

    std::ofstream ofs; ofs.open("spectrum.txt");
    ofs << "Ground state energy = " << std::fixed << std::setprecision(2) << energy[0] << " cm^-1\n"
        << "Max intensity = " << std::fixed << std::setprecision(6) << max_intensity << '\n'
        << "ΔE / cm^-1    convergence    amplitude    intensity    normalized intensity\n";
    for (size_t i = 0; i < N; i++) if (intensity[i] > threshold)
    ofs << std::fixed << std::setw(10) << std::setprecision(2) << E          [i] << "    "
        << std::fixed << std::setw(11) << std::setprecision(6) << convergence[i] << "    "
        << std::fixed << std::setw( 9) << std::setprecision(6) << amplitude  [i] << "    "
        << std::fixed << std::setw( 9) << std::setprecision(6) << intensity  [i] << "    "
        << std::fixed << std::setw( 9) << std::setprecision(6) << normalized_intensity[i] << '\n';
    ofs.close();
}