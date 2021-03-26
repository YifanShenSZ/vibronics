# Vibronics
Vibronic spectrum simulation package

In general, there are a several steps to simulate a vibronic spectrum:
1. Expand the diabatic electronic Hamiltonian in polynomials of normal modes (this is [*diabatz*](https://github.com/YifanShenSZ/diabatz)'s job)
2. Generate a vibrational basis with `vibration`
3. Generate a seed vector with `seed`
4. Run `Lanczos` to diagonalize the vibronic Hamiltonian
5. Process the tridiagonal matrix with `steqr`

For details and theories of each module, see into each directories

# Theory
The time-independent approach to simulate vibronic spectrum is to diagonalize the vibronic Hamiltonian

This program adopts diabatic electronic state and harmonic oscillator direct product vibrational basis, under which framework the vibronic Hamiltonian is written as
```
H^vibronic = H^harmonic + H^anharmonic
```
so that `H^harmonic` is diagonal with harmonic oscillator eigenvalues, while `H^anharmonic` can be very sparse in usual applications

For integral purpose, the dependency of `H^anharmonic` on molecular geometry is fitted as a polynomial of normal modes. In usual cases where the polynomial order is low (2nd order), `H^anharmonic` is almost a band matrix, making Lanczos algorithm extremely efficient

## Cite this work
> 1. [M. S. Schuurman, R. A. Young, D. Yarkony, Chem. Phys. 2008, 347, 57-64](https://doi.org/10.1016/j.chemphys.2007.09.040)
> 2. [M. S. Schuurman and D. R. Yarkony, J. Chem. Phys. 2008, 128, 044119](https://doi.org/10.1063/1.2826380)
> 3. [Y. Shen and D. R. Yarkony, J. Phys. Chem. Lett. 2020, 11, 17, 7245–7252](https://doi.org/10.1021/acs.jpclett.0c02199)