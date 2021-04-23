# Compute eigenvalues and eigenvectors
Lanczos algorithm produces a Krylov subspace spaned on which the original Hermite matrix becomes tridiagonal, whose eigenvalues approaches those of the original. This program diagonalizes this tridiagonal matrix to obtain its eigenvalues and eigenvectors, and optionally the eigenvectors of the original Hermite matrix

If the eigenvectors of the original Hermite matrix is desired, the Krylov vectors are necessary, otherwise only the tridiagonal matrix is requested

The eigenvalues tell us the energy of the vibronic levels.

The tridiagonal eigenvectors tell us the spectral intensity (when the seed vector is generated from `seed`)