# Lanczos seed vector
Generate the seed vector for Lanczos iteration

The normal coordinate origin can be placed on (`initial`) or off (`final`) the initial-state equilibirum geometry

This library depends on *vibron* to define vibronic wave function

## Theory
Let the initial vibronic state be `d`, the final vibronic state be `D`, the transition dipole be `miu`. The spectral amplitude `A` is:
```
A = d . miu . D    (1)
```

By Lanczos diagonalization, if we set the seed vector to the unit vector along `d . miu_x`, then `A_x` emerges from the seed vector's contribution to `D` times `||d . miu_x||`. Similarly `A_y` and `A_z`.

What we really want is the spectral strength:
```
A^2 = A_x^2 + A_y^2 + A_z^2    (2)
```
So under special circumstances we can replace the vector `miu` in Eq. (1) with a scalar `||miu||`, subsequently `A` also becomes `||A||`:
1. `<initial electron|miu_x|final electron> = <initial electron|miu_y|final electron> = <initial electron|miu_z|final electron>` holds for all final electronic states
2. Under Franck-Condon approximation, if all final electronic states share a same `<initial electron|miu|final electron>` value
3. If each final electronic state has a different `<initial electron|miu|final electron>` component, e.g. state 1 has only `<initial electron|miu_x|1>` & state 2 only has `<initial electron|miu_y|2>` & state 3 only has `<initial electron|miu_z|3>`
4. Under Franck-Condon approximation, the criterion in 3 can be loosen to "each final electronic state has different `<initial electron|miu|final electron>` components", e.g. state 1 has only `<initial electron|miu_z|1>` while state 2 only has `<initial electron|miu_x|2>` and `<initial electron|miu_y|2>`