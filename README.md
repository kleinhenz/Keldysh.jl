# Keldysh.jl

[![Build Status](https://github.com/kleinhenz/Keldysh.jl/workflows/CI/badge.svg?branch=master)](https://github.com/kleinhenz/Keldysh.jl/actions)

`Keldysh.jl` provides a set of tools for working with non-equilibrium Keldysh Green's functions.
It contains types to represent contours, grids defined on these contours, and two-time Green's functions defined on these grids.
Additionally, it provides functions for generating Green's functions, performing integration on a contour and hdf5 serialization.


Credit to Andrey Antipov and Igor Krivenko for designing a first version of the abstractions implemented here.

## Usage
[anderson_nca.jl](examples/anderson_nca.jl) implements a NCA solver for the anderson impurity model using `Keldysh.jl`.
