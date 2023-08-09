# pulmtln
[![Tests](https://github.com/lmdiazangulo/pulmtln/actions/workflows/builds-and-tests.yml/badge.svg)](https://github.com/lmdiazangulo/pulmtln/actions/workflows/builds-and-tests.yml)

Per Unit Length Multiconductor Transmission Line Network solver. Features:
- $C$ and $L$ matrix extraction.
- Third order isoparametric elements. 
- Support for dielectric materials.
- Open boundary problems.
- Domain decomposition.


## License
- ``` pulmtln ``` has the same license as MFEM.

## Compiling
Compilation needs vcpkg with the following packages:
- gtest
- nlohmann_json

Additionally needs:
- mfem

