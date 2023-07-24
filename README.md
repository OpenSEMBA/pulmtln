# pulmtln
Per Unit Length Multiconductor Transmission Line Network solver. Features:
- $C$ and $L$ matrix extraction.
- Third order isoparametric elements. 
- Support for dielectric materials.
- Open boundary problems.
- Domain decomposition.


## License
- Scripts in ```gmshWrapper``` folder are under GPL-2 license because they link to gmsh program.
- ``` pulmtln ``` has the same license as MFEM.

## Compiling
Compilation needs vcpkg with the following packages:
- gtest
- nlohmann_json

Additionally needs:
- mfem

