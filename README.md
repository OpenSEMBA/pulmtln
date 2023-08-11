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

### Compiling in windows (cmake)
Configure and build presets are available. To configure
``` 
    cmake 
        -DCMAKE_FIND_USE_PACKAGE_REGISTRY=FALSE
        --preset "msbuild-vcpkg"
        -S <project folder>
        -B <build folder>
```
which requires the following environment variables to be set (using ```export```)
```
    VCPKG_ROOT=<vcpkg root folder>
    MFEM_PKG=<mfem folder including cmake config package>
```

