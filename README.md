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
Compilation needs vcpkg with the packages stated in the ```vcpkg.json``` manifest. 

Additionally needs:
- mfem (with the version pointed by the external/mfem-geg submodule)

### Compiling in windows (cmake)

#### Manually (Windows/Linux)
Compile mfem in external/mfem-geg
```
    cmake -S external/mfem-geg -B mfem-build/rls
    cmake --build mfem-build/rls  --config Release
```
Launch cmake in root.
```
    cmake -S . -B pulmtln-build/rls -Dmfem_DIR=mfem-build/rls
    cmake --build pulmtln-build/rls --config Release
```

#### Using presets
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
    MFEM_PACKAGE=<mfem folder including cmake config package>
```
Using ```CMAKE_FIND_USE_PACKAGE_REGISTRY=FALSE``` warranties that no previously used package is used for compilation if any of the needed paths is not found (a questionable Windows _feature_). 
