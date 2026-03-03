# Documentation

## Overview

This program is based on the Seismic Unix (SU) and FFTW software packages. 
Before compiling, make sure that both SU and FFTW are installed on your system.

## Compilation

The program can be compiled using the provided Makefile. There are two compiler-specific Makefiles:

- `Makefile.gcc` — for GCC compiler
- `Makefile.icc` — for Intel ICC compiler

Before running `make`, create a symbolic link named `Makefile` pointing to the desired compiler Makefile. For example:

```bash
# Use ICC
ln -sf Makefile.gcc Makefile

# Or use ICC
# ln -sf Makefile.icc Makefile
```

After the symbolic link is set, you can use the following commands:

- `make` — Compile the program.
- `make remake` — Clean all previous build files and recompile the program.
- `make clean` — Remove all object files and executable files.

## Running the Programs

```bash
# Change to the working directory
cd ./demo/marm-ii/

# Run forward modeling
./modeling.sh

# Run full-waveform inversion
./fwi.sh
```

## Requirements

- Seismic Unix (SU)
- FFTW (Fastest Fourier Transform in the West)
- GCC or ICC compiler
- OpenMP support (for parallelization)

## Notes

- Ensure that the environment variables `CWPROOT` (for SU) and `FFTWROOT` (for FFTW) are correctly set to the installation paths of SU and FFTW, respectively.
- If using GCC, make sure to use the GCC-compatible version of FFTW to avoid Intel-specific library dependencies.