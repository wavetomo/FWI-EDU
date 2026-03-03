# FWI-EDU: Full Waveform Inversion Educational Tool

## Overview

FWI-EDU is an educational tool for Full Waveform Inversion (FWI) in seismic exploration, providing both C and MATLAB implementations. This project aims to offer an intuitive and extensible FWI solution for teaching and research purposes in the field of seismic exploration.

The project implements core algorithms for 2D acoustic full waveform inversion, including multi-scale frequency strategy, checkpoint wavefield reconstruction method (only c version), gradient computation, and other key technologies. It also supports forward modeling functionality for numerical simulation of seismic data and model inversion.

## Core Features

### 1. Full Waveform Inversion (FWI)
- **Multi-scale frequency strategy**: Start from low frequencies and gradually increase to high frequencies to improve inversion stability
- **Checkpoint method**: Used for wavefield reconstruction to reduce memory consumption
- **Gradient computation**: Calculate gradient based on the adjoint-state method
- **Step-length calculation**: Use line-search method to compute optimal step size
- **Model constraints**: Support model parameter range limits and water-layer fixing

### 2. Forward Modeling
- **2D acoustic wave equation solver**: Using finite-difference method
- **Boundary conditions**: Support free surface and MEAL (Modified EAL) absorbing boundary
- **Acquisition geometry**: Support reading acquisition configuration from JSON file or parameters
- **Data output**: Support generating seismic record data

### 3. Model Processing
- **Gardner density relationship**: Automatically compute density from velocity
- **Water-layer modeling**: Support setting water-layer parameters
- **Model smoothing**: Support smoothing of gradients and models
- **Model clipping**: Limit model parameters within reasonable ranges

### 4. Parallel Computing
- **C version**: Using OpenMP parallel computing
- **MATLAB version**: Support parallel computing


## Technical Implementation

### C Version
- **Dependencies**: Seismic Unix (SU), FFTW (Fastest Fourier Transform in the West)
- **Compilation**: Compile with Makefile, support GCC and Intel ICC compilers
- **Memory management**: Use memory allocation functions provided by SU
- **Data format**: Use Seismic Unix style binary files (32-bit float, column-major)

### MATLAB Version
- **Code structure**: Modular design, main function calls multiple helper functions
- **Visualization**: Provide visualization of models and results
- **Data format**: Support binary files and MAT files


## Installation Dependencies

### C Version
1. **Seismic Unix (SU)**: For seismic data processing and memory management
2. **FFTW**: For fast Fourier transform
3. **GCC or Intel ICC compiler**: For compiling C code
4. **OpenMP support**: For parallel computing

### MATLAB Version
1. **MATLAB R2019a or higher**: For running MATLAB code


## Directory Structure

```
FWI-EDU/
├── FWILab2d/
│   ├── FWILab2dCLangV1.2.6/        # C language implementation
│   │   ├── demo/                   # Example data and scripts
│   │   ├── include/                # Header files
│   │   ├── main/                   # Main function files
│   │   ├── src/                    # Source code files
│   │   ├── Makefile                # Compilation configuration
│   │   ├── Makefile.gcc            # GCC compiler configuration
│   │   ├── Makefile.icc            # Intel ICC compiler configuration
│   │   └── README.md               # Documentation
│   └── FWILab2dMATLAB.v1.2.2/      # MATLAB implementation
│       ├── model_geometry/         # Model geometry data
│       ├── acoustic2d_fwi.m        # Full waveform inversion main function
│       └── acoustic2d_modeling.m   # Forward modeling main function
│       └── README.md               # Documentation
├── MatlabAPP/                      # MATLAB applications
│   ├── Record2dView.mlappinstall   # 2D record data viewer
│   ├── SurveyDesign2D.mlappinstall # 2D survey design tool
│   └── README.md                   # Documentation
├── README.md                       # Project documentation
└── LICENSE                         # Open source license
```

## Usage

### Use C Version
#### 1. Compile C Version

```bash
cd FWILab2d/FWILab2dCLangV1.2.6
# Use GCC compiler
ln -sf Makefile.gcc Makefile
# Or use Intel ICC compiler
# ln -sf Makefile.icc Makefile
make
```

#### 2. Run Forward Modeling

```bash
./modeling.sh
```

#### 3. Run Full Waveform Inversion

```bash
./fwi.sh
```

### Use MATLAB Version

1. Open MATLAB software
2. Navigate to `FWILab2d/FWILab2dMATLAB.v1.2.2/` directory
3. Edit `acoustic2d_fwi.m` or `acoustic2d_modeling.m` files to set parameters
4. Run the script: type `acoustic2d_fwi` or `acoustic2d_modeling` in MATLAB command window

## Example Data

- **Marmousi II model**: Located in `FWILab2d/FWILab2dCLangV1.2.6/demo/marm-ii/` directory, including:
  - Velocity model file: `vp_00221_00601_12.5m.bin`
  - Density model file: `rho_00221_00601_12.5m.bin`
  - Source wavelet file: `ricker_6001_1ms_10Hz_delay0.15ms.bin`
  - Acquisition geometry file: `acquisition.json`
  - Run scripts: `fwi.sh` (full waveform inversion) and `modeling.sh` (forward modeling)

## License

This project is licensed under the [MIT License](LICENSE), copyright (c) 2026 WaveTomo.

## Contact

For questions or suggestions, please submit an Issue or contact the project maintainers.

## Citation

If you use this project in your research or teaching, please cite the following information:
- Project name: FWI-EDU: Full Waveform Inversion Educational Tool
- Version: v1.2.6 (C) / v1.2.2 (MATLAB)
- Copyright owner: WaveTomo
- Open source license: MIT License
