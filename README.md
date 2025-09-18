# HRS-f Enhanced: PhD-Level Harmonic Rayleigh Scattering Simulation

A professional-grade C++ framework for simulating harmonic electric field generation from dipole population

## 🎯 Project Overview

HRS-f Enhanced implements the full tensor formalism for hyperpolarizability-mediated nonlinear optical phenomena, including:
- Second harmonic generation (SHG)
- Electric-magnetic coupling effects (EME/MEE)
- Electric quadrupole contributions (QEE)
- Retardation and near-field effects

## 🏗️ Enhanced Project Structure

```
HRS-f/
├── src/                    # Source files
│   ├── main.cpp           # Legacy main program
│   ├── main_enhanced.cpp  # Enhanced main with CLI/config
│   ├── class.cpp          # Core physics implementation
│   └── config.cpp         # Configuration management
├── include/                # Header files
│   ├── class.h            # Core class declarations
│   └── config.h           # Configuration structures
├── tests/                  # Unit tests
│   ├── test_main.cpp      # Comprehensive test suite
│   └── catch2/            # Testing framework
├── validation/             # Analytical validation
│   ├── analytical/        # Benchmark implementations
│   └── reference_data/    # Reference datasets
├── config/                 # Configuration files
│   ├── experiments/       # Experiment configurations
│   └── sweeps/            # Parameter sweep definitions
├── build/                  # Build artifacts
├── output/                 # Simulation outputs
│   ├── data/              # Data files
│   ├── configs/           # Configuration snapshots
│   └── plots/             # Visualization outputs
├── docs/                   # Documentation
│   ├── THEORY.md          # Theoretical foundation
│   └── API.md             # API reference
├── profiling/              # Performance analysis
├── benchmarks/             # Benchmark results
└── external/               # External dependencies
```

## 🛠️ Building and Installation

### Prerequisites
- C++14 compatible compiler (clang++ or g++)
- Make build system
- OpenMP (optional, for parallelization)

### Quick Start
```bash
# Clone and build
git clone <repository>
cd HRS-f
make all

# Run with default configuration
make run

# Run high-precision simulation
make run-precision

# Run all tests
make test-all
```

### Build Targets
```bash
make all           # Build enhanced version (default)
make legacy        # Build original version
make debug         # Debug build with symbols
make release       # Optimized release build
make test          # Run unit tests
make validate      # Run validation tests
make benchmark     # Performance benchmarks
make docs          # Generate documentation
make clean         # Clean build files
```

## 📊 Usage Examples

### Basic Simulation
```bash
# Run with 1000 frames, 2.5° resolution
./build/HRS-f-enhanced --frames 1000 --dGamma 2.5 --radius 15

# Enable all physics contributions
./build/HRS-f-enhanced --enable-eme --enable-mee --enable-qee
```

### Configuration File
```bash
# Use predefined configuration
./build/HRS-f-enhanced --config config/experiments/high_precision.yaml

# Generate template configuration
./build/HRS-f-enhanced --generate-config my_config.yaml
```

### Validation and Testing
```bash
# Run validation tests
./build/HRS-f-enhanced --validation

# Uncertainty quantification with 5000 samples
./build/HRS-f-enhanced --uq 5000
```

## 📈 Performance Features

### Optimizations Implemented
- **Memory Management**: Pre-allocated arrays, minimal dynamic allocation
- **Computational**: Precomputed trigonometric functions, vectorized operations
- **Parallelization**: OpenMP support for frame and angle loops
- **I/O**: Efficient file operations with metadata tracking

### Benchmarking Results
```bash
# Run performance benchmarks
make benchmark

# Profile with detailed analysis
make profile

# Memory leak checking
make memcheck
```

## 🧪 Scientific Validation

### Analytical Benchmarks
1. **Single Dipole Rayleigh**: Validates against analytical sin²(θ) pattern
2. **Ring Symmetry**: Verifies rotational symmetry preservation
3. **Conservation Laws**: Checks energy and momentum conservation
4. **Multipole Contributions**: Validates EME/MEE/QEE implementations

### Validation Reports
```bash
# Run comprehensive validation
make validate

# View validation report
cat validation/validation_report.txt
```

## 📚 Theoretical Foundation

The simulation implements the full nonlinear optical response:

```
P_i^(2ω) = β_{ijk}^{EEE} E_j^(ω) E_k^(ω) 
         + β_{ijk}^{EME} E_j^(ω) B_k^(ω)
         + β_{ijk}^{MEE} B_j^(ω) E_k^(ω)
         + β_{ijkl}^{QEE} ∇_l E_j^(ω) E_k^(ω)
```

With proper:
- Frame transformations via Euler angles
- Phase matching and retardation effects
- Statistical averaging over orientations
- Far-field radiation pattern calculation

See [`docs/THEORY.md`](docs/THEORY.md) for complete theoretical derivation.

## 🔬 Research Applications

- **Chiral Plasmonics**: Helical nanostructure arrays
- **Metamaterial Design**: Nonlinear optical metamaterials  
- **Surface Science**: Interfacial nonlinear spectroscopy
- **Quantum Optics**: Collective dipole interactions
- **Materials Science**: Hyperpolarizability measurements

### Publication-Ready Features
- Rigorous error analysis and uncertainty quantification
- Comprehensive validation against analytical theory
- Performance benchmarking and scaling analysis
- Complete methodology documentation
- Reproducible computational protocols

## 📖 Documentation

- [`docs/THEORY.md`](docs/THEORY.md) - Complete theoretical foundation
- [`docs/API.md`](docs/API.md) - API reference and examples
- [`config/experiments/`](config/experiments/) - Example configurations
- [`validation/`](validation/) - Validation methodology



________

## 🏆 Academic Impact

- **Rigorous Validation**: Benchmarked against analytical theory
- **Professional Software Engineering**: Modern C++ practices, comprehensive testing
- **Scientific Methodology**: Uncertainty quantification, reproducible protocols
- **Performance Engineering**: Optimized algorithms, parallel computing
- **Complete Documentation**: Theory, implementation, and usage guides

**Research Level**: Graduate to PhD-level computational physics research

---

*HRS-f Enhanced: Professional-grade simulation framework for nonlinear optical phenomena*


