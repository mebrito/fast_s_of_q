# Fast S(q) Calculator

A high-performance C++ application for calculating the structure factor S(q) from particle trajectory data using the Debye formula. This tool is designed for analyzing molecular dynamics simulations.

## Overview

The structure factor S(q) is a fundamental quantity in condensed matter physics and materials science, providing insights into the spatial organization and correlations of particles in a system. This implementation uses the Debye formula to compute S(q) efficiently for systems with thousands of particles.

## Features

- **Fast Computation**: Optimized C++ implementation with OpenMP parallelization
- **VTF File Support**: Reads Visual Trajectory Format (.vtf) files containing particle positions
- **Physics-Based Q-range**: Wave vectors as multiples of fundamental q₀ = 2π/L for periodic boundary compatibility
- **Memory Efficient**: Streamlined algorithms optimized for large particle systems
- **Comprehensive Output**: Generates both S(q) data and system information files

## Requirements

- C++ compiler with C++23 support (GCC recommended)
- OpenMP support for parallel processing
- Make build system

## Installation

1. Clone or download the repository
2. Compile using the provided Makefile:

```bash
make
```

This will create the `fast_s_of_q` executable.

## Usage

### Basic Usage

```bash
./fast_s_of_q <vtf_file> <num_q> <q_min> [output_file] [info_file]
```

### Parameters

- `vtf_file`: Path to the .vtf trajectory file containing particle positions
- `num_q`: Number of q values to calculate (integer)
- `q_min`: Minimum q value threshold - used to determine starting multiple of q₀ = 2π/L
- `output_file`: (Optional) Output file for S(q) data (default: "s_of_q.txt")
- `info_file`: (Optional) System information file (default: "info.txt")

### Example

```bash
./fast_s_of_q trajectory.vtf 50 0.1 results.txt system_info.txt
```

This command:
- Reads particle positions from `trajectory.vtf`
- Calculates S(q) for 50 q values starting from the first multiple of q₀ ≥ 0.1
- Q values are generated as q₀ × (n_min, n_min+1, ..., n_min+49) where q₀ = 2π/box_length
- Outputs results to `results.txt` and system information to `system_info.txt`

## Input File Format

The program expects .vtf (Visual Trajectory Format) files with the following structure:

```
unitcell <box_x> <box_y> <box_z>
atom <id> [properties...]
bond <atom1> <atom2>
...
timestep <step>
<particle_id> <x> <y> <z>
<particle_id> <x> <y> <z>
...
```

### Key Requirements:
- Lines starting with 'u' contain unit cell information
- Lines starting with 'a' declare atoms
- Lines starting with 'b' declare bonds
- Lines starting with 't' indicate the start of position data
- Position lines contain particle ID followed by x, y, z coordinates

## Output Files

### S(q) Data File
Two-column format containing:
```
# q S(q)
1.256637 1.234567
1.884956 1.198765
2.513274 1.156432
...
```
*Note: Q values are multiples of q₀ = 2π/box_length*

### System Information File
Contains metadata about the calculation:
- Input file path
- Number of particles and bonds
- Box dimensions
- Number of q values and actual q-range (min/max)
- Fundamental wave vector q₀
- Calculation time

## Algorithm Details

The program implements the Debye formula for structure factor calculation:

S(q) = (1/N) ∑ᵢ ∑ⱼ sin(qrᵢⱼ)/(qrᵢⱼ)

Where:
- N is the number of particles
- q is the wave vector magnitude
- rᵢⱼ is the distance between particles i and j
- The sum runs over all particle pairs

### Q-Vector Generation

The program generates q values as integer multiples of the fundamental wave vector q₀ = 2π/L, where L is the box length:

- **Physical Basis**: This approach ensures compatibility with periodic boundary conditions
- **Starting Point**: The minimum multiplier n_min is determined from q_min: n_min = max(2, floor(q_min/q₀))
- **Sequence**: Q values are generated as q[i] = q₀ × (n_min + i) for i = 0, 1, ..., num_q-1
- **Advantages**: 
  - Ensures proper sampling of reciprocal space
  - Avoids artifacts from incompatible q vectors
  - Provides physically meaningful wave vectors for periodic systems

### Key Optimizations:
- **OpenMP Parallelization**: Utilizes multiple CPU cores for faster computation
- **Memory Optimization**: Efficient memory usage for large systems
- **Vectorized Operations**: Uses STL algorithms for optimal performance
- **Physics-Based Q-vectors**: Uses multiples of q₀ = 2π/L for consistency with periodic boundaries

## Performance

The implementation is optimized for:
- Large particle systems (tested with 10,000+ particles)
- Multiple q values (typically 50-100 points)
- Parallel execution on multi-core systems

Typical performance: ~1-10 seconds for 1000 particles with 50 q values on a modern CPU.

## Building from Source

### Prerequisites
- GCC with C++23 support
- OpenMP library
- Make

### Compilation Options
The Makefile uses the following flags:
- `-std=c++23`: C++23 standard for modern language features
- `-O3`: Maximum optimization for performance
- `-fopenmp`: OpenMP support for parallelization

### Manual Compilation
```bash
g++ -std=c++23 -O3 -fopenmp -o fast_s_of_q fast_s_of_q.cpp
```

## File Structure

```
├── fast_s_of_q.cpp    # Main program and VTF file parsing
├── s_of_q.hpp         # Core S(q) calculation algorithms
├── Makefile           # Build configuration
└── README.md          # This file
```

## Technical Details

### Core Functions
- `extractPartPositions()`: Parses VTF files and extracts particle coordinates
- `calculate_s_of_q()`: Computes structure factor using parallelized Debye formula
- `calculate_qs()`: Generates linearly spaced q values as multiples of q₀ = 2π/L

### Supported Features
- Cubic simulation boxes
- Periodic boundary conditions
- Isotropic averaging
- Parallel computation with OpenMP

## License

This project is provided as-is for academic and research purposes.

## Contributing

Contributions are welcome! Please ensure:
- Code follows the existing style
- New features include appropriate documentation
- Performance optimizations are tested

## Support

For questions or issues, please refer to the source code documentation or create an issue in the repository.
