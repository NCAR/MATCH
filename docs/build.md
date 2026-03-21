# Building MATCH

## Prerequisites

- **Fortran compiler**: `gfortran` (GCC Fortran)
- **C compiler**: `gcc`
- **NetCDF libraries**: netCDF-C and netCDF-Fortran with HDF5 support
- **Make**: GNU Make

## macOS (Apple Silicon) Setup

On macOS without admin/Homebrew access, compilers and libraries can
be installed entirely in user space via conda-forge:

```bash
# Install compilers
conda install -c conda-forge gfortran_osx-arm64 gcc_osx-arm64

# Install netCDF and HDF5
conda install -c conda-forge netcdf-fortran hdf5 libnetcdf
```

The conda compilers use target-prefixed names
(e.g. `arm64-apple-darwin20.0.0-gfortran`). Create short symlinks
in `~/bin` for convenience:

```bash
mkdir -p ~/bin
ln -sf ~/miniconda3/bin/arm64-apple-darwin20.0.0-gfortran ~/bin/gfortran
ln -sf ~/miniconda3/bin/arm64-apple-darwin20.0.0-gcc ~/bin/gcc
```

Ensure `~/bin` is on your `PATH`.

### macOS-specific Makefile settings

The Makefile shipped in `build/` is configured for macOS with
conda-forge toolchains. Key differences from a Linux build:

- **No `--static` linker flag** — macOS does not support fully
  static linking. Use `usrLDFLAGS = -g` instead.
- **GCC 15 compatibility flags** — Modern gfortran requires
  `-fallow-invalid-boz` (for octal BOZ literals in legacy code) and
  `-fallow-argument-mismatch` (for legacy Fortran calling conventions).
  These are included in `usrFFLAGS`.
- **`CONDA_PREFIX`** — The Makefile uses this variable to locate
  headers and libraries under `~/miniconda3`.

## Building the Model

### 1. Configure `build/Makefile`

Edit the top section of `build/Makefile` to set paths for your system:

```makefile
# Conda environment (adjust if not using default miniconda3 location)
CONDA_PREFIX = $(HOME)/miniconda3

# Compiler paths (use absolute paths or ensure ~/bin is on PATH)
FC = /Users/yourusername/bin/gfortran
CC = /Users/yourusername/bin/gcc

# NetCDF/HDF5 include and library paths
usrCPP_INCS = -I$(CONDA_PREFIX)/include
nchdflib = -L$(CONDA_PREFIX)/lib -lnetcdff -lnetcdf \
           -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lsz -lz
```

Other Makefile variables you may need to adjust:

| Variable | Default | Purpose |
|----------|---------|---------|
| `precision_mode` | `r8` | Double precision |
| `usrFFLAGS` | `-O3 -ffixed-line-length-132 -fallow-invalid-boz -fallow-argument-mismatch` | Fortran compile flags |
| `usrLDFLAGS` | `-g` | Linker flags (no `--static` on macOS) |
| `usrCPP_DEFS_AER` | `-DDST -DCAER -DSCYC -DWETDEP -DDRYDEP -DASSIM` | Aerosol module selection |

### 2. Configure `build/Filepath`

`Filepath` lists the source directories to search, one per line.
Use relative paths from the `build/` directory:

```
../readers/ncep
../readers
../src/spitfire
../dst
../src_scyc
../src_assim
../src
```

### 3. Build

```
cd build
make
```

This compiles all source listed in `Srcfiles`, placing object files
in `build/OBJ/`, and links the `match` executable in `build/`.

To rebuild dependencies after adding or removing source files:
```
make Depends
```

To clean and rebuild:
```
make clean
make
```

### Build System Files

| File | Purpose |
|------|---------|
| `Makefile` | Compiler flags, paths, pattern rules |
| `Filepath` | Directories to search for source (sets VPATH) |
| `Srcfiles` | List of all source files to compile |
| `Depends` | Auto-generated dependency rules (`#include` and `use` tracking) |
| `makdep` / `makdep_usrbinperl` | Perl script that scans source to generate `Depends` |
| `params.h` | Compile-time grid parameters |

### Special Compilation Rules

Some files require modified optimization or flags:
- `fft.F` — Cannot use bounds checking
- `advect.F`, `advfix.F` — Use `specFFLAGS2` (separate optimization level)

These are handled by explicit rules in the Makefile.

## Building the Utilities

The `utils/` directory builds two post-processing programs:

- **`ccm2nc`** — Converts CCM-format history files to netCDF
- **`hsum`** — Computes hash checksums of history files

```
cd utils
```

Edit `Makefile` to set the NetCDF path:
```makefile
CC = gcc
NETCDF_DIR = /path/to/netcdf
```

Then:
```
make
```

## Choosing a Reader

MATCH supports two meteorological input formats. The reader is
selected by which `dyninp.F` appears in the `Filepath` search order:

- **NCEP GRIB** (`readers/ncep/`) — reads GRIB-format reanalysis files via `grib.c`/`gribr.c`
- **NetCDF** (`readers/netcdf/`) — reads NetCDF-format reanalysis

The default `Filepath` lists `readers/ncep` first, selecting the GRIB reader.
