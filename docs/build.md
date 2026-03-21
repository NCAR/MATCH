# Building MATCH

## Prerequisites

- **Fortran compiler**: `gfortran` (GCC Fortran)
- **C compiler**: `gcc`
- **NetCDF libraries**: netCDF-C and netCDF-Fortran with HDF5 support
- **Make**: GNU Make

On macOS without admin access, compilers can be installed via conda:
```
conda install -c conda-forge gfortran_osx-arm64 gcc_osx-arm64
```

## Building the Model

### 1. Configure `build/Makefile`

Edit the top section of `build/Makefile` to set paths for your system:

```makefile
# Compiler paths
FC = gfortran
CC = gcc

# NetCDF/HDF5 include path
usrCPP_INCS = -I/path/to/netcdf/include

# NetCDF/HDF5 library path and link flags
nchdflib = -L/path/to/netcdf/lib -lnetcdff -lnetcdf \
           -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lsz -lz
```

Other Makefile variables you may need to adjust:

| Variable | Default | Purpose |
|----------|---------|---------|
| `precision_mode` | `r8` | Double precision |
| `usrFFLAGS` | `-O3 -ffixed-line-length-132` | Fortran compile flags |
| `usrLDFLAGS` | `--static -g` | Linker flags |
| `usrCPP_DEFS_AER` | `-DDST -DCAER -DSCYC -DWETDEP -DDRYDEP -DASSIM` | Aerosol module selection |

### 2. Configure `build/Filepath`

`Filepath` lists the source directories to search, one per line.
Update the paths to match your installation:

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
