# MATCH Architecture

## Overview

MATCH (Model for Atmospheric Transport and CHemistry) is a global
offline chemical transport model driven by meteorological reanalysis
data. It is used within the NASA CERES project to produce global
aerosol datasets via assimilation of satellite observations (MODIS).

The model reads NCEP reanalysis fields, advects trace species on a
Gaussian grid, and applies aerosol source/sink processes including
emission, dry and wet deposition, and chemistry. MODIS aerosol
optical depth observations are assimilated using a 3-D optimal
interpolation scheme.

## Source Directory Map

```
src/              Main model core: dynamics, physics, I/O, grid
  slt/            Semi-Lagrangian transport advection
  spitfire/       Spitfire spectral advection (FFT-based)
dst/              Mineral dust aerosol module (mobilization, transport, deposition)
src_scyc/         Sulfur cycle chemistry (DMS, SO2, sulfate)
src_assim/        3-D aerosol optical depth assimilation
readers/          Meteorological input readers
  ncep/           NCEP GRIB format reader
  netcdf/         NetCDF format reader
utils/            Post-processing utilities (ccm2nc, hsum)
build/            Build system (Makefile, dependency generation)
control/          Run scripts (csh) and run generation tools (Python)
grib/             GRIB field listing tables for NCEP reanalysis
reform/           SAGE aerosol data reformatting tools
```

## Preprocessor Flags

The build is configured via CPP defines that enable aerosol and
process modules. These are set in `build/Makefile` via `usrCPP_DEFS_AER`:

| Flag | Module | Purpose |
|------|--------|---------|
| `-DDST` | `dst/` | Mineral dust aerosol |
| `-DCAER` | `src/caer.F90`, `src/caerbnd.F90` | Carbonaceous aerosol (black/organic carbon) |
| `-DSCYC` | `src_scyc/` | Sulfur cycle (DMS, SO2, sulfate) |
| `-DWETDEP` | `src/wetdepdr.F` | Wet deposition |
| `-DDRYDEP` | `src/drydep.F`, `src/drydepdr.F` | Dry deposition |
| `-DASSIM` | `src_assim/` | AOD assimilation |
| `-DLINUX` | (always set) | Platform conventions |

Precision is controlled separately via `usrCPP_DEFS_PRC`:
- `-DRPREC=12` and `-DREALSIZE=8` select double precision (r8).

## Major Functional Areas

### Dynamics and Transport
- `main.F` — Program entry, time stepping loop
- `dyninp.F` — Reads and interpolates meteorological input
- `mass.F`, `drymass.F` — Mass diagnostics and conservation
- `etadot.F` — Vertical velocity computation
- `div.F` — Divergence
- `slt/advect.F` — Semi-Lagrangian advection
- `spitfire/` — Spectral advection (FFT, Lagrangian transport)
- `advfix.F` — Advection mass fixer

### Physics and Cloud
- `aphys.F` — Column physics driver
- `cldfrc.F` — Cloud fraction
- `cldwat.F` — Prognostic cloud water
- `cloud.F` — Cloud processes
- `conv_zhang.F`, `conv_hack.F`, `conv_pjr.F` — Convection schemes
- `vdiff_mod.F` — Vertical diffusion

### Aerosol Modules
- `srcsnk.F` — Tracer source/sink driver
- `dst/` — Dust: mobilization (`dstmbl`), chemistry (`dstchm`), deposition (`dstdpsdry`, `dstdpswet`), radiation (`dstrad`)
- `src_scyc/` — Sulfur: emissions (`sulemis`), chemistry (`sulchem`), boundary conditions (`sulbnd`, `soxbnd`, `dmsbnd`, `acbnd`)
- `caer.F90`, `caerbnd.F90` — Carbonaceous aerosol emission and transport
- `getsslt.F` — Sea salt aerosol
- `optDepth.F` — Aerosol optical depth computation

### Emission Tuning

Aerosol emissions can be scaled per source region or latitude band via
namelist parameters declared in `src/control.com` and read in `main.F`.
All defaults are 1.0 (identity) so unconfigured runs are unchanged.

| Species | Knob(s) | Applied in | Region count |
|---------|---------|------------|--------------|
| Dust    | `dst_rgn_scale`, `vwc_scale`/`vwc_offset` (legacy) | `dst/dstmbl.F` | 7 source regions |
| Sea salt | `sslt_scale` (global), `sslt_bands`/`n_sslt_bands` (latitude bands) | `src/getsslt.F` | up to 10 user-defined bands |
| Sulfate (SOx) | `so2_rgn_scale` | `src_scyc/sulemis.F90` | 4 anthropogenic regions |
| DMS | `dms_rgn_scale` | `src_scyc/sulemis.F90` | 6 ocean basins |
| Organic carbon | `oc_rgn_scale` | `src/caer.F90` | 9 regions (BB + FF + Boreal) |
| Black carbon | `bc_rgn_scale` | `src/caer.F90` | 9 regions (same as OC) |

Region boundaries are hardcoded in the application files and resolved
per gridpoint; volcanic SO2 (separate emission stream in `volcemist`)
is not affected by `so2_rgn_scale`. See [`RUN.md`](../RUN.md) for the
namelist reference and region tables, and [`AOD_TUNE.md`](../AOD_TUNE.md)
for the per-month tuning workflow against assimilated AOD climatologies.

### Data Assimilation (`src_assim/`)
- `aod_assim_3d.F90` — Main assimilation driver
- `analysis_3d.F90` — 3-D analysis update
- `covariances_3d.F90` — Background error covariances
- `obs_operator.F90` — Observation operator (model-to-obs mapping)
- `sat_obs.F90` — Satellite observation ingestion
- `avg_kernel.F90` — Averaging kernels

### I/O and Utilities
- `ccmhist.F`, `histout.F` — History file output (CCM format)
- `mass2nc.F`, `ncUtil.F` — NetCDF output
- `rec_io.F` — Record I/O
- `restart.F` — Checkpoint/restart
- `control.F` — Namelist and run control
- `grid.F` — Gaussian grid setup
- `calendar.F` — Date/time management

### Meteorological Readers (`readers/`)
- `readers/ncep/dyninp.F` + `grib.c`, `gribr.c` — Read NCEP GRIB reanalysis
- `readers/netcdf/dyninp.F` — Read NetCDF reanalysis
- `rgi.c`, `rgip.F` — Grid interpolation
- `utim.c` — Time utilities

### Post-processing Utilities (`utils/`)
- `ccm2nc` — Convert CCM history files to netCDF
- `hsum` — Compute hash checksums of history files

## Include Files

Two types of include files are used throughout the Fortran source:

- **`.h` headers** — Preprocessor parameter definitions and derived-type declarations (e.g., `pmgrid.h` defines grid dimensions, `dst.h` defines dust preprocessor tokens)
- **`.com` files** — Fortran COMMON block declarations shared across compilation units (e.g., `control.com`, `grid.com`, `calendar.com`)

Both are pulled in via `#include` and are required for compilation.

## Grid

MATCH uses a global Gaussian grid. Grid dimensions are set at compile
time via `pmgrid.h` and `params.h` in the build directory. The standard
CERES configuration uses the NCEP T62 resolution.
