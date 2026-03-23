# MATCH Run Setup — March 2026

## Goal

Run MATCH for March 2026 (spin-up: Feb 2026, target: Mar 2026).

## Completed

1. **Data files copied** from `~/Downloads/MATCH_Ed4/data/` to `~/EarthSystem/MATCH/data/`:
   - aer.nc, carbon_emissions.nc, DMS_emissions.nc, dst_bnd.nc, ginouxbasin,
     ncep_grid.nc, obs.source, oxid_3d.nc, SAGE_1020_OD_Filled.dat, SOx_emissions.nc

2. **Run script converted** from csh to bash:
   - Original: `control/run_sarb4_hourly_202602.csh`
   - New: `control/run_sarb4_hourly_202602.sh`
   - `doAODasm` set to `.false.` (no MODIS assimilation data available)

3. **Run directory created** with symlinks:
   - `~/Data/MATCH/2026/sarb4_hourly_202602/` → links to `data/*` and `grib/*`

4. **GRIB listings** (`grib/grbsanl.list`, `grib/grb2d.list`) already have Jan + Feb 2026 entries.
   GRIB files present in `~/Data/NCAR_A/`:
   - A27307-202601.grbsanl, A27326-202602.grbsanl
   - A27318-202601.grb2d, A27337-202602.grb2d

## Problem — Missing GRIB Data

Original plan was Feb 2026 output (spin-up from Jan 1). Failed because model needs
Dec 31, 2025 18Z meteorology for interpolation at startup — no Dec 2025 GRIB files exist.

Switched target to **March 2026** (spin-up Feb, ICDATE=20260201). This requires:
- Jan 2026 GRIB (have it — needed for Jan 31 18Z boundary)
- Feb 2026 GRIB (have it)
- **Mar 2026 GRIB (MISSING)**

## TODO — Before Running

### 1. Download March 2026 GRIB data from GDEX

Source: https://gdex.ucar.edu/datasets/d090000/dataaccess

Need both `grbsanl` and `grb2d` for March 2026. The GDEX site has data through
2026-03-01, so March data should be available. The web interface uses dynamic
JavaScript for downloads — could not extract direct URLs programmatically.

Existing GDEX archives are at `~/Data/GDEX/{grbsanl,grb2d}/2026/`.

**Alternative**: The CORe→GDEX translation pipeline at `~/EarthSystem/CDAS-CORe/`
can generate the GRIB files from CORe source data, BUT CORe data for March 2026
only goes through the 18th (`~/Data/CORe/{pgb,flx}/2026/03/` — incomplete).

### 2. Update grib list files

Add March 2026 entries to both:
- `grib/grbsanl.list`: `AXXXXX-202603.grbsanl <filesize> 202603010000`
- `grib/grb2d.list`: `AXXXXX-202603.grb2d <filesize> 202603010000`

### 3. Copy GRIB files to NCAR_A

Copy or link the new March archives into `~/Data/NCAR_A/`.

### 4. Update run script for March 2026

Update `control/run_sarb4_hourly_202602.sh` (or create new `_202603.sh`):
- `CASE=sarb4_hourly_202603`
- `ICDATE=20260201`
- `count_inc=59` (28 Feb + 31 Mar = 59 days) — same as current
- Update `wrkdir` path

### 5. Create March run directory with symlinks

```bash
mkdir -p ~/Data/MATCH/2026/sarb4_hourly_202603
cd ~/Data/MATCH/2026/sarb4_hourly_202603
ln -s ~/EarthSystem/MATCH/data/* .
ln -s ~/EarthSystem/MATCH/grib/* .
```

### 6. Run

```bash
cd ~/EarthSystem/MATCH/control
./run_sarb4_hourly_202603.sh
```

## Reference — CDAS-CORe Translation Pipeline

If downloading from GDEX directly isn't feasible, the translation pipeline can
generate MATCH-compatible GRIB from CORe reanalysis:

```bash
cd ~/EarthSystem/CDAS-CORe
conda activate sarb
python run_translation.py 2026 3 --all
```

But this requires complete CORe source data for the month (currently only through Mar 18).

## CORe Test Run — Feb 2026

### Completed

1. **CORe run script** created: `control/run_sarb4_hourly_202602_core.sh`
   - Uses CORe-translated GRIB via CDAS-CORe pipeline
   - ICDATE=20260201, hourly output (NHTFRQ=-1, MFILT=1)
   - 12 tracers: SO2, SO4, DMS, H2O2, DSTQ01-04, OCPHO, BCPHO, OCPHI, BCPHI
   - Run completed: 672 hourly history files in `~/Data/MATCH/2026/sarb4_hourly_202602_core/`

2. **Built `ccm2nc` converter** (`utils/ccm2nc`):
   - Fixed missing C standard library headers for modern GCC (`ccm_util.c`, `crpack.c`,
     `parse.c`, `strops.c`, `hsum.c`)
   - Makefile updated to use `$(CONDA_PREFIX)` for netCDF path — requires `conda activate`
     before `make`

3. **Converted all output to netCDF**:
   - `control/convert_core_nc.sh` — batch conversion script (skips existing `.nc` files)
   - 672 files converted successfully (h0001.nc–h0672.nc, ~53 MB each)

4. **AOD plot script** created: `control/plot_aod_core.py`
   - NCAR brand styling (DAVINCI-MONET) matching `~/EarthSystem/CDAS-CORe/plot_diagnostics.py`
   - Robinson projection, YlOrBr colormap, horizontal colorbar, 300 DPI PNG + PDF
   - Supports `--ref-dir` for comparison runs with RdBu_r difference plots
   - 5 monthly-mean AOD maps: total (AEROD), sulfate (SO4OD), dust (DSTODXC),
     sea salt (SSLTOD), carbonaceous (OC+BC sum)

5. **Conda environment** spec: `control/environment.yml` (`match-post`)

### Results — CORe Feb 2026

- Sulfate, sea salt, and carbonaceous AOD look physically reasonable
- **Dust AOD is very high** — Saharan/Arabian plume saturates above 0.5
  - Suspect: CORe surface winds or soil moisture over dust source regions differ
    from CDAS, amplifying dust emissions
  - Needs validation against MODIS/VIIRS and AERONET for Feb 2026

### CDAS Jan 2026 Run (existing)

- `~/Data/MATCH/2026/sarb4_hourly_202602/` — **actually covers Jan 2026** (20260101–20260201)
  despite directory name
- 747 history files, last file (h0747) truncated (run didn't finish cleanly)
- 746 files converted to netCDF; AOD plots generated with correct "Jan 2026" labels
- **Cannot diff CORe vs CDAS** — different months (Feb vs Jan)
- Proper comparison requires a CDAS Feb run or CORe Jan run

## Reference — Key Paths

| Item | Path |
|------|------|
| MATCH repo | `~/EarthSystem/MATCH/` |
| CDAS-CORe pipeline | `~/EarthSystem/CDAS-CORe/` |
| MATCH binary | `~/EarthSystem/MATCH/build/match` |
| Ancillary data | `~/EarthSystem/MATCH/data/` |
| GRIB list files | `~/EarthSystem/MATCH/grib/` |
| NCAR_A GRIB files | `~/Data/NCAR_A/` |
| GDEX archives | `~/Data/GDEX/{grbsanl,grb2d}/2026/` |
| CORe source data | `~/Data/CORe/{pgb,flx}/2026/` |
| Run directories | `~/Data/MATCH/2026/` |
| MODIS assim data | `~/Data/modis_assim_combo/2026/` (empty, not needed if doAODasm=.false.) |
| ccm2nc converter | `~/EarthSystem/MATCH/utils/ccm2nc` |
| Batch convert script | `~/EarthSystem/MATCH/control/convert_core_nc.sh` |
| AOD plot script | `~/EarthSystem/MATCH/control/plot_aod_core.py` |
| Conda env spec | `~/EarthSystem/MATCH/control/environment.yml` |
| CORe Feb plots | `~/Data/MATCH/2026/sarb4_hourly_202602_core/plots/` |
| CDAS Jan plots | `~/Data/MATCH/2026/sarb4_hourly_202602/plots/` |
