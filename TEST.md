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
