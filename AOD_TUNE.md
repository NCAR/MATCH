# AOD Tuning Plan

## Goal

Tune MATCH CORe aerosol emission parameters per calendar month against
AOD climatologies from MATCH assimilation runs, targeting two species:

- **Dust** — soil moisture rescaling (`vwc_scale`, `vwc_offset`)
- **Sea salt** — emission flux scaling (`sslt_scale`, TBD)

## Problem

CORe SOILW uses 18-bit packing (range 0.02-1.0) vs CDAS 9-bit (range
0.10-0.43).  A single linear rescaling cannot match CDAS dust across
months — Jan 2026 CORe dust is higher than CDAS, Feb 2026 is substantially
lower.  Sea salt AOD is systematically higher with CORe meteorology,
likely due to wind speed differences.  Rather than matching CDAS (which
itself has biases), we tune directly against assimilated AOD climatologies.

## Reference Climatologies

Monthly-mean total extinction AOD (`TOTEXTTAU`) on the MATCH T62 grid
(192 lon x 94 lat), produced by MATCH runs that assimilated MODIS or VIIRS
column AOD.  Located in `~/Data/MATCH_climo/`:

| File | Source | Record |
|------|--------|--------|
| `MATCH_Terra-Aqua-MODIS_Ed4.TOTEXTTAU_L192x94.nc4` | MATCH + MODIS Ed4 assimilation | 295 months |
| `MATCH_NOAA20-VIIRS_Ed1C.TOTEXTTAU_L192x94.nc4` | MATCH + VIIRS Ed1C assimilation | 83 months (May 2018 - Mar 2025) |

Both files have a `time` dimension (no coordinate values) concatenated via
`ncrcat` from individual monthly files.  The time-to-month mapping is
encoded in the `history` global attribute.

These are total AOD (all species combined).  Per-species climatologies will
be available later.

## Approach

### Step 1: Compute monthly climatologies

Average the assimilated TOTEXTTAU into 12 calendar-month means.  MODIS
provides the longer baseline; VIIRS serves as an independent check.

### Step 2: Tune sea salt first (clean-ocean regions)

Over remote ocean regions far from dust, pollution, and biomass burning
(Southern Ocean, central South Pacific, etc.), total AOD is dominated by
sea salt.  Compare free-running MATCH `AEROD` against the climatology in
these clean-ocean boxes to fit `sslt_scale` per calendar month.

Candidate clean-ocean regions:
- Southern Ocean: 50-65S, all longitudes
- Central South Pacific: 30-50S, 150W-90W
- Remote South Atlantic: 30-50S, 30W-10E

### Step 3: Tune dust (dust-dominated regions)

With sea salt scaling fixed, run MATCH and compare `AEROD` against the
climatology in dust-dominated regions to fit `vwc_scale` and `vwc_offset`
per calendar month.

Dust-dominated regions:
- Sahara/Sahel: 10-30N, 20W-30E
- Atlantic dust corridor: 5-25N, 60W-20W
- Arabian Peninsula: 15-30N, 35-60E
- Gobi/Taklamakan: 35-45N, 75-110E

### Step 4: Validate

Compare tuned MATCH output against:
- Held-out satellite data (train on MODIS, validate on VIIRS)
- CDAS baseline runs
- Per-species climatologies (when available)

## Namelist Integration

### Existing

```fortran
vwc_scale  = 0.72    ! dust soil moisture rescaling slope (per month)
vwc_offset = 0.09    ! dust soil moisture rescaling offset (per month)
```

### Needed

```fortran
sslt_scale = 1.0     ! sea salt emission flux multiplier (per month)
```

`sslt_scale` needs to be added to the namelist and applied in the sea salt
emission routine.  Default `1.0` = no change from current behavior.

## Output

A table of 12 monthly coefficient sets:

| Month | vwc_scale | vwc_offset | sslt_scale | Notes |
|-------|-----------|------------|------------|-------|
| Jan   | TBD | TBD | TBD | |
| Feb   | TBD | TBD | TBD | |
| Mar   | TBD | TBD | TBD | |
| Apr   | TBD | TBD | TBD | |
| May   | TBD | TBD | TBD | |
| Jun   | TBD | TBD | TBD | |
| Jul   | TBD | TBD | TBD | |
| Aug   | TBD | TBD | TBD | |
| Sep   | TBD | TBD | TBD | |
| Oct   | TBD | TBD | TBD | |
| Nov   | TBD | TBD | TBD | |
| Dec   | TBD | TBD | TBD | |
