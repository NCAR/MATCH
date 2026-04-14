# AOD Tuning Plan

## Goal

Tune MATCH CORe aerosol emission parameters per calendar month against
AOD climatologies from MATCH assimilation runs, targeting two species:

- **Dust** — per-region emission multiplier (`dst_rgn_scale`, 7 regions).
  An older global rescale of soil moisture (`vwc_scale`, `vwc_offset`)
  remains available but is superseded by the regional approach because
  a single linear rescale cannot match multiple source regions
  simultaneously.
- **Sea salt** — emission flux scaling (`sslt_scale`)

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

### Step 3: Tune dust per region

With sea salt scaling fixed, run MATCH and compare `AEROD` against the
climatology in each named dust source region.  Fit one multiplier per
region per calendar month via `dst_rgn_scale(7)`.

Regions (boundaries hardcoded in `dst/dstmbl.F`):

| Index | Region | Lat | Lon |
|-------|--------|-----|-----|
| 1 | Sahara/Sahel | 10-30N | 20W-30E |
| 2 | Arabia | 15-30N | 30-55E |
| 3 | Central Asia | 30-50N | 55-90E |
| 4 | Gobi/Taklamakan | 35-50N | 90-115E |
| 5 | Australia | 35-15S | 115-145E |
| 6 | SW N.America | 25-40N | 115W-100W |
| 7 | Patagonia | 55-35S | 75W-60W |

Iteration heuristic: per region, compute a target dust AOD as
`VIIRS_total - (baseline_AEROD - baseline_DST)`, i.e. VIIRS minus the
local non-dust background.  Initial scale ≈ `target_dust / baseline_dust`.
Short 10-day runs (5-day Feb spin-up + 5-day March sample) suffice for
iteration; dust lifetime (~1–5 d) is short relative to the sample.

Caveat: regions where baseline emission is near zero (e.g. SW N.America
with box as defined) have no leverage via a multiplier — the box would
need relocation to actual source zones.

The older global approach (`vwc_scale`, `vwc_offset`) remains wired in
for backwards compatibility but is superseded: a single linear rescale
cannot simultaneously suppress hyperarid Sahara/Arabia emission and
preserve moderate Gobi/Australia emission.

### Step 4: Validate

Compare tuned MATCH output against:
- Held-out satellite data (train on MODIS, validate on VIIRS)
- CDAS baseline runs
- Per-species climatologies (when available)

## Namelist Integration

All parameters are in namelist `/nlist/` with identity defaults.  Example:

```fortran
sslt_scale    = 1.2    ! sea salt emission flux multiplier (per month)
dst_rgn_scale = 0.28, 0.40, 1.2, 4.0, 3.0, 1.0, 0.5
                       ! Sahara, Arabia, C.Asia, Gobi, Australia, SW-NAm, Patagonia
! Legacy (superseded by dst_rgn_scale):
! vwc_scale  = 1.0     ! dust soil moisture rescaling slope
! vwc_offset = 0.0     ! dust soil moisture rescaling offset
```

- `dst_rgn_scale` — applied in `dst/dstmbl.F` to `flx_mss_hrz_slt_ttl`
  alongside `lnd_frc_mbl` and `flx_mss_fdg_fct`, per-gridpoint based on
  hardcoded region boxes.
- `vwc_scale` / `vwc_offset` — applied in `dst/dstmbl.F` after
  `vwc_sfc_get()`, clamped to `[0.10, 0.43]` when engaged.
- `sslt_scale` — applied in `src/getsslt.F` after concentration and unit
  conversion.

## Output

A table of 12 monthly coefficient sets.  Each row holds one
`dst_rgn_scale` array (7 values) plus `sslt_scale`:

| Month | Sahara | Arabia | C.Asia | Gobi | Aust | SW-NAm | Patag | sslt | Notes |
|-------|--------|--------|--------|------|------|--------|-------|------|-------|
| Jan   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Feb   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Mar   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | in progress |
| Apr   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| May   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Jun   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Jul   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Aug   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Sep   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Oct   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Nov   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Dec   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
