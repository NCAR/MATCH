# AOD Tuning Plan

## Goal

Tune MATCH CORe aerosol emission parameters per calendar month against
AOD climatologies from MATCH assimilation runs, targeting four species:

- **Dust** — per-region emission multiplier (`dst_rgn_scale`, 7 regions).
  An older global rescale of soil moisture (`vwc_scale`, `vwc_offset`)
  remains available but is superseded by the regional approach because
  a single linear rescale cannot match multiple source regions
  simultaneously.
- **Sea salt** — global emission flux multiplier (`sslt_scale`) plus
  optional per-latitude-band tuning (`sslt_bands`, `n_sslt_bands`) to
  address N/S asymmetry in wind-driven emission (strong Southern Ocean
  bias is the largest signal against MODIS/VIIRS).
- **Sulfate** — per-region SOx emission multiplier (`so2_rgn_scale`,
  4 anthropogenic source regions) and per-basin DMS multiplier
  (`dms_rgn_scale`, 6 ocean basins). Volcanic SO2 is unaffected.
- **Carbon** — per-region OC and BC emission multipliers
  (`oc_rgn_scale`, `bc_rgn_scale`, 9 regions each — same boxes for both,
  spanning major BB and FF source areas plus a Boreal NH catch-all).

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
these clean-ocean boxes to fit the sea salt scaling per calendar month.

Two knobs are available and compose multiplicatively:
- `sslt_scale` — a single global multiplier.
- `sslt_bands` + `n_sslt_bands` — per-latitude-band multipliers.  Each
  active band is a `(lat_min_deg, lat_max_deg, scale)` triple.  Useful
  when the bias has N/S asymmetry (typical: Southern Ocean stronger
  than Northern).

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

### Step 4: Tune sulfate and carbon per region

With dust and sea salt fixed, tune sulfate (SOx) and carbon (OC, BC)
against per-species climatologies (or AAOD for BC).  Each species has
its own predefined regional boxes; outside every box, the multiplier
is 1.0.  All scales default to 1.0 (identity).

**Sulfate — `so2_rgn_scale(4)`** (anthropogenic SOx; both 0m and >100m
emissions are scaled together; volcanic SO2 is unaffected):

| Index | Region | Lat | Lon |
|-------|--------|-----|-----|
| 1 | E. Asia    | 20-50N | 100-145E |
| 2 | S. Asia    |  5-35N |  65-95E  |
| 3 | Europe     | 35-65N |  10W-50E |
| 4 | N. America | 25-60N | 125W-65W |

**DMS — `dms_rgn_scale(6)`** (ocean basins):

| Index | Region | Lat | Lon |
|-------|--------|-----|-----|
| 1 | N. Pacific  |   0-60N | 120E-100W |
| 2 | S. Pacific  | 60S- 0  | 145E- 70W |
| 3 | N. Atlantic |   0-60N | 100W- 20E |
| 4 | S. Atlantic | 60S- 0  |  70W- 20E |
| 5 | Indian      | 30S-30N |  30E-120E |
| 6 | Southern    | 90S-60S | all       |

**Carbon — `oc_rgn_scale(9)` and `bc_rgn_scale(9)`** (same boxes,
tested in order; first match wins so anthropogenic FF regions
take precedence over the Boreal NH catch-all):

| Index | Region | Lat | Lon |
|-------|--------|-----|-----|
| 1 | Amazon          | 15S-10N |  80W- 50W |
| 2 | S. Africa       | 35S- 5N |  10E- 45E |
| 3 | SE Asia/Indo    | 10S-20N |  95E-150E |
| 4 | Australia       | 40S-10S | 110E-155E |
| 5 | E. Asia         | 20N-50N | 100E-145E |
| 6 | S. Asia         |  5N-35N |  65E- 95E |
| 7 | Europe          | 35N-65N |  10W- 50E |
| 8 | N. America      | 25N-60N | 125W- 65W |
| 9 | Boreal NH       | 50N-75N | all       |

### Step 5: Validate

Compare tuned MATCH output against:
- Held-out satellite data (train on MODIS, validate on VIIRS)
- CDAS baseline runs
- Per-species climatologies (when available)

## Namelist Integration

All parameters are in namelist `/nlist/` with identity defaults.  Example:

```fortran
sslt_scale    = 1.08   ! global sea salt multiplier (per month)
n_sslt_bands  = 3      ! N active latitude bands below (0 disables)
sslt_bands    = -90., -45., 1.15,  ! Southern Ocean
                -45.,  45., 1.00,  ! mid-latitudes
                 45.,  90., 0.95   ! Northern Ocean
dst_rgn_scale = 0.14, 0.21, 1.75, 3.75, 1.4, 1.0, 0.45
                       ! Sahara, Arabia, C.Asia, Gobi, Australia, SW-NAm, Patagonia
so2_rgn_scale = 1.0, 1.0, 1.0, 1.0
                       ! E.Asia, S.Asia, Europe, N.America
dms_rgn_scale = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                       ! N.Pac, S.Pac, N.Atl, S.Atl, Indian, Southern
oc_rgn_scale  = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
bc_rgn_scale  = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                       ! Amazon, S.Africa, SE.Asia, Aust, E.Asia, S.Asia,
                       ! Europe, N.America, Boreal NH
! Legacy (superseded by dst_rgn_scale):
! vwc_scale  = 1.0     ! dust soil moisture rescaling slope
! vwc_offset = 0.0     ! dust soil moisture rescaling offset
```

- `dst_rgn_scale` — applied in `dst/dstmbl.F` to `flx_mss_hrz_slt_ttl`
  alongside `lnd_frc_mbl` and `flx_mss_fdg_fct`, per-gridpoint based on
  hardcoded region boxes.
- `vwc_scale` / `vwc_offset` — applied in `dst/dstmbl.F` after
  `vwc_sfc_get()`, clamped to `[0.10, 0.43]` when engaged.
- `sslt_scale` and `sslt_bands` — applied in `src/getsslt.F` after
  concentration and unit conversion; band lookup uses the row latitude
  in radians passed down from `src/physlic.F`.  Validation runs at
  startup in `src/main.F` (latitude range, lat_min<lat_max, no overlap).
- `so2_rgn_scale` and `dms_rgn_scale` — applied in
  `src_scyc/sulemis.F90`.  SOx scaling multiplies all four
  SOx-derived `sflx` slots (98% SO2 + 2% SO4, both 0m and >100m).
  DMS scaling multiplies the DMS slot.  The latitude in radians is
  passed in from `src/surface.F` and longitude is recovered per
  gridpoint as `(i-1)*360/plon`.
- `oc_rgn_scale` and `bc_rgn_scale` — applied in `src/caer.F90`
  inside `caersf()` to the input components from `caerbnd`.
  OC scale multiplies BBOCSF, FFOCSF, NOCSF; BC scale multiplies
  BBBCSF, FFBCSF.  Latitude/longitude are resolved as for sulfate.

## Output

Per-species per-calendar-month coefficient tables follow.  Sea salt has
a single global multiplier (`sslt_scale`) and an optional set of latitude
bands (`sslt_bands`); use "—" in the bands column when no bands are active.
The other species use per-region arrays; "TBD" entries are unset (treat
as 1.0 at run time).

### Dust (`dst_rgn_scale`, 7 regions)

| Month | Sahara | Arabia | C.Asia | Gobi | Aust | SW-NAm | Patag | Notes |
|-------|--------|--------|--------|------|------|--------|-------|-------|
| Jan   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Feb   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Mar   | 0.14 | 0.21 | 1.75 | 3.75 | 1.4 | 1.0 | 0.45 | tune5 — VIIRS 2025 |
| Apr   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| May   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Jun   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Jul   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Aug   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Sep   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Oct   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Nov   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Dec   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |

### Sea salt (`sslt_scale` + `sslt_bands`)

`sslt_bands` column lists active `(lat_min, lat_max, scale)` triples, or
"—" when `n_sslt_bands = 0` (global `sslt_scale` only).

| Month | sslt_scale | sslt_bands | Notes |
|-------|------------|------------|-------|
| Jan   | TBD  | TBD | |
| Feb   | TBD  | TBD | |
| Mar   | 1.08 | —   | tune5 — VIIRS 2025 (pre-bands) |
| Apr   | TBD  | TBD | |
| May   | TBD  | TBD | |
| Jun   | TBD  | TBD | |
| Jul   | TBD  | TBD | |
| Aug   | TBD  | TBD | |
| Sep   | TBD  | TBD | |
| Oct   | TBD  | TBD | |
| Nov   | TBD  | TBD | |
| Dec   | TBD  | TBD | |

### Sulfate (`so2_rgn_scale`, 4 regions)

| Month | E.Asia | S.Asia | Europe | N.Am | Notes |
|-------|--------|--------|--------|------|-------|
| Jan   | TBD | TBD | TBD | TBD | |
| Feb   | TBD | TBD | TBD | TBD | |
| Mar   | TBD | TBD | TBD | TBD | |
| Apr   | TBD | TBD | TBD | TBD | |
| May   | TBD | TBD | TBD | TBD | |
| Jun   | TBD | TBD | TBD | TBD | |
| Jul   | TBD | TBD | TBD | TBD | |
| Aug   | TBD | TBD | TBD | TBD | |
| Sep   | TBD | TBD | TBD | TBD | |
| Oct   | TBD | TBD | TBD | TBD | |
| Nov   | TBD | TBD | TBD | TBD | |
| Dec   | TBD | TBD | TBD | TBD | |

### DMS (`dms_rgn_scale`, 6 ocean basins)

| Month | N.Pac | S.Pac | N.Atl | S.Atl | Indian | Southern | Notes |
|-------|-------|-------|-------|-------|--------|----------|-------|
| Jan   | TBD | TBD | TBD | TBD | TBD | TBD | |
| Feb   | TBD | TBD | TBD | TBD | TBD | TBD | |
| Mar   | TBD | TBD | TBD | TBD | TBD | TBD | |
| Apr   | TBD | TBD | TBD | TBD | TBD | TBD | |
| May   | TBD | TBD | TBD | TBD | TBD | TBD | |
| Jun   | TBD | TBD | TBD | TBD | TBD | TBD | |
| Jul   | TBD | TBD | TBD | TBD | TBD | TBD | |
| Aug   | TBD | TBD | TBD | TBD | TBD | TBD | |
| Sep   | TBD | TBD | TBD | TBD | TBD | TBD | |
| Oct   | TBD | TBD | TBD | TBD | TBD | TBD | |
| Nov   | TBD | TBD | TBD | TBD | TBD | TBD | |
| Dec   | TBD | TBD | TBD | TBD | TBD | TBD | |

### Organic carbon (`oc_rgn_scale`, 9 regions)

| Month | Amaz | S.Afr | SE.Asia | Aust | E.Asia | S.Asia | Europe | N.Am | Boreal | Notes |
|-------|------|-------|---------|------|--------|--------|--------|------|--------|-------|
| Jan   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Feb   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Mar   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Apr   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| May   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Jun   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Jul   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Aug   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Sep   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Oct   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Nov   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Dec   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |

### Black carbon (`bc_rgn_scale`, 9 regions)

| Month | Amaz | S.Afr | SE.Asia | Aust | E.Asia | S.Asia | Europe | N.Am | Boreal | Notes |
|-------|------|-------|---------|------|--------|--------|--------|------|--------|-------|
| Jan   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Feb   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Mar   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Apr   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| May   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Jun   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Jul   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Aug   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Sep   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Oct   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Nov   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |
| Dec   | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | |

## Evaluation — March 2026 (tune5 vs VIIRS Mar 2025)

Full 41-day run (10-day Feb spin-up + 31-day March), CORe meteorology,
compared against NOAA-20 VIIRS Ed1C `TOTEXTTAU` climatology.

Regional AEROD (area-weighted box mean):

| Region | Box (lat, lon) | Scale | TUNE5 | BASE | VIIRS | T5 − V |
|--------|----------------|-------|-------|------|-------|--------|
| Sahara/Sahel       | 10-30N, 20W-30E     | 0.14 | 0.41 | 2.03 | 0.38 | **+0.03** |
| Arabia             | 15-30N, 30-55E      | 0.21 | 0.65 | 3.12 | 0.42 | **+0.23** ⚠ |
| Central Asia       | 30-50N, 55-90E      | 1.75 | 0.24 | 0.41 | 0.25 | −0.01 |
| Gobi/Taklamakan    | 35-50N, 90-115E     | 3.75 | 0.25 | 0.16 | 0.25 | 0.00 |
| Australia          | 35-15S, 115-145E    | 1.4  | 0.12 | 0.09 | 0.11 | 0.00 |
| SW N.America       | 25-40N, 115W-100W   | 1.0  | 0.03 | 0.03 | 0.09 | −0.06 (box) |
| Patagonia          | 55-35S, 75W-60W     | 0.45 | 0.09 | 0.12 | 0.08 | 0.00 |
| Southern Ocean     | 50-65S, all         | —    | 0.10 | 0.10 | 0.09 | +0.02 |
| Central S.Pacific  | 30-50S, 150W-90W    | —    | 0.07 | 0.06 | 0.08 | −0.01 |
| Remote S.Atlantic  | 30-50S, 30W-10E     | —    | 0.10 | 0.10 | 0.09 | +0.01 |

**Global area-weighted mean:** TUNE5 = 0.122, BASE (identity) = 0.247, VIIRS = 0.157.

Outstanding issues:
- **Arabia** (+0.23 over): scale 0.21 still too large.  Nearby Sahara
  transport may set an effective emission floor independent of local scale.
- **SW N.America** (−0.06 under): baseline dust ≈ 0 in this box regardless
  of scale; the box lies outside the Great Basin / Chihuahua Desert source
  zones and should be relocated.
- **Global under-bias** (−0.035): the 7 source boxes cover only a small
  fraction of the globe; the tuneable knobs cannot raise AOD over
  non-source continents and most of the tropical/mid-latitude oceans,
  where VIIRS observes more AOD than MATCH produces.
