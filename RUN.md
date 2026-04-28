# Running MATCH

## Run directory setup

Every MATCH run needs its own work directory with the following files before
the binary can execute.  The run script (`control/run_*.sh`) creates the
directory and writes the `namelist`, but the ancillary files must already be
in place.

### Required symlinks (from repo)

These point to static data shipped with the MATCH source tree:

| Symlink | Target | Purpose |
|---------|--------|---------|
| `ginouxbasin` | `$MATCH/data/ginouxbasin` | Dust source basin map |
| `obs.source` | `$MATCH/data/obs.source` | Observing station metadata |
| `SAGE_1020_OD_Filled.dat` | `$MATCH/data/SAGE_1020_OD_Filled.dat` | Stratospheric aerosol climatology |
| `grib.table2.rean` | `$MATCH/grib/grib.table2.rean` | GRIB parameter table |

### GRIB file lists

The model reads meteorological input via `grb2d.list` and `grbsanl.list`.
Each line has the format:

    <filename> <byte_offset_or_0> <YYYYMMDDHHMM>

**CDAS runs** can symlink the default lists from `$MATCH/grib/`:

    ln -s $MATCH/grib/grb2d.list .
    ln -s $MATCH/grib/grbsanl.list .

**CORe runs** need custom list files because the GRIB filenames differ.
Example for Jan-Mar 2026 CORe:

    A27340-202601.grb2d 0 202601010000
    A27339-202602.grb2d 0 202602010000
    A27338-202603.grb2d 0 202603010000

The byte-offset field is `0` for CORe (one month per file); CDAS files from
GDEX use non-zero offsets because multiple months may share an archive file.

### Restart files (for continuation runs)

To restart from an earlier run, copy the restart binary and the `rstrt`
pointer file into the new work directory:

    cp <source_run>/r<NNNN> <new_run>/
    cp <source_run>/rstrt   <new_run>/

The `rstrt` file is a fixed-format text file.  The run script reads field 3
(day count) to compute the remaining integration length.

### Namelist

Written by the run script into `<wrkdir>/namelist`.  Key variables:

| Variable | Type | Description |
|----------|------|-------------|
| `DELT` | int | Timestep in seconds (typically 900) |
| `ICDATE` | int | Start date YYYYMMDD |
| `NESTEP` | int | Negative = total days; positive = timestep count |
| `RSTFLG` | int | 0 = cold start, 1 = restart |
| `RPTHDYN` / `LPTHDYN` | char | Path to met-data GRIB directory |
| `vwc_scale` | real | Dust soil moisture rescaling slope (default 1.0 = identity) |
| `vwc_offset` | real | Dust soil moisture rescaling offset (default 0.0) |
| `sslt_scale` | real | Global sea salt emission multiplier (default 1.0 = no change) |
| `n_sslt_bands` | int | Number of active sea salt latitude bands (0 disables). See below. |
| `sslt_bands(3,10)` | real array | Per-band sea salt tuning: `(lat_min, lat_max, scale)` triples. See below. |
| `dst_rgn_scale(7)` | real array | Per-region dust emission multiplier (default 1.0 each). See below. |
| `so2_rgn_scale(4)` | real array | Per-region SOx (SO2+SO4) emission multiplier (default 1.0). See below. |
| `dms_rgn_scale(6)` | real array | Per-basin DMS emission multiplier (default 1.0). See below. |
| `oc_rgn_scale(9)` | real array | Per-region organic carbon emission multiplier (default 1.0). See below. |
| `bc_rgn_scale(9)` | real array | Per-region black carbon emission multiplier (default 1.0). See below. |

#### Regional dust tuning (`dst_rgn_scale`)

A 7-element array that multiplies the horizontal saltation flux in
`dst/dstmbl.F` at each gridpoint, based on which named region the point
falls in. Gridpoints outside all regions use scale 1.0. The region boxes
are hardcoded in `dstmbl.F`:

| Index | Region | Lat bounds | Lon bounds |
|-------|--------|-----------|------------|
| 1 | Sahara/Sahel | 10N–30N | 20W–30E |
| 2 | Arabia | 15N–30N | 30E–55E |
| 3 | Central Asia | 30N–50N | 55E–90E |
| 4 | Gobi/Taklamakan | 35N–50N | 90E–115E |
| 5 | Australia | 35S–15S | 115E–145E |
| 6 | SW N.America | 25N–40N | 115W–100W |
| 7 | Patagonia | 55S–35S | 75W–60W |

Example:
```fortran
dst_rgn_scale = 0.28, 0.40, 1.2, 4.0, 3.0, 1.0, 0.5
```

Use when the global `vwc_scale`/`vwc_offset` rescale cannot match AOD
in multiple source regions simultaneously — Sahara/Arabia need strong
suppression (CORe SOILW ≈ 0.02 → max emission) while Gobi/Australia
are under-emitting. The scalar-per-region approach decouples these and
is easier to iterate against satellite climatologies than the
soil-moisture rescale.

#### Regional sea salt tuning (`sslt_bands`, `n_sslt_bands`)

User-defined latitude bands, each with its own multiplier, applied in
`src/getsslt.F`.  `sslt_bands` is a `(3, max_sslt_bands=10)` array.  Column
`ib` is the triple `(lat_min_deg, lat_max_deg, scale)`.  Only the first
`n_sslt_bands` columns are read; the rest are ignored.  Columns whose
latitude lies in `[lat_min, lat_max]` are multiplied by `scale` (first
matching band wins).  `sslt_scale` is applied multiplicatively on top of
the band scale (or alone, when `n_sslt_bands = 0`).

Validation at startup (run aborts if any fail):
- `0 <= n_sslt_bands <= 10`
- each band has `-90 <= lat_min < lat_max <= 90`
- no two bands overlap (touching at the boundary is allowed)

Example — three bands split at ±45°:
```fortran
n_sslt_bands = 3
sslt_bands   = -90., -45., 1.15,
               -45.,  45., 1.00,
                45.,  90., 0.95
```

#### Regional sulfate tuning (`so2_rgn_scale`, `dms_rgn_scale`)

SOx-derived emissions (98% SO2 + 2% SO4, both 0m and >100m) are scaled
together by `so2_rgn_scale(4)` — applied per gridpoint in
`src_scyc/sulemis.F90` based on which named anthropogenic source region
the point falls in. Volcanic SO2 is unaffected. Gridpoints outside all
regions use scale 1.0.

| Index | Region | Lat bounds | Lon bounds |
|-------|--------|-----------|------------|
| 1 | E. Asia | 20N–50N | 100E–145E |
| 2 | S. Asia | 5N–35N | 65E–95E |
| 3 | Europe | 35N–65N | 10W–50E |
| 4 | N. America | 25N–60N | 125W–65W |

DMS surface flux is scaled by `dms_rgn_scale(6)` — also applied in
`sulemis.F90`, by ocean basin:

| Index | Region | Lat bounds | Lon bounds |
|-------|--------|-----------|------------|
| 1 | N. Pacific | 0–60N | 120E–100W |
| 2 | S. Pacific | 60S–0 | 145E–70W |
| 3 | N. Atlantic | 0–60N | 100W–20E |
| 4 | S. Atlantic | 60S–0 | 70W–20E |
| 5 | Indian | 30S–30N | 30E–120E |
| 6 | Southern | 90S–60S | all |

Example:
```fortran
so2_rgn_scale = 0.85, 1.20, 0.90, 1.00
dms_rgn_scale = 1.00, 1.10, 1.00, 1.00, 1.05, 1.20
```

#### Regional carbon tuning (`oc_rgn_scale`, `bc_rgn_scale`)

Organic and black carbon are scaled independently by 9-element arrays,
applied in `src/caer.F90` `caersf()` to the input components from
`caerbnd` before they are composed into surface fluxes. OC scale
multiplies BBOCSF + FFOCSF + NOCSF; BC scale multiplies BBBCSF + FFBCSF.

Both arrays use the same 9 regions, tested in order; first match wins
(anthropogenic FF boxes precede the Boreal NH catch-all):

| Index | Region | Lat bounds | Lon bounds |
|-------|--------|-----------|------------|
| 1 | Amazon | 15S–10N | 80W–50W |
| 2 | S. Africa | 35S–5N | 10E–45E |
| 3 | SE Asia/Indonesia | 10S–20N | 95E–150E |
| 4 | Australia | 40S–10S | 110E–155E |
| 5 | E. Asia | 20N–50N | 100E–145E |
| 6 | S. Asia | 5N–35N | 65E–95E |
| 7 | Europe | 35N–65N | 10W–50E |
| 8 | N. America | 25N–60N | 125W–65W |
| 9 | Boreal NH | 50N–75N | all |

Example:
```fortran
oc_rgn_scale = 1.20, 1.10, 1.30, 0.95, 0.85, 1.05, 0.90, 1.00, 1.40
bc_rgn_scale = 1.20, 1.10, 1.30, 0.95, 0.80, 1.10, 0.85, 1.00, 1.40
```

### Quick checklist

```
wrkdir=/path/to/new/run
MATCH=$HOME/EarthSystem/MATCH

mkdir -p "$wrkdir"

# Static data
ln -s $MATCH/data/ginouxbasin        "$wrkdir/"
ln -s $MATCH/data/obs.source         "$wrkdir/"
ln -s $MATCH/data/SAGE_1020_OD_Filled.dat "$wrkdir/"
ln -s $MATCH/grib/grib.table2.rean   "$wrkdir/"

# GRIB lists (CDAS: symlink; CORe: create custom files)
ln -s $MATCH/grib/grb2d.list         "$wrkdir/"
ln -s $MATCH/grib/grbsanl.list       "$wrkdir/"

# Restart (if continuing from another run)
cp <source>/r<NNNN> "$wrkdir/"
cp <source>/rstrt   "$wrkdir/"
```

Then launch the run script from `control/`:

    cd $MATCH/control
    ./run_<case>.sh > ~/Data/MATCH/2026/<case>.log 2>&1 &

---

## Migrating from CDAS to CORe meteorology

NOAA CORe reanalysis can replace CDAS as the meteorological driver for MATCH.
The model binary and ancillary data are unchanged; only the meteorological
input configuration differs.

### What changes

| Item | CDAS | CORe |
|------|------|------|
| GRIB directory | `$HOME/Data/NCAR_A` (or similar) | `$HOME/Data/NCAR_A_CORe` |
| List byte offset | Non-zero (multi-month archive files from GDEX) | `0` (one standalone file per month) |
| GRIB filenames | `A20066-199801.grb2d` etc. | `A27340-202601.grb2d` etc. |
| SOILW encoding | 9-bit, range 0.10-0.43, floor artifact at 0.10 | 18-bit, full range 0.02-1.0 |

### Steps

1. **Download CORe GRIB files** into a dedicated directory (e.g.
   `$HOME/Data/NCAR_A_CORe`).  You need both `grb2d` (2-D surface fields)
   and `grbsanl` (3-D sigma-level fields) for each month.

2. **Create custom GRIB list files** in the run directory.  Use byte offset
   `0` for all entries.  Include the month before and after the target period
   for time interpolation at boundaries.  Example for a Feb 2026 run:

   `grb2d.list`:
   ```
   A27340-202601.grb2d 0 202601010000
   A27339-202602.grb2d 0 202602010000
   A27338-202603.grb2d 0 202603010000
   ```

   `grbsanl.list`:
   ```
   A27329-202601.grbsanl 0 202601010000
   A27328-202602.grbsanl 0 202602010000
   A27327-202603.grbsanl 0 202603010000
   ```

3. **Update the run script** namelist:
   - Set `RPTHDYN` and `LPTHDYN` to the CORe GRIB directory.
   - Optionally set `vwc_scale` and `vwc_offset` to tune dust emission
     globally via a linear rescale of the volumetric water content.
     Without rescaling (defaults 1.0 and 0.0), CORe's finer soil moisture
     resolution produces different dust AOD than CDAS.
   - Optionally set `dst_rgn_scale` (7-element array) to tune dust
     emission per named source region (Sahara, Arabia, C.Asia, Gobi,
     Australia, SW N.America, Patagonia) — more controllable than
     `vwc_scale`/`vwc_offset` when regional biases differ in sign.
   - Optionally set `sslt_scale` for a global sea salt emission multiplier,
     and/or `sslt_bands` + `n_sslt_bands` for per-latitude-band scaling
     (e.g. a stronger multiplier over the Southern Ocean where winds and
     sea-salt AOD are highest).  See the regional sea salt tuning section
     above.
   - Optionally set `so2_rgn_scale` (4-element) and `dms_rgn_scale`
     (6-element) for per-region sulfate tuning, and/or `oc_rgn_scale` /
     `bc_rgn_scale` (9-element each) for per-region carbon tuning.
     See the regional sulfate and carbon tuning sections above.
   - Coefficients should be tuned per month against assimilated AOD
     climatologies.  See `AOD_TUNE.md` for the plan.

4. **Everything else stays the same**: binary, ancillary symlinks, restart
   handling, output format.
