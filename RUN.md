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
| `sslt_scale` | real | Sea salt emission flux multiplier (default 1.0 = no change) |

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
   - Optionally set `vwc_scale` and `vwc_offset` to tune dust emission.
     Without rescaling (defaults 1.0 and 0.0), CORe's finer soil moisture
     resolution produces different dust AOD than CDAS.
   - Optionally set `sslt_scale` to tune sea salt emission (default 1.0).
     CORe wind speeds differ from CDAS, producing a systematic sea salt
     bias over oceans.
   - Both sets of coefficients should be tuned per month against
     assimilated AOD climatologies.  See `AOD_TUNE.md` for the plan.

4. **Everything else stays the same**: binary, ancillary symlinks, restart
   handling, output format.
