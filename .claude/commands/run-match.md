Set up and launch a new MATCH run directory.

Ask the user for (or infer from context):
1. **Case name** — e.g. `sarb4_hourly_202602_core_spinup`
2. **Met source** — CDAS or CORe (determines GRIB list handling)
3. **Run period** — start date, number of days
4. **Restart source** — path to an existing run to copy restart files from (optional)
5. **vwc_scale / vwc_offset** — soil moisture rescaling coefficients (optional, defaults to identity)

Then follow these steps:

## 1. Create the work directory

```
wrkdir=$HOME/Data/MATCH/<year>/<case>
mkdir -p "$wrkdir"
```

## 2. Symlink static ancillary data

```
MATCH=$HOME/EarthSystem/MATCH
ln -s $MATCH/data/ginouxbasin            "$wrkdir/"
ln -s $MATCH/data/obs.source             "$wrkdir/"
ln -s $MATCH/data/SAGE_1020_OD_Filled.dat "$wrkdir/"
ln -s $MATCH/grib/grib.table2.rean       "$wrkdir/"
```

## 3. Set up GRIB file lists

**CDAS:** symlink the default lists:
```
ln -s $MATCH/grib/grb2d.list    "$wrkdir/"
ln -s $MATCH/grib/grbsanl.list  "$wrkdir/"
```

**CORe:** create custom `grb2d.list` and `grbsanl.list` in the work directory.
Each line: `<filename> 0 <YYYYMMDDHHMM>`. List the CORe GRIB files in
`$HOME/Data/NCAR_A_CORe/` that cover the run period (need one month before
and after for interpolation).

## 4. Copy restart files (if restarting)

```
cp <source_run>/r<NNNN> "$wrkdir/"
cp <source_run>/rstrt   "$wrkdir/"
```

## 5. Verify the run directory

Before launching, confirm all required files are present:
- `ginouxbasin`, `obs.source`, `SAGE_1020_OD_Filled.dat`, `grib.table2.rean` (symlinks)
- `grb2d.list`, `grbsanl.list` (symlinks or files)
- `r<NNNN>` and `rstrt` (if restart run)
- That the GRIB files referenced in the lists actually exist in the dyndir

## 6. Launch

Either run an existing run script or create one. The run script writes the
namelist and invokes the binary. See `RUN.md` for namelist variable reference.

```
cd $MATCH/control
./run_<case>.sh > $HOME/Data/MATCH/<year>/<case>.log 2>&1 &
```

Report the PID and log path back to the user.
