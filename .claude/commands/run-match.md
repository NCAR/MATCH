Set up and launch a new MATCH run directory.

Ask the user for (or infer from context):
1. **Case name** — e.g. `sarb4_hourly_202602_core_spinup`
2. **Met source** — CDAS or CORe (determines GRIB list handling)
3. **Run period** — start date, number of days
4. **Restart source** — path to an existing run to copy restart files from (optional)
5. **vwc_scale / vwc_offset** — dust soil moisture rescaling coefficients (optional, default identity; superseded by dst_rgn_scale)
6. **dst_rgn_scale(7)** — per-region dust emission multiplier (optional, default 1.0 each; regions: Sahara, Arabia, C.Asia, Gobi, Aust, SW-NAm, Patagonia)
7. **sslt_scale** — global sea salt emission multiplier (optional, default 1.0)
8. **sslt_bands / n_sslt_bands** — optional per-latitude-band sea salt multipliers, `(lat_min, lat_max, scale)` triples (see RUN.md)
9. **so2_rgn_scale(4)** — per-region SOx emission multiplier (optional; regions: E.Asia, S.Asia, Europe, N.Am)
10. **dms_rgn_scale(6)** — per-basin DMS emission multiplier (optional; basins: N.Pac, S.Pac, N.Atl, S.Atl, Indian, Southern)
11. **oc_rgn_scale(9) / bc_rgn_scale(9)** — per-region organic / black carbon emission multipliers (optional; same 9 regions: Amazon, S.Africa, SE.Asia, Aust, E.Asia, S.Asia, Europe, N.Am, Boreal)

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
