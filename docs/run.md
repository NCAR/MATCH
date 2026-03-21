# Running MATCH

## Overview

CERES Edition 4 MATCH runs one calendar month at a time, with a
one-month spin-up. For example, to produce output for August 2017,
the run starts on 2017-07-01 and ends on 2017-08-31 (62 days total).

## External Data

### NCEP/NCAR Reanalysis

MATCH is driven by NCEP/NCAR reanalysis meteorological fields. These
datasets are now part of GEDEX on NCAR HPC machines.

The original RDA dataset references (ds090.0):
- **grbsanl** — Pressure-level analysis fields
- **grb2d** — 2-D (surface/flux) fields

Download to `$DATA_DIR/NCAR_A`.

After downloading new months, update the GRIB listing files:

**`grib/grbsanl.list`** — add a line:
```
A25426-201710.grbsanl 421724160 201710010000
```

**`grib/grb2d.list`** — add a line:
```
A25437-201710.grb2d 189038592 201710010000
```

Format: `filename filesize yyyymm010000`

### MODIS Assimilation Data

Satellite AOD observations for assimilation are produced in a
separate MODIS preprocessing step and stored at:
```
$DATA_DIR/modis_assim_combo/<year>/
```

## Directory Setup

### Environment

Set the top-level data directory:
```csh
setenv DATA_DIR /path/to/data
```

Pre/post processing scripts and run scripts read `$DATA_DIR`.

### Run Directory

Create the run directory following this naming convention:
```
$DATA_DIR/MATCH/<year>/<runname>
```

The run name follows the convention `sarb4_C6_hourly_<yyyymm>`, e.g.:
```
mkdir -p $DATA_DIR/MATCH/2017/sarb4_C6_hourly_201708
```

### Ancillary Data

MATCH reads ancillary datasets from the run directory. Create
symbolic links to the repository data and GRIB listing files:
```csh
cd $DATA_DIR/MATCH/2017/sarb4_C6_hourly_201708
ln -s $HOME/MATCH/data/* .
ln -s $HOME/MATCH/grib/* .
```

Also add `MATCH/build` to your `PATH` for the `list_nc_files` utility.

## Run Script Configuration

Example run scripts are in `control/`. The key variables to set:

| Variable | Description | Example |
|----------|-------------|---------|
| `count_inc` | Total run days including spin-up | `62` |
| `count_max` | Total run days including spin-up | `62` |
| `CASE` | Run name | `sarb4_C6_hourly_201708` |
| `binary` | Path to MATCH executable | `$HOME/MATCH/build/match` |
| `wrkdir` | Run/output directory | `$DATA_DIR/MATCH/$CASE` |
| `dyndir` | NCEP reanalysis GRIB files | `$DATA_DIR/NCAR_A` |
| `satdir` | MODIS assimilation files | `$DATA_DIR/modis_assim_combo/2017` |
| `ICDATE` | Run start date (yyyymmdd) | `20170701` |

`count_inc` and `count_max` should equal the total number of days
in the spin-up month plus the target month.

## Running

```csh
cd MATCH/control
./run_sarb4_hourly_yyyymm.csh
```

MATCH writes all output to `$wrkdir`.

## Run Scripts

| Script | Purpose |
|--------|---------|
| `run_match_yyyy01-06.csh` | Batch run for January–June |
| `run_match_yyyy07-12.csh` | Batch run for July–December |
| `run_sarb4_daily_yyyymm.csh` | Single-month run, daily output |
| `run_sarb4_hourly_yyyymm.csh` | Single-month run, hourly output |
| `run_sarb4_viirs_hourly_yyyymm.csh` | Single-month run with VIIRS, hourly |
| `generate_run.py` | Generate run scripts for a given period |
| `CalUtil.py` | Calendar utility (day counts, date math) |
