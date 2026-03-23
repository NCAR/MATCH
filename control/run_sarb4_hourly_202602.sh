#!/bin/bash
#
# MATCH run script for February 2026
# Spin-up: January 2026 (31 days) + Target: February 2026 (28 days) = 59 days
#

set -e

casedir="$(pwd)"
runfile="$(readlink -f "$0")"

date
hostname

test_mode=off
count_inc=59
count_max=59

debug=0

export CASE=sarb4_hourly_202602
vers='MATCH SARB4 C6'
binary="$HOME/EarthSystem/MATCH/build/match"

export DATA_DIR="$HOME/Data"
wrkdir="$DATA_DIR/MATCH/2026/$CASE"

dyndir="$DATA_DIR/NCAR_A"
dyn="$dyndir"

satdir="$DATA_DIR/modis_assim_combo/2026"

[[ -d "$dyndir" ]] || mkdir -p "$dyndir"
[[ -d "$wrkdir" ]] || mkdir -p "$wrkdir"

cd "$wrkdir"

if [[ "$test_mode" == "off" ]]; then
    # check that signal hasn't been given to quit
    if [[ -e quit ]]; then
        echo "found quit file, exiting"
        exit
    fi

    # get resubmit info
    rst=0
    ndays=$count_inc
    if [[ -e rstrt ]]; then
        fline=$(head -1 rstrt)
        if [[ -n "$fline" ]]; then
            field3=$(echo "$fline" | awk '{print $3}')
            ndays=$(( field3 + count_inc ))
            rst=1
        fi
    fi
    if (( ndays > count_max )); then ndays=$count_max; fi
    if [[ -z "$ndays" ]]; then exit 1; fi
    runtime="-$ndays"

    # record start time to catch model if quits too early
    start_hour=$(date +%H)
    rstfrq=-1
    histfrq=-24
else
    runtime=24
    histfrq=-6
    rst=0
    rstfrq=-1
    export OMP_NUM_THREADS=1
fi

histfrq=-1

cat > namelist << END1
 &NLIST
 TITLE    = '$vers '
 DELT     = 1800
 RSTFLG   = $rst
 NESTEP   = $runtime
 NHTFRQ   = $histfrq
 RSTFRQ   = $rstfrq
 MFILT    = 1
 NDENS    = 1
 IRT      = 0
 RSTRRT   = 0
 ICDATE   = 20260101
 RPTHDYN  = '$dyn'
 LPTHDYN  = '$dyn'
 tracfrht = .false.
 lpthtra  = '$wrkdir/init_file'
 readQIC  = .false.
 RMOUT    = .false.
 RMRST    = .false.
 wpasswd  = 'fillmore'
 doAODasm = .false.
 assimdir = '$satdir'
 DOWETDEP = .true.
 DODRYDEP = .true.
 VOLCEMIS = .false.
 QMIN     = 12*1.e-20, 1.e-14
 soxsfconst = .false.
 soxsfyr  = 1995
 TRACNAM  = 'SO2', 'SO4', 'DMS', 'H2O2',
            'DSTQ01', 'DSTQ02', 'DSTQ03', 'DSTQ04',
            'OCPHO', 'BCPHO', 'OCPHI', 'BCPHI'
 SFLXNAM  = 'SO2SF','SO4SF','DMSSF',,,,,,
            'OCPHOSF','BCPHOSF'
 outinst = 'ORO', 'PHIS', 'PS', 'Z', 'ZI',
           'SO4', 'SO4OD',
           'SSLT', 'SSLTOD',
           'OCPHO', 'OCPHOOD',
           'OCPHI', 'OCPHIOD',
           'BCPHO', 'BCPHOOD',
           'BCPHI', 'BCPHIOD',
           'DSTQ', 'DSTQ01', 'DSTQ02', 'DSTQ03', 'DSTQ04',
           'DSTODXC', 'DSTODX01', 'DSTODX02', 'DSTODX03', 'DSTODX04',
           'AEROD'
  /
END1

$binary < namelist

runstat=$?
if (( runstat != 0 )); then exit $runstat; fi

if [[ "$test_mode" == "on" ]]; then
    echo "run.sh: done with run with status $runstat"
    pwd; ls -l
    exit 0
fi

# check to see if have finished days
if (( ndays >= count_max )); then
    echo "ndays $ndays >= count_max $count_max"
    touch done
    exit 0
fi

# see if was too short of a run
end_hour=$(date +%H)
if [[ "$start_hour" == "$end_hour" ]]; then
    echo "run time less than an hour, quitting"
    exit
fi
# check to see if signal has been given to quit
if [[ -e quit ]]; then exit; fi

pwd; ls -l

cd "$casedir"

exec "$runfile"
