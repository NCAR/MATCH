#!/bin/csh -f
#$ -l arch=lx-amd64
#$ -q scf.q

set
source $HOME/.cshrc_ceres
printenv
set casedir = $cwd
set runfile = $0

set echo
date
hostname

set test = off
set count_inc = $nn
set count_max = $nn

set debug = 0

setenv CASE sarb4_viirs_hourly_$yyyy$mm
setenv BASECASE
set vers =  VERS
set binary = /homedir/dfillmor/MATCH/build/match_ed4

set wrkdir = /SCF/CERES/sarb/dfillmor/MATCH/$CASE

set rdyn = "/DSS"
set dyndir = "/SCF/CERES/sarb/dfillmor/NCAR_A"
set dyn = $dyndir
# for ncep, just need directory and grb{2d,sanl}.list

set satdir = /SCF/CERES/sarb/dfillmor/viirs_assim/$yyyy #assimilation data directory

if ( ! -d $dyndir ) mkdir -p $dyndir
if ( ! -d $wrkdir ) mkdir -p $wrkdir

cd $wrkdir

if ( $test == off ) then
     # check that signal hasn't been given to quit
    if ( -e quit ) then
        echo "found quit file, exiting"
        exit
    endif

    # get resubmit info
    set rst = 0
    set ndays = $count_inc
    if ( -e rstrt ) then
       set fline = `cat rstrt | head -1`
       if ( $#fline > 0 ) then
           @ ndays = $fline[3] + $count_inc
           set rst = 1
       endif
    endif
    if ( $ndays > $count_max ) set ndays = $count_max
    if ( x$ndays == x ) exit -1
    set runtime = -$ndays

    # record start time to catch model if quits too early
    set temp = `date`
    set start_hour = `echo "$temp[4]" | awk -F: '{print$1}'`
    # set write out times
    set rstfrq = -1
    set histfrq = -24  #averages written out every 1 hours
else
    set runtime = 24
    set histfrq = -6
    set rst = 0
    set rstfrq = -1
    setenv OMP_NUM_THREADS 1
endif

set histfrq = -1

cat >! namelist << END1
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
 ICDATE   = $yyyymmdd_begin
 RPTHDYN  = '$dyn'
 LPTHDYN  = '$dyn'
 tracfrht = .false.
 lpthtra  = '$wrkdir/init_file'
 readQIC  = .false.
 RMOUT    = .false.
 RMRST    = .false.
 wpasswd  = 'fillmore'
 doAODasm = .true.
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

set runstat = $status
if ($runstat != 0) exit

if ( $test == on ) then
    echo "run.csh: done with run with status $runstat"
    pwd;ls -l
    exit 0
endif

#check to see if have finished days
if ( $ndays >= $count_max ) then
    echo "ndays $ndays >= count_max $count_max "
    mail username@email.com <<EOF
$CASE has reached the $count_max days 
EOF
    touch done
    exit 0
endif

#see if was too short of a run
set temp = `date`
set end_hour = `echo "$temp[4]" | awk -F: '{print$1}'`
if ( $start_hour == $end_hour ) then
    echo "run time less than an hour, quitting"
    exit
endif
#check to see if signal has been given to quit
if ( -e quit ) exit

pwd;ls -l

cd $casedir

$runfile

exit 0

