#!/bin/csh
#$ -N MATCH_REFORMAT
#$ -o match_ed2_reformat.olog
#$ -e match_ed2_reformat.elog
#$ -cwd
# Code should work for both X86 and PPC
###########XX $ -l arch='lx26-ppc64'
###########XX $ -l arch=lx-amd64
echo CPU: $CPU

set YYYY = '2008'
set MM = '01'
set DD = '01'
echo current working directory: $cwd

setenv LD_LIBRARY_PATH /SPG_ops/utils/CERES/TOOLKIT/scf_toolkit5.2.18v1-x86_64/lib:/SPG_ops/utils/x86_64/gcc-4.5.3-ceres/lib64

set inpath  = $cwd/in/
set outpath = $cwd/out/

set YYYYMMDD = $YYYY$MM$DD
set file_ed2_nc = ${inpath}match_sarb_ed2_v1.00_20091113_$YYYYMMDD.nc

set aot_bin = ${outpath}match_ed2_aot.$YYYYMMDD
set prf_bin = ${outpath}match_ed2_prf.$YYYYMMDD
set aot_sage_bin = ${outpath}match_ed2_aot.sage.$YYYYMMDD
set prf_sage_bin = ${outpath}match_ed2_prf.sage.$YYYYMMDD

# Output Binaries are Big Endian set folloeing environment variable for x86 Gfortran
setenv GFORTRAN_CONVERT_UNIT BIG_ENDIAN


rm -f ${outpath}*.$YYYYMMDD


if ( ! -e $file_ed2_nc ) exit
echo MATCH_ED2_NETCDF_FILE :$file_ed2_nc


#getarg(1 : Ed2Match NETCDF File
#getarg(2 :  Ed2 Binary Reformated MATCH AOT file (Used for  ed2CRS)
#getarg(3 :  Ed2 Binary Reformated MATCH PRF file (Used for  ed2CRS)
#getarg(4 : MONTH 01 - 12
#getarg(6 :Ed2 Binary Reformated SAGED MATCH AOT file (Used for  ed3 SYN)
#getarg(7 :Ed2 Binary Reformated SAGED MATCH PRF file (Used for  ed3 SN)


reformat_ed2_match.$CPU.exe $file_ed2_nc $aot_bin $prf_bin $MM  $aot_sage_bin  $prf_sage_bin


# ARE THE RESILTS THE SAME AS EARLIER PROCEESSD FILES

diff -s $aot_bin /SCF/CERES/sarb/rose/match/processed/aots/MATCH_TERRA_AOTS_MODIS.C5.$YYYYMMDD 

diff -s $prf_bin /SCF/CERES/sarb/rose/match/processed/vertical/MATCH_TERRA_VERTICAL_MODIS.C5.$YYYYMMDD 

diff -s $aot_sage_bin /SCF/CERES/sarb/rose/match/destrat/saged/aot/C5/$YYYY/$MM/MATCH_TERRA_AOTS_MODIS.C5.$YYYYMMDD

diff -s $prf_sage_bin /SCF/CERES/sarb/rose/match/destrat/saged/profiles/C5/$YYYY/$MM/MATCH_TERRA_VERTICAL_MODIS.C5.$YYYYMMDD


exit
