c     $Header$ -*-fortran-*-

c     Purpose: Common block dst_chm_cmn stores chemical properties of dust

c     These variables are initialized in dstchm.F:dst_chm_cmn_ini()
c     dst_chm_cmn_ini() is called by dst_msc_cmn_ini()

c     dstchm.h MUST have access to dstgrd.h and to pmgrid.h to work

c     Usage: 
c     #include <dstchm.h> /* Chemical properties of dust */

      integer,parameter::idx_H2O2_gas=1 ! [idx] Index of gaseous H2O2
      integer,parameter::idx_HNO3_gas=2 ! [idx] Index of gaseous HNO3
      integer,parameter::idx_HO2_gas=3 ! [idx] Index of gaseous HO2
      integer,parameter::idx_N2O5_gas=4 ! [idx] Index of gaseous N2O5
      integer,parameter::idx_NO3_gas=5 ! [idx] Index of gaseous NO3
      integer,parameter::idx_O3_gas=6 ! [idx] Index of gaseous O3
      integer,parameter::idx_OH_gas=7 ! [idx] Index of gaseous OH
      integer,parameter::idx_SO2_gas=8 ! [idx] Index of gaseous SO2
      integer,parameter::idx_SO4_aer=9 ! [idx] Index of particulate SO4
      integer,parameter::idx_NO3_aer=10 ! [idx] Index of particulate NO3
      integer,parameter::idx_chm_end=10 ! [idx] Last chemical index

c     fxm: Eventually, uptake coefficients should be moved to dst_chm_rxr()
      real(r8) upt_cff_H2O2_dst     ! [frc] Uptake coefficient for H2O2 to dust
      real(r8) upt_cff_HNO3_dst     ! [frc] Uptake coefficient for HNO3 to dust
      real(r8) upt_cff_HO2_dst      ! [frc] Uptake coefficient for HO2 to dust
      real(r8) upt_cff_N2O5_dst     ! [frc] Uptake coefficient for N2O5 to dust
      real(r8) upt_cff_NO3_dst      ! [frc] Uptake coefficient for NO3 to dust
      real(r8) upt_cff_O3_dst       ! [frc] Uptake coefficient for O3 to dust
      real(r8) upt_cff_OH_dst       ! [frc] Uptake coefficient for OH to dust
      real(r8) upt_cff_SO2_dst      ! [frc] Uptake coefficient for SO2 to dust
      real(r8) vmr_HNO3_gas_mxm     ! [mlc mlc-1] Maximum HNO3 volume mixing ratio
      common /dst_chm_cmn/      ! dst_chm_cmn is initialized in dst_chm_cmn_ini()
     $     upt_cff_H2O2_dst,    ! [frc] Uptake coefficient for H2O2 to dust
     $     upt_cff_HNO3_dst,    ! [frc] Uptake coefficient for HNO3 to dust
     $     upt_cff_HO2_dst,     ! [frc] Uptake coefficient for HO2 to dust
     $     upt_cff_N2O5_dst,    ! [frc] Uptake coefficient for N2O5 to dust
     $     upt_cff_NO3_dst,     ! [frc] Uptake coefficient for NO3 to dust
     $     upt_cff_O3_dst,      ! [frc] Uptake coefficient for O3 to dust
     $     upt_cff_OH_dst,      ! [frc] Uptake coefficient for OH to dust
     $     upt_cff_SO2_dst,     ! [frc] Uptake coefficient for SO2 to dust
     $     vmr_HNO3_gas_mxm     ! [mlc mlc-1] Maximum HNO3 volume mixing ratio
      

