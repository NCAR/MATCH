c     $Header$ -*-fortran-*- 
      
c     Purpose: dstscv.h contains common blocks which store important 
c     microphysical properties of the mineral dust aerosol. 
c     Variables in dstscv are used each timestep by source and sink routines
c     Related, but strictly diagnostic, microphysical quantities are stored in dstszdst.h
      
c     Most of these variables are initialized in dst_szdst_ini()
c     dst_szdst_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
c     dstscv.h MUST have access to dstgrd.h and to pmgrid.h to work
      
c     Usage:
c     #include <dstscv.h> /* Aerosol microphysical properties */
      
c     dst_aer_cmn is computed at run time, initialized in dst_szdst_ini()
      real(r8) dmt_vwr(dst_nbr)     ! [m] Mass weighted diameter resolved
      real(r8) dns_aer(dst_nbr)     ! [kg m-3] Particle density
      real(r8) ovr_src_snk_mss(dst_src_nbr,dst_nbr) ! [frc] Overlap of src with snk
      real(r8) stk_crc(dst_nbr)     ! [frc] Correction to Stokes settling velocity
      real(r8) sfc_spc_rsl(dst_nbr) ! [m2 kg-1] Specific surface area resolved 
      common /dst_aer_cmn/      ! dst_aer_cmn is initialized in dst_szdst_ini()
     $     dmt_vwr,             ! [m] Mass weighted diameter resolved
     $     dns_aer,             ! [kg m-3] Particle density
     $     ovr_src_snk_mss,     ! [frc] Overlap of src with snk
     $     sfc_spc_rsl,         ! [m2 kg-1] Specific surface area resolved 
     $     stk_crc              ! [frc] Correction to Stokes settling velocity





