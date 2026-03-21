c     $Header$ -*-fortran-*- 
      
c     Purpose: dstszdst.h contains common blocks which store diagnostic 
c     microphysical properties of the mineral dust aerosol.
c     Variables in dstszdst are not needed every timestep by source and sink routines
c     Related microphysical quantities, needed every timestep, are stored in dstaer.h
      
c     Most of these variables are initialized in dst_szdst_ini()
c     dst_szdst_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
c     dstszdst.h MUST have access to dstgrd.h and to pmgrid.h to work
      
c     Usage:
c     #include <dstszdst.h> /* Diagnostic aerosol microphysical properties */
      
c     dst_src_bd_cmn is block data initialized in dstszdstbd()
      real(r8) dmt_vma_src(dst_src_nbr) ! [m] Mass median diameter
      real(r8) gsd_anl_src(dst_src_nbr) ! [frc] Geometric standard deviation
      real(r8) mss_frc_src(dst_src_nbr) ! [frc] Mass fraction
      common /dst_src_bd_cmn/   ! dst_aer_bd_cmn is initialized in dstszdstbd()
     $     dmt_vma_src,         ! [m] Mass median diameter
     $     gsd_anl_src,         ! [frc] Geometric standard deviation
     $     mss_frc_src          ! [frc] Mass fraction
      
c     dst_src_cmn is computed at run time in dst_szdst_ini()
      real(r8) ovr_src_snk_frc(dst_src_nbr,dst_nbr) ! [frc] Overlap of src with snk
      real(r8) ovr_src_snk_mss_ttl  ! [frc] Total transported mass fraction of dust flux
      common /dst_src_cmn/      ! dst_src_cmn is initialized in dst_szdst_ini()
     $     ovr_src_snk_frc,     ! [frc] Overlap of src with snk
     $     ovr_src_snk_mss_ttl  ! [frc] Total transported mass fraction of dust flux

c     dst_dgn_cmn is computed at run time in dst_szdst_ini()
      real(r8) dmt_grd(dst_nbr+1)   ! [m] Particle diameter grid
      real(r8) dmt_max(dst_nbr)     ! [m] Maximum diameter in bin
      real(r8) dmt_min(dst_nbr)     ! [m] Minimum diameter in bin
      real(r8) dmt_naa(dst_nbr)     ! [m] Number mean particle size
      real(r8) dmt_nma(dst_nbr)     ! [m] Number median particle diameter
      real(r8) dmt_nmr(dst_nbr)     ! [m] Number median diameter resolved
      real(r8) dmt_nwr(dst_nbr)     ! [m] Number mean diameter resolved
      real(r8) dmt_sma(dst_nbr)     ! [m] Surface median particle size
      real(r8) dmt_smr(dst_nbr)     ! [m] Surface area median diameter resolved
      real(r8) dmt_swa(dst_nbr)     ! [m] Surface area weighted mean particle size
      real(r8) dmt_swr(dst_nbr)     ! [m] Surface area weighted diameter resolved
      real(r8) dmt_vma(dst_nbr)     ! [m] Mass median diameter
      real(r8) dmt_vmr(dst_nbr)     ! [m] Mass median diameter resolved
      real(r8) dmt_vwa(dst_nbr)     ! [m] Mass weighted mean particle size
      real(r8) gsd_anl(dst_nbr)     ! [frc] Geometric standard deviation
      real(r8) nbr_spc_rsl(dst_nbr) ! [# kg-1] Specific concentration resolved 
      real(r8) vlm_spc_rsl(dst_nbr) ! [m3 kg-1] Specific volume resolved 
      common /dst_dgn_cmn/      ! dst_dgn_cmn is initialized in dst_szdst_ini()
     $     dmt_grd,             ! [m] Particle diameter grid
     $     dmt_max,             ! [m] Maximum diameter in bin
     $     dmt_min,             ! [m] Minimum diameter in bin
     $     dmt_naa,             ! [m] Number mean particle size
     $     dmt_nma,             ! [m] Number median particle diameter
     $     dmt_nmr,             ! [m] Number median diameter resolved
     $     dmt_nwr,             ! [m] Number mean diameter resolved
     $     dmt_sma,             ! [m] Surface median particle size
     $     dmt_smr,             ! [m] Surface area median diameter resolved
     $     dmt_swa,             ! [m] Surface area weighted mean particle size
     $     dmt_swr,             ! [m] Surface area weighted diameter resolved
     $     dmt_vma,             ! [m] Mass median diameter
     $     dmt_vmr,             ! [m] Mass median diameter resolved
     $     dmt_vwa,             ! [m] Mass weighted mean particle size
     $     gsd_anl,             ! [frc] Geometric standard deviation
     $     nbr_spc_rsl,         ! [# kg-1] Specific concentration resolved 
     $     vlm_spc_rsl          ! [m3 kg-1] Specific volume resolved 
      
