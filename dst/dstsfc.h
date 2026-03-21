c     $Header$ -*-fortran-*-

c     Purpose: Common block dstsfc stores fields needed by dust parameterization
c     It needs to work independently of the host model (MATCH, CCM, CSM, MOZART)

c     params.h needed for PLON and PLAT
#ifndef PARAMS_H
#include <params.h>
#endif /* not PARAMS_H */ 

c     Usage: 
c     #include <dstsfc.h> /* 2-D surface fields on PLON x PLAT grid */

c     Time-independent surface information supplied by boundary dataset
      integer sfc_typ(PLON,PLAT) ! [idx] LSM surface type (0..28)
      common / dst_sfc_i /      ! dst_sfc_i is initialized in CCM:dst/dstnc()
     $     sfc_typ              ! [idx] LSM surface type (0..28)

      real(r8) lnd_frc_dry(PLON,PLAT) ! [frc] Land surface that is not lake or wetland (by area)
      real(r8) mbl_bsn_fct(PLON,PLAT) ! [frc] Topographic source limitation
      real(r8) mss_frc_CaCO3(PLON,PLAT) ! [frc] Mass fraction of CaCO3
      real(r8) mss_frc_cly(PLON,PLAT) ! [frc] Mass fraction of clay
      real(r8) mss_frc_snd(PLON,PLAT) ! [frc] Mass fraction of sand
      common / dst_sfc_r /      ! dst_sfc_r is initialized in CCM:dst/dstnc()
     $     lnd_frc_dry,         ! [frc] Land surface that is not lake or wetland (by area)
     $     mbl_bsn_fct,         ! [frc] Topographic source limitation
     $     mss_frc_CaCO3,       ! [frc] Mass fraction of CaCO3
     $     mss_frc_cly,         ! [frc] Mass fraction of clay
     $     mss_frc_snd          ! [frc] Mass fraction of sand

c     Time-varying surface information supplied by CCM
c     fxm: Probably should store radiative fluxes elsewhere
      real(r8) flx_LW_dwn_sfc(PLON,PLAT) ! [W m-2] Longwave downwelling flux at surface
      real(r8) flx_SW_abs_sfc(PLON,PLAT) ! [W m-2] Solar flux absorbed by ground
      common / dst_rad_ccm /    ! dst_rad_lsm is initialized in CCM:physics/tphysbc()
     $     flx_LW_dwn_sfc,      ! [W m-2] Longwave downwelling flux at surface
     $     flx_SW_abs_sfc       ! [W m-2] Solar flux absorbed by ground

c     NB: Use PLON,PLAT tokens instead of plon,plat
c     params.h was dropped from LSM after CCM3.5, but we reinclude it for PLON,PLAT tokens
c     In LSM, plon and plat are passed as variables into lsmdrv()
c     In CCM, plon and plat are parameters in every routine with pmgrid.h or prgrid.h
c     dst_sfc needs to work in both LSM and CCM so use PLON,PLAT

c     Time-varying surface information supplied by LSM
      logical sfc_set_sgs_mbl(PLON,PLAT) ! [flg] Are surface properties already set?
      real(r8) tpt_gnd(PLON,PLAT)   ! [K] Ground temperature
      real(r8) tpt_soi(PLON,PLAT)   ! [K] Soil temperature
      real(r8) vwc_sfc(PLON,PLAT)   ! [m3 m-3] Volumetric water content
      common / dst_hyd_lsm /    ! dst_hyd_lsm is initialized in CCM:lsm/lsmdrv()
     $     sfc_set_sgs_mbl,     ! [flg] Are surface properties already set?
     $     tpt_gnd,             ! [K] Ground temperature
     $     tpt_soi,             ! [K] Soil temperature
     $     vwc_sfc              ! [m3 m-3] Volumetric water content








