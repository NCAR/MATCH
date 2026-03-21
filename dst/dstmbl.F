      module dstmbl
c     $Header$
      
c     Purpose: dstmbl.F controls mineral dust mobilization processes
      
c     Requires dst.h for DST_MSS_BDG
#include <dst.h> /* Dust preprocessor tokens */
      
      use precision
      use dstblm

c     $Header$ -*-fortran-*-
      
c     Purpose: Common block dst_blm_cmn stores fundamental physical constants
c     needed for boundary layer meteorology
      
c     These variables are initialized in dstcmnini:dst_cst_cmn_ini()
c     dst_cst_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun() 
      
c     Usage: 
c     #include <dstblm.h> /* Physical constants for boundary layer meteorology */
      
      real(r8)cp_vpr_rcp_cp_dry_m1 ! 0.83745 [frc] Constant for moist specific heat IrG81 p. 77
      real(r8)cst_Stefan_Boltzmann ! (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant
      real(r8)dns_H2O_lqd_std      ! [kg m-3] Density of liquid water
      real(r8)ltn_heat_fsn_H2O_std ! [J kg-1] Latent heat of fusion of H2O at 0 C, standard
      real(r8)ltn_heat_sbl_H2O_std ! [J kg-1] Latent heat of sublimation of H2O, standard
      real(r8)ltn_heat_vpr_H2O_std ! [J kg-1] Latent heat of vaporization of H2O, standard
      real(r8)one_mns_eps_H2O          ! (0.378) [frc] Constant for saturation specific humidity
      real(r8)spc_heat_H2O_ice_vlm ! [J m-3 K-1] Volumetric specific heat capacity of ice water
      real(r8)spc_heat_H2O_lqd_vlm ! [J m-3 K-1] Volumetric specific heat capacity of liquid water
      real(r8)spc_heat_dry_air     ! (1005.0) [J kg-1 K-1] IrG81 p. 25
      real(r8)tpt_frz_pnt          ! [K] Kelvin--Celsius scale offset
      common /dst_blm_cmn/      ! dst_blm_cmn is initialized in dst_cst_cmn_ini()
     $     cp_vpr_rcp_cp_dry_m1, ! 0.83745 [frc] Constant for moist specific heat IrG81 p. 77
     $     cst_Stefan_Boltzmann, ! (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant
     $     dns_H2O_lqd_std,     ! [kg m-3] Density of liquid water
     $     ltn_heat_fsn_H2O_std, ! [J kg-1] Latent heat of fusion of H2O at 0 C, standard
     $     ltn_heat_sbl_H2O_std, ! [J kg-1] Latent heat of sublimation of H2O, standard
     $     ltn_heat_vpr_H2O_std, ! [J kg-1] Latent heat of vaporization of H2O, standard
     $     one_mns_eps_H2O,         ! (0.378) [frc] Constant for saturation specific humidity
     $     spc_heat_H2O_ice_vlm, ! [J m-3 K-1] Volumetric specific heat capacity of ice water
     $     spc_heat_H2O_lqd_vlm, ! [J m-3 K-1] Volumetric specific heat capacity of liquid water
     $     spc_heat_dry_air,    ! (1005.0) [J kg-1 K-1] IrG81 p. 25
     $     tpt_frz_pnt          ! [K] Kelvin--Celsius scale offset
      
      

      contains
c      subroutine dst_mbl(hbuf,obuf,
      subroutine dst_mbl(obuf,
     $     doy,                 ! I [day] Day of year [1.0..366.0)
     $     hgt_mdp,             ! I [m] Midpoint height above surface
     $     lat_idx,             ! I [idx] Model latitude index
     $     lat_rdn,             ! I [rdn] Latitude
     $     oro,                 ! I [frc] Orography
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     prs_mdp,             ! I [Pa] Pressure
     $     prs_sfc,             ! I [Pa] Surface pressure
     $     q_H2O_vpr,           ! I [kg kg-1] Water vapor mixing ratio
     $     q_dst,               ! I/O [kg kg-1] Dust mixing ratio
     $     snw_hgt_lqd,         ! I [m] Equivalent liquid water snow depth
     $     tm_adj,              ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     tpt_mdp,             ! I [K] Temperature
     $     tpt_ptn_mdp,         ! I [K] Midlayer local potential temperature
     $     wnd_mrd_mdp,         ! I [m s-1] Meridional wind component
     $     wnd_znl_mdp)         ! I [m s-1] Zonal wind component
c     Purpose: Driver for aerosol mobilization
c     dst_mbl() is called by CCM:physics/tphysac(), MATCH:src/physlic()
c     NB: dst_mbl() is designed to require only single layer surface fields 
c     Eliminating multi-layer fields allows easier implementation in LSM
c++match
      use histout, only: outfld
c--match

      use dstmblutl

      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstaer.h> /* Aerosol microphysical properties */
#include <dstmsc.h> /* Miscellaneous variables for dust parameterization */
#include <dstnm.h> /* Nomenclature for outfld() */
!#include <dstblm.h> /* Physical constants for boundary layer meteorology */
c     Parameters
!      real(r8),parameter::flx_mss_fdg_fct=5.0e-4 ! [frc] Global mass flux tuning factor (empirical)
!pjr replace with convolution natalie and chao's fudge * 0.5
      real(r8),parameter::flx_mss_fdg_fct=7.0e-4 ! [frc] Global mass flux tuning factor (empirical)
      real(r8),parameter::hgt_rfr=10.0 ! [m] Reference height for mobilization processes
      real(r8),parameter::hgt_zpd_mbl=0.0 ! [m] Zero plane displacement for erodible surfaces
      real(r8),parameter::rgh_mmn_mbl=100.0e-6 ! [m] Roughness length momentum for erodible surfaces MaB95 p. 16420, GMB98 p. 6205
c     fxm: rgh_mmn_smt set to 33.3e-6 um, MaB95 p. 16426 recommend 10.0e-6
      real(r8),parameter::rgh_mmn_smt=33.3e-6 ! [m] Smooth roughness length MaB95 p. 16426, MaB97 p. 4392, GMB98 p. 6207
      real(r8),parameter::wnd_min_mbl=1.0 ! [m s-1] Minimum windspeed used for mobilization 
c     Input
      integer,intent(in)::lat_idx           ! I [idx] Model latitude index
      real(r8),intent(in)::doy                  ! I [day] Day of year [1.0..366.0)
c     real(r8),intent(in)::hbuf(*)              ! I [ptr] History buffer
      real(r8),intent(in)::hgt_mdp(plond)       ! I [m] Midlayer height above surface
      real(r8),intent(in)::lat_rdn              ! I [rdn] Latitude
      real(r8),intent(inout)::obuf(*)           ! I [ptr] Output buffer
      real(r8),intent(in)::oro(plond)           ! I [frc] Orography
      real(r8),intent(in)::prs_dlt(plond)       ! I [Pa] Pressure thickness
      real(r8),intent(in)::prs_mdp(plond)       ! I [Pa] Pressure
      real(r8),intent(in)::prs_sfc(plond)       ! I [Pa] Surface pressure
      real(r8),intent(in)::q_H2O_vpr(plond)     ! I [kg kg-1] Water vapor mixing ratio
      real(r8),intent(in)::snw_hgt_lqd(plond)   ! I [m] Equivalent liquid water snow depth
      real(r8),intent(in)::tm_adj               ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
      real(r8),intent(in)::tpt_mdp(plond)       ! I [K] Temperature
      real(r8),intent(in)::tpt_ptn_mdp(plond)   ! I [K] Midlayer local potential temperature
      real(r8),intent(in)::wnd_mrd_mdp(plond)   ! I [m s-1] Meridional wind component
      real(r8),intent(in)::wnd_znl_mdp(plond)   ! I [m s-1] Zonal wind component
c     Input/Output
      real(r8),intent(inout)::q_dst(plond,dst_nbr) ! I/O [kg kg-1] Dust mixing ratio
c     Local Output
      real(r8) dst_slt_flx_rat_ttl(plond) ! [m-1] Ratio of vertical dust flux to streamwise mass flux
      real(r8) flx_mss_hrz_slt_ttl(plond) ! [kg m-1 s-1] Vertically integrated streamwise mass flux
      real(r8) flx_mss_vrt_dst(plond,dst_nbr) ! [kg m-2 s-1] Vertical mass flux of dust
      real(r8) flx_mss_vrt_dst_ttl(plond) ! [kg m-2 s-1] Total vertical mass flux of dust
      real(r8) frc_thr_rgh_fct(plond) ! [frc] Threshold friction velocity increase from roughness
      real(r8) frc_thr_wet_fct(plond) ! [frc] Threshold friction velocity increase from moisture
      real(r8) hgt_zpd(plond)       ! [m] Zero plane displacement
      real(r8) lnd_frc_mbl(plond)   ! [frc] Bare ground fraction
      real(r8) mno_lng(plond)       ! [m] Monin-Obukhov length
      real(r8) snw_frc(plond)       ! [frc] Fraction of surface covered by snow
      real(r8) trn_fsh_vpr_soi_atm(plond) ! [frc] Transfer efficiency of vapor from soil to atmosphere
      real(r8) wnd_frc(plond)       ! [m s-1] Friction velocity
      real(r8) wnd_frc_slt(plond)   ! [m s-1] Saltating friction velocity
      real(r8) wnd_frc_thr_slt(plond) ! [m s-1] Threshold friction velocity for saltation
      real(r8) wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
      real(r8) wnd_rfr(plond)       ! [m s-1] Wind speed at reference height
      real(r8) wnd_rfr_thr_slt(plond) ! [m s-1] Threshold 10 m wind speed for saltation
c     Local
      logical flg_CaCO3         ! [flg] Activate CaCO3 tracer
      logical flg_mbl(plond)    ! [flg] Mobilization candidates
      character fl_out*80       ! [sng] Name of netCDF output file
      integer i                 ! [idx] Counting index
      integer m                 ! [idx] Counting index
      integer mbl_nbr           ! [nbr] Number of mobilization candidates
      integer sfc_typ(plond)    ! [idx] LSM surface type (0..28)
      real(r8) cnd_trm_soi(plond)   ! [W m-1 K-1] Soil thermal conductivity
      real(r8) dns_mdp(plond)       ! [kg m-3] Midlayer density
      real(r8) flx_LW_dwn_sfc(plond) ! [W m-2] Longwave downwelling flux at surface
      real(r8) flx_SW_abs_sfc(plond) ! [W m-2] Solar flux absorbed by ground
      real(r8) lnd_frc_dry(plond)   ! [frc] Dry land fraction
      real(r8) mbl_bsn_fct(plond)   ! [frc] Topographic source limitation
      real(r8) mss_frc_CaCO3(plond) ! [frc] Mass fraction of CaCO3
      real(r8) lvl_dlt(plond)       ! [m] Soil layer thickness
      real(r8) mpl_air(plond)       ! [kg m-2] Air mass path in layer
      real(r8) mss_frc_cly(plond)   ! [frc] Mass fraction of clay
      real(r8) mss_frc_snd(plond)   ! [frc] Mass fraction of sand
      real(r8) tm_dlt               ! [s] Mobilization timestep
      real(r8) tpt_gnd(plond)       ! [K] Ground temperature
      real(r8) tpt_soi(plond)       ! [K] Soil temperature
      real(r8) tpt_soi_frz          ! [K] Temperature of frozen soil
      real(r8) tpt_vrt_mdp          ! [K] Midlayer virtual temperature
      real(r8) vwc_dry(plond)       ! [m3 m-3] Dry volumetric water content (no E-T)
      real(r8) vwc_opt(plond)       ! [m3 m-3] E-T optimal volumetric water content 
      real(r8) vwc_sat(plond)       ! [m3 m-3] Saturated volumetric water content (sand-dependent)
      real(r8) vwc_sfc(plond)       ! [m3 m-3] Volumetric water content
c     GCM diagnostics
      real(r8) q_dst_tnd_mbl(plond,dst_nbr) ! [kg kg-1 s-1] Dust tendency due to gravitational settling
      real(r8) q_dst_tnd_mbl_ttl(plond) ! [kg kg-1 s-1] Total dust tendency due to gravitational settling
      
c     Main Code
c     Timesplit if desired
      tm_dlt=tm_adj             ! [s] (default CCM: 2*dt, MATCH: dt)
c     Assume water in soil freezes at 0 C
      tpt_soi_frz=tpt_frz_pnt   ! [K] Temperature of frozen soil
      
c     Initialize output fluxes and tendencies
      q_dst_tnd_mbl=0.0         ! [kg kg-1 s-1]
      q_dst_tnd_mbl_ttl=0.0     ! [kg kg-1 s-1]
      flx_mss_vrt_dst=0.0       ! [kg m-2 s-1]
      flx_mss_vrt_dst_ttl=0.0   ! [kg m-2 s-1]
      frc_thr_wet_fct=0.0       ! [frc]
      wnd_rfr=0.0               ! [m s-1]
      wnd_frc=0.0               ! [m s-1]
      wnd_frc_slt=0.0           ! [m s-1]
      wnd_frc_thr_slt=0.0       ! [m s-1]
      wnd_rfr_thr_slt=0.0       ! [m s-1]
      hgt_zpd=hgt_zpd_mbl       ! [m]
c++pjr nats mods
      lnd_frc_mbl =0.0
      vwc_sfc=0.
      tpt_gnd=0.
c--pjr
      
c     Compute necessary derived fields
      do i=1,plon
c     Stop occasional haywire model runs 
         if(tpt_mdp(i) > 350.0) stop 'ERROR: dst_mbl() reports tpt_mdp(i) > 350.0'
c     Midlayer virtual temperature
         tpt_vrt_mdp=tpt_mdp(i)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(i)) ! [K]
c     Density at center of gridbox
         dns_mdp(i)=prs_mdp(i)/(tpt_vrt_mdp*gas_cst_dry_air) ! [kg m-3]
c     Approximate surface virtual temperature (uses midlayer moisture)
c         tpt_vrt_sfc=tpt_sfc(i)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(i)) ! [K]
c     Surface density
c         dns_sfc(i)=prs_sfc(i)/(tpt_vrt_sfc*gas_cst_dry_air) ! [kg m-3]
c     Mass of air currently in gridbox
         mpl_air(i)=prs_dlt(i)*grv_sfc_rcp ! [kg m-2]
c     Mean surface layer horizontal wind speed 
         wnd_mdp(i)=sqrt(wnd_znl_mdp(i)*wnd_znl_mdp(i)+wnd_mrd_mdp(i)*wnd_mrd_mdp(i))
      end do                    ! end loop over lon
      
c     Gather input variables from common blocks
c     Surface type
      call sfc_typ_get(
     $     lat_idx,             ! I [idx] Latitude index
     $     plon,                ! I [frc] Dimension size (normally plon)
     $     sfc_typ)             ! O [idx] LSM surface type (0..28)
c     Soil texture and dry land fraction
      call soi_txt_get(
     $     lat_idx,             ! I [idx] Latitude index
     $     plon,                ! I [frc] Dimension size (normally plon)
     $     lnd_frc_dry,         ! O [frc] Dry land fraction
     $     mbl_bsn_fct,         ! O [frc] Topographic source limitation
     $     mss_frc_CaCO3,       ! O [frc] Mass fraction of CaCO3
     $     mss_frc_cly,         ! O [frc] Mass fraction of clay
     $     mss_frc_snd)         ! O [frc] Mass fraction of sand
c     Soil moisture
      call vwc_sfc_get(
     $     lat_idx,             ! I [idx] Latitude index
     $     plon,                ! I [frc] Dimension size (normally plon)
     $     vwc_sfc)             ! O [m3 m-3] Volumetric water content
c     Surface and soil temperature
      call tpt_gnd_soi_get(
     $     lat_idx,             ! I [idx] Latitude index
     $     plon,                ! I [frc] Dimension size (normally plon)
     $     tpt_gnd,             ! O [K] Ground temperature
     $     tpt_soi)             ! O [K] Soil temperature
c     Radiative fluxes at surface
      call flx_rad_sfc_get(
     $     lat_idx,             ! I [idx] Latitude index
     $     plon,                ! I [frc] Dimension size (normally plon)
     $     flx_LW_dwn_sfc,      ! O [W m-2] Longwave downwelling flux at surface
     $     flx_SW_abs_sfc)      ! O [W m-2] Solar flux absorbed by ground
      
c     All variables from common blocks have been assembled
c     Use these to compute derived fields
c     Fraction of surface covered by snow
      call snw_frc_get(
     $     snw_hgt_lqd,         ! I [m] Equivalent liquid water snow depth
     $     snw_frc)             ! O [frc] Fraction of surface covered by snow
      
c     Fraction of each gridcell suitable for dust mobilization
      call lnd_frc_mbl_get(
     $     doy,                 ! I [day] Day of year [1.0..366.0)
     $     flg_mbl,             ! O [flg] Mobilization candidate flag
     $     lat_rdn,             ! I [rdn] Latitude
     $     lnd_frc_dry,         ! I [frc] Dry land fraction
     $     lnd_frc_mbl,         ! O [frc] Bare ground fraction
     $     mbl_nbr,             ! O [flg] Number of mobilization candidates
     $     oro,                 ! I [frc] Orography
     $     sfc_typ,             ! I [idx] LSM surface type (0..28)
c-mahow     $     snw_frc)             ! I [frc] Fraction of surface covered by snow
c++mahow
c     $     snw_frc              ! I [frc] Fraction of surface covered by snow
c     ,   , lat_idx)
     $     snw_frc,              ! I [frc] Fraction of surface covered by snow
     $     lat_idx)
c--mahow
c     Much ado about nothing
      if (mbl_nbr == 0) goto 737

c     Hydrologic properties
      call hyd_prp_get(         ! NB: These properties are time-invariant
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     mss_frc_cly,         ! I [frc] Mass fraction clay 
     $     mss_frc_snd,         ! I [frc] Mass fraction sand
     $     vwc_dry,             ! O [m3 m-3] Dry volumetric water content (no E-T)
     $     vwc_opt,             ! O [m3 m-3] E-T optimal volumetric water content
     $     vwc_sat)             ! O [m3 m-3] Saturated volumetric water content
      
c     Soil thermal conductivity and layer thickness
      call cnd_trm_soi_get(
     $     cnd_trm_soi,         ! O [W m-1 K-1] Soil thermal conductivity
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     lvl_dlt,             ! O [m] Soil layer thickness
     $     mss_frc_cly,         ! I [frc] Mass fraction clay 
     $     mss_frc_snd,         ! I [frc] Mass fraction sand
     $     tpt_soi,             ! I [K] Soil temperature
     $     vwc_dry,             ! I [m3 m-3] Dry volumetric water content (no E-T)
     $     vwc_opt,             ! I [m3 m-3] E-T optimal volumetric water content
     $     vwc_sat,             ! I [m3 m-3] Saturated volumetric water content
     $     vwc_sfc)             ! I [m3 m-3] Volumetric water content
      
c     Transfer efficiency of vapor from soil to atmosphere
      call trn_fsh_vpr_soi_atm_get(
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     tpt_soi,             ! I [K] Soil temperature
     $     tpt_soi_frz,         ! I [K] Temperature of frozen soil
     $     trn_fsh_vpr_soi_atm, ! O [frc] Transfer efficiency of vapor from soil to atmosphere
     $     vwc_dry,             ! I [m3 m-3] Dry volumetric water content (no E-T)
     $     vwc_opt,             ! I [m3 m-3] E-T optimal volumetric water content
     $     vwc_sfc)             ! I [m3 m-3] Volumetric water content
      
c     Surface exchange properties over erodible surfaces
      call blm_mbl(
     $     cnd_trm_soi,         ! I [W m-1 K-1] Soil thermal conductivity
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     flx_LW_dwn_sfc,      ! I [W m-2] Longwave downwelling flux at surface
     $     flx_SW_abs_sfc,      ! I [W m-2] Solar flux absorbed by ground
     $     hgt_mdp,             ! I [m] Midlayer height above surface
     $     hgt_zpd_mbl,         ! I [m] Zero plane displacement
     $     lvl_dlt,             ! I [m] Soil layer thickness
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_vpr,           ! I [kg kg-1] Specific humidity
     $     rgh_mmn_mbl,         ! I [m] Roughness length momentum
     $     tpt_mdp,             ! I [K] Midlayer temperature
     $     tpt_ptn_mdp,         ! I [K] Midlayer local potential temperature
     $     tpt_soi,             ! I [K] Soil temperature
     $     trn_fsh_vpr_soi_atm, ! I [frc] Transfer efficiency of vapor from soil to atmosphere
     $     wnd_mrd_mdp,         ! I [m s-1] Surface layer meridional wind speed
     $     wnd_znl_mdp,         ! I [m s-1] Surface layer zonal wind speed
     $     mno_lng,             ! O [m] Monin-Obukhov length
     $     tpt_gnd,             ! I/O [K] Ground temperature
     $     wnd_frc)             ! O [m s-1] Surface friction velocity
      
c     Interpolate midlayer wind speed to 10 m
      call wnd_rfr_get(
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     hgt_mdp,             ! I [m] Midpoint height above surface
     $     hgt_rfr,             ! I [m] Reference height for mobilization processes
     $     hgt_zpd,             ! I [m] Zero plane displacement
     $     mno_lng,             ! I [m] Monin-Obukhov length
     $     wnd_frc,             ! I [m s-1] Surface friction velocity
     $     wnd_mdp,             ! I [m s-1] Surface layer mean wind speed
     $     wnd_min_mbl,         ! I [m s-1] Minimum windspeed
     $     wnd_rfr)             ! O [m s-1] Wind speed at reference height
      
c     Factor by which surface roughness increases threshold friction velocity 
      call frc_thr_rgh_fct_get(
     $     frc_thr_rgh_fct,     ! O [frc] Factor by which surface roughness increases threshold friction velocity
     $     rgh_mmn_mbl,         ! I [m] Roughness length momentum for erodible surfaces
     $     rgh_mmn_smt)         ! I [m] Smooth roughness length
      
c     Factor by which soil moisture increases threshold friction velocity 
      call frc_thr_wet_fct_get(
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     frc_thr_wet_fct,     ! O [frc] Factor by which moisture increases threshold friction velocity
     $     mss_frc_cly,         ! I [frc] Mass fraction of clay
     $     vwc_sfc)             ! I [frc] Volumetric water content
      
c     Threshold friction velocity for saltation over dry, bare, smooth ground
c     fxm: Use surface density not midlayer density
      call wnd_frc_thr_slt_get(
     $     dmt_vwr,             ! I [m] Mass weighted diameter resolved 
     $     dns_aer,             ! I [kg m-3] Particle density
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     wnd_frc_thr_slt)     ! O [m s-1] Threshold friction velocity for saltation
      
c     Adjust threshold friction velocity to account for moisture and roughness
      do i=1,plon
         wnd_frc_thr_slt(i)=    ! [m s-1] Threshold friction velocity for saltation
     $        wnd_frc_thr_slt(i)* ! [m s-1] Threshold for dry, flat ground
     $        frc_thr_wet_fct(i)* ! [frc] Adjustment for moisture
     $        frc_thr_rgh_fct(i) ! [frc] Adjustment for roughness
      end do                    ! end loop over lon
      
c     Threshold saltation wind speed
      do i=1,plon
         if (flg_mbl(i)) then
            wnd_rfr_thr_slt(i)= ! [m s-1] Threshold 10 m wind speed for saltation
     $           wnd_rfr(i)*wnd_frc_thr_slt(i)/wnd_frc(i)
         endif                  ! endif flg_mbl
      end do                    ! end loop over lon
      
c     Saltation increases friction speed by roughening surface
      call wnd_frc_slt_get(
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     wnd_frc,             ! I [m s-1] Surface friction velocity
     $     wnd_frc_slt,         ! O [m s-1] Saltating friction velocity
     $     wnd_rfr,             ! I [m s-1] Wind speed at reference height
     $     wnd_rfr_thr_slt)     ! I [m s-1] Threshold 10 m wind speed for saltation
      
c     Horizontal streamwise mass flux
      call flx_mss_hrz_slt_ttl_Whi79_get(
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     flx_mss_hrz_slt_ttl, ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
     $     wnd_frc_slt,         ! I [m s-1] Saltating friction velocity
     $     wnd_frc_thr_slt)     ! I [m s-1] Threshold friction speed for saltation
      
c     Apply land surface and vegetation limitations and global tuning factor
      do i=1,plon
         flx_mss_hrz_slt_ttl(i)=flx_mss_hrz_slt_ttl(i) ! [kg m-2 s-1]
     $        *lnd_frc_mbl(i)   ! [frc] Bare ground fraction
c++pjr nats mods for ginoux basins
c     $        *mbl_bsn_fct(i)   ! [frc] Topographic source limitation
c--pjr
     $        *flx_mss_fdg_fct  ! [frc] Global mass flux tuning factor (empirical)
c++mahow
!     &        *1.72          ! tuning tot source to 1500 mt/yr
c++pjr nats mods for ginoux basins
!     &        *2.02          ! tuning according to ginox basin factor
!     &        /1.2377        ! 5/16/01
c--pjr
c--mahow
      end do                    ! end loop over lon
      
c     Vertical dust mass flux
      call flx_mss_vrt_dst_ttl_MaB95_get(
     $     dst_slt_flx_rat_ttl, ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     flx_mss_hrz_slt_ttl, ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
     $     flx_mss_vrt_dst_ttl, ! O [kg m-2 s-1] Total vertical mass flux of dust
     $     mss_frc_cly)         ! I [frc] Mass fraction clay 
      
c     Partition vertical mass flux into transport bins
      call flx_mss_vrt_dst_prt(
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     flx_mss_vrt_dst,     ! O [kg m-2 s-1] Vertical mass flux of dust
     $     flx_mss_vrt_dst_ttl) ! I [kg m-2 s-1] Total vertical mass flux of dust
      
c     Mask dust mass flux by tracer mass fraction at source
      flg_CaCO3=.false.         ! [flg] Activate CaCO3 tracer
      if (flg_CaCO3) then
         call flx_mss_CaCO3_msk(
     $        dmt_vwr,          ! I [m] Mass weighted diameter resolved 
     $        flg_mbl,          ! I [flg] Mobilization candidate flag
     $        flx_mss_vrt_dst,  ! I/O [kg m-2 s-1] Vertical mass flux of dust
     $        mss_frc_CaCO3,    ! I [frc] Mass fraction of CaCO3
     $        mss_frc_cly,      ! I [frc] Mass fraction of clay
     $        mss_frc_snd)      ! I [frc] Mass fraction of sand
      endif                     ! endif flg_CaCO3

c     Fluxes are known, so adjust mixing ratios
      do i=1,plon               ! NB: Inefficient loop order
         if (flg_mbl(i)) then
            do m=1,dst_nbr
#ifdef DST_DBG
               if (q_dst(i,m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,1p,e8.1,0p,a)') 
     $              'dst_mbl:b lat= ',lat_idx,' q_dst(',i,',',m,') = ',q_dst(i,m),' kg kg-1'
#endif /* not DST_DBG */
               q_dst_tnd_mbl(i,m)= ! [kg kg-1 s-1]
     $              flx_mss_vrt_dst(i,m)/mpl_air(i)
               q_dst(i,m)=q_dst(i,m)+ ! [kg kg-1]
     $              tm_adj*q_dst_tnd_mbl(i,m)
#ifdef DST_DBG
               if (q_dst(i,m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,1p,e8.1,0p,a)') 
     $              'dst_mbl:a lat= ',lat_idx,' q_dst(',i,',',m,') = ',q_dst(i,m),' kg kg-1'
#endif /* not DST_DBG */
            end do              ! end loop over cst
#ifdef DST_DBG
            if (flx_mss_vrt_dst_ttl(i) > flx_mss_mxm) write(6,'(a,i2,a,i3,a,1p,e8.1,0p,a)') 
     $           'dst_mbl: lat = ',lat_idx,' flx_mss_vrt_dst_ttl(',i,') = ',flx_mss_vrt_dst_ttl(i),' kg m-2 s-1'
#endif /* not DST_DBG */
         endif                  ! endif flg_mbl
      end do                    ! end loop over lon
      
c     Integrate over all size categories for diagnostic output
      call dst_add_lon(q_dst_tnd_mbl,q_dst_tnd_mbl_ttl)
c     Recompute flx_mss_vrt_dst_ttl as sum of transported dust
      call dst_add_lon(flx_mss_vrt_dst, ! I [kg m-2 s-1] Vertical mass flux of dust
     $     flx_mss_vrt_dst_ttl) ! O [kg m-2 s-1] Total vertical mass flux of dust
      
c     Jump to here when no points are mobilization candidates
 737  continue

#ifdef DST_DBG
      if (lat_idx == lat_dbg.and..false.) then
c      if (lat_idx == lat_dbg) then
         write (6,'(a,1(a,f9.5,a),17(a,1p,e9.2,0p,a))') 
     $        'dst: Diagnostics at (lat_dbg,lon_dbg) set in dstcmnini() :',
     $        'doy = ',doy,' day, ',
     $        'flx_LW_dwn_sfc = ',flx_LW_dwn_sfc(lon_dbg),' W m-2, ',
     $        'flx_SW_abs_sfc = ',flx_SW_abs_sfc(lon_dbg),' W m-2, ',
     $        'dst_slt_flx_rat_ttl = ',dst_slt_flx_rat_ttl(lon_dbg),' kg m-1 s-1, ',
     $        'flx_mss_hrz_slt_ttl = ',flx_mss_hrz_slt_ttl(lon_dbg),' kg m-1 s-1, ',
     $        'frc_thr_rgh_fct = ',frc_thr_rgh_fct(lon_dbg),' frc, ',
     $        'frc_thr_wet_fct = ',frc_thr_wet_fct(lon_dbg),' frc, ',
     $        'lnd_frc_mbl = ',lnd_frc_mbl(lon_dbg),' frc, ',
     $        'mbl_bsn_fct = ',mbl_bsn_fct(lon_dbg),' frc, ',
     $        'mss_frc_CaCO3 = ',mss_frc_CaCO3(lon_dbg),' frc, ',
     $        'mno_lng = ',mno_lng(lon_dbg),' m, ',
     $        'rgh_mmn = ',rgh_mmn_mbl,' m, ',
     $        'tpt_gnd = ',tpt_gnd(lon_dbg),' K, ',
     $        'tpt_ptn_mdp = ',tpt_ptn_mdp(lon_dbg),' K, ',
     $        'tpt_soi = ',tpt_soi(lon_dbg),' K, ',
     $        'vwc_sfc = ',vwc_sfc(lon_dbg),' m3 m-3, ',
     $        'wnd_frc = ',wnd_frc(lon_dbg),' m s-1, ',
     $        'wnd_frc_slt = ',wnd_frc_slt(lon_dbg),' m s-1, ',
     $        'wnd_frc_thr_slt = ',wnd_frc_thr_slt(lon_dbg),' m s-1, '
      endif                     ! endif dbg
#endif /* not DST_DBG */ 

#ifdef BXM
c     Update netCDF file
      fl_out='aer.nc'           ! [sng] Name of netCDF output file
      call strnul(fl_out)
      call mbl2nc(            
     $     dst_slt_flx_rat_ttl, ! I [m-1] Ratio of vertical dust flux to streamwise mass flux
     $     fl_out,              ! I [sng] Name of netCDF output file
     $     flx_mss_hrz_slt_ttl, ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
     $     flx_mss_vrt_dst,     ! I [kg m-2 s-1] Vertical mass flux of dust
     $     flx_mss_vrt_dst_ttl, ! I [kg m-2 s-1] Total vertical mass flux of dust
     $     frc_thr_rgh_fct,     ! I [frc] Threshold friction velocity increase from roughness
     $     frc_thr_wet_fct,     ! I [frc] Threshold friction velocity increase from moisture
     $     hgt_zpd_mbl,         ! I [m] Zero plane displacement height
     $     lat_idx,             ! I [idx] Model latitude index
     $     lnd_frc_mbl,         ! I [frc] Bare ground fraction
     $     mno_lng,             ! I [m] Monin-Obukhov length
     $     rgh_mmn_mbl,         ! I [m] Roughness length momentum
     $     snw_frc,             ! I [frc] Fraction of surface covered by snow
     $     wnd_frc,             ! I [m s-1] Friction velocity
     $     wnd_frc_slt,         ! I [m s-1] Saltating friction velocity
     $     wnd_frc_thr_slt,     ! I [m s-1] Threshold friction velocity for saltation
     $     wnd_mdp,             ! I [m s-1] Surface layer mean wind speed
     $     wnd_rfr,             ! I [m s-1] Wind speed at reference height
     $     wnd_rfr_thr_slt)     ! I [m s-1] Threshold 10 m wind speed for saltation
#endif /* not BXM */
      
#ifdef CCM
      call outfld('FRC_WET',frc_thr_wet_fct,plond,lat_idx,hbuf)
      call outfld('TPT_GND',tpt_gnd,plond,lat_idx,hbuf)
      call outfld('VWC_SFC',vwc_sfc,plond,lat_idx,hbuf)
      call outfld('WND_FRC',wnd_frc,plond,lat_idx,hbuf)
      call outfld('WND_FRCS',wnd_frc_slt,plond,lat_idx,hbuf)
      call outfld('WND_FRCT',wnd_frc_thr_slt,plond,lat_idx,hbuf)
      call outfld('WND_RFR',wnd_rfr,plond,lat_idx,hbuf)
      call outfld('WND_RFRT',wnd_rfr_thr_slt,plond,lat_idx,hbuf)
      call outfld('DSTSFMBL',flx_mss_vrt_dst_ttl,plond,lat_idx,hbuf)
      do m=1,dst_nbr
         call outfld(flx_mss_mbl_sfc_nm(m),flx_mss_vrt_dst(1,m),plond,lat_idx,hbuf)
      end do                    ! end loop over cst
#else /* not CCM */
      call outfld('FRC_WET',frc_thr_wet_fct,plond,lat_idx,obuf)
      call outfld('TPT_GND',tpt_gnd,plond,lat_idx,obuf)
      call outfld('VWC_SFC',vwc_sfc,plond,lat_idx,obuf)
      call outfld('WND_FRC',wnd_frc,plond,lat_idx,obuf)
      call outfld('WND_FRCS',wnd_frc_slt,plond,lat_idx,obuf)
      call outfld('WND_FRCT',wnd_frc_thr_slt,plond,lat_idx,obuf)
      call outfld('WND_RFR',wnd_rfr,plond,lat_idx,obuf)
      call outfld('WND_RFRT',wnd_rfr_thr_slt,plond,lat_idx,obuf)
      call outfld('DSTSFMBL',flx_mss_vrt_dst_ttl,plond,lat_idx,obuf)
      do m=1,dst_nbr
         call outfld(flx_mss_mbl_sfc_nm(m),flx_mss_vrt_dst(1,m),plond,lat_idx,obuf)
      end do                    ! end loop over cst
#endif /* not CCM */
#ifdef DST_MSS_BDG
      call bdg_aal('dst_sf_mbl',lat_idx,flx_mss_vrt_dst_ttl)
c     NB: Only call to bdg_gam_wet_2D is made here
      call bdg_gam_wet_2D('dst_ss_mbl',lat_idx,prs_dlt,q_dst_tnd_mbl_ttl)
#endif /* not DST_MSS_BDG */
      return
      end subroutine dst_mbl                      ! end dst_mbl()
      
      subroutine mbl2nc(            
     $     dst_slt_flx_rat_ttl, ! I [m-1] Ratio of vertical dust flux to streamwise mass flux
     $     fl_out,              ! I [sng] Name of netCDF output file
     $     flx_mss_hrz_slt_ttl, ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
     $     flx_mss_vrt_dst,     ! I [kg m-2 s-1] Vertical mass flux of dust
     $     flx_mss_vrt_dst_ttl, ! I [kg m-2 s-1] Total vertical mass flux of dust
     $     frc_thr_rgh_fct,     ! I [frc] Threshold friction velocity increase from roughness
     $     frc_thr_wet_fct,     ! I [frc] Threshold friction velocity increase from moisture
     $     hgt_zpd,             ! I [m] Zero plane displacement height
     $     lat_idx,             ! I [idx] Model latitude index
     $     lnd_frc_mbl,         ! I [frc] Bare ground fraction
     $     mno_lng,             ! I [m] Monin-Obukhov length
     $     rgh_mmn,             ! I [m] Roughness length momentum
     $     snw_frc,             ! I [frc] Fraction of surface covered by snow
     $     wnd_frc,             ! I [m s-1] Friction velocity
     $     wnd_frc_slt,         ! I [m s-1] Saltating friction velocity
     $     wnd_frc_thr_slt,     ! I [m s-1] Threshold friction velocity for saltation
     $     wnd_mdp,             ! I [m s-1] Surface layer mean wind speed
     $     wnd_rfr,             ! I [m s-1] Wind speed at reference height
     $     wnd_rfr_thr_slt)     ! I [m s-1] Threshold 10 m wind speed for saltation
c     Purpose: Output aerosol mobilization properties to netCDF file
      use netcdf
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstctl.h> /* Model control variables */
#include <dbg.com> /* Debugging constants */
#include <netcdf.inc>
c     Input
      character*(*) fl_out      ! [sng] Name of netCDF output file
      integer    lat_idx        ! I [idx] Model latitude index
      real(r8) dst_slt_flx_rat_ttl(plond) ! [m-1] Ratio of vertical dust flux to streamwise mass flux
      real(r8) flx_mss_hrz_slt_ttl(plond) ! [kg m-1 s-1] Vertically integrated streamwise mass flux
      real(r8) flx_mss_vrt_dst(plond,dst_nbr) ! [kg m-2 s-1] Vertical mass flux of dust
      real(r8) flx_mss_vrt_dst_ttl(plond) ! [kg m-2 s-1] Total vertical mass flux of dust
      real(r8) frc_thr_rgh_fct(plond) ! [frc] Threshold friction velocity increase from roughness
      real(r8) frc_thr_wet_fct(plond) ! [frc] Threshold friction velocity increase from moisture
      real(r8) hgt_zpd              ! [m] Zero plane displacement
      real(r8) lnd_frc_mbl(plond)   ! [frc] Bare ground fraction
      real(r8) mno_lng(plond)       ! [m] Monin-Obukhov length
      real(r8) rgh_mmn              ! [m] Roughness length momentum
      real(r8) snw_frc(plond)       ! [frc] Fraction of surface covered by snow
      real(r8) wnd_frc(plond)       ! [m s-1] Friction velocity
      real(r8) wnd_frc_slt(plond)   ! [m s-1] Saltating friction velocity
      real(r8) wnd_frc_thr_slt(plond) ! [s m-1] Threshold friction velocity for saltation
      real(r8) wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
      real(r8) wnd_rfr(plond)       ! [m s-1] Wind speed at reference height
      real(r8) wnd_rfr_thr_slt(plond) ! [m s-1] Threshold 10 m wind speed for saltation
c     Output
c     Local
c     File metadata and dimension IDs
      integer dim_lon_sz_time(3) ! [enm] Dimension IDs
      integer dim_lon_time(2)   ! [enm] Dimension IDs
      integer srt_lon_time(2)   ! Starting index array
      integer cnt_lon_time(2)   ! Count array
      integer srt_lon_sz_time(3) ! Starting index array
      integer cnt_lon_sz_time(3) ! Count array
      integer fll_mode_old      ! Old fill mode
      integer lon_dim_id        ! [enm] Dimension ID for lon
      integer time_dim_id       ! [enm] Dimension ID for time
      integer nc_id             ! File handle
      integer rcd               ! [rcd] Return success code
      integer sz_dim_id         ! [enm] Dimension ID for sz
      integer sz_grd_dim_id     ! [enm] Dimension ID for sz grid
c     Variable IDs
      integer dst_slt_flx_rat_ttl_id ! [enm] Variable ID
      integer flx_mss_hrz_slt_ttl_id ! [enm] Variable ID
      integer flx_mss_vrt_dst_id ! [enm] Variable ID
      integer flx_mss_vrt_dst_ttl_id ! [enm] Variable ID
      integer frc_thr_rgh_fct_id ! [enm] Variable ID
      integer frc_thr_wet_fct_id ! [enm] Variable ID
      integer hgt_zpd_id        ! [enm] Variable ID
      integer lnd_frc_mbl_id    ! [enm] Variable ID
      integer mno_lng_id        ! [enm] Variable ID
      integer rgh_mmn_id        ! [enm] Variable ID
      integer snw_frc_id        ! [enm] Variable ID
      integer wnd_frc_id        ! [enm] Variable ID
      integer wnd_frc_slt_id    ! [enm] Variable ID
      integer wnd_frc_thr_slt_id ! [enm] Variable ID
      integer wnd_mdp_id        ! [enm] Variable ID
      integer wnd_rfr_id        ! [enm] Variable ID
      integer wnd_rfr_thr_slt_id ! [enm] Variable ID
c     Variable data
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
      
c     Begin netCDF output routines
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf90_open(fl_out,NF_WRITE,nc_id)
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR opening ',fl_out(1:f_strlen(fl_out)),
     $     ' in mbl2nc(), rcd = ',rcd
      rcd=rcd+nf90_redef(nc_id)
      rcd=rcd+nf90_inq_dimid(nc_id,'lon',lon_dim_id)
      rcd=rcd+nf90_inq_dimid(nc_id,'time',time_dim_id)
      rcd=rcd+nf90_inq_dimid(nc_id,'sz',sz_dim_id)
c     Add global attributes
c     Define dimension IDs
c     Assemble ID and count vectors for each multidimensional combination of dimensions
      dim_lon_time=(/lon_dim_id,time_dim_id/)
      srt_lon_time=(/1,nstep/)   ! Starting index array
      cnt_lon_time=(/plon,1/)    ! Count array

      dim_lon_sz_time=(/lon_dim_id,sz_dim_id,time_dim_id/)
      srt_lon_sz_time=(/1,1,nstep/)
      cnt_lon_sz_time=(/plon,dst_nbr,1/)

      if (nstep == 1) then
c     Variable definitions
         rcd=rcd+nf90_def_var(nc_id,'frc_thr_rgh_fct',NF_DOUBLE,dim_lon_time,frc_thr_rgh_fct_id)
         rcd=rcd+nf90_def_var(nc_id,'dst_slt_flx_rat_ttl',NF_DOUBLE,dim_lon_time,dst_slt_flx_rat_ttl_id)
         rcd=rcd+nf90_def_var(nc_id,'flx_mss_hrz_slt_ttl',NF_DOUBLE,dim_lon_time,flx_mss_hrz_slt_ttl_id)
         rcd=rcd+nf90_def_var(nc_id,'flx_mss_vrt_dst',NF_DOUBLE,dim_lon_sz_time,flx_mss_vrt_dst_id)
         rcd=rcd+nf90_def_var(nc_id,'flx_mss_vrt_dst_ttl',NF_DOUBLE,dim_lon_time,flx_mss_vrt_dst_ttl_id)
         rcd=rcd+nf90_def_var(nc_id,'frc_thr_wet_fct',NF_DOUBLE,dim_lon_time,frc_thr_wet_fct_id)
         rcd=rcd+nf90_def_var(nc_id,'hgt_zpd_mbl',NF_DOUBLE,0,hgt_zpd_id)
         rcd=rcd+nf90_def_var(nc_id,'lnd_frc_mbl',NF_DOUBLE,dim_lon_time,lnd_frc_mbl_id)
         rcd=rcd+nf90_def_var(nc_id,'mno_lng_mbl',NF_DOUBLE,dim_lon_time,mno_lng_id)
         rcd=rcd+nf90_def_var(nc_id,'rgh_mmn_mbl',NF_DOUBLE,0,rgh_mmn_id)
         rcd=rcd+nf90_def_var(nc_id,'snw_frc',NF_DOUBLE,dim_lon_time,snw_frc_id)
         rcd=rcd+nf90_def_var(nc_id,'wnd_frc_mbl',NF_DOUBLE,dim_lon_time,wnd_frc_id)
         rcd=rcd+nf90_def_var(nc_id,'wnd_frc_slt',NF_DOUBLE,dim_lon_time,wnd_frc_slt_id)
         rcd=rcd+nf90_def_var(nc_id,'wnd_frc_thr_slt',NF_DOUBLE,dim_lon_time,wnd_frc_thr_slt_id)
         rcd=rcd+nf90_def_var(nc_id,'wnd_mdp',NF_DOUBLE,dim_lon_time,wnd_mdp_id)
         rcd=rcd+nf90_def_var(nc_id,'wnd_rfr_mbl',NF_DOUBLE,dim_lon_time,wnd_rfr_id)
         rcd=rcd+nf90_def_var(nc_id,'wnd_rfr_thr_slt',NF_DOUBLE,dim_lon_time,wnd_rfr_thr_slt_id)
c     Add english text descriptions
         rcd=rcd+nf90_put_att(nc_id,dst_slt_flx_rat_ttl_id,'long_name','Ratio of vertical dust flux to streamwise mass flux')
         rcd=rcd+nf90_put_att(nc_id,flx_mss_hrz_slt_ttl_id,'long_name','Vertically integrated streamwise mass flux')
         rcd=rcd+nf90_put_att(nc_id,flx_mss_vrt_dst_id,'long_name','Vertical mass flux of dust')
         rcd=rcd+nf90_put_att(nc_id,flx_mss_vrt_dst_ttl_id,'long_name','Total vertical mass flux of dust')
         rcd=rcd+nf90_put_att(nc_id,frc_thr_rgh_fct_id,'long_name','Threshold friction velocity increase from roughness')
         rcd=rcd+nf90_put_att(nc_id,frc_thr_wet_fct_id,'long_name','Threshold friction velocity increase from moisture')
         rcd=rcd+nf90_put_att(nc_id,hgt_zpd_id,'long_name','Zero plane displacement height')
         rcd=rcd+nf90_put_att(nc_id,lnd_frc_mbl_id,'long_name','Bare ground fraction')
         rcd=rcd+nf90_put_att(nc_id,mno_lng_id,'long_name','Monin-Obukhov length')
         rcd=rcd+nf90_put_att(nc_id,rgh_mmn_id,'long_name','Roughness length momentum')
         rcd=rcd+nf90_put_att(nc_id,snw_frc_id,'long_name','Fraction of surface covered by snow')
         rcd=rcd+nf90_put_att(nc_id,wnd_frc_id,'long_name','Friction velocity')
         rcd=rcd+nf90_put_att(nc_id,wnd_frc_slt_id,'long_name','Saltating friction velocity')
         rcd=rcd+nf90_put_att(nc_id,wnd_frc_thr_slt_id,'long_name','Threshold friction velocity for saltation')
         rcd=rcd+nf90_put_att(nc_id,wnd_mdp_id,'long_name','Surface layer mean wind speed')
         rcd=rcd+nf90_put_att(nc_id,wnd_rfr_id,'long_name','Wind speed at reference height')
         rcd=rcd+nf90_put_att(nc_id,wnd_rfr_thr_slt_id,'long_name','Threshold 10 m wind speed for saltation')
c     Add units
         rcd=rcd+nf90_put_att(nc_id,dst_slt_flx_rat_ttl_id,'units','meter-1')
         rcd=rcd+nf90_put_att(nc_id,flx_mss_hrz_slt_ttl_id,'units','kilogram meter-1 second-1')
         rcd=rcd+nf90_put_att(nc_id,flx_mss_vrt_dst_id,'units','kilogram meter-2 second-1')
         rcd=rcd+nf90_put_att(nc_id,flx_mss_vrt_dst_ttl_id,'units','kilogram meter-2 second-1')
         rcd=rcd+nf90_put_att(nc_id,frc_thr_rgh_fct_id,'units','fraction')
         rcd=rcd+nf90_put_att(nc_id,frc_thr_wet_fct_id,'units','fraction')
         rcd=rcd+nf90_put_att(nc_id,hgt_zpd_id,'units','meter')
         rcd=rcd+nf90_put_att(nc_id,lnd_frc_mbl_id,'units','fraction')
         rcd=rcd+nf90_put_att(nc_id,mno_lng_id,'units','meter')
         rcd=rcd+nf90_put_att(nc_id,rgh_mmn_id,'units','meter')
         rcd=rcd+nf90_put_att(nc_id,snw_frc_id,'units','fraction')
         rcd=rcd+nf90_put_att(nc_id,wnd_frc_id,'units','meter second-1')
         rcd=rcd+nf90_put_att(nc_id,wnd_frc_slt_id,'units','meter second-1')
         rcd=rcd+nf90_put_att(nc_id,wnd_frc_thr_slt_id,'units','meter second-1')
         rcd=rcd+nf90_put_att(nc_id,wnd_mdp_id,'units','meter second-1')
         rcd=rcd+nf90_put_att(nc_id,wnd_rfr_id,'units','meter second-1')
         rcd=rcd+nf90_put_att(nc_id,wnd_rfr_thr_slt_id,'units','meter second-1')
      else                      ! endif nstep == 1
c     Get variable IDs
         rcd=rcd+nf90_inq_varid(nc_id,'dst_slt_flx_rat_ttl',dst_slt_flx_rat_ttl_id)
         rcd=rcd+nf90_inq_varid(nc_id,'flx_mss_hrz_slt_ttl',flx_mss_hrz_slt_ttl_id)
         rcd=rcd+nf90_inq_varid(nc_id,'flx_mss_vrt_dst',flx_mss_vrt_dst_id)
         rcd=rcd+nf90_inq_varid(nc_id,'flx_mss_vrt_dst_ttl',flx_mss_vrt_dst_ttl_id)
         rcd=rcd+nf90_inq_varid(nc_id,'frc_thr_rgh_fct',frc_thr_rgh_fct_id)
         rcd=rcd+nf90_inq_varid(nc_id,'frc_thr_wet_fct',frc_thr_wet_fct_id)
         rcd=rcd+nf90_inq_varid(nc_id,'hgt_zpd_mbl',hgt_zpd_id)
         rcd=rcd+nf90_inq_varid(nc_id,'lnd_frc_mbl',lnd_frc_mbl_id)
         rcd=rcd+nf90_inq_varid(nc_id,'mno_lng_mbl',mno_lng_id)
         rcd=rcd+nf90_inq_varid(nc_id,'rgh_mmn_mbl',rgh_mmn_id)
         rcd=rcd+nf90_inq_varid(nc_id,'snw_frc',snw_frc_id)
         rcd=rcd+nf90_inq_varid(nc_id,'wnd_frc_mbl',wnd_frc_id)
         rcd=rcd+nf90_inq_varid(nc_id,'wnd_frc_slt',wnd_frc_slt_id)
         rcd=rcd+nf90_inq_varid(nc_id,'wnd_frc_thr_slt',wnd_frc_thr_slt_id)
         rcd=rcd+nf90_inq_varid(nc_id,'wnd_mdp',wnd_mdp_id)
         rcd=rcd+nf90_inq_varid(nc_id,'wnd_rfr_mbl',wnd_rfr_id)
         rcd=rcd+nf90_inq_varid(nc_id,'wnd_rfr_thr_slt',wnd_rfr_thr_slt_id)
      endif                     ! endif nstep /= 1
c     End define mode
      rcd=rcd+nf90_enddef(nc_id)
c     Write data
      if (nstep == 1) then
         rcd=rcd+nf90_put_var(nc_id,hgt_zpd_id,hgt_zpd)
         rcd=rcd+nf90_put_var(nc_id,rgh_mmn_id,rgh_mmn)
      endif                     ! endif nstep == 1
      rcd=rcd+nf90_put_var(nc_id, dst_slt_flx_rat_ttl_id, dst_slt_flx_rat_ttl, srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, flx_mss_hrz_slt_ttl_id, flx_mss_hrz_slt_ttl, srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, flx_mss_vrt_dst_id,     flx_mss_vrt_dst,     srt_lon_sz_time, cnt_lon_sz_time)
      rcd=rcd+nf90_put_var(nc_id, flx_mss_vrt_dst_ttl_id, flx_mss_vrt_dst_ttl, srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, frc_thr_rgh_fct_id,     frc_thr_rgh_fct,     srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, frc_thr_wet_fct_id,     frc_thr_wet_fct,     srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, lnd_frc_mbl_id,         lnd_frc_mbl,         srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, mno_lng_id,             mno_lng,             srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, snw_frc_id,             snw_frc,             srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, wnd_frc_id,             wnd_frc,             srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, wnd_frc_slt_id,         wnd_frc_slt,         srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, wnd_frc_thr_slt_id,     wnd_frc_thr_slt,     srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, wnd_mdp_id,             wnd_mdp,             srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, wnd_rfr_id,             wnd_rfr,             srt_lon_time,    cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id, wnd_rfr_thr_slt_id,     wnd_rfr_thr_slt,     srt_lon_time,    cnt_lon_time)
c     Close output file
      rcd=rcd+nf90_close(nc_id)
      if (nstep == 1) then
         write (6,'(a,a43,a)') prg_nm(1:f_strlen(prg_nm)),': Initialized mobilization data archive in ',fl_out(1:f_strlen(fl_out))
      endif                     ! endif
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR writing ',fl_out(1:f_strlen(fl_out)),
     $     ' in mbl2nc(), rcd = ',rcd
      return 
      end subroutine mbl2nc                      ! end mbl2nc()
      
      end module dstmbl
