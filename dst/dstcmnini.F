c     $Header$
      
c     Purpose: dstcmnini.F contains routines to initialize the following
c     common blocks for the dust model: 
c     dstblm.h, dstchm.h, dstcst.h, dstmsc.h, dstodx.h

c     dst.h needed for DST_CHM token
#include <dst.h> /* Dust preprocessor tokens */

      subroutine dst_cst_cmn_ini(gas_cst_dry_airx,grv_sfcx,rds_Earthx)
c     Initialize time-invariant physical constants common block dstcstcmn
c     dst_cst_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
c     CCM:physics/inti() passes values from CCM:physics/comcon.h which were set by CCM:eul/initcom()
      use precision
      implicit none
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstblm.h> /* Physical constants for boundary layer meteorology */
c     Parameters for source and sink processes
      real(r8),parameter::gas_cst_unvx=8.31441 ! [J mol-1 K-1] Universal gas constant
      real(r8),parameter::mmw_H2Ox=1.8015259e-02 ! [kg mol-1] Mean molecular weight of water HITRAN
      real(r8),parameter::mmw_dry_airx=28.9644e-3 ! [kg mol-1] Mean molecular weight of dry air CCM2/3:physics/radcsw()
      real(r8),parameter::cst_von_krmx=0.4 ! [fraction] Von Karman constant

c     Parameters for boundary layer meteorology
      real(r8),parameter::cst_Stefan_Boltzmannx=5.67032e-8 ! [W m-2 K-4] Stefan-Boltzmann constant
      real(r8),parameter::dns_H2O_lqd_stdx=1000.0 ! [kg m-3] Density of liquid water
      real(r8),parameter::ltn_heat_fsn_H2O_stdx=0.3336e06 ! [J kg-1] Latent heat of fusion of H2O at 0 C, standard CCM:lsm/phyconi.F 
      real(r8),parameter::ltn_heat_sbl_H2O_stdx=2.8440e06 ! [J kg-1] Latent heat of sublimation of H2O, standard CCM:lsm/phyconi.F
      real(r8),parameter::ltn_heat_vpr_H2O_stdx=2.5104e06 ! [J kg-1] Latent heat of vaporization of H2O, standard CCM:lsm/phyconi.F
      real(r8),parameter::spc_heat_H2O_ice_vlmx=2.094e06 ! [J m-3 K-1] Volumetric specific heat capacity of ice water CCM:lsm/phyconi.F
      real(r8),parameter::spc_heat_H2O_lqd_vlmx=4.188e06 ! [J m-3 K-1] Volumetric specific heat capacity of liquid water CCM:lsm/phyconi.F
      real(r8),parameter::tpt_frz_pntx=273.15 ! [K] Kelvin--Celsius scale offset Bol80
c     Input
      real(r8) grv_sfcx             ! [m s-2] (9.80616) Mean gravitational acceleration at Earth's surface
      real(r8) gas_cst_dry_airx     ! [J kg-1 K-1] (287.05) Gas constant of dry air
      real(r8) rds_Earthx           ! [m] (6.37122e+6) Radius of sphere of same volume as Earth
c     Local
      real(r8) gas_cst_H2O          ! (461.65) [J kg-1 K-1] Gas constant of H2O
      real(r8) spc_heat_H2O_vpr     ! (1850.0) [J kg-1 K-1] IrG81 pp. 77, 245
c     Main code

c     Parameters for source and sink processes
      grv_sfc=grv_sfcx          ! [m s-2] (9.80616) Mean gravitational acceleration at Earth's surface
      grv_sfc_rcp=1.0/grv_sfcx  ! [s2 m-1] (0.101977) Reciprocal of gravity
      gas_cst_dry_air=gas_cst_dry_airx ! [J kg-1 K-1] (287.05) Gas constant of dry air
      rds_Earth=rds_Earthx      ! [m] (6.37122e+6) Radius of sphere of same volume as Earth
      mmw_H2O=mmw_H2Ox          ! [kg mol-1] Mean molecular weight of dry air
      mmw_dry_air=mmw_dry_airx  ! [kg mol-1] Mean molecular weight of dry air
      gas_cst_unv=gas_cst_unvx  ! [J mol-1 K-1] Universal gas constant
      cst_von_krm=cst_von_krmx  ! [frc] Von Karman's constant
      cst_von_krm_rcp=1.0/cst_von_krmx ! [frc] Reciprocal of Von Karman's constant
      eps_H2O=mmw_H2Ox/mmw_dry_airx ! [frc] (0.622) mmw(H2O)/mmw(dry air)
      eps_H2O_rcp_m1=-1.0+mmw_dry_airx/mmw_H2Ox ! [frc] (0.60777) Constant for virtual temperature

c     Parameters for boundary layer meteorology
      cst_Stefan_Boltzmann=cst_Stefan_Boltzmannx ! (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant
      dns_H2O_lqd_std=dns_H2O_lqd_stdx ! [kg m-3] Density of liquid water
      ltn_heat_fsn_H2O_std=ltn_heat_fsn_H2O_stdx ! [J kg-1] Latent heat of fusion of H2O at 0 C, standard
      ltn_heat_sbl_H2O_std=ltn_heat_sbl_H2O_stdx ! [J kg-1] Latent heat of sublimation of H2O, standard
      ltn_heat_vpr_H2O_std=ltn_heat_vpr_H2O_stdx ! [J kg-1] Latent heat of vaporization of H2O, standard
      one_mns_eps_H2O=1.0-eps_H2O   ! (0.378) [frc] Constant for saturation specific humidity
      spc_heat_H2O_ice_vlm=spc_heat_H2O_ice_vlmx ! [J m-3 K-1] Volumetric specific heat capacity of ice water
      spc_heat_H2O_lqd_vlm=spc_heat_H2O_lqd_vlmx ! [J m-3 K-1] Volumetric specific heat capacity of liquid water
      tpt_frz_pnt=tpt_frz_pntx  ! [K] Kelvin--Celsius scale offset

c     Derived 
      spc_heat_dry_air=7.0*gas_cst_dry_airx/2.0 ! (1005.0) [J kg-1 K-1]
      kappa_dry_air=gas_cst_dry_airx/spc_heat_dry_air ! (0.286 = 2/7) [frc] Constant in potential temperature (IrG81 p. 25, Tre922 p. 72) 
      gas_cst_H2O=gas_cst_unv/mmw_H2O ! (461.65) [J kg-1 K-1] Gas constant of H2O
      spc_heat_H2O_vpr=4.0*gas_cst_H2O ! (1850.0) [J kg-1 K-1] Specific heat of H2O vapor IrG81 pp. 77, 245
      cp_vpr_rcp_cp_dry_m1=spc_heat_H2O_vpr/spc_heat_dry_air-1.0 ! (0.83745) [frc] Constant for moist specific heat IrG81 p. 77

c     Sanity checks on some derived constants
      if (abs(461.52-gas_cst_H2O)/461.52 > 1.0e-4) stop 'dst_cst_cmn_ini(): gas_cst_H2O error'
      if (abs(1846.08-spc_heat_H2O_vpr)/1846.08 > 1.0e-4) stop 'dst_cst_cmn_ini(): spc_heat_H2O_vpr error'
      if (abs(0.83745-cp_vpr_rcp_cp_dry_m1)/0.83745 > 1.0e-3) then
         write (6,'(a,f11.6)') 'cp_vpr_rcp_cp_dry_m1 = ',cp_vpr_rcp_cp_dry_m1
         stop 'dst_cst_cmn_ini(): cp_vpr_rcp_cp_dry_m1 error'
      endif
      if (abs(0.28571-kappa_dry_air)/0.28571 > 1.0e-4) stop 'dst_cst_cmn_ini(): kappa_dry_air error'
      if (abs(273.15-tpt_frz_pnt)/273.15 > 1.0e-4) stop 'dst_cst_cmn_ini(): tpt_frz_pnt error'
      if (abs(1004.697-spc_heat_dry_air)/1004.697 > 1.0e-4) stop 'dst_cst_cmn_ini(): spc_heat_dry_air error'
      if (abs(0.60777-eps_H2O_rcp_m1)/0.60777 > 1.0e-4) stop 'dst_cst_cmn_ini(): eps_H2O_rcp_m1 error'

      return
      end                       ! end dst_cst_cmn_ini()
      
      subroutine dst_msc_cmn_ini()
c     Purpose: Initialize parameter-dependent fields in miscellaneous common blocks
c     dst_msc_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
      use precision
      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstmsc.h> /* Miscellaneous variables for dust parameterization */
c     Input
c     Output
c     Input/Output
c     Local
      integer idx               ! [idx] Counting index 
c     Main code

c     Initialize optical depth information
      call dst_odx_cmn_ini()

c     Initialize wet scavenging information
c      call dst_scv_cmn_ini()

#ifdef DST_CHM
c     Initialize dust chemistry
      call dst_chm_cmn_ini()
#endif /* not DST_CHM */

c     All source/sink routines compare new dust mixing ratio to q_dst_mxm
c     q_dst > q_dst_mxm indicates unbelievably high dust emissions
c     Model should abort until veracity can be established
c     q_dst_mxm should exceed greatest recorded atmospheric dust concentration
      q_dst_mxm=1.0e-4          ! [kg kg-1] Maximum plausible mass mixing ratio of dust
      if(q_dst_mxm < 1.0e-7) stop 'dst_msc_cmn_ini() reports q_dst_mxm too small, will trip for reasonable dust concentrations'
      
c     Source/sink routines compare new dust mixing ratio to q_dst_sgn
c     q_dst < q_dst_sgn indicates unbelievably low dust emissions
c     probably due to numerical artifacts like rounding, underflow, order of operations
c     q_dst_sgn = 0.0 is fine, but numeric artifacts will not be eliminated
c     This can result in underflows from operations involving very small values
c     q_dst_sgn should not exceed about 1.0e-20 kg kg-1
c     WARNING: Model is free to (and does) eliminate concentrations lower than q_dst_sgn
      q_dst_sgn=1.0e-30         ! [kg kg-1] Minimum plausible mass mixing ratio of dust
      if(q_dst_sgn > 1.0e-20) stop 'dst_msc_cmn_ini() reports q_dst_sgn too large, fixers may cause significant mass imbalances'
      
c     Threshold flux indicating presence of "significant" dust
c     Evaporative source routine prints warning when evaporative dust flux 
c     exceeds both the available precipitating flux and flx_mss_sgn
      flx_mss_sgn=1.0e-30       ! [kg m-2 s-1] Minimum significant flux of dust
      
c     Source routine compares new dust flux to flx_mss_mxm
c     flx_mss_vrt_dst_ttl > flx_mss_mxm indicates suspiciously high dust fluxes
      flx_mss_mxm=100.0e-9      ! [kg m-2 s-1] Maximum plausible flux of dust
      if(flx_mss_mxm < 1.0e-9) stop 'dst_msc_cmn_ini() reports flx_mss_mxm too small, will trip for reasonable mass fluxes'
      
c     Scavenging ratios
c     fxm: z_scv is only used in old wet deposition scheme, remove once PJR/MKB scheme fully adopted
c-mahow      z_scv(1)=1000.0           ! [(kg dust/kg rain)/(kg dust/kg air)] Uematsu
c     z_scv(idx)=500.0         ! [(kg dust/kg rain)/(kg dust/kg air)] dstmch19
c-mahow      do idx=2,dst_nbr          ! NB: loop starts at 2
c-mahow   z_scv(idx)=300.0       ! [(kg dust/kg rain)/(kg dust/kg air)] Duce
c     z_scv(idx)=750.0         ! [(kg dust/kg rain)/(kg dust/kg air)] TeF94 p. 22901

c++mahow
      do idx=1,dst_nbr          ! NB: loop starts at 2
         z_scv(idx)=750.0               ! [(kg dust/kg rain)/(kg dust/kg air)] Duce
      end do                    ! end loop over cst
c--mahow

      
c     Initialize common block used for debugging
#ifdef BXM
      lat_dbg=1                 ! [idx] Latitude  index for diagnostics
      lev_dbg=1                 ! [idx] Level     index for diagnostics
      lon_dbg=1                 ! [idx] Longitude index for diagnostics
      time_dbg=0                ! [idx] Time      index for diagnostics
#endif /* not BXM */
#ifdef T5
c     T5 lon = i = 13 = 112.5W and lat = j = 3 = 31.7S (land): Chile
      lat_dbg=3                 ! [idx] Latitude  index for diagnostics
      lev_dbg=11                ! [idx] Level     index for diagnostics
      lon_dbg=13                ! [idx] Longitude index for diagnostics
      time_dbg=0                ! [idx] Time      index for diagnostics
#endif /* not T5 */
#ifdef T42
c     T42 lon = i = 99 = 84.3W and lat = j = 37 = 12.5N (land): Panama isthmus
c     T42 lon = i = 98 = 87.1W and lat = j = 37 = 12.5N (ocean): Pacific
c     T42 lon = i = 115 = 39.8W and lat = j = 61 = 79.5N (ocean): Greenland
c     T42 lon = i = 125 = 11.25W and lat = j = 40 = 20.9N (land): Desert at center of Mauritania
c     fxm: pass in coordinate arrays to simplify index identification
c     lat_dbg=vec_val2idx(lat,lat_nbr,20.9) ! [idx] Latitude  index for diagnostics
c     lon_dbg=vec_val2idx(lon,lon_nbr,348.75) ! [idx] Longitude index for diagnostics
      lat_dbg=40                ! [idx] Latitude  index for diagnostics
      lev_dbg=9                 ! [idx] Level     index for diagnostics
      lon_dbg=125               ! [idx] Longitude index for diagnostics
      time_dbg=0                ! [idx] Time      index for diagnostics
#endif /* not T42 */ 
#ifdef T62
c     T62 lon = i = 187 = 11.25W and lat = j = 58 = 20.0N (land): Desert at center of Mauritania
c     T62 lon = i = 57 = 105E and lat = j = 66 = 35.23N (land): Desert near Takla Makan
      lat_dbg=66                ! [idx] Latitude  index for diagnostics
      lev_dbg=9                 ! [idx] Level     index for diagnostics
      lon_dbg=57                ! [idx] Longitude index for diagnostics
      time_dbg=0                ! [idx] Time      index for diagnostics
#endif /* not T62 */ 
      
      return
      end                       ! end dst_msc_cmn_ini()

      subroutine dst_nm_cmn_ini()
c     Purpose: Initialize names of dust fields NOT stored in comtracnm
c     dst_nm_cmn_ini() is called by CCM:physics/inti(), MATCH:inirun()
c     NB: The call to dst_nm_cmn_ini() could be lumped into dst_msc_cmn_ini()
c     except that MATCH requires dst_nm_cmn_ini() be called "earlier" than CCM
      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstnm.h> /* Nomenclature for outfld() */
c     Local
      character*3 trc_nbr         ! Advected species number (Character)
      integer idx               ! cst index 

c     Constituent-dependent arrays in dst_nm_cmn
      do idx=1,dst_nbr
         write (unit=trc_nbr,fmt='(i3)') idx+100
         flx_mss_mbl_sfc_nm(idx)='DSTSFM'//trc_nbr(2:3)
         flx_mss_pcp_sfc_nm(idx)='DSTSFP'//trc_nbr(2:3)
         flx_mss_trb_sfc_nm(idx)='DSTSFT'//trc_nbr(2:3)
         flx_mss_grv_sfc_nm(idx)='DSTSFG'//trc_nbr(2:3)
         flx_mss_dry_sfc_nm(idx)='DSTSFD'//trc_nbr(2:3)
         mpc_dst_nm(idx)='DSTMPC'//trc_nbr(2:3)
         odxc_dst_nm(idx)='DSTODX'//trc_nbr(2:3)
      end do                    ! end loop over cst

c     Nomenclature for standard radiative fields (currently also in dst_nm_cmn)
c     NB: These names are used by CCM, but not MATCH or MOZART
      flat_nm='FLAT    '
      flatc_nm='FLATC   '
      flds_nm='FLDS    '
      fsat_nm='FSAT    '
      fsatc_nm='FSATC   '
      ftat_nm='FTAT    '
      ftns_nm='FTNS    '
      ftnt_nm='FTNT    '
      qrt_nm='QRT     '

c     Nomenclature for radiative forcing fields in dst_frc_nm_cmn
c     NB: These names are used by CCM, but not MATCH or MOZART
      flat_dst_nm='FLATDST '
      flat_frc_nm='FLATFRC '
      flatc_dst_nm='FLATCDST'
      flatc_frc_nm='FLATCFRC'
      flds_dst_nm='FLDSDST '
      flds_frc_nm='FLDSFRC '
      flns_dst_nm='FLNSDST '
      flns_frc_nm='FLNSFRC '
      flnsc_dst_nm='FLNSCDST'
      flnsc_frc_nm='FLNSCFRC'
      flnt_dst_nm='FLNTDST '
      flnt_frc_nm='FLNTFRC '
      flntc_dst_nm='FLNTCDST'
      flntc_frc_nm='FLNTCFRC'
      fsat_dst_nm='FSATDST '
      fsat_frc_nm='FSATFRC '
      fsatc_dst_nm='FSATCDST'
      fsatc_frc_nm='FSATCFRC'
      fsds_dst_nm='FSDSDST '
      fsds_frc_nm='FSDSFRC '
      fsns_dst_nm='FSNSDST '
      fsns_frc_nm='FSNSFRC '
      fsnsc_dst_nm='FSNSCDST'
      fsnsc_frc_nm='FSNSCFRC'
      fsnt_dst_nm='FSNTDST '
      fsnt_frc_nm='FSNTFRC '
      fsntc_dst_nm='FSNTCDST'
      fsntc_frc_nm='FSNTCFRC'
      ftat_dst_nm='FTATDST '
      ftat_frc_nm='FTATFRC '
      ftns_dst_nm='FTNSDST '
      ftns_frc_nm='FTNSFRC '
      ftnt_dst_nm='FTNTDST '
      ftnt_frc_nm='FTNTFRC '
      qrl_dst_nm='QRLDST  '
      qrl_frc_nm='QRLFRC  '
      qrs_dst_nm='QRSDST  '
      qrs_frc_nm='QRSFRC  '
      qrt_dst_nm='QRTDST  '
      qrt_frc_nm='QRTFRC  '

      return
      end                       ! end dst_nm_cmn_ini()

      subroutine dst_odx_cmn_ini()
c     Purpose: Initialize optical depth information
c     dst_odx_cmn_ini() is called by dst_msc_cmn_ini()
c     Dependencies: Routine requires dmt_swr, must be called after dst_szdst_ini()
c     NB: Nomenclature encodes wavelength in microns as last three digits
c     Thus ..._063 means "evaluated at 0.63 microns" (AVHRR Channel 1)
c     Thus ..._050 means "evaluated at 0.50 microns" (WMO std wavelength)
      use precision
      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstszdst.h> /* Diagnostic aerosol microphysical properties */
#include <dstodx.h> /* Optical depth information */
#include <xtr.com> /* Extrapolation constants */
c     Parameters
      integer,parameter::xtr_typ_LHS=xtr_vrb ! [enm] LHS Extrapolation type
      integer,parameter::xtr_typ_RHS=xtr_vrb ! [enm] RHS Extrapolation type
c     High resolution specific extinction
c     Data generated by mie on Sun May 21 12:39:02 2000
c     mie version Unknown dated Unknown
c     Saharan dust aerosol mass extinction coefficient in m2 kg-1
c     Command line: mie --no_wrn_ntp --hrz --dbg=1 --dst_lbl=01 --aer_typ=saharan_dust --dist=lognormal --wvl_grd=reg --wvl_nbr=1 --wvl_mnm=0.625 --wvl_mxm=0.635 --bnd_nbr=2 --sz_grd=log --sz_mnm=0.005 --sz_mxm=50.0 --sz_nbr=64 --rds_nma=0.2986 --gsd_anl=2.0 /data/zender/dst/mie/aer_saharan_dust_hrz_063.nc
      integer,parameter::dst_nbr_hrz=64 ! [nbr] Number of high resolution bins
      real(r8),dimension(dst_nbr_hrz),parameter::ext_cff_mss_dst_hrz=(/ 
     $     2.0708963e+01, 2.0787266e+01, 2.0905321e+01, 2.1083765e+01, 
     $     2.1354128e+01, 2.1764294e+01, 2.2388514e+01, 2.3339769e+01, 
     $     2.4791828e+01, 2.7011854e+01, 3.0411173e+01, 3.5623936e+01, 
     $     4.3629124e+01, 5.5939003e+01, 7.4890823e+01, 1.0408587e+02, 
     $     1.4903189e+02, 2.1799266e+02, 3.2283185e+02, 4.7891080e+02, 
     $     7.0142871e+02, 9.9449951e+02, 1.3420798e+03, 1.7650575e+03, 
     $     2.4588586e+03, 3.1149944e+03, 3.1692219e+03, 3.7535420e+03, 
     $     3.4849006e+03, 3.5819019e+03, 2.9263354e+03, 2.1190962e+03, 
     $     1.3966730e+03, 8.7480927e+02, 7.7586859e+02, 8.3421100e+02, 
     $     8.8776550e+02, 6.5199640e+02, 4.4449576e+02, 5.0862491e+02, 
     $     3.6941647e+02, 3.4139664e+02, 2.5660623e+02, 2.5565089e+02, 
     $     1.9882104e+02, 1.7131413e+02, 1.5897310e+02, 1.3040485e+02, 
     $     1.0873315e+02, 9.4183754e+01, 8.2329193e+01, 7.0274239e+01, 
     $     6.2494713e+01, 5.3757774e+01, 4.6389805e+01, 3.9296574e+01, 
     $     3.4052509e+01, 2.9806177e+01, 2.5394520e+01, 2.2001064e+01, 
     $     1.9136105e+01, 1.6531433e+01, 1.4268397e+01, 1.2337605e+01 /)

      real(r8),dimension(dst_nbr_hrz),parameter::dmt_ctr_hrz=(/ 
     $     1.0773910e-08, 1.2441516e-08, 1.4367239e-08, 1.6591027e-08, 
     $     1.9159019e-08, 2.2124489e-08, 2.5548960e-08, 2.9503479e-08, 
     $     3.4070084e-08, 3.9343519e-08, 4.5433190e-08, 5.2465424e-08, 
     $     6.0586132e-08, 6.9963768e-08, 8.0792901e-08, 9.3298183e-08, 
     $     1.0773906e-07, 1.2441511e-07, 1.4367232e-07, 1.6591021e-07, 
     $     1.9159012e-07, 2.2124479e-07, 2.5548951e-07, 2.9503465e-07, 
     $     3.4070069e-07, 3.9343499e-07, 4.5433163e-07, 5.2465396e-07, 
     $     6.0586092e-07, 6.9963721e-07, 8.0792842e-07, 9.3298115e-07, 
     $     1.0773897e-06, 1.2441503e-06, 1.4367223e-06, 1.6591009e-06, 
     $     1.9158997e-06, 2.2124464e-06, 2.5548932e-06, 2.9503444e-06, 
     $     3.4070044e-06, 3.9343472e-06, 4.5433130e-06, 5.2465357e-06, 
     $     6.0586044e-06, 6.9963671e-06, 8.0792788e-06, 9.3298058e-06, 
     $     1.0773891e-05, 1.2441495e-05, 1.4367213e-05, 1.6590999e-05, 
     $     1.9158986e-05, 2.2124450e-05, 2.5548916e-05, 2.9503426e-05, 
     $     3.4070024e-05, 3.9343449e-05, 4.5433106e-05, 5.2465330e-05, 
     $     6.0586011e-05, 6.9963629e-05, 8.0792743e-05, 9.3298113e-05 /)
c     fxm: dmt_grd only needed for rbn_vec(), not ntp_vec()
      real(r8),dimension(dst_nbr_hrz+1),parameter::dmt_grd_hrz=(/ 
     $     9.9999999e-09, 1.1547820e-08, 1.3335214e-08, 1.5399264e-08, 
     $     1.7782792e-08, 2.0535246e-08, 2.3713731e-08, 2.7384189e-08, 
     $     3.1622768e-08, 3.6517402e-08, 4.2169638e-08, 4.8696737e-08, 
     $     5.6234114e-08, 6.4938142e-08, 7.4989394e-08, 8.6596401e-08, 
     $     9.9999959e-08, 1.1547814e-07, 1.3335207e-07, 1.5399257e-07, 
     $     1.7782784e-07, 2.0535238e-07, 2.3713721e-07, 2.7384178e-07, 
     $     3.1622753e-07, 3.6517383e-07, 4.2169614e-07, 4.8696711e-07, 
     $     5.6234080e-07, 6.4938098e-07, 7.4989345e-07, 8.6596339e-07, 
     $     9.9999886e-07, 1.1547806e-06, 1.3335199e-06, 1.5399247e-06, 
     $     1.7782772e-06, 2.0535224e-06, 2.3713706e-06, 2.7384158e-06, 
     $     3.1622731e-06, 3.6517358e-06, 4.2169586e-06, 4.8696675e-06, 
     $     5.6234039e-06, 6.4938054e-06, 7.4989293e-06, 8.6596283e-06, 
     $     9.9999825e-06, 1.1547799e-05, 1.3335190e-05, 1.5399237e-05, 
     $     1.7782761e-05, 2.0535210e-05, 2.3713690e-05, 2.7384142e-05, 
     $     3.1622712e-05, 3.6517336e-05, 4.2169562e-05, 4.8696649e-05, 
     $     5.6234010e-05, 6.4938016e-05, 7.4989250e-05, 8.6596228e-05,
     $     9.9999997e-05 /)

c     Precomputed values for common configurations
c     Data generated by $HOME/dst/szdst.pl calling mie:
c     $HOME/dst/szdst.pl --dbg=1 --ODX --aer_nbr=4 --dst_lbl={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --wvl_odx_mnm=0.625 --wvl_odx_mxm=0.635 > & ! $DATA/dst/mie/szdst_ODX.txt
c     mie --no_wrn_ntp --dbg=1 --dst_lbl=01 --aer_typ=saharan_dust --dist=lognormal --wvl_grd=reg --wvl_nbr=1 --wvl_mnm=0.625 --wvl_mxm=0.635 --bnd_nbr=2 --sz_grd=log --sz_mnm=0.05 --sz_mxm=0.5 --sz_nbr=200 --rds_nma=0.2986 --gsd_anl=2.0 /data/zender/dst/mie/aer_saharan_dust_01_063.nc
      real(r8),dimension(4),parameter::ext_cff_mss_dst_063_04bn=(/ 2892.77,835.046,382.509,196.051 /) ! [m2 kg-1] Mass extinction coefficient of dust for diagnostic optical depth
      real(r8),dimension(4),parameter::ext_cff_mss_dst_0532_04bn=(/ 2834.69, 714.585, 330.916, 156.794 /)
      real(r8),dimension(4),parameter::ext_cff_mss_dst_0550_04bn=(/ 2834.69, 714.585, 330.916, 156.794 /)
      real(r8),dimension(4),parameter::ext_cff_mss_dst_0558_04bn=(/ 2834.69, 714.585, 330.916, 156.794 /)
      real(r8),dimension(4),parameter::ext_cff_mss_dst_0630_04bn=(/ 2872.42, 753.6, 330.37, 158.369 /)
c     Input
c     Output
c     Input/Output
c     Local
      integer sz_idx            ! [idx] Counting index 
c     Main code
c     Sanity check
      if(dmt_swr(1) < 0.0001e-6.or.dmt_swr(1) > 100.0e-6) 
     $     stop 'ERROR: dst_odx_cmn_ini() probably called before dst_szdst_ini()'
c     Initialize common block used for diagnostic optical depth output
c     NB: dmt_swr is only variable that requires dstszdst.h in this routine
      call ntp_vec(dst_nbr_hrz,dmt_ctr_hrz,ext_cff_mss_dst_hrz,dst_nbr,dmt_swr,ext_cff_mss_dst_dgn,xtr_typ_LHS,xtr_typ_RHS)
c      call rbn_vec(dst_nbr_hrz,dmt_grd_hrz,ext_cff_mss_dst_hrz,
c     $     dst_nbr,dmt_grd,ext_cff_mss_dst_dgn,
c     $     xtr_typ_LHS,xtr_typ_RHS)

c     Override automatic grid with preset grid if available
      if (dst_nbr == 4) then
         do sz_idx=1,dst_nbr             
#ifdef LITE
            ext_cff_mss_dst_dgn(sz_idx)=ext_cff_mss_dst_0532_04bn(sz_idx)
#endif
#ifdef MODIS
            ext_cff_mss_dst_dgn(sz_idx)=ext_cff_mss_dst_0550_04bn(sz_idx)
#endif
#ifdef MISR
            ext_cff_mss_dst_dgn(sz_idx)=ext_cff_mss_dst_0558_04bn(sz_idx)
#endif
#ifdef AVHRR
            ext_cff_mss_dst_dgn(sz_idx)=ext_cff_mss_dst_0630_04bn(sz_idx)
#endif
#ifdef OLDOPTICS
            ext_cff_mss_dst_dgn(sz_idx)=ext_cff_mss_dst_063_04bn(sz_idx) ! [m2 kg-1] Mass extinction coefficient of dust for diagnostic optical depth
#endif
         end do                 ! end loop over cst
      endif                     ! endif dst_nbr == 4

      write (6,'(a)') 'Optical depth information:'
      write (6,'(a3,2(a10,1x))') 'idx','dmt_swr','ext_spc'
      write (6,'(a3,2(a10,1x))') '','um','m2 kg-1'
      do sz_idx=1,dst_nbr
         write (6,'(1p,i3,2(e10.3,1x)),0p') sz_idx,
     $        dmt_swr(sz_idx)*1.0e6,ext_cff_mss_dst_dgn(sz_idx)
      enddo                     ! end loop over cst

      return
      end                       ! end dst_odx_cmn_ini()

      subroutine dst_scv_cmn_ini()
c     Purpose: Initialize wet scavenging parameters
c     dst_scv_cmn_ini() is called by dst_msc_cmn_ini()
      use precision
      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstaer.h> /* Aerosol microphysical properties */
#include <xtr.com> /* Extrapolation constants */
c     Parameters
      integer,parameter::xtr_typ_LHS=xtr_vrb ! [enm] LHS Extrapolation type
      integer,parameter::xtr_typ_RHS=xtr_vrb ! [enm] RHS Extrapolation type
c     High resolution scavenging coefficients
c     Data generated by mie on Tue Jul 11 10:44:24 2000
c     mie version Unknown dated Unknown
c     Saharan dust aerosol mass extinction coefficient in m2 kg-1
c     Saharan dust aerosol scavenging coefficient in s-1
c     Command line: mie --no_wrn_ntp -no_mie --hrz --dbg=1 --dst_lbl=01 --dmt_pcp=1000.0e-6 --flx_hgt_pcp=1.0e-6 --aer_typ=saharan_dust --dist=lognormal --sz_grd=log --sz_mnm=0.005 --sz_mxm=50.0 --sz_nbr=64 --rds_nma=0.2986 --gsd_anl=2.0 /data/zender/mie/out.nc
      integer,parameter::dst_nbr_hrz=64 ! [nbr] Number of high resolution bins
      real(r8),dimension(dst_nbr_hrz),parameter::scv_cff_dst_hrz=(/ 
     $     6.5642057e-06, 5.5623586e-06, 4.7186900e-06, 4.0076829e-06, 
     $     3.4080399e-06, 2.9019666e-06, 2.4745857e-06, 2.1134440e-06, 
     $     1.8081125e-06, 1.5498485e-06, 1.3313188e-06, 1.1463688e-06, 
     $     9.8982969e-07, 8.5735991e-07, 7.4531175e-07, 6.5062244e-07, 
     $     5.7072310e-07, 5.0346443e-07, 4.4705709e-07, 4.0002405e-07, 
     $     3.6116458e-07, 3.2952812e-07, 3.0439827e-07, 2.8528660e-07, 
     $     2.7193553e-07, 2.6433227e-07, 2.6273474e-07, 2.6771221e-07, 
     $     2.8020460e-07, 3.0160646e-07, 3.3388361e-07, 3.7973254e-07, 
     $     4.4279730e-07, 9.1151329e-04, 9.1287307e-04, 9.1413822e-04, 
     $     9.1533037e-04, 9.1647520e-04, 9.1760309e-04, 9.1875100e-04, 
     $     9.1996416e-04, 9.2129898e-04, 9.2282752e-04, 9.2464121e-04, 
     $     9.2685974e-04, 9.2963898e-04, 9.3318376e-04, 9.3776471e-04, 
     $     9.4373996e-04, 9.5158460e-04, 9.6193026e-04, 9.7561633e-04, 
     $     9.9376054e-04, 1.0178511e-03, 1.0498704e-03, 1.0924604e-03, 
     $     1.1491410e-03, 1.2246047e-03, 1.3251060e-03, 1.4589831e-03, 
     $     1.5000001e-03, 1.5000001e-03, 1.5000001e-03, 1.3411095e-03 /)

      real(r8),dimension(dst_nbr_hrz),parameter::dmt_ctr_hrz=(/ 
     $     1.0773910e-08, 1.2441516e-08, 1.4367239e-08, 1.6591027e-08, 
     $     1.9159019e-08, 2.2124489e-08, 2.5548960e-08, 2.9503479e-08, 
     $     3.4070084e-08, 3.9343519e-08, 4.5433190e-08, 5.2465424e-08, 
     $     6.0586132e-08, 6.9963768e-08, 8.0792901e-08, 9.3298183e-08, 
     $     1.0773906e-07, 1.2441511e-07, 1.4367232e-07, 1.6591021e-07, 
     $     1.9159012e-07, 2.2124479e-07, 2.5548951e-07, 2.9503465e-07, 
     $     3.4070069e-07, 3.9343499e-07, 4.5433163e-07, 5.2465396e-07, 
     $     6.0586092e-07, 6.9963721e-07, 8.0792842e-07, 9.3298115e-07, 
     $     1.0773897e-06, 1.2441503e-06, 1.4367223e-06, 1.6591009e-06, 
     $     1.9158997e-06, 2.2124464e-06, 2.5548932e-06, 2.9503444e-06, 
     $     3.4070044e-06, 3.9343472e-06, 4.5433130e-06, 5.2465357e-06, 
     $     6.0586044e-06, 6.9963671e-06, 8.0792788e-06, 9.3298058e-06, 
     $     1.0773891e-05, 1.2441495e-05, 1.4367213e-05, 1.6590999e-05, 
     $     1.9158986e-05, 2.2124450e-05, 2.5548916e-05, 2.9503426e-05, 
     $     3.4070024e-05, 3.9343449e-05, 4.5433106e-05, 5.2465330e-05, 
     $     6.0586011e-05, 6.9963629e-05, 8.0792743e-05, 9.3298113e-05 /)

c     Precomputed values for common configurations
      real(r8),dimension(4),parameter::scv_cff_dst_04bn=(/ 3.27722e-07,0.000771208,0.000920177,0.000928821 /) ! [s-1] Mass mean scavenging coefficient
c     Input
c     Output
c     Input/Output
c     Local
      integer sz_idx            ! [idx] Counting index 
      real(r8) scv_cff_dst(dst_nbr) ! [s-1] Scavenging coefficient
c     Main code
c     Sanity check
      if(dmt_vwr(1) < 0.0001e-6.or.dmt_vwr(1) > 100.0e-6) 
     $     stop 'ERROR: dst_scv_cmn_ini() probably called before dst_szdst_ini()'
c     Initialize common block used for diagnostic optical depth output
      call ntp_vec(dst_nbr_hrz,dmt_ctr_hrz,scv_cff_dst_hrz,dst_nbr,dmt_vwr,scv_cff_dst,xtr_typ_LHS,xtr_typ_RHS)
c      call rbn_vec(dst_nbr_hrz,dmt_grd_hrz,scv_cff_dst_hrz,
c     $     dst_nbr,dmt_grd,scv_cff_dst,
c     $     xtr_typ_LHS,xtr_typ_RHS)

c     Override automatic grid with preset grid if available
      if (dst_nbr == 4) then
         do sz_idx=1,dst_nbr             
            scv_cff_dst(sz_idx)=scv_cff_dst_04bn(sz_idx) ! [s-1] Scavenging coefficient
         end do                 ! end loop over cst
      endif                     ! endif dst_nbr == 4

      write (6,'(a)') 'Scavenging coefficient information:'
      write (6,'(a3,2(a10,1x))') 'idx','dmt_vwr','scv_cff'
      write (6,'(a3,2(a10,1x))') '','um','s-1'
      do sz_idx=1,dst_nbr
         write (6,'(1p,i3,2(e10.3,1x)),0p') sz_idx,
     $        dmt_vwr(sz_idx)*1.0e6,scv_cff_dst(sz_idx)
      enddo                     ! end loop over cst

      return
      end                       ! end dst_scv_cmn_ini()
