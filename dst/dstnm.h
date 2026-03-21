c     $Header$ -*-fortran-*-

c     Purpose: Common block dstnm stores variables names used in outfld()
c     dst_nm_cmn is initialized in dst_nm_cmn_ini() which is called by dst_msc_cmn_ini()

c     Usage: 
c     #include <dstnm.h> /* Nomenclature for outfld() */
c     params.h needed for DST_NBR
#ifndef PARAMS_H
#include <params.h>
#endif /* not PARAMS_H */ 

c     Nomenclature info for outfld()
      common / dst_nm_cmn /     ! dst_nm_cmn is initialized in dst_nm_cmn_ini()
     $     flx_mss_dry_sfc_nm(dst_nbr), ! [sng] Name of surface dry deposition flux
     $     flx_mss_grv_sfc_nm(dst_nbr), ! [sng] Name of surface gravitational deposition flux
     $     flx_mss_mbl_sfc_nm(dst_nbr), ! [sng] Name of surface mobilization flux
     $     flx_mss_pcp_sfc_nm(dst_nbr), ! [sng] Name of surface precipitation flux
     $     flx_mss_trb_sfc_nm(dst_nbr), ! [sng] Name of surface turbulent deposition flux
     $     mpc_dst_nm(dst_nbr), ! [sng] Name of column mass path of dust
     $     odxc_dst_nm(dst_nbr), ! [sng] Name of column extinction optical depth of dust
     $     flat_nm,             ! Flux Longwave Absorbed Atmosphere
     $     flatc_nm,            ! Flux Longwave Absorbed Atmosphere, clear sky
     $     flds_nm,             ! Flux Longwave Downwelling Surface
     $     fsat_nm,             ! Flux Shortwave Absorbed Atmosphere
     $     fsatc_nm,            ! Flux Shortwave Absorbed Atmosphere, clear sky
     $     ftat_nm,             ! Flux Total Absorbed Atmosphere
     $     ftns_nm,             ! Flux Total Net Surface
     $     ftnt_nm,             ! Flux Total Net TOA
     $     qrt_nm               ! Total radiative heating rate
      character*8 flx_mss_dry_sfc_nm ! [sng] Name of surface dry deposition flux
      character*8 flx_mss_grv_sfc_nm ! [sng] Name of surface gravitational deposition flux
      character*8 flx_mss_mbl_sfc_nm ! [sng] Name of surface mobilization flux
      character*8 flx_mss_pcp_sfc_nm ! [sng] Name of surface precipitation flux
      character*8 flx_mss_trb_sfc_nm ! [sng] Name of surface turbulent deposition flux
      character*8 mpc_dst_nm    ! [sng] Name of column mass path of dust
      character*8 odxc_dst_nm   ! [sng] Name of column extinction optical depth of dust
      character*8 flat_nm       ! Flux Longwave Absorbed Atmosphere
      character*8 flatc_nm      ! Flux Longwave Absorbed Atmosphere, clear sky
      character*8 flds_nm       ! Flux Longwave Downwelling Surface
      character*8 fsat_nm       ! Flux Shortwave Absorbed Atmosphere
      character*8 fsatc_nm      ! Flux Shortwave Absorbed Atmosphere, clear sky
      character*8 ftat_nm       ! Flux Total Absorbed Atmosphere
      character*8 ftns_nm       ! Flux Total Net Surface
      character*8 ftnt_nm       ! Flux Total Net TOA
      character*8 qrt_nm        ! Total radiative heating rate

c     Nomenclature info for outfld()
      common / dst_frc_nm_cmn / ! dst_frc_nm_cmn is initialized in dst_nm_cmn_ini()
     $     flat_dst_nm,
     $     flat_frc_nm,
     $     flatc_dst_nm,
     $     flatc_frc_nm,
     $     flds_dst_nm,
     $     flds_frc_nm,
     $     flns_dst_nm,
     $     flns_frc_nm,
     $     flnsc_dst_nm,
     $     flnsc_frc_nm,
     $     flnt_dst_nm,
     $     flnt_frc_nm,
     $     flntc_dst_nm,
     $     flntc_frc_nm,
     $     fsat_dst_nm,
     $     fsat_frc_nm,
     $     fsatc_dst_nm,
     $     fsatc_frc_nm,
     $     fsds_dst_nm,
     $     fsds_frc_nm,
     $     fsns_dst_nm,
     $     fsns_frc_nm,
     $     fsnsc_dst_nm,
     $     fsnsc_frc_nm,
     $     fsnt_dst_nm,
     $     fsnt_frc_nm,
     $     fsntc_dst_nm,
     $     fsntc_frc_nm,
     $     ftat_dst_nm,
     $     ftat_frc_nm,
     $     ftns_dst_nm,
     $     ftns_frc_nm,
     $     ftnt_dst_nm,
     $     ftnt_frc_nm,
     $     qrl_dst_nm,
     $     qrl_frc_nm,
     $     qrs_dst_nm,
     $     qrs_frc_nm,
     $     qrt_dst_nm,
     $     qrt_frc_nm
      character*8 flat_dst_nm
      character*8 flat_frc_nm
      character*8 flatc_dst_nm
      character*8 flatc_frc_nm
      character*8 flds_dst_nm
      character*8 flds_frc_nm
      character*8 flns_dst_nm
      character*8 flns_frc_nm
      character*8 flnsc_dst_nm
      character*8 flnsc_frc_nm
      character*8 flnt_dst_nm
      character*8 flnt_frc_nm
      character*8 flntc_dst_nm
      character*8 flntc_frc_nm
      character*8 fsat_dst_nm
      character*8 fsat_frc_nm
      character*8 fsatc_dst_nm
      character*8 fsatc_frc_nm
      character*8 fsds_dst_nm
      character*8 fsds_frc_nm
      character*8 fsns_dst_nm
      character*8 fsns_frc_nm
      character*8 fsnsc_dst_nm
      character*8 fsnsc_frc_nm
      character*8 fsnt_dst_nm
      character*8 fsnt_frc_nm
      character*8 fsntc_dst_nm
      character*8 fsntc_frc_nm
      character*8 ftat_dst_nm
      character*8 ftat_frc_nm
      character*8 ftns_dst_nm
      character*8 ftns_frc_nm
      character*8 ftnt_dst_nm
      character*8 ftnt_frc_nm
      character*8 qrl_dst_nm
      character*8 qrl_frc_nm
      character*8 qrs_dst_nm
      character*8 qrs_frc_nm
      character*8 qrt_dst_nm
      character*8 qrt_frc_nm
