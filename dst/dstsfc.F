c     $Header$
      
c     Purpose: dstsfc.F contains interface routines to access 
c     2-D surface variables stored in common blocks in dstsfc.h
      
      subroutine soi_txt_get(
     $     lat_idx,             ! I [idx] Latitude index
     $     lon_nbr_out,         ! I [nbr] Dimension size (normally plon)
     $     lnd_frc_dry_out,     ! O [frc] Dry land fraction
     $     mbl_bsn_fct_out,     ! O [frc] Topographic source limitation
     $     mss_frc_CaCO3_out,   ! O [frc] Mass fraction of CaCO3
     $     mss_frc_cly_out,     ! O [frc] Mass fraction of clay
     $     mss_frc_snd_out)     ! O [frc] Mass fraction of sand
c     Return latitude slice of soil texture
c     soi_txt_get() is called by dst_mbl()

      use precision

      implicit none
#include <dstsfc.h> /* Surface fields for dust parameterization */
c     Parameters
      integer,parameter::lon_nbr=PLON ! [nbr] Number of longitudes
c     Input
      integer,intent(in)::lat_idx ! I [idx] Latitude index
      integer,intent(in)::lon_nbr_out ! I [nbr] Dimension size (normally plon)
c     Output
      real(r8),intent(out)::lnd_frc_dry_out(lon_nbr_out) ! O [frc] Dry land fraction
      real(r8),intent(out)::mbl_bsn_fct_out(lon_nbr_out) ! O [frc] Topographic source limitation
      real(r8),intent(out)::mss_frc_CaCO3_out(lon_nbr_out) ! O [frc] Mass fraction of CaCO3
      real(r8),intent(out)::mss_frc_cly_out(lon_nbr_out) ! O [frc] Mass fraction of clay
      real(r8),intent(out)::mss_frc_snd_out(lon_nbr_out) ! O [frc] Mass fraction of sand
c     Local
      integer lon_idx           ! [idx] Counting index for lon
c     Main code
      if (lon_nbr_out > lon_nbr) then
         write (6,'(a)') 'dst: ERROR soi_txt_get() reports lon_nbr_out > lon_nbr'
         stop
      endif                     ! endif err
      
      do lon_idx=1,lon_nbr
         lnd_frc_dry_out(lon_idx)=lnd_frc_dry(lon_idx,lat_idx)
         mbl_bsn_fct_out(lon_idx)=mbl_bsn_fct(lon_idx,lat_idx)
c     fxm: CaCO3 currently has missing value of 1.0e36 which causes problems
         if (mss_frc_CaCO3(lon_idx,lat_idx) <= 1.0) then
            mss_frc_CaCO3_out(lon_idx)=mss_frc_CaCO3(lon_idx,lat_idx)
         else
            mss_frc_CaCO3_out(lon_idx)=0.0
         endif
c     fxm: Temporarily set mss_frc_cly used in mobilization to globally uniform SGS value of 0.20
c     mss_frc_cly_out(lon_idx)=mss_frc_cly(lon_idx,lat_idx)
         mss_frc_cly_out(lon_idx)=0.20
         mss_frc_snd_out(lon_idx)=mss_frc_snd(lon_idx,lat_idx)
      end do                    ! end loop over lon
      
      return
      end                       ! end soi_txt_get()
      
      subroutine sfc_typ_get(
     $     lat_idx,             ! I [idx] Latitude index
     $     lon_nbr_out,         ! I [nbr] Dimension size (normally plon)
     $     sfc_typ_out)         ! O [idx] LSM surface type (0..28)
c     Return latitude slice of soil texture
c     sfc_typ_get() is called by dst_mbl(), dst_dps_dry()

      use precision

      implicit none
#include <dstsfc.h> /* Surface fields for dust parameterization */
c     Parameters
      integer,parameter::lon_nbr=PLON ! [nbr] Number of longitudes
c     Input
      integer,intent(in)::lat_idx ! I [idx] Latitude index
      integer,intent(in)::lon_nbr_out ! I [nbr] Dimension size (normally plon)
c     Output
      integer,intent(out)::sfc_typ_out(lon_nbr_out) ! O [idx] LSM surface type (0..28)
c     Local
      integer lon_idx           ! [idx] Counting index for lon
c     Main code
      if (lon_nbr_out > lon_nbr) then
         write (6,'(a)') 'dst: ERROR sfc_typ_get() reports lon_nbr_out > lon_nbr'
         stop
      endif                     ! endif err
      do lon_idx=1,lon_nbr
         sfc_typ_out(lon_idx)=sfc_typ(lon_idx,lat_idx)
      end do                    ! end loop over lon
      return
      end                       ! end sfc_typ_get()
      
      subroutine vwc_sfc_get(
     $     lat_idx,             ! I [idx] Latitude index
     $     lon_nbr_out,         ! I [nbr] Dimension size (normally plon)
     $     vwc_sfc_out)         ! O [m3 m-3] Volumetric water content
c     Return latitude slice of soil texture
c     vwc_sfc_get() is called by dst_mbl()

      use precision

      implicit none
#include <dstsfc.h> /* Surface fields for dust parameterization */
c     Parameters
      integer,parameter::lon_nbr=PLON ! [nbr] Number of longitudes
c     Input
      integer,intent(in)::lat_idx ! I [idx] Latitude index
      integer,intent(in)::lon_nbr_out ! I [nbr] Dimension size (normally plon)
c     Output
      real(r8),intent(out)::vwc_sfc_out(lon_nbr_out) ! O [m3 m-3] Volumetric water content
c     Local
      integer lon_idx           ! [idx] Counting index for lon
c     Main code
      if (lon_nbr_out > lon_nbr) then
         write (6,'(a)') 'dst: ERROR vwc_sfc_get() reports lon_nbr_out > lon_nbr'
         stop
      endif                     ! endif err
      do lon_idx=1,lon_nbr
         vwc_sfc_out(lon_idx)=vwc_sfc(lon_idx,lat_idx)
      end do                    ! end loop over lon
      return
      end                       ! end vwc_sfc_get()
      
      subroutine tpt_gnd_soi_get(
     $     lat_idx,             ! I [idx] Latitude index
     $     lon_nbr_out,         ! I [nbr] Dimension size (normally plon)
     $     tpt_gnd_out,         ! O [K] Ground temperature
     $     tpt_soi_out)         ! O [K] Soil temperature
c     Return latitude slice of soil temperature and ground temperature
c     tpt_gnd_soi_get() is called by dst_mbl()

      use precision

      implicit none
#include <dstsfc.h> /* Surface fields for dust parameterization */
c     Parameters
      integer,parameter::lon_nbr=PLON ! [nbr] Number of longitudes
c     Input
      integer,intent(in)::lat_idx ! I [idx] Latitude index
      integer,intent(in)::lon_nbr_out ! I [nbr] Dimension size (normally plon)
c     Output
      real(r8),intent(out)::tpt_gnd_out(lon_nbr_out) ! O [K] Ground temperature
      real(r8),intent(out)::tpt_soi_out(lon_nbr_out) ! O [K] Soil temperature
c     Local
      integer lon_idx           ! [idx] Counting index for lon
c     Main code
      if (lon_nbr_out > lon_nbr) then
         write (6,'(a)') 'dst: ERROR tpt_gnd_soi_get() reports lon_nbr_out > lon_nbr'
         stop
      endif                     ! endif err
      do lon_idx=1,lon_nbr
         tpt_gnd_out(lon_idx)=tpt_gnd(lon_idx,lat_idx)
         tpt_soi_out(lon_idx)=tpt_soi(lon_idx,lat_idx)
      end do                    ! end loop over lon
      return
      end                       ! end tpt_gnd_soi_get()
      
      subroutine flx_rad_sfc_get(
     $     lat_idx,             ! I [idx] Latitude index
     $     lon_nbr_out,         ! I [nbr] Dimension size (normally plon)
     $     flx_LW_dwn_sfc_out,  ! O [W m-2] Longwave downwelling flux at surface
     $     flx_SW_abs_sfc_out)  ! O [W m-2] Solar flux absorbed by ground
c     Return latitude slice of surface radiative fluxes
c     flx_rad_sfc_get() is called by dst_mbl()

      use precision

      implicit none
#include <dstsfc.h> /* Surface fields for dust parameterization */
c     Parameters
      integer,parameter::lon_nbr=PLON ! [nbr] Number of longitudes
c     Input
      integer,intent(in)::lat_idx ! I [idx] Latitude index
      integer,intent(in)::lon_nbr_out ! I [nbr] Dimension size (normally plon)
c     Output
      real(r8),intent(out)::flx_LW_dwn_sfc_out(lon_nbr_out) ! O [W m-2] Longwave downwelling flux at surface
      real(r8),intent(out)::flx_SW_abs_sfc_out(lon_nbr_out) ! O [W m-2] Solar flux absorbed by ground
c     Local
      integer lon_idx           ! [idx] Counting index for lon
c     Main code
      if (lon_nbr_out > lon_nbr) then
         write (6,'(a)') 'dst: ERROR flx_rad_sfc_get() reports lon_nbr_out > lon_nbr'
         stop
      endif                     ! endif err
      do lon_idx=1,lon_nbr
         flx_LW_dwn_sfc_out(lon_idx)=flx_LW_dwn_sfc(lon_idx,lat_idx)
         flx_SW_abs_sfc_out(lon_idx)=flx_SW_abs_sfc(lon_idx,lat_idx)
      end do                    ! end loop over lon
      return
      end                       ! end flx_rad_sfc_get()
      
      subroutine flx_rad_sfc_set(
     $     lat_idx,             ! I [idx] Latitude index
     $     lon_nbr_in,          ! I [nbr] Dimension size (normally plon)
     $     flx_LW_dwn_sfc_in,   ! I [W m-2] Longwave downwelling flux at surface
     $     flx_SW_abs_sfc_in)   ! I [W m-2] Solar flux absorbed by ground
c     Set latitude slice of surface radiative fluxes
c     flx_rad_sfc_set() is called by CCM:physics/tphysbc(), MATCH:physlic()

      use precision

      implicit none
#include <dstsfc.h> /* Surface fields for dust parameterization */
c     Parameters
      integer,parameter::lon_nbr=PLON ! [nbr] Number of longitudes
c     Input
      integer,intent(in)::lat_idx ! I [idx] Latitude index
      integer,intent(in)::lon_nbr_in ! I [nbr] Dimension size (normally plond)
      real(r8),intent(in)::flx_LW_dwn_sfc_in(lon_nbr_in) ! I [W m-2] Longwave downwelling flux at surface
      real(r8),intent(in)::flx_SW_abs_sfc_in(lon_nbr_in) ! I [W m-2] Solar flux absorbed by ground
c     Local
      integer lon_idx           ! [idx] Counting index for lon
c     Main Code
      if (lon_nbr_in < lon_nbr) then
         write (6,'(a)') 'dst: ERROR flx_rad_sfc_set() reports lon_nbr_in < lon_nbr'
         stop
      endif                     ! endif err
      do lon_idx=1,lon_nbr
         flx_LW_dwn_sfc(lon_idx,lat_idx)=flx_LW_dwn_sfc_in(lon_idx)
         flx_SW_abs_sfc(lon_idx,lat_idx)=flx_SW_abs_sfc_in(lon_idx)
      end do                    ! end loop over lon
      return
      end                       ! end flx_rad_sfc_set()
      
      subroutine vwc_sfc_set(
     $     lat_idx,             ! I [idx] Latitude index
     $     lon_nbr_in,          ! I [nbr] Dimension size (normally plon)
     $     vwc_sfc_in)          ! I [m3 m-3] Volumetric water content
c     Set latitude slice of surface volumetric water content
c     vwc_sfc_set() is called by MATCH:physlic()

      use precision

      implicit none
#include <dstsfc.h> /* Surface fields for dust parameterization */
c     Parameters
      integer,parameter::lon_nbr=PLON ! [nbr] Number of longitudes
c     Input
      integer,intent(in)::lat_idx ! I [idx] Latitude index
      integer,intent(in)::lon_nbr_in ! I [nbr] Dimension size (normally plond)
      real(r8),intent(in)::vwc_sfc_in(lon_nbr_in) ! I [m3 m-3] Volumetric water content
c     Local
      integer lon_idx           ! [idx] Counting index for lon
c     Main Code
      if (lon_nbr_in < lon_nbr) then
         write (6,'(a)') 'dst: ERROR vwc_sfc_set() reports lon_nbr_in < lon_nbr'
         stop
      endif                     ! endif err
      do lon_idx=1,lon_nbr
         vwc_sfc(lon_idx,lat_idx)=vwc_sfc_in(lon_idx) ! [m3 m-3]
      end do                    ! end loop over lon
      return
      end                       ! end vwc_sfc_set()
      
      subroutine tpt_gnd_soi_set(
     $     lat_idx,             ! I [idx] Latitude index
     $     lon_nbr_in,          ! I [nbr] Dimension size (normally plon)
     $     tpt_gnd_in,          ! I [K] Ground temperature
     $     tpt_soi_in)          ! I [K] Soil temperature
c     Set latitude slice of surface and soil temperatures
c     tpt_gnd_soi_set() is called by MATCH:physlic()

      use precision

      implicit none
#include <dstsfc.h> /* Surface fields for dust parameterization */
c     Parameters
      integer,parameter::lon_nbr=PLON ! [nbr] Number of longitudes
c     Input
      integer,intent(in)::lat_idx ! I [idx] Latitude index
      integer,intent(in)::lon_nbr_in ! I [nbr] Dimension size (normally plond)
      real(r8),intent(in)::tpt_gnd_in(lon_nbr_in) ! I [K] Ground temperature
      real(r8),intent(in)::tpt_soi_in(lon_nbr_in) ! I [K] Soil temperature
c     Local
      integer lon_idx           ! [idx] Counting index for lon
c     Main Code
      if (lon_nbr_in < lon_nbr) then
         write (6,'(a)') 'dst: ERROR tpt_gnd_soi_set() reports lon_nbr_in < lon_nbr'
         stop
      endif                     ! endif err
      do lon_idx=1,lon_nbr
         tpt_gnd(lon_idx,lat_idx)=tpt_gnd_in(lon_idx)
         tpt_soi(lon_idx,lat_idx)=tpt_soi_in(lon_idx)
      end do                    ! end loop over lon
      return
      end                       ! end tpt_gnd_soi_set()
