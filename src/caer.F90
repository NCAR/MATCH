#include <params.h>

module caer

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Simple carbon aerosol model for MATCH.
   !
   ! Author: B. Eaton
   !-----------------------------------------------------------------------

   use pmgrid
   use precision

   implicit none
   save
   private
   public :: &
      caerini,    &! initialization
      caersf,     &! return surface fluxes
      caercv,     &! 1st order chemistry
      docaer,     &! returns .true. when carbon aerosol species are present
      caer_ncomp, &! returns number of carbon aerosol components
      caer_idx1    ! returns constituent index of 1st carbon aerosol component

   real(r8), parameter ::&
      freq = 1./(1.2*24.*3600.),  &! 1 / (e-folding time in seconds)
      mwrat = 1.                   ! ratio of mol. wts of cphil to cphob

   integer, parameter :: &
      ncomp = CAER_NBR,    &! number of carbon aerosol components
      idx1 = CAER_IDX_SRT   ! start constituent index for first carbon aerosol component

!##############################################################################
contains
!##############################################################################

   subroutine caerini()

      !----------------------------------------------------------------------- 
      ! Purpose: Initialize carbon aerosol module.
      ! 
      ! Author: B. Eaton
      !-----------------------------------------------------------------------

      implicit none

      !-----------------------------------------------------------------------

      if ( docaer() ) then
         write(6,*)'Carbon aerosol configured with ',ncomp,' components.'
         write(6,*)'  Hydrophobic components converted to hydrophilic with 1.2 day e-folding time.'
      else
         write(6,*)'Carbon aerosols disabled.'
      end if

   end subroutine caerini

!#######################################################################

   subroutine caersf( nemis, lat, rlat, oc_rgn_scale, bc_rgn_scale, sflx )

      ! Set carbon aerosol surface fluxes.

      use caerbnd, only: caerbndget

      implicit none

      integer, intent(in) ::&
         nemis,             &! number of carbon aerosol species with emissions
         lat                 ! model latitude index

      real(r8), intent(in) ::&
         rlat,              &! [rdn] latitude of this row
         oc_rgn_scale(9),   &! per-region OC multiplier
         bc_rgn_scale(9)     ! per-region BC multiplier

      real(r8), intent(out) ::&
         sflx(plond,nemis)   ! surface flux in kg/m^2/s

      ! Local variables:
      integer i, m
      real(r8) ::&
         tmp(plond,5)
      real(r8) :: lat_deg, lon_deg, oc_scl, bc_scl
      integer :: rgn
      !------------------------------------------------------------------------------

      ! Get latitude slice of carbon aerosol emissions data.
      ! tmp(:,1) BBOCSF  tmp(:,2) BBBCSF  tmp(:,3) FFOCSF
      ! tmp(:,4) FFBCSF  tmp(:,5) NOCSF
      call caerbndget( lat, plond, tmp )

      ! Apply per-region OC and BC tuning to the input components.
      ! Regions tested in order; first match wins (FF boxes precede the
      ! Boreal NH catch-all so they take priority in overlap zones).
      lat_deg = rlat * 57.29577951308232_r8
      do i = 1, plon
         lon_deg = (i-1) * 360.0_r8 / plon
         rgn = 0
         if      ( lat_deg >= -15.0 .and. lat_deg <= 10.0 .and. &
                   lon_deg >= 280.0 .and. lon_deg <= 310.0 ) then
            rgn = 1                            ! Amazon
         else if ( lat_deg >= -35.0 .and. lat_deg <=  5.0 .and. &
                   lon_deg >=  10.0 .and. lon_deg <= 45.0 ) then
            rgn = 2                            ! S. Africa
         else if ( lat_deg >= -10.0 .and. lat_deg <= 20.0 .and. &
                   lon_deg >=  95.0 .and. lon_deg <= 150.0 ) then
            rgn = 3                            ! SE Asia/Indonesia
         else if ( lat_deg >= -40.0 .and. lat_deg <= -10.0 .and. &
                   lon_deg >= 110.0 .and. lon_deg <= 155.0 ) then
            rgn = 4                            ! Australia
         else if ( lat_deg >=  20.0 .and. lat_deg <= 50.0 .and. &
                   lon_deg >= 100.0 .and. lon_deg <= 145.0 ) then
            rgn = 5                            ! E. Asia
         else if ( lat_deg >=   5.0 .and. lat_deg <= 35.0 .and. &
                   lon_deg >=  65.0 .and. lon_deg <=  95.0 ) then
            rgn = 6                            ! S. Asia
         else if ( lat_deg >=  35.0 .and. lat_deg <= 65.0 .and. &
                   ( lon_deg >= 350.0 .or. lon_deg <= 50.0 ) ) then
            rgn = 7                            ! Europe
         else if ( lat_deg >=  25.0 .and. lat_deg <= 60.0 .and. &
                   lon_deg >= 235.0 .and. lon_deg <= 295.0 ) then
            rgn = 8                            ! N. America
         else if ( lat_deg >=  50.0 .and. lat_deg <= 75.0 ) then
            rgn = 9                            ! Boreal NH (catch-all)
         end if

         oc_scl = 1.0_r8
         bc_scl = 1.0_r8
         if ( rgn > 0 ) then
            oc_scl = oc_rgn_scale(rgn)
            bc_scl = bc_rgn_scale(rgn)
         end if
         tmp(i,1) = oc_scl * tmp(i,1)          ! BBOCSF
         tmp(i,3) = oc_scl * tmp(i,3)          ! FFOCSF
         tmp(i,5) = oc_scl * tmp(i,5)          ! NOCSF
         tmp(i,2) = bc_scl * tmp(i,2)          ! BBBCSF
         tmp(i,4) = bc_scl * tmp(i,4)          ! FFBCSF
      end do

      if ( nemis .eq. 1 ) then
         do i = 1, plon
            sflx(i,1) = tmp(i,1)+tmp(i,2)+tmp(i,3)+tmp(i,4)+tmp(i,5)  ! total carbon
         end do
      else if ( nemis .eq. 2 ) then
         do i = 1, plon
            sflx(i,1) = tmp(i,1) + tmp(i,3) + tmp(i,5)  ! organic carbon
            sflx(i,2) = tmp(i,2) + tmp(i,4)             ! black carbon
         end do
      else if ( nemis .eq. 5 ) then
         do m = 1, 5
            do i = 1, plon
               sflx(i,m) = tmp(i,m)
            end do
         end do
      end if

   end subroutine caersf

!##############################################################################

   subroutine caercv( A, deltat, Asrc, Bsrc )

      ! A undergoes 1st order decay to B.

      implicit none

      real(r8), intent(in) ::&
         A(plond,plev),         &! mixing ratio (kg A/(kg moist air))
         deltat                  ! time step in seconds

      real(r8), intent(out) ::&
         Asrc(plond,plev),      &! conversion rate (kg A /(s kg moist air))
         Bsrc(plond,plev)        ! conversion rate (kg B /(s kg moist air))

      ! Local variables:
      integer i, k
      !-----------------------------------------------------------------------

      ! calculate tendencies using Backward Euler

      do k = 1,plev
         do i = 1,plon
            Asrc(i,k) = -A(i,k)*freq / (1. + freq*deltat)
            Bsrc(i,k) = -Asrc(i,k)*mwrat
         end do
      end do

   end subroutine caercv

!#######################################################################

   logical function docaer ()

      !----------------------------------------------------------------------- 
      ! Purpose: Return .true. when carbon aerosols are active.
      ! 
      ! Author: B. Eaton
      !-----------------------------------------------------------------------

      implicit none

      docaer = ( ncomp > 0 )

   end function docaer

!#######################################################################

   integer function caer_ncomp()
      implicit none
      caer_ncomp = ncomp
   end function caer_ncomp

!#######################################################################

   integer function caer_idx1()
      implicit none
      caer_idx1 = idx1
   end function caer_idx1

!#######################################################################

end module caer
