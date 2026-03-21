module sulemis

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Add the sulfur emissions to the tracer array.  This is done here rather
   ! than in the vertical diffusions routines because the emissions are added
   ! to more than just the surface layer (stack emissions above 100 m).
   ! 
   ! Author: module coded by B. Eaton.
   !-----------------------------------------------------------------------

   use precision
   use pmgrid

   implicit none
   save
   private
   public :: &
      addsulemis     ! Add emissions for SO2, SO4, and DMS to the tracer array.

!##############################################################################
contains
!##############################################################################

   subroutine addsulemis( lat, dtime, gravit, rpdel, zi, &
                          chtr, obuf )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Add emissions for SO2, SO4, and DMS to the tracer array.
      !
      ! The emissions are specified from data provided by Benkovitz.
      ! SO2 and SO4 are emitted at 2 levels (below and above 100m).
      ! DMS is emitted from the surface only. 
      ! 
      ! Author: M. Barth
      !-----------------------------------------------------------------------

      use scyc, only: useGEIA
      use dmsbnd, only: dmsbndget
      use soxbnd, only: soxbndget
      use sulbnd, only: sulbndget
      use histout, only: outfld
#ifdef SCYC_MASSBGT
      use massbgt, only: aal
#endif

      implicit none

      integer, intent(in) ::&
         lat                    ! model latitude index
      real(r8), intent(in) ::&
         dtime,                &
         gravit,               &
         rpdel(plond,plev),    &! reciprocal of pdel
         zi(plond,plev)

      real(r8), intent(inout) ::&
         chtr(plond,plev,3),   &! chemical tracer array (SO2,SO4,DMS)
         obuf(*)                ! history output buffer

      ! Local variables:
      integer :: i, m, il, kp
      real(r8) ::&
         scl,                  &! scale factor to put flux in correct units
         sflx(plond,6),        &! sulfur emissions for a given latitude 
                                !  1,4 SO2 emissions at 0m, >100m
                                !  2,5 SO4 emissions at 0m, >100m
                                !  3,6 DMS emissions at 0m, >100m
         tmp(plond,2)
      !-----------------------------------------------------------------------

      ! Get latitude slice of emissions data.
      if ( useGEIA() ) then
         call sulbndget( plond, 6, lat, sflx )
      else
         call soxbndget( lat, plond, tmp )
         call dmsbndget( lat, plond, sflx(1,3) )
         do i = 1, plon
            ! Split SOx emissions, 98% --> SO2, 2% --> SO4
            sflx(i,1) = .98*tmp(i,1) * 2.  ! *2. converts kg S --> kg SO2
            sflx(i,4) = .98*tmp(i,2) * 2.
            sflx(i,2) = .02*tmp(i,1) * 3.  ! *3. converts kg S --> kg SO4
            sflx(i,5) = .02*tmp(i,2) * 3.
         end do
      end if

      ! Put fluxes on history file.
      do i = 1, plon
         tmp(i,1) = sflx(i,1) + sflx(i,4)
         tmp(i,2) = sflx(i,2) + sflx(i,5)
      end do
      call outfld( 'SO2SF   ', tmp(1,1), plond, lat, obuf )
      call outfld( 'SO4SF   ', tmp(1,2), plond, lat, obuf )
      call outfld( 'DMSSF   ', sflx(1,3), plond, lat, obuf )
#ifdef SCYC_MASSBGT
      call aal( 'so2sf', lat, tmp(:,1) )
      call aal( 'so4sf', lat, tmp(:,2) )
      call aal( 'dmssf', lat, sflx(:,3) )
#endif

      ! Add surface flux to lowest tracer level.
      do il = 1, 3
         if( il .eq. 1 ) m = 1
         if( il .eq. 2 ) m = 2
         if( il .eq. 3 ) m = 3
         do i = 1, plon
            scl = gravit*rpdel(i,plev)  ! scales emissions from kg/m2-s to kg/kg-s
            chtr(i,plev,m) = chtr(i,plev,m) + dtime*sflx(i,il)*scl
         end do
      end do

      ! Add "above 100 m emissions" to bottom and next to bottom levels as
      ! appropriate.  (SO2 and SO4 only)
      do il = 4, 5
         if( il .eq. 4 ) m = 1
         if( il .eq. 5 ) m = 2
         do i = 1, plon
            if( 100. .lt. zi(i,plev) ) then
               scl = gravit*rpdel(i,plev)    
               chtr(i,plev,m) = chtr(i,plev,m) + &
                                dtime*sflx(i,il)*scl * (zi(i,plev)-100.)/ &
                                (zi(i,plev-1)-100.)

               scl = gravit*rpdel(i,plev-1)    
               chtr(i,plev-1,m) = chtr(i,plev-1,m) + &
                                  dtime*sflx(i,il)*scl * (zi(i,plev-1)-zi(i,plev))/ &
                                  (zi(i,plev-1)-100.)
            else
               scl = gravit*rpdel(i,plev-1)    
               chtr(i,plev-1,m) = chtr(i,plev-1,m) +dtime*sflx(i,il)*scl 
            endif
         end do
      end do

   end subroutine addsulemis

!#######################################################################

end module sulemis
