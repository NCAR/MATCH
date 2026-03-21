module dmsbnd

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! This code does time interpolation for DMS boundary data in a netCDF
   ! file.  Assumptions on the data in the netCDF file are:
   ! 1. Coordinates are ordered (lon,lat,time)
   ! 2. The time coordinate is in days, and the data is assumed to be periodic
   !    annually.
   !
   ! It is assumed that the model calling this interface has been
   ! compiled so that 8 byte real data are being used.  On non-CRAY
   ! machines this implies compiling with a "-r8" flag.
   ! 
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 

   use precision
   use pmgrid
   use netcdf

   implicit none
   save
   private
   public :: &
      dmsbndini,      &! initialize dmsbnd module
      dmsbndint,      &! interpolate dmsbnd data to requested date/time
      dmsbndget        ! return latitude slice data at current date/time

   ! private module data

   real(r8), allocatable, dimension(:) :: &
      time            ! time coordinate (calander days + frac)
   real(r8), dimension(plon,plat,2) :: &
      dmsin           ! input data
   real(r8), dimension(plon,plat) :: &
      dms             ! interpolated data

   integer :: &
      ncid,          &! ID for netCDF file
      nrec,          &! number of records (time samples)
      lotim,         &! time(lotim) .le. current time
      hitim,         &! current time .lt. time(hitim)
      loin,          &! index into input data array containing time(lotim) data
      hiin,          &! index into input data array containing time(hitim) data
      start(3),      &! start vector for netCDF hyperslabs
      count(3)        ! count vector for netCDF hyperslabs

!##############################################################################
contains
!##############################################################################

   subroutine dmsbndini( calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Open netCDF file containing DMS emissions data.  Initialize arrays
      ! with the data to be interpolated to the current time.
      !
      ! It is assumed that the time coordinate is increasing and represents
      ! calendar days (range = [1.,366.)).
      ! 
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: alloc_err, handle_ncerr

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         calday  ! current time in calendar days + fraction.

      ! Local variables:
      integer ::&
         did,   &
         istat, &
         recid, &! record ID
         nlon,  &
         vid
      !-----------------------------------------------------------------------

      start(1) = 1
      start(2) = 1
      start(3) = 1
      count(1) = plon
      count(2) = plat
      count(3) = 1

      ! Open file.
      call handle_ncerr( nf90_open( 'DMS_emissions.nc', NF_NOWRITE, ncid ), &
         'dmsbndini: error opening file '//'DMS_emissions.nc' )

      ! get the record id
      call handle_ncerr( nf90_inquire( ncid, unlimiteddimid=recid), &
         'dmsbndini: no record variables ' )

      ! Check that input data is a right resolution.
      call handle_ncerr( nf90_inq_dimid( ncid, 'lon', did ), 'dmsbndini: ' )
      call handle_ncerr( nf90_inquire_dimension( ncid, did, len=nlon ), 'dmsbndini: ' )
      if ( nlon .ne. plon ) then
         write(6,*)'dmsbndini: model plon = ',plon,', dataset nlon = ',nlon
         call endrun
      end if

      ! Get size of unlimited dimension.
      call handle_ncerr( nf90_inquire_dimension( ncid, recid, len=nrec ), 'dmsbndini: ' )

      ! Allocate space for time coordinate data.
      allocate( time(nrec), stat=istat )
      call alloc_err( istat, 'dmsbndini', 'time', nrec )

      ! Get time coordinate.
      call handle_ncerr( nf90_inq_varid( ncid, 'time', vid ), &
         'dmsbndini: cannot find time coordinate variable' )
      call handle_ncerr( nf90_get_var( ncid, vid, time ), &
         'dmsbndini: error getting time coordinate data' )

      ! Make sure the time coordinate looks like calander day, and is
      ! increasing.
      call chktime( time, nrec )

      ! Find indices for time samples that bound the current time.
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      ! Read data.
      loin = 1
      hiin = 2

      call handle_ncerr( nf90_inq_varid( ncid, 'DMS', vid ), &
         'dmsbndini: cannot find variable '//'DMS' )

      start(3) = lotim
      call handle_ncerr( nf90_get_var( ncid, vid,  dmsin(:,:,loin), start, count ), &
         'dmsbndini: cannot read data for '//'DMS' )

      start(3) = hitim
      call handle_ncerr( nf90_get_var( ncid, vid, dmsin(:,:,hiin), start, count ), &
         'dmsbndini: cannot read data for '//'DMS' )

      write(6,*)'dmsbndini: calendar day = ',calday, ' : read data for days ', &
         time(lotim), ' and ',time(hitim)

   end subroutine dmsbndini

!#######################################################################

   subroutine dmsbndint( calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Interpolate DMS data to the current time.  Update the input data
      ! as necessary.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: handle_ncerr

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         calday  ! current time in calendar days + fraction.

      ! Local variables:
      integer ::&
         oldhitim,        &
         vid
      real(r8) ::&
         dt, dt1, tint
      !-----------------------------------------------------------------------

      ! Check to see if model time is still bounded by dataset times.
      oldhitim = hitim
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      if ( hitim .ne. oldhitim ) then
         ! Read in new hitim data.  Replace old lotim data.
         loin = hiin
         hiin = mod( loin, 2 ) + 1
         start(3) = hitim
         call handle_ncerr( nf90_inq_varid( ncid, 'DMS', vid ), &
            'dmsbndint: cannot find variable '//'DMS' )
         call handle_ncerr( nf90_get_var( ncid, vid, dmsin(:,:,hiin), start, count ), &
            'dmsbndint: cannot read data for '//'DMS' )
         write(6,*)'dmsbndint: read data for day ',time(hitim)

         if ( lotim .ne. oldhitim ) then
            ! Read in new lotim data.  Replace old hitim data.
            start(3) = lotim
            call handle_ncerr( nf90_inq_varid( ncid, 'DMS', vid ), &
               'dmsbndint: cannot find variable '//'DMS' )
            call handle_ncerr( nf90_get_var( ncid, vid, dmsin(:,:,loin), start, count), &
               'dmsbndint: cannot read data for '//'DMS' )
            write(6,*)'dmsbndint: read data for day ',time(lotim)
         end if

      end if


      ! Linear interpolation...  Start by computing the number of days between
      !                          the lower and upper bounds, and days between
      !                          the model time and lower bound.

      if( time(hitim) .lt. time(lotim) )then
         dt = 365. - time(lotim) + time(hitim)
         if( calday .le. time(hitim) )then
            dt1 = 365. - time(lotim) + calday
         else
            dt1 = calday - time(lotim)
         end if
      else
         dt = time(hitim) - time(lotim)
         dt1 = calday - time(lotim)
      end if

      tint = dt1/dt
      call linintp( plon*plat, 0._r8, 1._r8, tint, dmsin(1,1,loin), &
                    dmsin(1,1,hiin), dms )

   end subroutine dmsbndint

!#######################################################################

   subroutine dmsbndget( lat, nd1, x )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Return DMS emission data for the requested latitude.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) ::&
         lat,        &! requested latitude index
         nd1          ! dimension for output array

      real(r8), intent(out) ::&
         x(nd1)       ! DMS emissions (kg DMS/m2/s)

      ! Local variables:
      integer ::&
         i
      !-----------------------------------------------------------------------

      do i = 1, plon
         x(i) = dms(i,lat)
      end do

   end subroutine dmsbndget

!#######################################################################

end module dmsbnd
