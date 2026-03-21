c $Id$
c $Id$
      subroutine nfPutVarReal( ncid, vid, var )

      use precision
      use netcdf

      implicit none
c-----------------------------------------------------------------------
#include <netcdf.inc>
c-----------------------------------------------------------------------

c Input
      integer ncid, vid

      real(r8) var

c Local
      integer status
c-----------------------------------------------------------------------

      status = NF90_PUT_VAR( ncid, vid, var )
      call handle_ncerr( status, 'nfPutVarReal: error return:' )
      return
      end

c#######################################################################

      subroutine handle_err(status)

c     NetCDF library error handling routine.

      implicit none
#include <netcdf.inc>
      integer status

      if (status .ne. NF_NOERR) then
         write(*,*)'handle_err:', nf_strerror(status)
         stop
      endif

      end

c#######################################################################

      subroutine handle_ncerr( ret, mes )

c Check netCDF library function return code.  If error detected 
c issue error message then abort.

      implicit none
c-----------------------------------------------------------------------
#include <netcdf.inc>
c-----------------------------------------------------------------------

c Input arguments:
      integer ret            ! return code from netCDF library routine
      character*(*) mes      ! message to be printed if error detected
c-----------------------------------------------------------------------

      if ( ret .ne. NF_NOERR ) then
         write(*,*) mes
         write(*,*) nf_strerror( ret )
         call abort
      endif

      return
      end
