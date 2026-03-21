module error_messages

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! General purpose routines for issuing error messages.
   ! 
   ! Author: B. Eaton
   ! 
   !----------------------------------------------------------------------- 

   use netcdf

   implicit none
   save
   private
   public :: &
      alloc_err,      &! Issue error message after non-zero return from an allocate statement.
      dealloc_err,    &! Issue error message after non-zero return from a deallocate statement.
      handle_ncerr     ! Handle error returns from netCDF library procedures.

!##############################################################################
contains
!##############################################################################

   subroutine alloc_err( istat, routine, name, nelem )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from an allocate statement.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) ::&
         istat           ! status from allocate statement
      character(len=*), intent(in) ::&
         routine,       &! routine that called allocate
         name            ! name of array
      integer, intent(in) ::&
         nelem           ! number of elements attempted to allocate
      !-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         write(*,*)'error_message.F90: ERROR trying to allocate memory in routine: ' &
                   //trim(routine)
         write(*,*)'  Variable name: '//trim(name)
         write(*,*)'  Number of elements: ',nelem
         call endrun
      end if

   end subroutine alloc_err

!##############################################################################

   subroutine dealloc_err( istat, routine, name )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from an deallocate statement.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) ::&
         istat           ! status from allocate statement
      character(len=*), intent(in) ::&
         routine,       &! routine that called allocate
         name            ! name of array
      !-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         write(*,*)'ERROR trying to deallocate memory in routine: ' &
                   //trim(routine)
         write(*,*)'  Variable name: '//trim(name)
         call endrun
      end if

   end subroutine dealloc_err

!##############################################################################

   subroutine handle_ncerr( ret, mes )
      
      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Check netCDF library function return code.  If error detected 
      ! issue error message then abort.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      integer, intent(in) ::&
         ret                 ! return code from netCDF library routine
      character(len=*), intent(in) ::&
         mes                 ! message to be printed if error detected
      !-----------------------------------------------------------------------

      if ( ret .ne. NF90_NOERR ) then
         write(6,*) 'error detected by handle_ncerr'
         write(6,*) mes
         write(6,*) nf90_strerror( ret )
         call endrun
      endif

!      if ( ret .ne. NF_NOERR ) then
!         write(6,*) mes
!         write(6,*) nf_strerror( ret )
!         call endrun
!      endif

   end subroutine handle_ncerr

!##############################################################################

end module error_messages
