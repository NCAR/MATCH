c $Id$
      subroutine setr( n, val, f )

      use precision
      implicit none

      integer n
      real(r8) f(n), val
      integer i

      do i = 1, n
         f(i) = val
      end do

      return
      end
c##############################################################################

      subroutine seti( n, val, f )
      integer n
      integer f(n), val
      integer i

      do i = 1, n
         f(i) = val
      end do

      return
      end
c##############################################################################

      subroutine copyr( n, from, to )

      use precision
      implicit none

      integer n
      real(r8) from(n), to(n)
      integer i

      do i = 1, n
         to(i) = from(i)
      end do

      return
      end
c##############################################################################

      subroutine scale( sc, nd1, nx, ny, f )
      use precision
      implicit none
      integer nd1, nx, ny
      integer i, j
      real(r8) sc, f(nd1,ny)
      do j = 1, ny
         do i = 1, nx
            f(i,j) = sc*f(i,j)
         end do
      end do
      return
      end
c##############################################################################

      subroutine calctend( f0, f1, dt, tend )

c Calculate a process tendency.

      use precision
      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  f0(plond,plev)       ! before process values
     $, f1(plond,plev)       ! after process values
     $, dt                   ! timestep size (s)

c Output argument:
      real(r8)
     $  tend(plond,plev)     ! process tendency

c Local variables:
      integer i, k
c-----------------------------------------------------------------------

      do k = 1, plev
         do i = 1, plon
            tend(i,k) = (f1(i,k) - f0(i,k))/dt
         end do
      end do

      return
      end

c##############################################################################

      subroutine cpslic( n, from, to )

c copy latitude slices

      use precision
      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  n                    ! number of slices
      real(r8)
     $  from(plond,plev,n)

c Output argument:
      real(r8)
     $  to(plond,plev,n)

c Local variables:
      integer i, k, m
c-----------------------------------------------------------------------

      do m = 1, n
         do k = 1, plev
            do i = 1, plon
               to(i,k,m) = from(i,k,m)
            end do
         end do
      end do
      return
      end

c##############################################################################

      subroutine alloc_err( istat, routine, name, nelem )

c Issue error message after non-zero return from an allocate statement.

      implicit none

c Input arguments:
      integer :: istat   ! status from allocate statement
      character(len=*) ::
     $  routine          ! routine that called allocate
     $, name             ! name of array
      integer :: nelem   ! number of elements attempted to allocate
c-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         write(*,*)'utils.F: ERROR trying to allocate memory in routine: '
     $          //trim(routine)
         write(*,*)'  Variable name: '//trim(name) 
         write(*,*)'  Number of elements: ',nelem
         call endrun
      end if
      return
      end

c##############################################################################

      subroutine endrun

c Abort the model.

      implicit none
c ----------------------------------------------------------------------

      call abort
      end

c##############################################################################

