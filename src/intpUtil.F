c $Id$
      subroutine chktime( time, nrec )

c Make sure the time coordinate looks like calander day, and is increasing.
c Calendar day can either start with 1 Jan 0Z = day 1.0  or
c  1 Jan 0Z = day 0.0


      use precision

      implicit none

c Input arguments:

      integer
     $  nrec
      real(r8)
     $  time(nrec)

c Local varibles:

      integer
     $  i
c-----------------------------------------------------------------------

      if ( time(1) .lt. 0. .or. time(1) .ge. 366. ) then
         write(*,*)'chktime: illegal time coordinate ',time(1)
         call endrun
      end if
      do i = 2, nrec
         if ( time(i) .lt. 0. .or. time(i) .ge. 366. ) then
            write(*,*)'chktime: illegal time coordinate ', time(i)
            call endrun
         end if
         if ( time(i) .le. time(i-1) ) then
            write(*,*)'chktime: time not increasing ', time(i-1)
     $,               time(i)
            write(*,*)'DRBDBG ',time
            call endrun
         end if
      end do

      return
      end

c#######################################################################

      subroutine findplb( x, nx, xval, index )

c "find periodic lower bound"
c Search the input array for the lower bound of the interval that
c contains the input value.  The returned index satifies:
c x(index) .le. xval .lt. x(index+1)
c Assume the array represents values in one cycle of a periodic coordinate.
c So, if xval .lt. x(1), then index=nx.


      use precision

      implicit none

c Input arguments:
      integer
     $  nx
      real(r8)
     $  x(nx)  ! strictly increasing array
     $, xval

c Output arguments:
      integer
     $  index

c Local variables:
      integer i
c-----------------------------------------------------------------------

      if ( xval .lt. x(1) .or. xval .ge. x(nx) ) then
         index = nx
         return
      end if

      do i = 2, nx
         if ( xval .lt. x(i) ) then
            index = i-1
            return
         end if
      end do

      end

c#######################################################################

      subroutine linintp( npts, t1, t2, tint, f1, f2, fint )

c Linearly interpolate between f1(t1) and f2(t2) to fint(tint).


      use precision

      implicit none

c Input arguments:
      integer
     $  npts
      real(r8)
     $  t1               ! time level of f1
     $, t2               ! time level of f2
     $, tint             ! interpolant time
     $, f1(npts)         ! field at time t1
     $, f2(npts)         ! field at time t2

c Output arguments:
      real(r8)
     $  fint(npts)       ! field at time tint

c Local variables:
      integer i
      real(r8) factor
c------------------------------------------------------------------------------

      factor = ( tint - t1 )/( t2 - t1)

      do i = 1, npts
         fint(i) = f1(i) + ( f2(i) - f1(i) )*factor
      end do

      return
      end

c#######################################################################
