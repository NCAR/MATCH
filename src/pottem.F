c $Id$
      subroutine inipottem( xcappa )

c Initialize /pottem0/.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pottem.com>
c-----------------------------------------------------------------------

c Input arguments:

      real(r8)
     $  xcappa      ! R / Cp for dry air
c-----------------------------------------------------------------------

      cappa = xcappa
      p0 = 1.e5

      return
      end

c##############################################################################

      subroutine t2pt( t, pmid, th )

c Compute potential temperature.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <pottem.com>
c-----------------------------------------------------------------------

c Input arguments

      real(r8)
     $  t(plond,plev)     ! temperature (K)
     $, pmid(plond,plev)  ! midpoint pressures (Pa)

c Output arguments

      real(r8)
     $  th(plond,plev)    ! potential temperature (K)

c Local variables

      integer i, k
c----------------------------------------------------------------------

      do k = 1, plev
         do i = 1, plon
            th(i,k) = t(i,k)*(p0/pmid(i,k))**cappa
         end do
      end do

      return
      end

c##############################################################################

      subroutine pt2t( th, pmid, t )

c Convert potential temperature back to temperature.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <pottem.com>
c-----------------------------------------------------------------------

c Input arguments

      real(r8)
     $  th(plond,plev)    ! potential temperature (K)
     $, pmid(plond,plev)  ! midpoint pressures (Pa)

c Output arguments

      real(r8)
     $  t(plond,plev)    ! temperature (K)

c Local variables

      integer i, k
c----------------------------------------------------------------------

      do k = 1, plev
         do i = 1, plon
            t(i,k) = th(i,k)*(pmid(i,k)/p0)**cappa
         end do
      end do

      return
      end

