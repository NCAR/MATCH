c $Id$
      subroutine inivtem( xzvir )

c Initialize /virtem0/.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <virtem.com>
c-----------------------------------------------------------------------

c Input arguments:

      real(r8)
     $  xzvir      ! rh2o/rair - 1.
c-----------------------------------------------------------------------

      zvir = xzvir

      return
      end

c##############################################################################

      subroutine virtem(t       ,q       ,tv      )
C-----------------------------------------------------------------------
C
C Compute the virtual temperature.
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Boville
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, J. Hack, August 1992
C
C-----------------------------------------------------------------------

      use precision

      implicit none
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <virtem.com>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) t(plond,plev),      ! temperature
     $     q(plond,plev)       ! specific humidity
C
C Output arguments
C
      real(r8) tv(plond,plev)      ! virtual temperature
C
C---------------------------Local storage-------------------------------
C
      integer
     $     i, k                ! longitude and level indexes
C
C-----------------------------------------------------------------------
C
      do k=1,plev
         do i=1,plon
            tv(i,k) = t(i,k)*(1.0 + zvir*q(i,k))
         end do
      end do
C
      return
      end
