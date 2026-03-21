c $Id$
      subroutine iniplevs( xhyai, xhyam, xhybi, xhybm )

c Initialize /plevs/.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <plevs.com>
c-----------------------------------------------------------------------

c Input arguments:

      real(r8)
     $  xhyai(plevp)   ! p0 component of hybrid coordinate - interfaces
     $, xhyam(plev)    ! p0 component of hybrid coordinate - midpoints
     $, xhybi(plevp)   ! ps component of hybrid coordinate - interfaces
     $, xhybm(plev)    ! ps component of hybrid coordinate - midpoints

c Local variables:

      integer k                 ! Level index
c-----------------------------------------------------------------------

      ps0 = 1.0e5               ! Base state surface pressure (pascals)

      do k = 1, plev
         hyai(k) = xhyai(k)
         hybi(k) = xhybi(k)
         hyam(k) = xhyam(k)
         hybm(k) = xhybm(k)
      end do
      hyai(plevp) = xhyai(plevp)
      hybi(plevp) = xhybi(plevp)

      do k = 1, plev
         dhya(k) = hyai(k+1) - hyai(k)
         dhyb(k) = hybi(k+1) - hybi(k)
      end do

      ptop = hyai(1)*ps0

      return
      end

c==============================================================================
c==============================================================================

      subroutine plevs0( ps, pint, pmid, pdel )
C-----------------------------------------------------------------------
C
C Define the pressures of the interfaces and midpoints from the
C coordinate definitions and the surface pressure.
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Boville
C Standardized:      J. Rosinski, June 1992
C Reviewed:          B. Boville, P. Rasch, August 1992
C
C-----------------------------------------------------------------------

      use precision

      implicit none
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <plevs.com>
C-----------------------------------------------------------------------
C
C Input arguments
C
      real(r8) ps(plond)          ! Surface pressure (pascals)
C
C Output arguments
C
      real(r8) pint(plond,plevp), ! Pressure at model interfaces
     $     pmid(plond,plev),  ! Pressure at model levels
     $     pdel(plond,plev)   ! Layer thickness (pint(k+1) - pint(k))
C
C---------------------------Local workspace-----------------------------
C
      integer
     $     i,k                ! Longitude, level indices
C-----------------------------------------------------------------------
C
C Set interface pressures
C
      do k=1,plev+1
         do i=1,plon
            pint(i,k) = hyai(k)*ps0 + hybi(k)*ps(i)
         end do
      end do
C
C Set midpoint pressures and layer thicknesses
C
      do k=1,plev
         do i=1,plon
            pmid(i,k) = hyam(k)*ps0 + hybm(k)*ps(i)
            pdel(i,k) = pint(i,k+1) - pint(i,k)
         end do
      end do
C
      return
      end

c==============================================================================
c==============================================================================

      subroutine pdelab( ps, pdela, pdelb )

c Contribution to pressure differences across layers from the A and B
c components of the hybrid coordinate.


      use precision

      implicit none
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <plevs.com>
C-----------------------------------------------------------------------

c Input arguments

      real(r8) ps(plond)          ! Surface pressure (pascals)

c Output arguments

      real(r8)
     $  pdela(plond,plev)   ! Pressure diff due to A component
     $, pdelb(plond,plev)   ! Pressure diff due to B component

c Local varibles

      integer
     $     i,k                ! Longitude, level indices
c-----------------------------------------------------------------------

      do k = 1, plev
         do i = 1, plon
            pdela(i,k) = dhya(k)*ps0
            pdelb(i,k) = dhyb(k)*ps(i)
         end do
      end do

      return
      end

c##############################################################################

      function p_top()

c Return pressure at top of model.


      use precision

      implicit none

      real(r8) p_top

c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <plevs.com>
c-----------------------------------------------------------------------

      p_top = ptop
      return
      end

c##############################################################################

      subroutine getPressure( ps, p , dp )

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments
      real(r8)
     $  ps(plond,plat)      ! surface pressure in Pa

c Output arguments
      real(r8)
     $  p(plev,plat,plon)   ! pressure in Pa
     $, dp(plev,plat,plon)  ! pressure change across layer in Pa

c Local
      integer i, j, k
      real(r8)
     $  pint(plond,plevp)      ! pressure at layer interfaces
     $, pmid(plond,plev)       ! pressure at layer midpoints
     $, pdel(plond,plev)       ! pressure difference across layers
c-----------------------------------------------------------------------

      do j = 1, plat
         call plevs0( ps(1,j), pint, pmid, pdel )
         do i = 1, plon
            do k = 1, plev
               p(k,j,i) = pmid(i,k)
               dp(k,j,i) = pdel(i,k)
            end do
         end do
      end do

      return
      end

c#######################################################################
