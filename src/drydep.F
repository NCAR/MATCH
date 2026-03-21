c $Id$
#include <params.h>

      subroutine inidrydep( xrair, xgravit, xphi )

c Initialize dry deposition parameterization.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <drydep.h>
c-----------------------------------------------------------------------
#include <netcdf.inc>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  xrair                ! Gas constant for dry air
     $, xgravit              ! Gravitational acceleration
     $, xphi(plat)           ! grid latitudes (radians)

c Local variables:
      integer
     $  i, j
     $, ncid, vid
c-----------------------------------------------------------------------

      rair = xrair
      gravit = xgravit
      do j = 1, plat
         phi(j) = xphi(j)
      end do

      return
      end

c##############################################################################

      subroutine setdvel( ioro, vgl, vgo, vgsi, vg )

c Set the deposition velocity depending on whether we are over
c land, ocean, and snow/ice


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments:

      integer ioro(plon)        ! flag for land type

      real(r8) vgl                  ! dry deposition velocity in m/s (land)
      real(r8) vgo                  ! dry deposition velocity in m/s (ocean)
      real(r8) vgsi                 ! dry deposition velocity in m/s (snow/ice)

c Output arguments:
      real(r8)
     $  vg(plon) ! dry deposition velocity in m/s

c Local variables:

      integer i
      real(r8) a


      do i = 1, plon
         if (ioro(i).eq.0) then
            vg(i) = vgo
         else if (ioro(i).eq.1) then
            vg(i) = vgl
         else
            vg(i) = vgsi
         endif
      end do

      return
      end

c##############################################################################

      subroutine ddflux( vg, q, p, tv, flux )

c Compute surface flux due to dry deposition processes.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <drydep.h>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  vg(plon)  ! dry deposition velocity in m/s
     $, q(plon)   ! tracer conc. in surface layer (kg tracer/kg moist air)
     $, p(plon)   ! midpoint pressure in surface layer (Pa)
     $, tv(plon)  ! midpoint virtual temperature in surface layer (K)

c Output arguments:

      real(r8)
     $  flux(plon) ! flux due to dry deposition in kg/m^s/sec

c Local variables:

      integer i

      do i = 1, plon
         flux(i) = -vg(i) * q(i) * p(i) /(tv(i) * rair)
      end do

      return
      end

c##############################################################################
