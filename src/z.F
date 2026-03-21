c $Id$
      subroutine iniz(  xrair, xgravit )

c Initialize /z0/.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <z.com>
c-----------------------------------------------------------------------

c Input arguments:

      real(r8)
     $  xrair       ! Gas constant for dry air
     $, xgravit     ! acceleration due to gravity (m/s^2)
c-----------------------------------------------------------------------

      rair = xrair
      gravit = xgravit

      return
      end

c##############################################################################

      subroutine zrel( pint, pmln, tv, zm, zi )

c Compute the geopotential height above the surface at layer midpoints and
c interfaces.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments

      real(r8)
     $  pint(plond,plevp)    ! interface pressures
     $, pmln(plond,plev)     ! Log midpoint pressures
     $, tv(plond,plev)       ! Virtual temperature

c Output arguments

      real(r8)
     $  zm(plond,plev)       ! Height above surface at midpoints
     $, zi(plond,plevp)      ! Height above surface at interfaces

c Local variables.

      integer i, k
      real(r8) piln(plond,plevp)
c-----------------------------------------------------------------------

      do k = 2, plevp
         do i = 1, plon
            piln(i,k) = log( pint(i,k) )
         end do
      end do
      do i = 1, plon
         if ( pint(i,1) .ne. 0.0 ) then
            piln(i,1) = log( pint(i,1) )
         else
            piln(i,1) = log( pint(i,2)/10. )
         end if
      end do

      call zmid( piln, pmln, tv, zm )
      call zint( piln, pmln, tv, zm, zi )
      end

c##############################################################################

      subroutine zint( piln, pmln, tv, zm, zi )
C-----------------------------------------------------------------------
C
C Compute the geopotential height at the interface points from *HEIGHT
C ABOVE THE SURFACE* at midlayer midpoints, using the supplied virtual
C temperatures and pressures.
C
C---------------------------Code history--------------------------------
C
C Original version:  D. Williamson, J. Hack
C Standardized:      L. Buja, June 1992
C Reviewed:          J. Hack, B. Boville, August 1992
C
C-----------------------------------------------------------------------

      use precision

      implicit none
C------------------------------Commons----------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <z.com>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) piln(plond,plevp),   ! Log interface pressures
     $     pmln(plond,plev),    ! Log midpoint pressures
     $     tv(plond,plev),      ! Virtual temperature
     $     zm(plond,plev)       ! Height above surface at midpoints
C
C Output arguments
C
      real(r8) zi(plond,plevp)      ! Height above surface at interfaces
C
C---------------------------Local variables-----------------------------
C
      integer i,                ! Longitude index
     $        k                 ! Level index
      real(r8) rog                  ! Rair / gravit
C
C Add increment to midlayer height
C
      rog = rair/gravit
      do i=1,plon
         zi(i,1) = zm(i,1) + rog*(pmln(i,1) - piln(i,1))*tv(i,1)
      end do
C
      do k=2,plev
         do i=1,plon
            zi(i,k) = zm(i,k) + rog*(pmln(i,k) - piln(i,k))*0.5*
     $                (2.0*tv(i,k) - (tv(i,k) - tv(i,k-1))/
     $                 (pmln(i,k) - pmln(i,k-1))*
     $                 (pmln(i,k) - piln(i,k)))
         end do
      end do
C
C The surface height is zero by definition.
C
      do i=1,plon
         zi(i,plevp) = 0.0
      end do
C
      return
      end

c##############################################################################

      subroutine zmid(piln    ,pmln    ,tv      ,z       )
C-----------------------------------------------------------------------
C
C Compute the geopotential height *ABOVE THE SURFACE* at layer
C midpoints from the virtual temperatures and pressures.
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Boville, D. Williamson, January 1990
C Standardized:      L. Buja, June 1992
C Reviewed:          J. Hack, B. Boville, August 1992
C
C-----------------------------------------------------------------------

      use precision

      implicit none
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <z.com>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) piln(plond,plevp),    ! Log interface pressures
     $     pmln(plond,plev),     ! Log midpoint pressures
     $     tv(plond,plev)        ! Virtual temperature
C
C Output arguments
C
      real(r8) z(plond,plev)         ! Height above surface at midpoints
C
C---------------------------Local variables-----------------------------
C
      integer i,                 ! Longitude index
     $        k,                 ! Level index
     $        l                  ! Level index
      real(r8) rog                   ! Rair / gravit
C
C Diagonal term of hydrostatic equation
C
      rog = rair/gravit
      do k=1,plev-1
         do i=1,plon
            z(i,k) = rog*tv(i,k)*0.5*(pmln(i,k+1) - pmln(i,k))
         end do
      end do
      do i=1,plon
         z(i,plev) = rog*tv(i,plev)*(piln(i,plevp) - pmln(i,plev))
      end do
C
C Bottom level term of hydrostatic equation
C
      do  k=1,plev-1
         do i=1,plon
            z(i,k) = z(i,k) + rog*tv(i,plev)*(piln(i,plevp) -
     $                        0.5*(pmln(i,plev-1) + pmln(i,plev)))
         end do
      end do
C
C Interior terms of hydrostatic equation
C
      do k=1,plev-2
         do l=k+1,plev-1
            do i=1,plon
               z(i,k) = z(i,k) + rog*(tv(i,l))*
     $                           0.5*(pmln(i,l+1) - pmln(i,l-1))
            end do
         end do
      end do
C
      return
      end

c##############################################################################

      subroutine zr2a( zr, nlev, phis, za )

c Convert from relative to absolute geopotential height.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <z.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  nlev            ! number of levels

      real(r8)
     $  zr(plond,nlev)  ! geopotential height relative to the surface
     $, phis(plond)     ! surface geopotential

c Output argument:

      real(r8)
     $  za(plond,nlev)  ! absolute geopotential height

c Local variables:

      integer i, k
      real(r8) rgrav
c-----------------------------------------------------------------------
	
      rgrav = 1./gravit
      do k = 1, nlev
         do i = 1, plon
            za(i,k) = zr(i,k) + phis(i)*rgrav
         end do
      end do

      return
      end

c##############################################################################

      subroutine z2phi( z, nlev, phi )

c Convert from geopotential height to geopotential.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <z.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  nlev            ! number of levels

      real(r8)
     $  z(plond,nlev)   ! geopotential height

c Output argument:

      real(r8)
     $  phi(plond,nlev)      ! geopotential

c Local variables:

      integer i, k
c-----------------------------------------------------------------------
	
      do k = 1, nlev
         do i = 1, plon
            phi(i,k) = z(i,k)*gravit
         end do
      end do

      return
      end

c##############################################################################

      subroutine zr2phi( z, nlev, phis, phi )

c Convert from geopotential height ABOVE SURFACE to geopotential.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <z.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  nlev            ! number of levels

      real(r8)
     $  z(plond,nlev)   ! geopotential height
     $, phis(plond)     ! surface geopotential

c Output argument:

      real(r8)
     $  phi(plond,nlev)      ! geopotential

c Local variables:

      integer i, k
c-----------------------------------------------------------------------
	
      do k = 1, nlev
         do i = 1, plon
            phi(i,k) = z(i,k)*gravit + phis(i)
         end do
      end do

      return
      end
