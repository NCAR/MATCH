c $Id$
#include <params.h>
      subroutine grdxy( pi, dlam, lam, phi, latwts )

c Define the "extended" grid used in the semi-Lagrangian transport
c scheme.  The longitudes are equally spaced and the latitudes are
c Gaussian.  The global grid is extended to include "wraparound" points
c on all sides.

      use precision

      implicit none
#if ( defined RS6000 )
      implicit automatic (a-z)
#endif
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <grid.com>
c-----------------------------------------------------------------------

c Input argument:

      real(r8)
     $  pi                  ! 3.14...

c Input/Output argument:
      real(r8)
     $  phi(platd)          ! lat. coords. in non-extended grid on input (radians),
c                             extended on output.  If latitudes are Gaussian then
c                             ensure that values are double precision on output.
     $, latwts(plat)        ! latitude weights for area integrals.  If lats are Gaussian
                            ! then ensure that values are double precision on output.

c Output arguments:

      real(r8)
     $  dlam                ! longitudinal increment
     $, lam(plond)          ! long. coords. in extended grid

c Local variables:

      integer
     $  istart              ! index for first model long.
     $, jstart              ! index for first model lat.
     $, jstop               ! index for last  model lat.
      parameter(istart = nxpt+1,
     $          jstart = nxpt+jintmx+1,
     $          jstop  = jstart-1+plat)

      integer
     $  i, j                ! indices

      real(r8)
     $  lam0             ! lamda = 0
     $, wrk(platd)       ! work space
     $, gl(plat)         ! Gaussian latitudes
     $, gw(plat)         ! Gaussian weights
     $, tmp, rms
     $, rad2deg
     $, dlamDeg
c-----------------------------------------------------------------------

      rad2deg = 180./pi

c     Compute Gauss latitudes and weights.  On return; phi contains the
c     sine of the latitudes starting closest to the north pole and going
c     toward the south; w contains the corresponding Gauss weights.
      call gauaw( gl, gw, plat )
c     Reorder and compute latitude values.  Do not need to reorder weights
c     since they are symmetric.
      do j = 1,plat/2
         tmp = asin( gl(j) )
         gl(j) = asin( gl(plat-j+1) )
         gl(plat-j+1) = tmp
      end do

c     Check to see if the latitudes from the dynamics data are Gaussian.
c     If they are replace them with the ones computed by gauaw to ensure
c     that we have full double precision values.
      rms = 0.0
      do j = 1, plat
         rms = rms + ( gl(j) - phi(jstart-1+j) )**2
      end do
      rms = sqrt( rms/plat )
      if ( rms .lt. 1.e-4 ) then      ! latitudes are Gaussian
         do j = 1, plat
            phi(jstart-1+j) = gl(j)
            latwts(j) = gw(j)
         end do
         write(*,*)'grdxy: Gaussian latitudes and weights have been '
     $,            'calculated internally.'
      end if

      do j = 1, plat
         lat(j) = phi(jstart-1+j)
         latDeg(j) = lat(j) * rad2deg
      end do

      lam0 = 0.0

c     Interval length in equally spaced longitude grid.
      dlam = 2.*pi/float(plon)
      dlamDeg = 360./float(plon)

c     Longitude values on extended grid.
      do i = 1,plond
         lam(i) = float(i-istart)*dlam + lam0
      end do
      do i = 1, plon
         lon(i) = lam(istart+i-1)
         lonDeg(i) = float(i-1) * dlamDeg
      end do

c     Extend latitude grid to north and south poles.
      phi(jstart-1) = -pi/2.0
      phi(jstop +1) =  pi/2.0

c     Extend Gauss latitudes below south pole so that the spacing about
c     the pole is symmetric, and phi is decreasing, i.e., phi < -pi/2
      if( jstart .gt. 2 )then
         do j = 1,jstart-2
            phi(j) = -pi - phi(2*jstart-2-j)
         end do
      end if

c     Analogously for Northern Hemisphere
      if( platd .gt. jstop+1 )then
         do j = jstop+2,platd
            phi(j) = pi - phi(2*jstop+2-j)
         end do
      end if

      return
      end

c##############################################################################

      function clat( index )

      use precision

c Return the latitude at the requested index.

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <grid.com>
c-----------------------------------------------------------------------

c Input argument.
      integer
     $  index   ! latitude index

c Return value is latitude at requested index (radians)
      real(r8) clat

c-----------------------------------------------------------------------
      clat = lat(index)

      return
      end

c##############################################################################

      subroutine getLatLon( units, grid_lat, grid_lon )

c Return the latitude and longitude arrays.
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <grid.com>
c-----------------------------------------------------------------------

c Input
      character*(*) units   ! 'radians' or 'degrees'

c Output
      real(r8)
     $  grid_lat(plat)  ! grid latitudes in radians
     $, grid_lon(plon)  ! grid longitudes in radians

c Local
      integer i
c-----------------------------------------------------------------------
      if ( units .eq. 'radians' ) then
         do i = 1, plat
            grid_lat(i) = lat(i)
         end do
         do i = 1, plon
            grid_lon(i) = lon(i)
         end do
      else if ( units .eq. 'degrees' ) then
         do i = 1, plat
            grid_lat(i) = latDeg(i)
         end do
         do i = 1, plon
            grid_lon(i) = lonDeg(i)
         end do
      else
         write(6,*)'getLatLon: unrecognised units: ',units
         stop
      end if

      return
      end

c##############################################################################

      subroutine hycoef( hyai, hyam, hybi, hybm, hybd, nprlev
     $,                  etamid, etaint )

c Initialize hybrid coordinates.
c Hybrid level definitions: p = a*p0 + b*ps
c        interfaces   p(k) = hyai(k)*p0 + hybi(k)*ps
c        midpoints    p(k) = hyam(k)*p0 + hybm(k)*ps

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments:

      real(r8)
     $  hyam(plev)    ! p0 component of hybrid coordinate - midpoints
     $, hybm(plev)    ! ps component of hybrid coordinate - midpoints

c Input/Output arguments:

      real(r8)
     $  hyai(plevp)   ! p0 component of hybrid coordinate - interfaces
     $, hybi(plevp)   ! ps component of hybrid coordinate - interfaces
c                       hyai(1) and hybi(1) are set to hyai(2)/10. and
c                       hybi(2)/10. resp. if they are BOTH zero on input.

c Output arguments:

      real(r8)
     $  hybd(plev)    ! difference of the hybi's
     $, etamid(plev)  ! hyam + hybm
     $, etaint(plevp) ! hyai + hybi

      integer
     $  nprlev        ! index of lowest pure pressure interface

c Local variables:

      integer
     $  k             ! indices
c-----------------------------------------------------------------------

      do k = 1, plev
         etamid(k) = hyam(k) + hybm(k)
      end do
      do k = 1, plevp
         etaint(k) = hyai(k) + hybi(k)
      end do

c     Because of the scheme used in the SLT for estimating vertical departure
c     points, the top of the model can't be at eta=0.
c      if ( etaint(1) .eq. 0. ) then
c         hyai(1) = hyai(2) / 10.
c         hybi(1) = hybi(2) / 10.
c         etaint(1) = hyai(1) + hybi(1)
c      end if

C     Set delta sigma part of layer thickness
      do k=1,plev
         hybd(k) = hybi(k+1) - hybi(k)
      end do

c     Set nprlev to the interface above, the first time a 
c     nonzero surface pressure contribution is found. "nprlev" 
c     identifies the lowest pure pressure interface.
      nprlev = 0
      do k = 1, plev
         if ( nprlev .eq. 0 .and. hybi(k) .ne. 0.0 ) nprlev = k - 1
      end do
C     Set nprlev if no nonzero b's have been found. All interfaces are 
C     pure pressure. A pure pressure model requires other changes as well. 
      if (nprlev.eq.0) nprlev = plev + 2

      return
      end
