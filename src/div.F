c $Id$
      subroutine inidiv( pi, phi, dmu, dlam, rearth )

c Initialize /div0/.

      use precision

      implicit none
c---------------------------------------------------------------------
#include <pmgrid.h>
c---------------------------------------------------------------------
#include <div.com>
c---------------------------------------------------------------------

c Input arguments:

      real(r8)
     $  pi         ! 3.14159...
     $, phi(plat)  ! Gaussian latitudes (radians)
     $, dmu(plat)  ! delta mu for cells containing Gaussian latitudes
     $, dlam       ! delta lambda (radians)
     $, rearth     ! radius of earth (m)

c Local variables:

      integer
     $  j          ! index
c---------------------------------------------------------------------

c     Cell wall positions in mu coords.

      wall(1) = -1.
      do j = 1,plat
         wall(j+1) = dmu(j) + wall(j)
         if (wall(j+1) .gt. 1.0) wall(j+1) = 1.0
      end do
      do j = 1,plat+1
         wall(j) = asin( wall(j) )
      end do

c     Calculate weights for interpolation to cell walls.

      do j = 1, plat

         if ( j .eq. 1 ) then
            w1s(j) = ( wall(j) + pi/2. ) / ( phi(j) + pi/2. )
         else
            w1s(j) = ( wall(j) - phi(j-1) ) / ( phi(j) - phi(j-1) )
         end if
         w2s(j) = 1 - w1s(j)

         if ( j .eq. plat ) then
            w1n(j) = ( wall(j+1) - phi(j) ) / ( pi/2. - phi(j) )
         else
            w1n(j) = ( wall(j+1) - phi(j) ) / ( phi(j+1) - phi(j) )
         end if
         w2n(j) = 1 - w1n(j)

      end do

c     Reciprocals of coordinate differences for derivative calcs.

      do j = 1, plat
         ddlam(j) = 1./( 2.*dlam*rearth*cos( phi(j) ) )
         ddmu(j) = 1./( rearth*dmu(j) )
      end do

c     Cosine of cell wall latitudes.

      do j = 1, plat+1
         cw(j) = cos( wall(j) )
      end do

      return
      end

c##############################################################################

      subroutine dcal( lat, u, v, d )

c Finite difference estimate of the divergence.  It is assumed that the 
c velocity fields have been extended to the slt grid.

      use precision

      implicit none
c---------------------------------------------------------------------
#include <pmgrid.h>
c---------------------------------------------------------------------
#include <div.com>
c---------------------------------------------------------------------

c Input arguments:

      integer
     $  lat                  ! current latitude

      real(r8)
     $  u(plond,plevd,3)     ! zonal velocity
     $, v(plond,plevd,3)     ! meridional velocity

c Output argument:

      real(r8)
     $  d(plon,plev)         ! estimate of divergence field

c Local variables:

      integer
     $  i, ix, k             ! indices

      real(r8)
     $  vn                   ! estimate of v at north cell wall
     $, vs                   ! estimate of v at south cell wall
c---------------------------------------------------------------------

      do k = 1,plev
         do i = 1, plon
            ix = i + nxpt

c           Interpolate velocities to north and south cell walls.
            vn = (w1n(lat)*v(ix,k,3) + w2n(lat)*v(ix,k,2))*cw(lat+1)
            vs = (w1s(lat)*v(ix,k,2) + w2s(lat)*v(ix,k,1))*cw(lat)

c           Divergence.
            d(i,k) = ( u(ix+1,k,2) - u(ix-1,k,2) )*ddlam(lat)
     $               + ( vn - vs )*ddmu(lat)

         end do
      end do

      return
      end
