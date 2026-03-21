c $Id$
      subroutine setqth( ps, q, t, x0, x1 )

c Initialize the array of misc advected variables with q and theta

      use precision

      implicit none
c---------------------------------------------------------------------
#include <pmgrid.h>
c---------------------------------------------------------------------
#include <control.com>
c---------------------------------------------------------------------

c Input arguments:

      real(r8)
     $  ps(plond,plat)           ! surface pressure
     $, q(plond,plev,plat)       ! specific humidity
     $, t(plond,plev,plat)       ! temperature
     $, x1(plond,plev,padv,plat) ! misc advected variables at end of
c                                ! previous timestep

c Output arguments:

      real(r8)
     $  x0(plond,plev,padv,plat) ! misc advected variables at beginning of
c                                ! current timestep

c Local variables:

      integer
     $  i, j, k
      real(r8)
     $  pint(plond,plevp)      ! pressure at layer interfaces
     $, pmid(plond,plev)       ! pressure at layer midpoints
     $, pdel(plond,plev)       ! pressure difference across layers
c---------------------------------------------------------------------

      do j = 1, plat
c        Update specific humidity.
         if ( qrelax .eq. 1. .or. (begstep .and. .not. readQIC) ) then
            do k = 1, plev
               do i = 1, plon
                  x0(i,k,1,j) = q(i,k,j)
               end do
            end do
         else
            do k = 1, plev
               do i = 1, plon
                  x0(i,k,1,j) = x1(i,k,1,j) + 
     $                          qrelax*( q(i,k,j) - x1(i,k,1,j) )
               end do
            end do
         end if
c        Update potential temperature.
         call plevs0( ps(1,j), pint, pmid, pdel )
         call t2pt( t(1,1,j), pmid, x0(1,1,2,j) )
      end do

      return
      end
