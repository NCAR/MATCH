c $Id$
      subroutine inimass( gravit, xgw, xetamid )

c Initialize /mass/.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:

      real(r8)
     $  gravit         ! Gravitational acceleration (m/s^2)
     $, xgw(plat)      ! Gauss weights
     $, xetamid(plev)  ! hybrid coord at layer midpoints

c Local variables:

      integer i, k
c-----------------------------------------------------------------------

      rgrav = 1./gravit
      cwava = 1./(plon*gravit)

      do i = 1, plat
         gw(i) = xgw(i)
      end do

      do k = 1, plev
         etamid(k) = xetamid(k)
      end do

      return
      end

c##############################################################################

      subroutine qmassa( lat, ncnst, q3, pdel, hw1lat )

c Calculate area average mass of tracer fields.
c Input concentration units are wet mixing ratio (kg tracer / kg moist air).
c Output masses are in kg/m^2.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  lat                  ! latitude index
     $, ncnst                ! number of fields

      real(r8)
     $  q3(plond,plev,ncnst) ! advected species (wet mixing ratio)
     $, pdel(plond,plev)     ! pressure difference across model layer

c Output argument:

      real(r8)
     $  hw1lat(ncnst)        ! latitude contribution to mass of each specie

c Local variables:

      integer i,k,m
c-----------------------------------------------------------------------

      do m = 1, ncnst
         hw1lat(m) = 0.0
         do k = 1, plev
            do i = 1, plon
               hw1lat(m) = hw1lat(m) + q3(i,k,m)*pdel(i,k)
            end do
         end do
         hw1lat(m) = cwava * gw(lat) * hw1lat(m) * 0.5
      end do

      return
      end

c##############################################################################

      subroutine qmassd( lat, ncnst, q1, q2, pdel, hwn )

c Compute latitude contribution to integral of q2*etamid*|q2 - q1|
c Input concentration units are wet mixing ratio (kg tracer / kg moist air).
c Output units are in kg/m^2.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  lat                  ! latitude index
     $, ncnst                ! number of fields

      real(r8)
     $  q1(plond,plev,ncnst) ! tracer before slt
     $, q2(plond,plev,ncnst) ! tracer after slt
     $, pdel(plond,plev)     ! pressure difference across model layers

c Output argument:

      real(r8)
     $  hwn(ncnst) ! latitude contribution to integral

c Local variables:

      integer i,k,m
      real(r8) hwava
c-----------------------------------------------------------------------

      do m = 1, ncnst
         hwava = 0.0
         do k = 1, plev
            do i = 1, plon
               hwava = hwava + q2(i,k,m)*etamid(k)
     $                 * abs( q1(i,k,m) - q2(i,k,m) ) * pdel(i,k)
            end do
         end do
         hwn(m) = cwava * gw(lat) * hwava * 0.5
      end do

      return
      end

c##############################################################################

      subroutine gatmass( lat, ps, tmassj )

c Compute the contribution to the global average of the total atmospheric
c mass from a latitude slice.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  lat             ! latitude index

      real(r8)
     $  ps(plond)       ! pressure in pascals

c Output arguments:

      real(r8)
     $  tmassj          ! total mass in kg/m^2

c Local variables:

      integer i
      real(r8) a
c-----------------------------------------------------------------------

      tmassj = 0.0

      do i = 1, plon
         tmassj = tmassj + ps(i)
      end do

      tmassj = tmassj * cwava * gw(lat) * .5

      return
      end

c##############################################################################

      subroutine gamdry( lat, ncnst, as, q, pdel, massj )

c Calculate global average mass of tracer fields.
c Input concentration units are dry mixing ratio (kg tracer / kg dry air).
c Output masses are in kg/m^2.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  lat                  ! latitude index
     $, ncnst                ! number of fields

      real(r8)
     $  as(plond,plev,ncnst) ! advected species (dry mixing ratio)
     $, q(plond,plev)        ! specific humidity
     $, pdel(plond,plev)     ! pressure difference across model layer

c Output argument:

      real(r8)
     $  massj(ncnst)         ! latitude contribution to mass of each specie

c Local variables:

      integer i,k,m
c-----------------------------------------------------------------------

      do m = 1, ncnst
         massj(m) = 0.0
         do k = 1, plev
            do i = 1, plon
               massj(m) = massj(m) + as(i,k,m)*(1.-q(i,k))*pdel(i,k)
            end do
         end do
         massj(m) = cwava * gw(lat) * massj(m) * 0.5
      end do

      return
      end

c##############################################################################

      subroutine gamdry2( lat, ncnst, as, pdel, massj )

c Calculate global average mass of tracer fields.
c Input concentration units are dry mixing ratio (kg tracer / kg dry air).
c Output masses are in kg/m^2.

c This routine is the same as gamdry except that the pdel argument assumes 
c a dry atmosphere.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  lat                  ! latitude index
     $, ncnst                ! number of fields

      real(r8)
     $  as(plond,plev,ncnst) ! advected species (dry mixing ratio)
     $, pdel(plond,plev)     ! pressure difference across model layer in
c                            ! a dry atmosphere.

c Output argument:

      real(r8)
     $  massj(ncnst)         ! latitude contribution to mass of each specie

c Local variables:

      integer i,k,m
c-----------------------------------------------------------------------

      do m = 1, ncnst
         massj(m) = 0.0
         do k = 1, plev
            do i = 1, plon
               massj(m) = massj(m) + as(i,k,m)*pdel(i,k)
            end do
         end do
         massj(m) = cwava * gw(lat) * massj(m) * 0.5
      end do

      return
      end

c##############################################################################

      subroutine ciwet( mr, pdel, ci )

c Calculate column integral of tracer fields for a latitude slice.
c Input concentration units are wet mixing ratio (kg tracer / kg moist air).
c Output masses are in kg/m^2.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  mr(plond,plev)       ! constituent mixing ratio (wet)
     $, pdel(plond,plev)     ! pressure difference across model layer

c Output argument:
      real(r8)
     $  ci(plond)            ! column integral, kg/m^2

c Local variables:

      integer i,k,m
c-----------------------------------------------------------------------

      do i = 1, plon
         ci(i) = 0.0
      end do
      do k = 1, plev
         do i = 1, plon
            ci(i) = ci(i) + mr(i,k)*pdel(i,k)*rgrav
         end do
      end do

      return
      end

c##############################################################################

      subroutine cidry( mr, pdel, q, ci )

c Calculate column integral of tracer fields for a latitude slice.
c Input concentration units are dry mixing ratio (kg tracer / kg dry air).
c Output masses are in kg/m^2.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  mr(plond,plev)       ! constituent mixing ratio (dry)
     $, pdel(plond,plev)     ! pressure difference across model layer
     $, q(plond,plev)        ! specific humidity (kg/kg)

c Output argument:
      real(r8)
     $  ci(plond)            ! column integral, kg/m^2

c Local variables:

      integer i,k,m
c-----------------------------------------------------------------------

      do i = 1, plon
         ci(i) = 0.0
      end do
      do k = 1, plev
         do i = 1, plon
            ci(i) = ci(i) + mr(i,k)*pdel(i,k)*(1.-q(i,k))*rgrav
         end do
      end do



      return
      end

c##############################################################################

      subroutine fixint( lat, ncnst, as1, as2, q, pdel, hwn )

c Compute latitude contribution to integral of as2*etamid*|as2 - as1|
c Tracer concentration units are dry mixing ratio (kg tracer / kg dry air).


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  lat                  ! latitude index
     $, ncnst                ! number of fields

      real(r8)
     $  as1(plond,plev,ncnst) ! tracer before slt
     $, as2(plond,plev,ncnst) ! tracer after slt
     $, q(plond,plev)         ! specific humidity
     $, pdel(plond,plev)      ! pressure difference across model layers

c Output argument:

      real(r8)
     $  hwn(ncnst) ! latitude contribution to integral

c Local variables:

      integer i,k,m
      real(r8) hwava
c-----------------------------------------------------------------------

      do m = 1, ncnst
         hwava = 0.0
         do k = 1, plev
            do i = 1, plon
               hwava = hwava + as2(i,k,m)*etamid(k)*(1.-q(i,k))
     $                 * abs( as1(i,k,m) - as2(i,k,m) ) * pdel(i,k)
            end do
         end do
         hwn(m) = cwava * gw(lat) * hwava * 0.5
      end do

      return
      end

c##############################################################################

      subroutine fixint2( lat, ncnst, as1, as2, pdel, hwn )

c Compute latitude contribution to integral of as2*etamid*|as2 - as1|
c Tracer concentration units are dry mixing ratio (kg tracer / kg dry air).
c Use dry pdel.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  lat                  ! latitude index
     $, ncnst                ! number of fields

      real(r8)
     $  as1(plond,plev,ncnst) ! tracer before slt
     $, as2(plond,plev,ncnst) ! tracer after slt
     $, pdel(plond,plev)      ! dry pressure difference across model layers

c Output argument:

      real(r8)
     $  hwn(ncnst) ! latitude contribution to integral

c Local variables:

      integer i,k,m
      real(r8) hwava
c-----------------------------------------------------------------------

      do m = 1, ncnst
         hwava = 0.0
         do k = 1, plev
            do i = 1, plon
               hwava = hwava + as2(i,k,m)*etamid(k)
     $                 * abs( as1(i,k,m) - as2(i,k,m) ) * pdel(i,k)
            end do
         end do
         hwn(m) = cwava * gw(lat) * hwava * 0.5
      end do

      return
      end

c##############################################################################

      subroutine gamk( as, q, ps, trmk, mixk )

c Calculate global average mass and global average mixing ratio
c of tracer fields at each level.
c Tracer concentration units are dry mixing ratio (kg tracer / kg dry air).


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  as(plond,plev,pcnst,plat) ! advected species (dry mixing ratio)
     $, q(plond,plev,plat)        ! specific humidity
     $, ps(plond,plat)            ! surface pressure

c Output arguments:
      real(r8)
     $  trmk(plev,pcnst)          ! global average mass (kg/m2)
     $, mixk(plev,pcnst)          ! global average m.r. (kg/kg)

c Local variables:
      integer
     $  i, k, m, lat
      real(r8)
     $  trmjk
     $, mixjk
     $, pint(plond,plevp) ! pressure at layer interfaces
     $, pmid(plond,plev)  ! pressure at layer midpoints
     $, pdel(plond,plev)  ! pressure difference across layers
     $, drym(plond,plev)  ! dry mass
c-----------------------------------------------------------------------

      do m = 1, pcnst
         do k = 1, plev
            trmk(k,m) = 0.0
            mixk(k,m) = 0.0
         end do
      end do

      do lat = 1, plat
         call plevs0( ps(1,lat), pint, pmid, pdel )
         do k = 1, plev
            do i = 1, plon
               drym(i,k) = (1.-q(i,k,lat))*pdel(i,k)
            end do
         end do

         do m = 1, pcnst
            do k = 1, plev
               trmjk = 0.0
               mixjk = 0.0
               do i = 1, plon
                  trmjk = trmjk + as(i,k,m,lat)*drym(i,k)
                  mixjk = mixjk + as(i,k,m,lat)
               end do
               trmk(k,m) = trmk(k,m) + trmjk*gw(lat)
               mixk(k,m) = mixk(k,m) + mixjk*gw(lat)
            end do
         end do
      end do

      do m = 1, pcnst
         do k = 1, plev
            trmk(k,m) = trmk(k,m)*cwava*0.5
            mixk(k,m) = mixk(k,m)*0.5/plon
         end do
      end do

      return
      end

c##############################################################################

      subroutine gamk2( as, ps, trmk, mixk )

c Calculate global average mass and global average mixing ratio
c of tracer fields at each level.
c Tracer concentration units are dry mixing ratio (kg tracer / kg dry air).

c This is gamk() but uses a dry ps.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  as(plond,plev,pcnst,plat) ! advected species (dry mixing ratio)
     $, ps(plond,plat)            ! dry surface pressure

c Output arguments:
      real(r8)
     $  trmk(plev,pcnst)          ! global average mass (kg/m2)
     $, mixk(plev,pcnst)          ! global average m.r. (kg/kg)

c Local variables:
      integer
     $  i, k, m, lat
      real(r8)
     $  trmjk
     $, mixjk
     $, pint(plond,plevp) ! pressure at layer interfaces
     $, pmid(plond,plev)  ! pressure at layer midpoints
     $, drym(plond,plev)  ! dry mass (in pressure units)
c-----------------------------------------------------------------------

      do m = 1, pcnst
         do k = 1, plev
            trmk(k,m) = 0.0
            mixk(k,m) = 0.0
         end do
      end do

      do lat = 1, plat
         call plevs0( ps(1,lat), pint, pmid, drym )

         do m = 1, pcnst
            do k = 1, plev
               trmjk = 0.0
               mixjk = 0.0
               do i = 1, plon
                  trmjk = trmjk + as(i,k,m,lat)*drym(i,k)
                  mixjk = mixjk + as(i,k,m,lat)
               end do
               trmk(k,m) = trmk(k,m) + trmjk*gw(lat)
               mixk(k,m) = mixk(k,m) + mixjk*gw(lat)
            end do
         end do
      end do

      do m = 1, pcnst
         do k = 1, plev
            trmk(k,m) = trmk(k,m)*cwava*0.5
            mixk(k,m) = mixk(k,m)*0.5/plon
         end do
      end do

      return
      end

c##############################################################################

      subroutine gamk3( as, drym, trmk, mixk )

c Calculate global average mass and global average mixing ratio
c of tracer fields at each level.
c Tracer concentration units are dry mixing ratio (kg tracer / kg dry air).

c This is gamk() but take dry pressure diffs as input


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <mass.com>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  as(plond,plev,pcnst,plat) ! advected species (dry mixing ratio)
     $, drym(plond,plev,plat)     ! dry delta pressure across layers

c Output arguments:
      real(r8)
     $  trmk(plev,pcnst)          ! global average mass (kg/m2)
     $, mixk(plev,pcnst)          ! global average m.r. (kg/kg)

c Local variables:
      integer
     $  i, j, k, m
      real(r8)
     $  trmjk
     $, mixjk
c-----------------------------------------------------------------------

      do m = 1, pcnst
         do k = 1, plev
            trmk(k,m) = 0.0
            mixk(k,m) = 0.0
         end do
      end do

      do j = 1, plat
         do m = 1, pcnst
            do k = 1, plev
               trmjk = 0.0
               mixjk = 0.0
               do i = 1, plon
                  trmjk = trmjk + as(i,k,m,j)*drym(i,k,j)
                  mixjk = mixjk + as(i,k,m,j)
               end do
               trmk(k,m) = trmk(k,m) + trmjk*gw(j)
               mixk(k,m) = mixk(k,m) + mixjk*gw(j)
            end do
         end do
      end do

      do m = 1, pcnst
         do k = 1, plev
            trmk(k,m) = trmk(k,m)*cwava*0.5
            mixk(k,m) = mixk(k,m)*0.5/plon
         end do
      end do

      return
      end

c##############################################################################

