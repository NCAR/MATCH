c $Id$
      subroutine inidry( ps, q, psdry )

c Initialize /drymass/.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <drymass.com>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  ps(plond,plat)         ! surface pressure (pascals)
     $, q(plond,plev,plat)     ! specific humidity

c Output arguments:
      real(r8)
     $  psdry(plond,plat)      ! dry surface pressure (pascals)

c Local variables:
      integer i, j, k
      real(r8)
     $  tmass, tmassj(plat)    ! total atm mass
     $, qmass, qmassj(plat)    ! total mass of moisture
     $, pint(plond,plevp)      ! pressure at layer interfaces
     $, pmid(plond,plev)       ! pressure at layer midpoints
     $, pdel(plond,plev)       ! pressure difference across layers
c-----------------------------------------------------------------------

      tmass = 0.0
      qmass = 0.0
      psdry = 0.0

      do j = 1, plat
         call plevs0( ps(1,j), pint, pmid, pdel )
         call gatmass( j, ps(1,j), tmassj(j) )
         call qmassa( j, 1, q(1,1,j), pdel, qmassj(j) )
         do k = 1, plev
            do i = 1, plon
               psdry(i,j) = psdry(i,j) + pdel(i,k)*(1.-q(i,k,j))
            end do
         end do
      end do

      do j = 1, plat
         tmass = tmass + tmassj(j)
         qmass = qmass + qmassj(j)
      end do

      mdryatm = tmass - qmass

      write(*,*)'inidry: dry mass of atmosphere will be held at '
     $,         mdryatm, ' kg/m^2'

      return
      end

c##############################################################################

      subroutine adjdry( psin, q, psout, diagprt )

c Adjust ps so that the total dry mass remains constant.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <drymass.com>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  psin(plond,plat)       ! surface pressure (pascals)
     $, q(plond,plev,plat)     ! specific humidity
      logical
     $  diagprt

c Output arguments:
      real(r8)
     $  psout(plond,plat)      ! adjusted surface pressure (pascals)

c Local variables:
      integer i, j, k
      real(r8)
     $  tmass, tmassj(plat)    ! total atm mass
     $, aqmass, aqmassj(plat)  ! A contribution to mass of moisture
     $, bqmass, bqmassj(plat)  ! B contribution to mass of moisture
     $, tmp1(plond,plev)
     $, tmp2(plond,plev)
     $, pint(plond,plevp)      ! pressure at interfaces
     $, fac
c-----------------------------------------------------------------------

      tmass = 0.0
      aqmass = 0.0
      bqmass = 0.0

!      call pgavg('adjdry: psin ',psin)

      do j = 1, plat
         call pdelab( psin(1,j), tmp1, tmp2 )
         call gatmass( j, psin(1,j), tmassj(j) )
         call qmassa( j, 1, q(1,1,j), tmp1, aqmassj(j) )
         call qmassa( j, 1, q(1,1,j), tmp2, bqmassj(j) )
      end do

      do j = 1, plat
         tmass = tmass + tmassj(j)
         aqmass = aqmass + aqmassj(j)
         bqmass = bqmass + bqmassj(j)
      end do

      fac = (mdryatm + aqmass) / (tmass - bqmass)

      do j = 1, plat
         do i = 1, plon
            psout(i,j) = psin(i,j) * fac
         end do
      end do
!      call pgavg('adjdry: psout wet ',psout)

c     Convert psout to a dry ps.
      do j = 1, plat
!        tmp1 contains pmid, tmp2 contains pdel
         call plevs0( psout(1,j), pint, tmp1, tmp2 )
         psout(:,j) = 0.0
         do k = 1, plev
            do i = 1, plon
               psout(i,j) = psout(i,j) + tmp2(i,k)*(1.-q(i,k,j))
            end do
         end do
      end do
!     call pgavg('adjdry: psout,dry ',psout)

      if ( diagprt ) then
         write(*,*)'adjdry: dry mass of atmosphere will be held at '
     $,             mdryatm, ' kg/m^2'
      end if

      return
      end

c##############################################################################

      subroutine rdryrst( lun )

c Read restart data for dry mass module.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <drymass.com>
c-----------------------------------------------------------------------

c Input argument:

      integer
     $  lun    ! logical unit for fortran i/o.
c-----------------------------------------------------------------------

c     Dry mass of atmosphere.

      read( lun ) mdryatm
      write (6,*) ' rdryrst: mdryatm is ', mdryatm

      return
      end

c##############################################################################

      subroutine wdryrst( lun )

c Write restart data for dry mass module.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <drymass.com>
c-----------------------------------------------------------------------

c Input argument:

      integer
     $  lun    ! logical unit for fortran i/o.
c-----------------------------------------------------------------------

c     Dry mass of atmosphere.

      write( lun ) mdryatm

      return
      end

c##############################################################################

      subroutine mr2wet( plond, plon, plev, ntrac, xdry, q, xwet )

c Convert tracer mixing ratios from a dry to a wet basis.

      use precision

      implicit none

c Input arguments:

      integer
     $  plond  ! 1st dim of arrays
     $, plon   ! no. of distinct longitudes
     $, plev   ! no. of levels
     $, ntrac  ! number of tracers
      real(r8)
     $  xdry(plond,plev,ntrac) ! Tracer concentration (kg tracer/kg dry air)
     $, q(plond,plev)         ! Specific humidity (kg water vapor/kg moist air)

c Output arguments:

      real(r8)
     $ xwet(plond,plev,ntrac)  ! Tracer concentration (kg tracer/kg moist air)

c Local variables:

      integer i, k, m
c-----------------------------------------------------------------------

      do m = 1, ntrac
         do k = 1, plev
            do i = 1, plon
               xwet(i,k,m) = xdry(i,k,m) * (1. - q(i,k))
            end do
         end do
      end do

      return
      end

c##############################################################################

      subroutine mr2dry( plond, plon, plev, ntrac, xwet, q, xdry )

c Convert tracer mixing ratios from a wet to a dry basis.

      use precision

      implicit none

c Input arguments:

      integer
     $  plond  ! 1st dim of arrays
     $, plon   ! no. of distinct longitudes
     $, plev   ! no. of levels
     $, ntrac  ! number of tracers
      real(r8)
     $  xwet(plond,plev,ntrac) ! Tracer concentration (kg tracer/kg moist air)
     $, q(plond,plev)         ! Specific humidity (kg water vapor/kg moist air)

c Output arguments:

      real(r8)
     $ xdry(plond,plev,ntrac)  ! Tracer concentration (kg tracer/kg dry air)

c Local variables:

      integer i, k, m
c-----------------------------------------------------------------------

      do m = 1, ntrac
         do k = 1, plev
            do i = 1, plon
               xdry(i,k,m) = xwet(i,k,m) / (1. - q(i,k))
            end do
         end do
      end do

      return
      end

c##############################################################################

      subroutine pgavg(title,ps)
      use precision
      implicit none
#include <pmgrid.h>
#include <mass.com>
      real(r8) gavg, ps(plond,plat), pi, dxlon
      character*(*) title
      integer i, j
      return
      pi = 4*atan(1.0)
      dxlon = 2*pi/plon
      gavg = 0.
      do j = 1,plat
         do i = 1,plon
            gavg = gavg+ps(i,j)*gw(j)
         end do
      end do
      gavg = gavg*dxlon/(4.*pi)
      write (6,*) 'PGAVG: ', title, gavg, 'or in kg/m2 ',gavg*rgrav
      return
      end 

      subroutine pdnavg(title,dn,deta)
      use precision
      implicit none
#include <pmgrid.h>
#include <mass.com>
      real(r8) gavg, dn(plond,plat,plev), pi, dxlon, ps(plond,plat), deta(plev)
      character*(*) title
      integer i, j, k

      return

      pi = 4*atan(1.0)
      dxlon = 2*pi/plon

      ps = 0.
      gavg = 0.
      do k = 1,plev
         gavg = gavg + deta(k)
         do j = 1,plat
            do i = 1,plon
               ps(i,j) = ps(i,j) + dn(i,j,k)*deta(k)
            end do
         end do
      end do
!     write (6,*) ' sum of detas are ', gavg

      gavg = 0.
      do j = 1,plat
         do i = 1,plon
            gavg = gavg+ps(i,j)*gw(j)
         end do
      end do
      gavg = gavg*dxlon/(4.*pi)
      write (6,*) 'PDNAVG: ', title, gavg, 'or in kg/m2 ',gavg*rgrav
      return
      end 
