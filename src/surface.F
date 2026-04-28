      module surface

      use precision

      real(r8)
     $  gravit       ! Gravitational acceleration
     $, rair         ! Gas constant for dry air

      real(r8) :: zero = 0.

      contains

c $Id$
      subroutine inisflx( xrair, xgravit, xrearth, xlam, xphi
     $,                   xgw, xoro )

      use precision

      implicit none
c----------------------------------------------------------------------------
#include <pmgrid.h>
c----------------------------------------------------------------------------

c Input arguments:
      real(r8) 
     $  xrair             ! Gas constant for dry air
     $, xgravit           ! Gravitational acceleration
     $, xrearth           ! radius of earth (m)
     $, xlam(plon)        ! longitudes in grid (radians)
     $, xphi(plat)        ! Gauss latitudes in grid (radians)
     $, xgw(plat)         ! Gauss weights (sum to 2.0)
     $, xoro(plond,plat)  ! ocean=0, land=1, sea ice=2
c----------------------------------------------------------------------------

c     Set constants for surface flux routines.
      rair = xrair
      gravit = xgravit

      return
      end subroutine inisflx

c##############################################################################

      subroutine sflxdr( lat, dtime, calday, ncdate, ioro
     $,                  as, pmid, rpdel, tv, zi
     $,                  obuf, sflx
     $)

c Set surface fluxes or make adjustment to field due to surface flux.
c Dry deposition routines are called from here because they work by setting
c surface fluxes.

      use histout, only: outfld
      use acsf, only: getsf
#ifdef SCYC
c++scyc
      use scyc, only: doscyc, scyc_idx1, add_volc_emis
      use sulemis, only: addsulemis
      use volcemission, only: volcemist
c--scyc
#endif
#ifdef CAER
c++caer
      use caer, only: docaer, caer_ncomp, caer_idx1, caersf
c--caer
#endif

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <tracnam.h>
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  lat                ! latitude index
     $, ncdate             ! date at begin of current timestep in yymmdd format
     $, ioro(plond)        ! orography 
      real(r8)
     $  dtime             ! timestep size (s)
     $, calday            ! Julian day + fraction (Greenwich)
     $, pmid(plond,plev)  ! pressure at layer midpoints
     $, rpdel(plond,plev) ! 1/pdel
     $, tv(plond,plev)    ! virtual temperature
     $, zi(plond,plevp)   ! potential height above surface at interfaces

c Input/Output arguments:
      real(r8)
     $  as(plond,plev,pcnst)   ! advected species
     $, obuf(*)                ! output buffer

c Output argument:
      real(r8)
     $  sflx(plond,pcnst)      ! surface flux for advected species

c Local variables:
      integer i, k, m
      real(r8) rlat                ! [rdn] latitude of this row
      real(r8) clat                ! external function (grid.F)
#ifdef SCYC
c++scyc
      real(r8) so2tend(plond,plev)
c--scyc
#endif
#ifdef CAER
c++caer
      integer mm, nphob
c--caer
#endif

c-----------------------------------------------------------------------

      call setr( pcnst*plond, zero, sflx )

      rlat = clat(lat)

c     Interface for surface fluxes from netCDF annual cycle data.
c     (do nothing if no annual cycle surface fluxes have been specified)
      call getsf( lat, sflx )

#ifdef SCYC
c++scyc
      if ( doscyc() ) then
         m = scyc_idx1()
c        Emissions are added to tracer array inside addsulemis.
         call addsulemis( lat, rlat, so2_rgn_scale, dms_rgn_scale
     $,                   dtime, gravit, rpdel, zi, as(1,1,m), obuf )
c        Volcano emissions.
         if ( add_volc_emis() ) then
            m = scyc_idx1()
            call volcemist( lat, gravit, rpdel, so2tend )
            call outfld( 'SO2VOLC', so2tend, plond, lat, obuf )
            do k = 1, plev
               do i = 1, plon
                  as(i,k,m) = as(i,k,m) + dtime*so2tend(i,k)
               end do
            end do
         end if
      end if
c--scyc
#endif
#ifdef SCYC
c++caer
      if ( docaer() ) then
c        Carbon aerosol surface flux(es).
         mm = caer_idx1()
         nphob = caer_ncomp()/2
         call caersf( nphob, lat, rlat, oc_rgn_scale, bc_rgn_scale
     $,               sflx(1,mm) )
         do m = 1, nphob
            mm = caer_idx1() + m - 1
            call outfld( sflxnam(mm), sflx(1,mm), plond, lat, obuf )
         end do
      end if
c--caer
#endif

      return
      end subroutine sflxdr

c##############################################################################

      end module surface
