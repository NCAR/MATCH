c $Id$
#include <params.h>

      subroutine drydepdr( lat, month, dtime, oro, sn, prect
     $,                    t2m, wvflx, shflx, pblh, ustar, obklen
     $,                    pmbot, rpdel, tvmbot, zi, um1, vm1, tm1
     $,                    swradsf, tr, obuf )

c Driver for dry deposition parameterization.  The default version of this
c routine does nothing.  It is provided only as a template.
c N.B. The input arguments may not all be required for some parameterizations
c      but they are included in an attempt to create a general interface that
c      is meets the needs of a range of parameterizations.

      use histout, only: outfld
c++scyc
      use scyc, only: doscyc, scyc_idx1
c--scyc
c++caer
      use caer, only: docaer, caer_ncomp, caer_idx1
c--caer
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <drydep.h>
c-----------------------------------------------------------------------
#include <tracnam.h>
c-----------------------------------------------------------------------

c Input arguments
      integer
     $  lat                   ! latitude index
     $, month                 ! month of year [1-12]
      real(r8)
     $  dtime                 ! timestep size
     $, oro(plond)            ! orography
     $, sn(plond)             ! snow; equivalent water amount (m)
     $, prect(plond)          ! total precipitation (units arbitrary as long
c                             ! the value is > 1.e-12 when it's raining)
     $, t2m(plond)            ! temperature at 2 m
     $, wvflx(plond)          ! water vapor flux (kg/m^2/s)
     $, shflx(plond)          ! sensible heat flux (w/m^2)
     $, pblh(plond)           ! planetary boundary layer height
     $, ustar(plond)          ! friction velocity
     $, obklen(plond)         ! Obukhov length
     $, pmbot(plond)          ! midpoint pressure in bottom layer
     $, rpdel(plond,plev)     ! 1./(pint(k+1)-pint(k))
     $, tvmbot(plond)         ! midpoint virtual temperature in bottom layer
     $, zi(plond,plevp)       ! height at interface of levels
     $, um1(plond,plev)       ! u horizontal wind component
     $, vm1(plond,plev)       ! v horizontal wind component
     $, tm1(plond,plev)       ! temperature
     $, swradsf(plond)        ! solar radiation at the surface (W/m2)

c Input/Output arguments.
      real(r8)
     $  tr(plond,plev,pcnst)  ! transported species (kg tracer/kg moist air)

c Output arguments.
      real(r8)
     $  obuf(*)               ! output buffer

c Local variables
      integer
     $  i, m, mm
     $, ioro(plond)            ! season indice found in getseas
      real(r8)
     $  dvel(plond)            ! deposition velocity
     $, flx(plond,pcnst)       ! dry deposition fluxes
c----------------------------------------------------------------------c


      do i = 1, plon
         ioro(i) = nint( oro(i) )
      end do
      do m = 1, pcnst
         do i = 1, plon
            flx(i,m) = 0.0
         end do
      end do

c++scyc
      if ( doScyc() ) then
c        Dry Deposition Velocities:
c        Here are the values used by Hans Feichter for SO2:
c        .6cm/s over land
c        .8cm/s over ocean
c        .1cm/s over ice or snow
c        For SO4 he uses .2cm/s everywhere

c        Dry deposition of SO2
c        #####################
         mm = scyc_idx1()
         call setdvel( ioro, .006_r8, .008_r8, .001_r8, dvel )
         call ddflux( dvel, tr(1,plev,mm), pmbot, tvmbot, flx(1,mm) )
         call outfld( trim(tracnam(mm))//'DRY', flx(1,mm), plond
     $,               lat, obuf )
#ifdef SCYC_MASSBGT
         call aal( 'so2dry', lat, flx(1,mm) )
#endif

c        Dry deposition of SO4
c        #####################
         mm = scyc_idx1()+1
         call setdvel( ioro, .002_r8, .001_r8, .002_r8, dvel )
         call ddflux( dvel, tr(1,plev,mm), pmbot, tvmbot, flx(1,mm) )
         call outfld( trim(tracnam(mm))//'DRY', flx(1,mm), plond
     $,               lat, obuf )
#ifdef SCYC_MASSBGT
         call aal( 'so4dry', lat, flx(1,mm) )
#endif

c        Dry deposition of H2O2
c        #####################
         mm = scyc_idx1()+3
         call setdvel( ioro, .006_r8, .008_r8, .001_r8, dvel )
         call ddflux( dvel, tr(1,plev,mm), pmbot, tvmbot, flx(1,mm) )
         call outfld( trim(tracnam(mm))//'DRY', flx(1,mm), plond
     $,               lat, obuf )
      end if
c--scyc

c++caer
      if ( docaer() ) then

c        Dry deposition of Carbon Aerosols
c        #################################
         call setdvel( ioro, .001_r8, .001_r8, .001_r8, dvel )
         do m = 1, caer_ncomp()
            mm = caer_idx1() + m - 1
            call ddflux( dvel, tr(1,plev,mm), pmbot, tvmbot, flx(1,mm) )
            call outfld( trim(tracnam(mm))//'DRY', flx(1,mm), plond
     $,                  lat, obuf )
         end do
      end if
c--caer

c     Update tracers with dry deposition contribution.
      do m = 1, pcnst
         do i = 1, plon
            tr(i,plev,m) = tr(i,plev,m) + 
     $                     flx(i,m)*dtime*gravit*rpdel(i,plev)
         end do
      end do

      return
      end

c##############################################################################

