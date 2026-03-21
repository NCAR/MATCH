      module cloud

      use precision

      contains
c $Id$

      subroutine inicld( xphi, xcappa, xgravit, xrair, xlatvap, xcpair
     $,                  xrhow, xtfh2o, etamid )

      implicit none
c----------------------------------------------------------------------------
#include <pmgrid.h>
c----------------------------------------------------------------------------
#include <cloud.com>
c----------------------------------------------------------------------------

c Input arguments:

      real(r8) 
     $  xphi(plat)  ! Gaussian latitudes (radians)
     $, xcappa      ! R/Cp
     $, xgravit     ! Gravitational acceleration
     $, xrair       ! Gas constant for dry air
     $, xlatvap     ! latent heat of vaporization for liq. water
     $, xcpair      ! Specific heat of dry air
     $, xrhow       ! water density (g/cm^3)
     $, xtfh2o      ! freezing point of water at STP (K)
     $, etamid(plev)! eta coordinate at layer midpoints

c Local variables:

      integer j, k, i
      real(r8) detamn, tmp
c----------------------------------------------------------------------------

c  set the physical constants for cloud routines

      do j = 1,plat
         phi(j) = xphi(j)
         coslat(j) = cos(phi(j))
      end do

      cappa = xcappa
      gravit = xgravit
      rair = xrair
      latvap = xlatvap
      cpair = xcpair
      rhow = xrhow
      tfh2o = xtfh2o

c     find level closest to 700mb over ocean
      detamn = abs( .7 - etamid(1) )
      do k = 2, plev
         tmp = abs( .7 - etamid(k) )
         if ( tmp .lt. detamn ) then
            k700 = k
            detamn = tmp
         end if
      end do

      return
      end subroutine inicld

c##############################################################################

      subroutine clouddr( lat, dtime, oro, pint, pmid, pdel, rpdeli
     $,                   ts, t, sh, omga, cmfmc, zdu, cldtop, cldbot
     $,                   dlf, tm1, shm1, cmfdqr, cwat
     $,                   cldt, cldc, cme, evapr, prain, precl
     $,                   obuf )

c Cloud physics calcs.

      use histout, only: outfld
      use cldwat, only: pcond

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <cloud.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  lat         ! latitude index

      real(r8)
     $  dtime                ! timestep size (s)
     $, oro(plond)           ! orography
     $, pint(plond,plevp)    ! pressure at layer interfaces
     $, pmid(plond,plev)     ! pressure at layer midpoints
     $, pdel(plond,plev)     ! pressure difference across layers
     $, rpdeli(plond,plev)   ! 1./(pmid(k+1)-pmid(k))
     $, ts(plond)            ! surface temperature
     $, omga(plond,plev)     ! vertical pressure velocity
     $, cmfmc(plond,plev)    ! convective mass flux (kg/m^2/s)
     $, zdu(plond,plev)      ! du2 from conv_ccm, 1/s
     $, cldtop(plond)        ! top level of convection
     $, cldbot(plond)        ! bottom level of convection
     $, dlf(plond,plev)      ! detraining cloud water from convection
     $, tm1(plond,plev)      ! temperature (K), tdyn(n-1)
     $, shm1(plond,plev)     ! specific humidity (kg/kg), shdyn(n-1)
     $, cmfdqr(plond,plev)   ! dq/dt due to convective rainout 

c Input/output arguments:

      real(r8)
     $  t(plond,plev)        ! temperature (K)
c                              tdyn(n-1)+adv+vdiff+conv on input,
c                              tendency due to cloud added on output
     $, sh(plond,plev)       ! specific humidity (kg/kg)
c                              shdyn(n-1)+adv+vdiff+conv on input
c                              tendency due to cloud added on output
     $, cwat(plond,plev)     ! cloud water (kg/kg)

c Output arguments:

      real(r8)
     $  cldt(plond,plevp)    ! total cloud fraction
     $, cldc(plond,plev)     ! convective cloud fraction
     $, cme(plond,plev)      ! rate of cond-evap within the cloud
     $, evapr(plond,plev)    ! rate of evaporation of falling precipitation (kg/kg/s)
     $, prain(plond,plev)    ! rate of conversion of condensate to precipitation (kg/kg/s)
     $, precl(plond)         ! large scale precip (m/s)
      real(r8) obuf(*)         ! output buffer

c Local variables:

      integer i, k

c     cloud fraction
      real(r8)
     $  clc(plond)           ! column convective cloud amount
     $, cldst(plond,plev)    ! cloud fraction (strat?)
     $, snowh(plond)         ! snow depth (liquid water equivalent)
      real(r8) cltot(plond)      ! diagnosed tot random overlap cld cover
      real(r8) cllow(plond)      ! diagnosed low random overlap cld cover
      real(r8) clmed(plond)      ! diagnosed med random overlap cld cover
      real(r8) clhgh(plond)      ! diagnosed hgh random overlap cld cover

c     moist
      real(r8)
     $  qtend(plond,plev)    ! moisture tendencies (kg/kg/s)
     $, ttend(plond,plev)    ! temp tendencies (K/s)
     $, rmelt(plond,plev)    ! heating rate due to precip phase change (K/s) (DISABLED)
     $, pcflx(plond,plevp)   ! convective precip level by level (kg/m2/s) (DISABLED)
     $, fracw(plond,plev)    ! relative importance of collection of liquid by rain
     $, fsaci(plond,plev)    ! relative importance of collection of ice by snow
     $, fsacw(plond,plev)    ! relative importance of collection of liquid by snow
     $, fsaut(plond,plev)    ! relative importance of ice auto conversion
     $, fwaut(plond,plev)    ! relative importance of warm cloud auto conversion

c-----------------------------------------------------------------------

c     ******************
c     * Cloud fraction *
c     ******************

      do i = 1, plon
         snowh(i) = 0.
      end do

c     Compute cloud amount.
      call cldfrc( pmid, rpdeli, t, sh, omga
     $,            cldtop, cldbot, cldt, clc, pdel
     $,            cmfmc, oro, snowh, cldc, cldst
     $,            ts, pint(1,plevp), zdu, lat )

c     Cloud cover diagnostics.
      call cldsav( cldt, pmid, cltot, cllow, clmed, clhgh )
      call outfld( 'CLOUD   ', cldt,   plond, lat, obuf )
      call outfld( 'CONCLD  ', cldc,   plond, lat, obuf )
      call outfld( 'CLDTOT  ', cltot,  plond, lat, obuf )
      call outfld( 'CLDLOW  ', cllow,  plond, lat, obuf )
      call outfld( 'CLDMED  ', clmed,  plond, lat, obuf )
      call outfld( 'CLDHGH  ', clhgh,  plond, lat, obuf )

c     Put the detraining cloud water into the cloud and environement
c     in proportion to the cloud fraction.
      do k = 1, plev
         do i = 1, plon
            sh(i,k) = sh(i,k) + dtime*dlf(i,k)*(1.-cldt(i,k))
            t(i,k) = t(i,k) 
     $              - latvap/cpair*dtime*dlf(i,k)*(1.-cldt(i,k))
            cwat(i,k) = cwat(i,k) + dtime*dlf(i,K)*cldt(i,k)
            if ( cwat(i,k) .lt. 1.e-12 ) cwat(i,k) = 0.0
            qtend(i,k) = (sh(i,k) - shm1(i,k))/dtime
            ttend(i,k) = (t(i,k) - tm1(i,k))/dtime
         end do
      end do

c     ******************************************************
c     * Stratiform condensation via prognostic cloud water *
c     ******************************************************

      do i = 1, plon
         precl(i) = 0.
      end do
      do k = 1, plevp
         do i = 1, plon
            pcflx(i,k) = 0.
         end do
      end do

      call pcond( qtend, ttend, omga, lat
     $,           shm1, cwat, tm1, pmid, pdel, oro
     $,           cldt, cldt, cme, evapr, prain, rmelt, dtime, pcflx
     $,           fwaut, fsaut, fracw, fsacw, fsaci )

      call outfld( 'CME     ', cme,    plond, lat, obuf )
      call outfld( 'PRAIN   ', prain,  plond, lat, obuf )
      call outfld( 'EVAPR   ', evapr,  plond, lat, obuf )


      do k = 1,plev
         do i = 1,plon
            t(i,k) = t(i,k) + latvap/cpair*(cme(i,k)-evapr(i,k))*dtime
     $                      + rmelt(i,k)*dtime
            sh(i,k) = sh(i,k) - (cme(i,k) - evapr(i,k))*dtime
            cwat(i,k) = cwat(i,k) + (cme(i,k) - prain(i,k))*dtime
            precl(i) = precl(i) 
     $               + (prain(i,k)-evapr(i,k))*pdel(i,k)/gravit
            if ( cwat(i,k) .lt. 1.e-12 ) cwat(i,k) = 0.0
         end do
      end do

c     Convert precl from kg/m2/s to m/s
      do i = 1, plon
         precl(i) = precl(i)*.001
      end do
      call outfld( 'PRECL   ', precl,  plond, lat, obuf )

      return
      end subroutine clouddr

c##############################################################################

      subroutine clddiag( t, pmid, pdel, cmfdqr, cldt, cme, evapr, prain
     $,                   cldv, rain )

c Diagnostic cloud physics calcs.

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <cloud.com>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  t(plond,plev)        ! temperature (K)
     $, pmid(plond,plev)     ! pressure at layer midpoints
     $, pdel(plond,plev)     ! pressure difference across layers
     $, cmfdqr(plond,plev)   ! dq/dt due to convective rainout 
     $, cldt(plond,plevp)    ! total cloud fraction
     $, cme(plond,plev)      ! rate of cond-evap within the cloud
     $, evapr(plond,plev)    ! rate of evaporation of falling precipitation (kg/kg/s)
     $, prain(plond,plev)    ! rate of conversion of condensate to precipitation (kg/kg/s)

c Output arguments:
      real(r8)
     $  cldv(plond,plev)     ! fraction occupied by rain or cloud water 
     $, rain(plond,plev)     ! mixing ratio of rain (kg/kg)

c Local variables:
      integer
     $  i, k
      real(r8)
     $  convfw         ! used in fallspeed calculation; taken from findmcnew
     $, sumppr(plond)        ! precipitation rate (kg/m2-s)
     $, sumpppr(plond)       ! sum of positive precips from above
     $, cldv1(plond)         ! precip weighted cloud fraction from above
     $, lprec                ! local production rate of precip (kg/m2/s)
     $, lprecp               ! local production rate of precip (kg/m2/s) if positive
     $, rho                  ! air density
     $, vfall
c-----------------------------------------------------------------------

c estimate the cloudy volume which is occupied by rain or cloud water as
c the max between the local cloud amount or the
c
c sum above of (cloud*positive precip production)      sum total precip from above
c              ----------------------------------   x ------------------------
c sum above of     (positive precip           )        sum positive precip from above

      convfw = 1.94*2.13*sqrt(rhow*1000.*gravit*2.7e-4)
      do i=1,plon
         sumppr(i) = 0.
         cldv1(i) = 0.
         sumpppr(i) = 1.e-36
      end do

      do k = 1,plev
         do i = 1,plon
            cldv(i,k) =
     $           max(min(1._r8,
     $           cldv1(i)/sumpppr(i)
     $           )*sumppr(i)/sumpppr(i),
     $           cldt(i,k)
     $           )
            lprec = pdel(i,k)/gravit
     $           *(prain(i,k)+cmfdqr(i,k)-evapr(i,k))
            lprecp = max(lprec,1.e-30_r8)
            cldv1(i) = cldv1(i)  + cldt(i,k)*lprecp
            sumppr(i) = sumppr(i) + lprec
            sumpppr(i) = sumpppr(i) + lprecp

            rain(i,k) = 0.
            if(t(i,k) .gt. tfh2o) then
c  qr = ppr/(rho*Vr)
               rho = pmid(i,k)/(rair*t(i,k))
               vfall = convfw/sqrt(rho)
               rain(i,k) = sumppr(i)/(rho*vfall)
               if (rain(i,k).lt.1.e-14) rain(i,k) = 0.
            endif
         end do
      end do

      return
      end subroutine clddiag

      end module cloud
