c $Id$
      subroutine physcan( dtime, calday, ncdate, ncsec, nstep
     $,                   ps, oro, phis, x0, t1, sh1
     $,                   uv, ts, taux, tauy, hflx, shflx
c     $,                   fsds, fsus, flds, vwcsfc, snowh, qrad
     $,                   fsds, fsus, flds, vwcsfc, snowh
     $,                   etadot, omga, cgs, kvh
     $,                   zmu, zmd, zeu, hketa, hkbeta, cmfdqr, conicw
     $,                   cldfr, concld, cme, prain, evapr, cwat
     $,                   as0, as1, x1
     $,                   wr, mcar, mca 
     $)

c Latitude loop over physics parameterizations.

      use precision
      implicit none
c------------------------ parameters -----------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  ncdate      ! date at beginning of current timestep in yymmdd format
     $, ncsec       ! seconds relative to ncdate
     $, nstep       ! timestep counter
      real(r8)
     $  dtime                   ! timestep
     $, calday                  ! calandar day in range [1., 366.)
     $, ps(plond,plat)          ! surface pressure (Pa)
     $, oro(plond,plat)         ! orography
     $, phis(plond,plat)        ! surface geopotential (m^2/s^2)
     $, x0(plond,plev,padv,plat)! misc advected variables before slt
     $, t1(plond,plev,plat)     ! temperature (K)
     $, sh1(plond,plev,plat)    ! specific humidity (kg/kg)
     $, uv(plond,plev,2,plat)   ! uv wind
     $, ts(plond,plat)          ! surface temperature
     $, taux(plond,plat)        ! x surface stress (n)
     $, tauy(plond,plat)        ! y surface stress (n)
     $, hflx(plond,plat)        ! surface sensible heat flux (w/m2)
     $, shflx(plond,plat)       ! surface water vapor flux (kg/m2/s)
     $, fsds(plond,plat)        ! downward solar flux at surface (W/m2)
     $, fsus(plond,plat)        ! upward solar flux at surface (W/m2)
     $, flds(plond,plat)        ! downward longwave flux at surface (W/m2)
     $, vwcsfc(plond,plat)      ! volumetric soil water content at surface (frac)
     $, snowh(plond,plat)       ! water equivalent accumulated snow depth (m)
     $, qrad(plond,plev,plat)   ! radiative heating tendency (K/s)
     $, etadot(plon,plevp,plat) ! vertical velocity in hybrid coord
     $, wr(plon,plevp,plat)     ! raw vertical velocity in hybrid coord
     $, mca(plon,plevp,plat)    ! mass convergence above
     $, mcar(plon,plevp,plat)   ! raw mass convergence above
     $, omga(plond,plev,plat)   ! vertical velocity in pressure coord
     $, cgs(plond,plevp,plat)   ! counter-gradient coefficient
     $, kvh(plond,plevp,plat)   ! vertical diffusion coefficient
     $, zmu(plond,plev,plat)    ! mu2 from conv_ccm, kg/m2/s
     $, zmd(plond,plev,plat)    ! md2 from conv_ccm, kg/m2/s
     $, zeu(plond,plev,plat)    ! eu2 from conv_ccm, 1/s
     $, hketa(plond,plev,plat)  ! convective mass flux, Hack (kg/m^2/s)
     $, hkbeta(plond,plev,plat) ! overshoot parameter, Hack (fraction)
     $, cmfdqr(plond,plev,plat) ! dq/dt due to convective rainout 
     $, cldfr(plond,plevp,plat) ! cloud fraction
     $, concld(plond,plev,plat) ! convective cloud fraction
     $, cme(plond,plev,plat)    ! rate of cond-evap within the cloud
     $, prain(plond,plev,plat)  ! rate of release of stratiform precip (1/s)
     $, evapr(plond,plev,plat)  ! rate of evaporation of precipitation (1/s)
     $, cwat(plond,plev,plat)   ! cloud water (kg/kg)
     $, as0(plond,plev,pcnst,plat) ! advected species before slt

c Input/Output argument:
      real(r8)
     $  as1(plond,plev,pcnst,plat) ! advected species after slt on input,
c                                    adjustment and tendency physics applied
c                                    on output.
     $, x1(plond,plev,padv,plat) ! misc advected variables after slt on input.
c                                  Contrib from physics added on output.
     $, conicw(plond,plev,plat)  ! in cloud water mixing ratio (convective)

c Local variables:
      integer j

c-----------------------------------------------------------------------

! OpenMP directive
!$OMP PARALLEL DO
      do j = 1, plat
         call physlic( j, dtime, calday, ncdate, ncsec, nstep
     $,                ps(1,j), oro(1,j), phis(1,j)
     $,                x0(1,1,2,j), x0(1,1,1,j), t1(1,1,j), sh1(1,1,j)
     $,                uv(1,1,1,j), ts(1,j), taux(1,j), tauy(1,j)
     $,                hflx(1,j), shflx(1,j), fsds(1,j), fsus(1,j)
     $,                flds(1,j), vwcsfc(1,j), snowh(1,j)
     $,                qrad(1,1,j), etadot(1,1,j), omga(1,1,j)
     $,                cgs(1,1,j), kvh(1,1,j)
     $,                zmu(1,1,j), zmd(1,1,j), zeu(1,1,j)
     $,                hketa(1,1,j), hkbeta(1,1,j)
     $,                cmfdqr(1,1,j), conicw(1,1,j), cldfr(1,1,j)
     $,                concld(1,1,j), cme(1,1,j)
     $,                prain(1,1,j), evapr(1,1,j), cwat(1,1,j)
     $,                as0(1,1,1,j), as1(1,1,1,j), x1(1,1,2,j)
     $,                x1(1,1,1,j)
     $,                wr(1,1,j), mca(1,1,j), mcar(1,1,j)
     $)

      end do

      return
      end


