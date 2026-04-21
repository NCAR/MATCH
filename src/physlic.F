c $Id$
#undef CALC_ADV_FLUX
#ifdef DST
c++csz
#include <dst.h> /* Dust preprocessor tokens */
c--csz
#endif
      subroutine physlic( lat, dtime, calday, ncdate, ncsec, nstep
     $,                   ps, oro, phis
     $,                   pt0, sh0, t1, sh1
     $,                   uv, ts, taux, tauy, hflx, wvflx
     $,                   fsds, fsus, flds, vwcsfc, snowh, qrad
     $,                   etadot, omga
     $,                   cgs, kvh, zmu, zmd, zeu, hketa, hkbeta
     $,                   cmfdqr, conicw, cldt, cldc
     $,                   cme, prain, evapr, cwat
     $,                   as0, as1, ptadv, shadv
     $,                   wr, mca, mcar
     $)

c Physics parameterizations and output fields.

      use histout, only: outfld, qhsti

#ifdef CALC_ADV_FLUX
        use advect, only: advfluxout
#endif           


#ifdef SCYC
c++scyc
      use scyc, only: doscyc, scyc_idx1
      use sulchem, only: chemwdepdr
c--scyc
#endif
#ifdef CAER
c++caer
      use caer, only: docaer, caer_idx1, caer_ncomp
c--caer
#endif

      use cloud, only: clddiag
      use surface, only: sflxdr
      use vdiff_mod, only: vdiffdr
      use cloud, only: clouddr
#ifdef WETDEP
      use wetdep
#endif
#ifdef DST
      use dstdpswet
      use dstdpsdry, only: dst_dps_dry
      use dstmbl, only: dst_mbl
#endif
      use precision
      use rec_io, only: wrarr
#ifdef ASSIM
      use aod_assim_3d, only: asmfluxout
#endif

      implicit none
c------------------------ parameters -----------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------
#include <tracnam.h>
c-----------------------------------------------------------------------
C+WDC 10/26/00 common for sea-salt
#ifdef SCYC
#include <ssltcom.h>
#endif
c-----------------------------------------------------------------------

c++csz
#ifdef DST
#include <dstgrd.h> /* Dust grid sizes */
#include <dstnm.h> /* Nomenclature for outfld() */
c     Input
c     Local
      real(r8) lat_rdn              ! [rdn] Latitude
      real(r8) fsns(plond)          ! [W m-2] Solar flux absorbed by ground
      real(r8) mno_lng(plond)       ! [m] Monin-Obukhov length
      real(r8) wnd_frc(plond)       ! [m s-1] Surface friction velocity
      real(r8) wnd_rfr(plond)       ! [m s-1] Wind speed at reference height
c     Diagnostic totals for outfld
      real(r8) mpc_dst(plond,dst_nbr) ! [kg m-2] Column mass path of dust
      real(r8) mpc_dst_ttl(plond)   ! [kg m-2] Total column mass path of dust
      real(r8) odxc_dst(plond,dst_nbr) ! [frc] Column dust optical depth
      real(r8) odxc_dst_ttl(plond)  ! [frc] Total column dust optical depth
      real(r8) q_dst_ttl(plond,plev) ! [kg kg-1] Total dust mixing ratio
      real(r8) hbuf
#endif /* not DST */ 
c--csz

c Input arguments:
      integer
     $  lat         ! latitude index for this slice
     $, ncdate      ! date at beginning of current timestep in yyyymmdd format
     $, ncsec       ! seconds relative to ncdate
     $, nstep       ! timestep counter
      real(r8)
     $  dtime              ! timestep
     $, calday             ! calendar day in range [1., 366.)
     $, ps(plond)          ! surface pressure (Pa)
     $, oro(plond)         ! orography
     $, phis(plond)        ! surface geopotential (m^2/s^2)
     $, pt0(plond,plev)    ! potential temperature from input data
                           ! interpolated to beginning of current timestep, (K)
     $, sh0(plond,plev)    ! advected specific humidity (possibly nudged by input data)
                           ! at beginning of current timestep, (kg/kg)
     $, t1(plond,plev)     ! temperature from input data interpolated to
                           !  end of current timestep, (K)
     $, sh1(plond,plev)    ! specific humidity from input data interpolated
                           ! to end of current timestep, (kg/kg)
     $, uv(plond,plev,2)   ! uv wind
     $, ts(plond)          ! surface temperature
     $, taux(plond)        ! x surface stress (n)
     $, tauy(plond)        ! y surface stress (n)
     $, hflx(plond)        ! surface sensible heat flux (w/m2)
     $, wvflx(plond)       ! surface water vapor flux (kg/m2/s)
     $, fsds(plond)        ! downward solar flux at surface (W/m2)
     $, fsus(plond)        ! upward solar flux at surface (W/m2)
     $, flds(plond)        ! downward longwave flux at surface (W/m2)
     $, vwcsfc(plond)      ! volumetric soil water content at surface (frac)
     $, snowh(plond)       ! water equivalent accumulated snow depth (m)
     $, qrad(plond,plev)   ! radiative heating tendency (K/s)
     $, etadot(plon,plevp) !vertical velocity (etadot)
     $, wr(plon,plevp)     ! raw vertical velocity (etadot)
     $, mca(plon,plevp)    ! mass convergence above
     $, mcar(plon,plevp)   ! raw mass convergence above
     $, omga(plond,plev)   !vertical velocity (dp/dt)
     $, cgs(plond,plevp)   ! counter-gradient coefficient
     $, kvh(plond,plevp)   ! vertical diffusion coefficient
     $, zmu(plond,plev)    ! mu2 from conv_ccm, kg/m2/s
     $, zmd(plond,plev)    ! md2 from conv_ccm, kg/m2/s
     $, zeu(plond,plev)    ! eu2 from conv_ccm, 1/s
     $, hketa(plond,plev)  ! convective mass flux, Hack convection (kg/m^2/s)
     $, hkbeta(plond,plev) ! overshoot parameter, Hack convection (fraction)
     $, cmfdqr(plond,plev) ! dq/dt due to convective rainout 
     $, cldt(plond,plevp)  ! total cloud fraction
     $, cldc(plond,plev)   ! convective cloud fraction
     $, prain(plond,plev)  ! rate of release of stratiform precip (1/s)
     $, evapr(plond,plev)  ! rate of evaporation of precipitation (1/s)
     $, cwat(plond,plev)   ! cloud water (kg/kg)
     $, as0(plond,plev,pcnst) ! advected species at beginning of current
                              ! timestep.

c Input/Output arguments:
      real(r8)
     $  as1(plond,plev,pcnst) ! transported species after advection on input.
                              ! vdiff and convection tends added on output.
     $, ptadv(plond,plev)     ! potential temperature from input data t0
                              ! advected one timestep on input.  vdiff, qrad,
                              ! convection and cloud tendencies added on output.
     $, shadv(plond,plev)     ! internal specific humidity after
                              ! advection on input.  vdiff,
                              ! convection and cloud tendencies added on output.
     $, conicw(plond,plev)    ! in cloud water mixing ratio (convective)
     $, cme(plond,plev)       ! rate of cond-evap within the cloud

c Local variables:

      real(r8), dimension(:), allocatable ::
     $  obuf                ! output buffer for CCM history files

      integer
     $  i, k, m  ! indices
     $, istat    ! return status
     $, maxsiz   ! output buffer length
     $, lunhst   ! fortran unit for output history file
     $, lunnas   ! fortran unit for out of core non-advected species
     $, ioro(plond) ! nint(oro)
     $, month    ! current month

      real(r8)
     $  pint(plond,plevp)      ! pressure at layer interfaces
     $, pmid(plond,plev)       ! pressure at layer midpoints
     $, pdel(plond,plev)       ! pressure difference across layers
     $, pmln(plond,plev)       ! log of midpoint pressure
     $, rpdel(plond,plev)      ! 1/pdel
     $, rpdeli(plond,plev)     ! 1./(pmid(k+1)-pmid(k))
     $, t0(plond,plev)         ! temperature corresponding to pt0
     $, tadv(plond,plev)       ! temperature corresponding to ptadv
     $, tv(plond,plev)         ! virtual temperature
     $, zm(plond,plev)         ! potential height above surface at midpoints
     $, zi(plond,plevp)        ! potential height above surface at interfaces
     $, zma(plond,plev)        ! absolute potential height at midpoints
     $, zia(plond,plevp)       ! absolute potential height at interfaces
     $, phima(plond,plev)      ! absolute geopotential at midpoints
     $, xtmp(plond,plev)       ! temp storage for output tracers

c     process tendency calc.
      integer mlast
      real(r8)
     $  mrsav(plond,plev,pcnst)! save mixing ratio before process
     $, shsav(plond,plev)      ! save sh before process
     $, cwsav(plond,plev)      ! save cwat before process
     $, tend(plond,plev)       ! process tendency
     $, flux(plond,plev)       ! advective mass flux

c     vertical diffusion
      real(r8)
     $  sflx(plond,pcnst)     ! surface flux for advected species
     $, pblh(plond)           ! planetary boundary layer height
     $, tpert(plond)          ! convective temperature excess
     $, qpert(plond)          ! convective humidity excess
     $, wpert3d(plond,plev)
     $, ustar(plond)          ! surface friction velocity
     $, obklen(plond)         ! Obukhov length

c     convective mass flux
      integer
     $  ctf(pcnst)            ! flag indicating species to transport
      real(r8)
     $  cmfmc(plond,plev)      ! moist convection cloud mass flux
     $, zdu(plond,plev)        ! du2 from conv_ccm, 1/s
     $, cldtop(plond)          ! top level of convective activity   
     $, cldbot(plond)          ! bottom level of convective activity
     $, dlf(plond,plev)        ! detraining cloud water from convection
     $, pflx(plond,plevp)      ! conv rain flux thru btm of that lev
     $, precc(plond)           ! convective precip (m/s)
     $, fracis(plond,plev,pcnst)  ! fraction of tracer in interstitial gas
c     gathered arrays of driving data for zhang convective transport
      real(r8)
     $  zmug(plond,plev)
     $, zmdg(plond,plev)
     $, zdug(plond,plev)
     $, zeug(plond,plev)
     $, zedg(plond,plev)
     $, zdpg(plond,plev)
     $, dsubcld(plond)
      integer
     $  zjtg(plond)
     $, zjbg(plond)
     $, ideep(plond)
     $, lengath

c     prognostic cloud water
      real(r8)
     $  cldv(plond,plev)     ! fraction occupied by rain or cloud water 
     $, rain(plond,plev)     ! mixing ratio of rain (kg/kg)
     $, precl(plond)         ! large scale precip (m/s)

c     Dry Deposition
      real(r8)
     $  prect(plond)         ! total precip (m/s)

c     Diagnostics
      real(r8)
     $  ci(plond)            ! column integral
     $, dqr(plond,plev)      ! dq/dt due to convective + large scale rainout 
     $, rh(plond,plev)       ! relative humidity

c++sslt
      real(r8) u10(plond)        ! wind speed at 10m for neutral BL
      real(r8) tausslt(plond)    ! estimate of sea salt aer opt depth
      real(r8) sslt(plond,plev)  ! sea salt mixing ratio
c--sslt
c++optdep
      integer mm
      real(r8) aerod(plond,pcnst), totaerod(plond)
c--optdep

c Externals:
      logical dryDep, outputTS
c++mahow
      logical doWetDepNew
c--mahow
c++scyc
      real(r8) clat
c--scyc

      real(r8) debdst(plond,plev,pcnst), fracisdeb(plond,plev,pcnst)
     $     , debscavt(plond,plev), debiscavt(plond,plev)
      real(r8) unused(plond,plev), qdstmx
      integer idstmx, kdstmx

      real(r8) cldc2(plond,plev), nomvel, sdiag1(plond), sdiag2(plond), fudge
      real(r8) mdiag1(plond,plev)
      integer iscat

#undef DEBUG
#ifdef DEBUG
      integer ilook,jlook,klook
      ilook = 99
      jlook = 72
      klook = 23
#endif

c-----------------------------------------------------------------------

c++mahow
      doWetDepNew = .false.
c--mahow
c     *************************
c     * Set up output buffer. *
c     *************************

#ifdef DEBUG
         if ( lat.eq.jlook ) write(6,*)'physlic 0.1 ',as1(ilook,klook,pcnst)
     $     ,as0(ilook,klook,pcnst)
#endif 

      if ( outputTS() ) then
         istat = qhsti( 'MAXSIZ', maxsiz )
         if ( istat .ne. 0 ) then
            write(*,*)'physlic: qhsti returned error for MAXSIZ'
            stop
         end if
         allocate( obuf(maxsiz), stat=istat )
         call alloc_err( istat, 'physlic', 'obuf', maxsiz )
         obuf = 0.0
      end if

c     *****************************
c     * Initialize some variables *
c     *****************************

      month = mod( ncdate, 10000 ) / 100

      do i = 1, plon
         ioro(i) = nint( oro(i) )
      end do

      call plevs0( ps, pint, pmid, pdel )

      do k = 1, plev
         do i = 1, plon
            pmln(i,k) = log( pmid(i,k) )
            rpdel(i,k) = 1. / pdel(i,k)
         end do
      end do
      do k = 1, plev-1
         do i = 1, plon
            rpdeli(i,k) = 1. / ( pmid(i,k+1) - pmid(i,k) )
         end do
      end do

      if ( physics ) then
         if ( advqth ) then
            call pt2t( pt0, pmid, t0 )
            call pt2t( ptadv, pmid, tadv )

            if ( minval(tadv(:plon,:plev)).lt.140. ) then
               write(6,*) 'bad tadv 1',minval(tadv(:plon,:plev)),lat
               do k = 1,plev
                  do i = 1,plon
                     if ( tadv(i,k).lt.140. ) write(6,*)' tadv lt 140 at i,k',i,k
                     tadv(i,k)= max(tadv(i,k),150._r8)
                  end do
               end do
               call endrun()
            endif

            call virtem( tadv, shadv, tv )
         else
            call virtem( t1, sh1, tv )
         endif
         call zrel( pint, pmln, tv, zm, zi )
         call zr2a( zm, plev, phis, zma )
         call zr2a( zi, plevp, phis, zia )
         call z2phi( zma, plev, phima )
      end if

      if ( .not. rdomega ) then
         do k = 1, plev
            do i = 1, plon
               omga(i,k) = etadot(i,k)*ps(i)
            end do
         end do
      end if

      if ( dicldphys ) then
c        This is a hack so that the cloud water is not advected.
         do k = 1, plev
            do i = 1, plon
               as1(i,k,pcnst) = as0(i,k,pcnst)
            end do
         end do
      end if

#ifdef DEBUG
         if ( lat.eq.jlook ) write(6,*)'physlic 0.3 ',as1(ilook,klook,pcnst)
     $     ,as0(ilook,klook,pcnst)
#endif 
c     role of fixer
      call qneg3( 'befix', lat, as1 )

      if ( fixmas ) then
c        Apply global mass fixer.
         call advfix( pcnst, dtime, lat, as0, as1 )
      end if

c     Modify negative values created by advection.
      call qneg3( 'afterSPF', lat, as1 )

c     Calculate total advection tendency.
      do m = 1, pcnst
         call calctend( as0(1,1,m), as1(1,1,m), dtime, tend )
         call outfld( tanam(m), tend, plond, lat, obuf )
      end do
      if ( advqth ) then
         call calctend( sh0, shadv, dtime, tend )
         call outfld( 'TAQ', tend, plond, lat, obuf )
         call calctend( pt0, ptadv, dtime, tend )
         call outfld( 'TATHETA', tend, plond, lat, obuf )
      end if

#ifdef DEBUG
         if ( lat.eq.jlook ) write(6,*)'physlic 0.4 ',as1(ilook,klook,pcnst)
#endif 
c     ******************************
c     * Surface fluxes / Emissions *
c     ******************************


      call sflxdr( lat, dtime, calday, ncdate, ioro, as1, pmid, rpdel
     $,            tv, zi, obuf, sflx )


c     **********************
c     * Vertical Diffusion *
c     **********************

      if ( divdiff .or. arvdiff ) then
         call cpslic( pcnst, as1, mrsav )
         if ( advqth ) then
            call cpslic( 1, shadv, shsav )
         end if
      end if

      call vdiffdr( divdiff, arvdiff
     $,             lat, uv, uv(1,1,2), tadv, t1, as1
     $,             pmid, pint, rpdel, rpdeli, dtime
     $,             ptadv, zm, taux, tauy, hflx
     $,             sflx, shadv, wvflx, cgs, kvh
     $,             obuf, ptadv, as1, shadv
     $,             pblh, tpert, qpert, wpert3d
     $,             ustar, obklen
     $)


c++csz
#ifdef DST
c     Driver for aerosol dry deposition
c     dst_dps_dry() is called by CCM:physics/tphysac(), MATCH:src/physlic()
c     call dst_dps_dry(hbuf,obuf,

      call dst_dps_dry(obuf,
     $     zm(1,plev),          ! I [m] Midpoint height above surface
     $     lat,                 ! I [idx] Model latitude index
     $     mno_lng,             ! O [m] Monin-Obukhov length
     $     oro,                 ! I [frc] Orography
     $     pdel,                ! I [Pa] Pressure thickness
     $     pmid,                ! I [Pa] Midlayer pressure
     $     shadv,               ! I [kg kg-1] Water vapor mixing ratio
     $     as1(1,1,dst_idx_srt), ! I/O [kg kg-1] Dust mixing ratio
     $     snowh,               ! I [m] Equivalent liquid water snow depth
     $     dtime,               ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     t1,                  ! I [K] Temperature
     $     pt0(1,plev),         ! I [K] Potential temperature
     $     ts,                  ! I [K] Surface temperature
     $     wnd_frc,             ! O [m s-1] Surface friction velocity
     $     uv(1,plev,2),           ! I [m s-1] Meridional wind component
     $     wnd_rfr,             ! O [m s-1] Wind speed at reference height
     $     uv(1,plev,1))            ! I [m s-1] Zonal wind component


#ifdef DST_DBG
      call qneg3('dst_dry ',lat,as1)
#endif /* not DST_DBG */


#ifdef DST_MSS_BDG
c     Archive mass path of dust after dry deposition
      call dst_add_lon_lev(as1(1,1,dst_idx_srt),q_dst_ttl)
      call bdg_gam_wet('mpc_dst_a_dry',lat,pdel,q_dst_ttl) ! Tracer is moist mass mixing ratio
#endif /* not DST_MSS_BDG */ 


c     Set fields used by mobilization routines
      fsns(1:plon)=fsds(1:plon)-fsus(1:plon)               ! I [W m-2] Solar flux absorbed by ground
      call flx_rad_sfc_set(
     $     lat,                 ! I [idx] Latitude index
     $     plond,               ! I [nbr] Dimension size (normally plon)
     $     flds,                ! I [W m-2] Longwave downwelling flux at surface
     $     fsns)                ! I [W m-2] Solar flux absorbed by ground
      call tpt_gnd_soi_set(
     $     lat,                 ! I [idx] Latitude index
     $     plond,               ! I [nbr] Dimension size (normally plon)
     $     ts,                  ! I [K] Ground temperature
c     fxm: Assume soil temperature is same as surface temperature
     $     ts)                  ! I [K] Soil temperature
      call vwc_sfc_set(
     $     lat,                 ! I [idx] Latitude index
     $     plond,               ! I [nbr] Dimension size (normally plon)
     $     vwcsfc)              ! I [m3 m-3] Volumetric water content

c     Driver for dust mobilization
c     dst_mbl() is called by CCM:physics/tphysbc(), MATCH:src/physlic()
      lat_rdn=clat(lat)         ! [rdn] Latitude

c      call dst_mbl(hbuf,obuf,
      call dst_mbl(obuf,
     $     calday,              ! I [day] Day of year [1.0..366.0)
     $     zm(1,plev),          ! I [m] Midpoint height above surface
     $     lat,                 ! I [idx] Model latitude index
     $     lat_rdn,             ! I [rdn] Latitude
     $     oro,                 ! I [frc] Orography
     $     pdel(1,plev),        ! I [Pa] Pressure thickness
     $     pmid(1,plev),        ! I [Pa] Pressure
     $     pint(1,plevp),       ! I [Pa] Surface pressure NB: plevp
     $     shadv(1,plev),       ! I [kg kg-1] Water vapor mixing ratio
     $     as1(:,plev,dst_idx_srt:dst_idx_end), ! I/O [kg kg-1] Dust mixing ratio
     $     snowh,               ! I [m] Equivalent liquid water snow depth
     $     dtime,               ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     t1(1,plev),          ! I [K] Temperature
     $     pt0(1,plev),         ! I [K] Potential temperature
     $     uv(1,plev,2),        ! I [m s-1] Meridional wind component
     $     uv(1,plev,1))        ! I [m s-1] Zonal wind component

#ifdef DST_DBG
      call qneg3('dst_mbl ',lat,as1)
#endif /* not DST_DBG */


#endif /* not DST */
c--csz

#ifdef DST_MSS_BDG
c     Archive mass path of dust after mobilization
      call dst_add_lon_lev(as1(1,1,dst_idx_srt),q_dst_ttl)
      call bdg_gam_wet('mpc_dst_a_mbl',lat,pdel,q_dst_ttl) ! Tracer is moist mass mixing ratio
#endif /* not DST_MSS_BDG */

c     Calculate vertical diffusion tendency.
      if ( divdiff .or. arvdiff ) then
         do m = 1, pcnst
            call calctend( mrsav(1,1,m), as1(1,1,m), dtime, tend )
            call outfld( vdnam(m), tend, plond, lat, obuf )
         end do
         if ( advqth ) then
            call calctend( shsav, shadv, dtime, tend )
            call outfld( 'VDQ', tend, plond, lat, obuf )
         end if
      end if

c     Update the temperature field to be consistent with the potential
c     temperature field which has just been modified by vertical diffusion.
      if ( advqth ) then
         call pt2t( ptadv, pmid, tadv )

         if ( minval(tadv(:plon,:plev)).lt.140. ) then
            write(6,*) 'bad tadv 2',minval(tadv(:plon,:plev)),lat
            do k = 1,plev
               do i = 1,plon
                  if ( tadv(i,k).lt.140. ) write(6,*)' tadv lt 140 at i,k',i,k
                  tadv(i,k)= max(tadv(i,k),150._r8)
               end do
            end do
               call endrun()
         endif

         if ( addqrad ) then
c           Add radiative heating tendencies to advected/diffused T
            do k = 1, plev
               do i = 1, plon
                  tadv(i,k) = tadv(i,k) + qrad(i,k)*dtime
               end do
            end do
            call t2pt( tadv, pmid, ptadv )
         end if
      end if


c     ***********************************
c     * Diagnose convective mass fluxes *
c     ***********************************

c     CCM Scheme (Zhang/Hack).
      if ( diconvccm ) then

         call cpslic( 1, shadv, shsav )

c        Do not do deep convective transport yet (except for cloud water).
         if ( dicldphys ) call cpslic( 1, as1(1,1,pcnst), cwsav )

#ifdef DEBUG
         if ( lat.eq.jlook ) write(6,*)'physlic 0.5 before aphys',as1(ilook,klook,pcnst)
#endif 

         call aphys( lat, nstep, dtime, pmid, pint
     $,              pdel, rpdel, zm, zi, tpert
     $,              qpert, phis, ts
     $,              pblh, tadv, shadv, as1, cmfmc
     $,              zdu, cmfdqr, conicw, precc, cldtop
     $,              cldbot, dlf, pflx
     $,              omga, wpert3d
     $,              obuf
     $,              zmug, zmdg, zdug, zeug
     $,              zedg, zdpg, dsubcld, zjtg, zjbg
     $,              ideep, lengath )
         call calctend( shsav, shadv, dtime, tend )
         call outfld( 'DCQ', tend, plond, lat, obuf )
         if ( dicldphys ) then
c           Calculate moist convection tendency for cloud water.
            call calctend( cwsav, as1(1,1,pcnst), dtime, tend )
            call outfld( dcnam(pcnst), tend, plond, lat, obuf )
         end if


c        update theta (advected q is updated in aphys)
         call t2pt( tadv, pmid, ptadv )
      end if

c     *****************
c     * Cloud Physics *
c     *****************

      if ( dicldphys ) then
         call clouddr( lat, dtime, oro, pint, pmid, pdel, rpdeli, ts
     $,                tadv, shadv, omga, cmfmc, zdu, cldtop, cldbot
     $,                dlf, t0, sh0, cmfdqr, as1(1,1,pcnst)
     $,                cldt, cldc, cme, evapr, prain, precl
     $,                obuf )
c        update theta
         call t2pt( tadv, pmid, ptadv )


c        Diagnostics
         do i = 1, plon
            prect(i) = precc(i) + precl(i)
         end do
         do k = 1, plev
            do i = 1, plon
               dqr(i,k) = cmfdqr(i,k) + prain(i,k) - evapr(i,k)
            end do
         end do
         call outfld( 'PRECT',  prect, plond, lat, obuf )
         call outfld( 'DQR',      dqr, plond, lat, obuf )
      end if
c     ******************
c     * Dry Deposition *
c     ******************

#ifdef DRYDEP
c     For now substitute ts for t2m.
      if ( dryDep() ) then
         call drydepdr( lat, month, dtime, oro, snowh, prect
     $,              ts, wvflx, hflx, pblh, ustar, obklen
     $,              pmid(1,plev), rpdel, tv(1,plev), zi, uv, uv(1,1,2)
     $,              tadv, fsds, as1, obuf )
         call qneg3( 'DryDep', lat, as1 )
      end if
#endif

c     ******************
c     * Wet Deposition *
c     ******************

c     Default fraction of tracer in interstitial gas is one.  This may be
c     modified for individual constituents by the appropriate call to one
c     of the wet deposition routines.
      do m = 1, pcnst
         do k = 1, plev
            do i = 1, plon
               fracis(i,k,m) = 1.0
            end do
         end do
      end do
#ifdef WETDEP
      if ( dowetdep ) then
         if ( dicldphys ) then
!           calculate cldv and rain
            call clddiag( tadv, pmid, pdel, cmfdqr, cldt, cme, evapr
     $,                   prain, cldv, rain )
            call wetdepdr( lat, dtime, tadv, shadv, pmid, pdel, as1
     $,                    cmfdqr, conicw, as1(:,:,pcnst), cldt, cldc
     $,                    cldv, cme, evapr, prain, rain, fracis
     $,                    obuf )

#ifdef SCYC
c++scyc, chemistry and wetdep routines are called by the sulfur chemistry driver.
            if ( doscyc() ) then
               m = scyc_idx1()
               call chemwdepdr( nstep, lat, clat(lat), calday, dtime
     $,                         pmid, pdel, zm, tadv, shadv
     $,                         as1(:,:,pcnst), cldt, cldc, cldv, cmfdqr
     $,                         prain, evapr, cme, rain, conicw
     $,                         as1(:,:,m:m+3), obuf, fracis(:,:,m:m+3))
            end if
c--scyc
#endif
         else if ( arcldphys ) then
            call clddiag( t1, pmid, pdel, cmfdqr, cldt, cme, evapr
     $,                   prain, cldv, rain )
            call wetdepdr( lat, dtime, t1, sh1, pmid, pdel, as1
     $,                    cmfdqr, conicw, cwat, cldt, cldc
     $,                    cldv, cme, evapr, prain, rain, fracis
     $,                    obuf )
         end if
         call qneg3( 'WetDep', lat, as1 )
      else
      end if
#endif
c++csz
#ifdef DST
c     Cloud diagnostics needed for wet deposition (regardless of dowetdep)
c     Currently cld_dgn() is identical to roundoff to clddiag()
c     cld_dgn() called by MATCH:physlic() and CCM:aphys()
      call cld_dgn( 
     $     cldt,                ! I [frc] Local total cloud fraction
     $     cldv,                ! O [frc] Cloud volume
     $     pdel,                ! I [Pa] Pressure thickness
     $     pmid,                ! I [Pa] Pressure
     $     prain,               ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
     $     cme,                 ! I [kg kg-1 s-1] Net H2O condensate formation tendency
     $     rain,                ! O [kg kg-1] Rainwater mixing ratio
     $     evapr,               ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
     $     cmfdqr,              ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
     $     t1)                  ! I [K] Temperature

c     Wet deposition


c     dst_dps_wet() called by CCM:physics/aphys(), MATCH:physlic()

      if ( doWetDepNew ) then ! use new wet deposition scheme
         call dst_dps_wet(obuf,
     $        cldt,                     ! I [frc] Local total cloud fraction
     $        cldc,                     ! I/O [frc] Convective cloud fraction
     $        cldv,                     ! I [frc] Cloud volume
     $        lat,                      ! I [idx] Model latitude index
     $        pdel,                     ! I [Pa] Pressure thickness
     $        pmid,                     ! I [Pa] Pressure
c    $     cwat,                ! I [kg kg-1] Condensed H2O mixing ratio
     $        as1(1,1,pcnst),           ! I [kg kg-1] Condensed H2O mixing ratio
     $        prain,                    ! I/O [kg kg-1 s-1] Condensed H2O to precipitation tendency
     $        evapr,                    ! I/O [kg kg-1 s-1] H2O precipitation to vapor tendency
     $        cmfdqr,                   ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
     $        as1(1,1,dst_idx_srt),     ! I/O [kg kg-1] Dust mixing ratio
     $        dtime,                    ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $        t1)                       ! I [K] Temperature
      else
         debdst = as1
         fracisdeb(:plon,:,:) = fracis(:plon,:,:)
!         call dst_dps_wet_old(obuf,
!     $        cldt,                     ! I [frc] Local total cloud fraction
!     $        lat,                      ! I [idx] Model latitude index
!     $        pdel,                     ! I [Pa] Pressure thickness
!     $        pmid,                     ! I [Pa] Pressure
!     $        prain,                    ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
!     $        evapr,                    ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
!     $        cmfdqr,                   ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
!     $        as1(1,1,dst_idx_srt),     ! I/O [kg kg-1] Dust mixing ratio
!     $        dtime,                    ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
!     $        t1)                       ! I [K] Temperature
!         call dst_add_lon_lev(as1(1,1,dst_idx_srt),q_dst_ttl)
         sdiag2 = 0.
         mdiag1 = 0.
         do mm = dst_idx_srt,dst_idx_end
            call wetdepa( lat, t1, pmid, unused, pdel, cldt, cldc
     $,                cmfdqr, conicw, prain, unused, evapr
     $,                as1(1,1,pcnst), debdst(1,1,mm), dtime
     $,                debscavt, debiscavt, cldv, fracisdeb(1,1,mm), 0.3_r8  )
            sdiag1 = 0.
            do k = 1,plev
               do i = 1,plon
                  sdiag1(i) = sdiag1(i) + debscavt(i,k)*pdel(i,k)/9.8
                  debdst(i,k,mm) = debdst(i,k,mm) + debscavt(i,k)*dtime
               end do
            end do
            call outfld(flx_mss_pcp_sfc_nm(mm-dst_idx_srt+1), sdiag1, plond, lat, obuf)
            sdiag2 = sdiag2 + sdiag1
            mdiag1 = mdiag1 + debscavt
         end do
         call outfld('DSTSFPCP', sdiag2, plond, lat, obuf)
         call outfld('DSTSSPCP', mdiag1, plond, lat, obuf)

!        overwrite the species array with the new values of dust
         as1(:plon,:plev,dst_idx_srt:dst_idx_end) =  debdst(:plon,:plev,dst_idx_srt:dst_idx_end)
         fracis(:plon,:plev,dst_idx_srt:dst_idx_end) =  fracisdeb(:plon,:plev,dst_idx_srt:dst_idx_end)

      endif                             ! endif do new wet deposition

c
        
#ifdef DST_DBG
      call qneg3('dst_wet ',lat,as1)
#endif /* not DST_DBG */
#endif /* not DST */
c--csz

#ifdef DST_MSS_BDG
c     Archive mass path of dust after wet deposition
      call dst_add_lon_lev(as1(1,1,dst_idx_srt),q_dst_ttl)
      call bdg_gam_wet('mpc_dst_a_wet',lat,pdel,q_dst_ttl) ! Tracer is moist mass mixing ratio
#endif /* not DST_MSS_BDG */

c     ************************************
c     * Convective mass flux adjustments *
c     ************************************

      if ( diconvccm .or. arconvccm ) then
         call cpslic( pcnst, as1, mrsav )
c        CCM Scheme (Zhang/Hack).
         call cmfadj( arconvccm, lat, nstep, dtime, pdel, rpdel
     $,            zmu, zmd, zeu, hketa, hkbeta, as1
     $,            fracis
     $,            zmug, zmdg, zdug, zeug
     $,            zedg, zdpg, dsubcld, zjtg, zjbg
     $,            ideep, lengath )

c        Calculate moist convection tendency.
         mlast = pcnst
         if ( dicldphys ) mlast = mlast - 1
         do m = 1, mlast
            call calctend( mrsav(1,1,m), as1(1,1,m), dtime, tend )
            call outfld( dcnam(m), tend, plond, lat, obuf )
         end do
      end if

c     ******************************************************************
c     * Atmospheric source/sink terms (simple chemistry/decay schemes) *
c     ******************************************************************

      if ( advqth ) then
         call srcsnkdr( lat, dtime, calday, ioro, shadv, as1, obuf )
      else
         call srcsnkdr( lat, dtime, calday, ioro, sh1, as1, obuf )
      end if
#ifdef DEBUG
c         if ( lat.eq.jlook ) write(6,*)'physlic 0.9 ',as1(ilook,klook,pcnst)
#endif 

c     ***************************************************
c     * Output fields corresponding to time index nstep *
c     ***************************************************

c++csz
#ifdef DST
c     Outfield total dust mixing ratio 
c     Called by CCM:dynamics/eul/linemsbc(), MATCH:physlic()
      call dst_add_lon_lev(as1(1,1,dst_idx_srt),q_dst_ttl) ! q3m1
      call outfld('DSTQ    ',q_dst_ttl,plond,lat,obuf)

#ifdef DST_MSS_BDG
c     Archive column mass path of dust at end of timestep 
      call bdg_gam_dry('mpc_dst_mdl',lat,pdel,sh1,q_dst_ttl) ! q3m1
      call bdg_gam_wet('mpc_H2O',lat,pdel,sh1) ! q3m1 NB: Only H2O uses bdg_gam_wet() outside of physlic()
#endif /* not DST_MSS_BDG */

      call mmrd2mpc(sh1,as1(1,1,dst_idx_srt),pdel,mpc_dst,mpc_dst_ttl)
      call outfld('DSTMPC  ',mpc_dst_ttl,plond,lat,obuf)
      do m=1,dst_nbr
         call outfld(mpc_dst_nm(m),mpc_dst(1,m),plond,lat,obuf)
      end do                    ! end loop over cst
c     Diagnostic optical depths
      call mpc2odxc(mpc_dst,odxc_dst,odxc_dst_ttl)
      call outfld('DSTODXC',odxc_dst_ttl,plond,lat,obuf)
      do m=1,dst_nbr
         call outfld(odxc_dst_nm(m),odxc_dst(1,m),plond,lat,obuf)
      end do                    ! end loop over cst
#endif /* not DST */
c--csz

      call outfld( 'PS',    ps,           plond, lat, obuf )
      call outfld( 'U',     uv,           plond, lat, obuf )
      call outfld( 'V', uv(1,1,2),        plond, lat, obuf )
      call outfld( 'QIN',   sh1,          plond, lat, obuf )
      call outfld( 'T',     t1,           plond, lat, obuf )
      call outfld( 'Z',     zm,           plond, lat, obuf )
      call outfld( 'ZI',    zi,           plond, lat, obuf )
      call outfld( 'ORO',   oro,          plond, lat, obuf )
      call outfld( 'PHIS', phis,          plond, lat, obuf )
      call outfld( 'TS',   ts,            plond, lat, obuf )
      call outfld( 'ETADOT', etadot(1,2), plon,  lat, obuf )
      call outfld( 'ETDOTR', wr(1,2),     plon,  lat, obuf )
      call outfld( 'MCA'   , mca(1,2),    plon,  lat, obuf )
      call outfld( 'MCAR'  , mcar(1,2),   plon,  lat, obuf )
      call outfld( 'OMEGA', omga,         plond, lat, obuf )
      call outfld( 'TAUX',  taux,         plond, lat, obuf )
      call outfld( 'TAUY',  tauy,         plond, lat, obuf )
      call outfld( 'SHFLX', hflx,         plond, lat, obuf )
      call outfld( 'QFLX', wvflx,         plond, lat, obuf )
      call outfld( 'FSDS',  fsds,         plond, lat, obuf )
      call outfld( 'FSUS',  fsus,         plond, lat, obuf )
      call outfld( 'FLDS',  flds,         plond, lat, obuf )
      call outfld( 'VWCSFC',vwcsfc,       plond, lat, obuf )
      call outfld( 'SNOWH',snowh,         plond, lat, obuf )


      do m = 1, pcnst
         call outfld( tracnam(m), as1(1,1,m), plond, lat, obuf )
      end do

#ifdef CALC_ADV_FLUX
      call advfluxout( lat, obuf )
#endif

      if ( advqth ) then
         call ciwet( shadv, pdel, ci )
         call outfld( 'TCWV',      ci, plond, lat, obuf )
         call outfld( 'Q',      shadv, plond, lat, obuf )
         call outfld( 'XTHETA', ptadv, plond, lat, obuf )
         call relhum( ps, tadv, shadv, rh )
         call outfld( 'RELHUM',    rh, plond, lat, obuf )
      end if


#ifdef SCYC
c++sslt
c     calculate the 10m wind over the ocean
      call getu10 (ioro, pmid, uv(1,1,1), uv(1,1,2), t1, shadv, zm,
     $     ts, u10)
      call outfld('U10',         u10, plond, lat, obuf )

c     calculate the optical depth of a sea salt aerosol over the ocean
c+drb  changes made to fit WDC 8/31/2K changes
      call getsslt (lat, clat(lat), ioro, u10, zm, zi, rpdel, sslt)
      call outfld('SSLT',       sslt, plond, lat, obuf )

!WDC 8/31/2K -- computed later using optDep routines
!      call outfld('SSLTOD', tausslt, plond, lat, obuf )
c-drb
c--sslt



c++optdep
c     Output non-dust optical depths.
      aerod = 0.0
      totaerod = 0.0
      if ( doscyc() ) then
c        SO4
         mm = scyc_idx1() + 1
         call optdep( tracnam(mm), as1(1,1,mm), pdel, shadv, rh
     $,               aerod(1,mm) )
         call outfld( trim(tracnam(mm))//'OD', aerod(1,mm), plond
     $,               lat, obuf )



#ifdef JUNK
c        if ( lat .le. 3 )  then
         if ( aerod(70,mm).gt.3)  then
            print *,'physlic.F:aerod(lon=70,SO4)',aerod(70,mm), lat
            print *,'as1 ', as1(70,:,mm)
            print *,'rh ', rh(70,:)
            call abort()
         end if
         if (nstep .ge. 43) then
            write (6,*) ' debug quit'
            call abort()
         endif
#endif  /* JUNK */
      end if     !doscyc
#endif  /* SCYC */
#ifdef DST
c     Dust
      do m = 1, dst_nbr
         mm = dst_idx_srt + m - 1
         call optdep( tracnam(mm), as1(1,1,mm), pdel, shadv, rh
     $,               aerod(1,mm) )
      end do
#endif /* not DST */ 
#ifdef CAER
      if ( docaer() ) then
c        Carbon
         do m = 1, caer_ncomp()
            mm = caer_idx1() + m - 1
            call optdep( tracnam(mm), as1(1,1,mm), pdel, shadv, rh
     $,               aerod(1,mm) )
            call outfld( trim(tracnam(mm))//'OD', aerod(1,mm), plond
     $,                  lat, obuf )
         end do
      end if

c-drb: added WDC changes 20020104
!+WDC 8/31/2K

      call optdep( 'SSLT', sslt, pdel, shadv, rh, tausslt )
      call outfld('SSLTOD', tausslt, plond, lat, obuf )


C+WDC 10/26/00 -- store sea-salt optical depth for use in assimilation (called after physlic)
      do i = 1, plon
         tsslt2d(i,lat) = tausslt(i)
      end do

!
!-WDC
!
c-drb: added WDC changes

#ifdef ASSIM
      !output tendencies from assimilation to history file
      call asmfluxout( lat, obuf )
#endif

c     Total aerosol optical depth (including sea salt).
      do m = 1, pcnst
         do i = 1, plon
            totaerod(i) = totaerod(i) + aerod(i,m)
         end do
      end do
      do i = 1, plon
         totaerod(i) = totaerod(i) + tausslt(i)
      end do
      call outfld('AEROD',totaerod, plond, lat, obuf )
c--optdep
#endif

c     ***************************************
c     * Write output buffer to history file *
c     ***************************************
      if ( outputTS() ) then
         istat = qhsti( 'LUNHST', lunhst )
         if ( istat .ne. 0 ) then
            write(*,*)'physlic: qhsti returned error for LUNHST'
            stop
         end if
         obuf(1) = lat
         obuf(2) = plon
         call wrarr( lunhst, maxsiz, obuf )
         deallocate( obuf )
      end if

#ifdef DEBUG
         if ( lat.eq.jlook ) write(6,*)'physlic 1.0 ',as1(ilook,klook,pcnst)
#endif 

      return
      end
      

