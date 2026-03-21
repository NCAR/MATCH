c $Id$
      subroutine aphys(lat     ,nstep   ,tdt     ,pmid    ,pint    ,
     $                 pdel    ,rpdel   ,zm      ,zi      ,tpert   ,
     $                 qpert   ,phis    ,ts      ,
     $                 pblh    ,t       ,q1      ,q       ,cmfmc   ,
     $                 zdu     ,cmfdqr  ,conicw  ,precc   ,cnt     ,
     $                 cnb     ,dlf     ,pflx    ,
     $                 omga    ,wpert3d ,
     $                 obuf    ,
     $               zmug    ,zmdg     ,zdug      ,zeug    ,
     $               zedg    ,zdpg     ,dsubcld   ,zjtg    ,zjbg     ,
     $               ideep   ,lengath )

c Convective adjustments using CCM3 parameterizations (Zhang + Hack).

      use constants, only: rgravit
      use histout, only: outfld
      use conv_zhang, only: conv_ccm,convtran,arconvtran
      use conv_hack, only: cmfmca,hkcmfadj
      use conv_pjr, only: conv_ccm_pjr
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments
      integer
     $     lat,                    ! latitude index (S->N)
     $     nstep                   ! current time step index
      real(r8) tdt,                    ! 2 delta-t (seconds)
     $     pmid(plond,plev),       ! pressure at midpoints
     $     pint(plond,plevp),      ! pressure at interfaces
     $     pdel(plond,plev),       ! delta-p across layers
     $     rpdel(plond,plev),      ! 1./pdel
     $     zm(plond,plev),         ! height above sfc at layer midpoints
     $     zi(plond,plevp)         ! interface height above sfc
      real(r8) tpert(plond),           ! PBL perturbation theta
     $     qpert(plond),           ! PBL perturbation specific humidity 
     $     phis(plond),            ! surface geopotential
     $     ts(plond),              ! surface temperature
     $     pblh(plond)             ! PBL height (provided by PBL routine)
      real(r8) omga(plond,plev)        ! vertical p velocity
      real(r8) wpert3d(plond,plev)     ! perturbation vertical velocities

c Input/output arguments
      real(r8) t(plond,plev),          ! temperature (t bar)
     $     q1(plond,plev),         ! specific humidity (sh bar)
     $     q(plond,plev,pcnst)     ! passive tracers

c Output arguments
      real(r8)
     $     cmfmc(plond,plev),      ! moist convection cloud mass flux
     $     zdu(plond,plev),        ! du2 from conv_ccm, 1/s
     $     cmfdqr(plond,plev),     ! dq/dt due to convective rainout 
     $     conicw(plond,plev),     ! in cloud water mixing ratio (convective)
     $     precc(plond),           ! convective precipitation rate
     $     cnt(plond),             ! top level of convective activity   
     $     cnb(plond),             ! bottom level of convective activity
     $     dlf(plond,plev),        ! detraining cloud water from convection
     $     pflx(plond,plevp)       ! conv rain flux thru btm of that lev
      real(r8)
     $  obuf(*)                    ! output buffer
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

c Local variables:

      integer
     $  i, k

      real(r8)
     $     zmdt(plond,plev),       ! zhang convective temperature tendency
     $     zmdq(plond,plev),       ! zhang convective moisture tendency
     $     cme(plond,plev),        ! zhang condensation - evaporation
     $     zmu(plond,plev),        ! mu2 from conv_ccm, kg/m2/s
     $     zmd(plond,plev),        ! md2 from conv_ccm, kg/m2/s
     $     zeu(plond,plev)         ! eu2 from conv_ccm, 1/s
     $,    fracis(plond,plev)      ! fraction of constituent in interstitial gas

      real(r8)
     $     precc2(plond),           ! convective-scale preciptn rate
     $     cnt2(plond),             ! top level of convective activity
     $     cnb2(plond),             ! bottom level of convective activity
     $     hketa(plond,plev),       ! eta from cmfmca
     $     hkbeta(plond,plev),      ! beta from cmfmca
     $     cmfdt(plond,plev),       ! dt/dt due to moist convection
     $     cmfdq(plond,plev),       ! dq/dt due to moist convection
     $     qc2(plond,plev),         ! dq/dt due to rainout terms
     $     cmfdqr2(plond,plev),     ! dq/dt due to moist convective rainout 
     $     conicw2(plond,plev),     ! in cloud water mixing ratio (convective)
     $     cmfmc2(plond,plev),      ! moist convection cloud mass flux
     $     cmfsl2(plond,plev),      ! moist convection lw stat energy flux
     $     cmflq2(plond,plev),      ! moist convection total water flux
     $     tpert2(plond),           ! perturbation T
     $     qpert2(plond)            ! perturbation q
      real(r8) mup(plond,plevp)
      real(r8) mdn(plond,plevp)
      real(r8) wpls(plond,plev)
      real(r8) wprime3d(plond,plev)     ! perturbation vertical velocities
      real(r8) capblfc(plond)           ! negative cape below level of free conv
      real(r8) pblt(plond)

      real(r8) :: zero = 0. 
      real(r8) :: one = 1. 

      integer ireset(plond)

      logical progCwat, altConv
c-----------------------------------------------------------------------

#undef DEBUG
#ifdef DEBUG
      integer jlook
      ilook = 99
      jlook = 72
      klook = 23
#endif

#ifdef DEBUG
      if ( lat == jlook ) then
         write(6,*)'aphys 0.0,i,k,m,q',ilook,klook,pcnst,q(ilook,klook,pcnst)
      endif
#endif

      


c     Zhang convection.

      call setr( plond*plev,zero, zmdt   )
      call setr( plond*plev,zero, zmdq   )
      call setr( plond*plev,zero, cmfmc  )
      call setr( plond*plev,zero, cme    )
      call setr( plond*plev,zero, zmu    )
      call setr( plond*plev,zero, zmd    )
      call setr( plond*plev,zero, zeu    )
      call setr( plond*plev,zero, cmfdqr )
      call setr( plond*plev, one, fracis )
      call setr( plond,zero, precc )


      call qneg3( 'APHYS_bz', lat, q )

      if ( altConv() ) then
c        add the large scale vert vel (in m/s) to the
c        pert velocities from flow over orography and pbl processes
         do k = 1,plev
            do i = 1,plon
               wpls(i,k) = - omga(i,k)*(zi(i,k)-zi(i,k+1))/pdel(i,k)
               wprime3d(i,k) = wpert3d(i,k) +wpls(i,k)
            end do
         end do
         call conv_ccm_pjr
     $                (t       ,q1      ,precc   ,cnt     ,cnb     ,
     $                 pblh    ,zm      ,phis    ,zi      ,zmdq    ,
     $                 zmdt    ,pmid    ,pint    ,pdel    ,ts      ,
     $                 .5*tdt  ,cmfmc   ,cme     ,nstep   ,lat     ,
     $                 tpert   ,qpert   ,dlf     ,pflx    ,zdu     ,
     $                 cmfdqr  ,mup     ,mdn     ,wprime3d,capblfc ,
     $                 pblt    ,ireset  ,
     $                 zmug    ,zmdg    ,zdug    ,zeug    ,zedg    ,
     $                 zdpg    ,dsubcld ,zjtg    ,zjbg    ,ideep   ,
     $                 lengath ,conicw  )
      else
         call conv_ccm(t       ,q1      ,precc   ,cnt     ,cnb     ,
     $                 pblh    ,zm      ,phis    ,zi      ,zmdq    ,
     $                 zmdt    ,pmid    ,pint    ,pdel    ,ts      ,
     $                 .5*tdt  ,cmfmc   ,cme     ,nstep   ,lat     ,
     $                 tpert   ,qpert   ,dlf     ,pflx    ,zdu     ,
     $                 cmfdqr  ,
     $                 zmug    ,zmdg    ,zdug    ,zeug    ,zedg    ,
     $                 zdpg    ,dsubcld ,zjtg    ,zjbg    ,ideep   ,
     $                 lengath ,conicw  )
      end if

      if ( progCwat() ) then
c        transport cloud water only
         call convtran
     $        (q(1,1,pcnst),1  ,zmug     ,zmdg    ,
     $         zdug    ,zeug    ,zedg    ,zdpg    ,
     $         dsubcld ,zjtg    ,zjbg    ,ideep   ,
     $         1       ,lengath ,nstep   ,lat     ,
     $         .5*tdt  ,fracis  )
      end if


      call qneg3( 'APHYS_az', lat, q )

      call outfld('ZPRECC  ', precc,   plond, lat, obuf )
      call outfld('ZCMFDQR ', cmfdqr,  plond, lat, obuf )
      call outfld('ZCMFMC  ', cmfmc,   plond, lat, obuf )

#ifdef DEBUG
      if ( lat == jlook ) then
         !         write(6,*)'aphys 0.5,i,k,m,q',ilook,klook,pcnst,q(ilook,klook,pcnst)
      endif
#endif

c     Hack convection.

      call setr( plond*plev,zero, cmfdt   )
      call setr( plond*plev,zero, cmfdq   )
      call setr( plond*plev,zero, cmfmc2  )
      call setr( plond*plev,zero, cmfdqr2 )
      call setr( plond*plev,zero, cmfsl2  )
      call setr( plond*plev,zero, cmflq2  )
      do i=1,plon
         tpert2(i) = 0.
         qpert2(i) = qpert(i)
      end do


      call cmfmca(lat     ,nstep   ,tdt     ,pmid    ,pdel    ,
     $            rpdel   ,zm      ,tpert2  ,qpert2  ,phis    ,
     $            pblh    ,t       ,q1      ,cmfdt   ,cmfdq   ,
     $            cmfmc2  ,cmfdqr2 ,cmfsl2  ,cmflq2  ,precc2  ,
     $            qc2     ,cnt2    ,cnb2    ,
     $            conicw2 ,
     $            q, pcnst, hketa, hkbeta   )


      call qneg3( 'APHYS_ah', lat, q )

c     Merge Zhang and Hack diagnostics.

      do i=1,plon
         precc(i) = precc(i) + precc2(i)
         if (cnt2(i) .lt. cnt(i)) cnt(i) = cnt2(i)
         if (cnb2(i) .gt. cnb(i)) cnb(i) = cnb2(i)
      end do
c
      do k=1,plev
         do i=1,plon
            cmfmc(i,k)  = cmfmc(i,k)  + cmfmc2(i,k)
            cmfdqr(i,k) = cmfdqr(i,k) + cmfdqr2(i,k)
            conicw(i,k)  = conicw(i,k) + conicw2(i,k)
         end do 
      end do

      do k = 1, plev
         do i = 1, lengath
c           Convert mass flux from mb/s to kg/m^2/s
            zmu(ideep(i),k) = zmug(i,k) * 100. * rgravit
            zmd(ideep(i),k) = zmdg(i,k) * 100. * rgravit
            zeu(ideep(i),k) = zeug(i,k)
         end do
      end do
      call outfld('ZMU     ', zmu,     plond, lat, obuf )
      call outfld('ZMD     ', zmd,     plond, lat, obuf )
      call outfld('ZEU     ', zeu,     plond, lat, obuf )
      call outfld('ZDU     ', zdu,     plond, lat, obuf )
      call outfld('CONBETA ', hkbeta,  plond, lat, obuf )
      call outfld('CONETA  ', hketa,   plond, lat, obuf )
      call outfld('PRECC   ', precc,   plond, lat, obuf )
      call outfld('CMFDQR  ', cmfdqr,  plond, lat, obuf )
      call outfld('CMFMC   ', cmfmc,   plond, lat, obuf )
      call outfld('CONICW  ', conicw,  plond, lat, obuf )


#ifdef DEBUG
      if ( lat == jlook ) then
         write(6,*)'aphys 1.0,i,k,m,q',ilook,klook,pcnst,q(ilook,klook,pcnst)
      endif
#endif


      return
      end

c###############################################################################

      subroutine cmfadj( arconvccm, lat, nstep, tdt, pdel, rpdel
     $,                  zmu, zmd, zeu, hketa, hkbeta, q
     $,                  fracis
     $,                  zmug, zmdg, zdug, zeug
     $,                  zedg, zdpg, dsubcld, zjtg, zjbg
     $,                  ideep, lengath )

c Convective adjustments using CCM3 parameterizations (Zhang + Hack).

      use conv_zhang
      use conv_hack
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
C
C Input arguments
C
      logical
     $  arconvccm              ! T => use archived mass fluxes
      integer
     $  lat                    ! latitude index (S->N)
     $, nstep                  ! current time step index
      real(r8)
     $  tdt                    ! 2 delta-t (seconds)
     $, pdel(plond,plev)       ! delta-p across layers
     $, rpdel(plond,plev)      ! 1./pdel
     $, zmu(plond,plev)        ! mu2 from conv_ccm, kg/m2/s
     $, zmd(plond,plev)        ! md2 from conv_ccm, kg/m2/s
     $, zeu(plond,plev)        ! eu2 from conv_ccm, 1/s
     $, hketa(plond,plev)      ! eta from cmfmca
     $, hkbeta(plond,plev)     ! beta from cmfmca
     $, fracis(plond,plev,pcnst)  ! fraction of tracer in interstitial gas
C
C Input/output arguments
C
      real(r8)
     $  q(plond,plev,pcnst)    ! passive tracers

c     gathered arrays of driving data for zhang convective transport.
c     When diagnosing mass fluxes these arrays are input.  When using archived 
c     mass fluxes these arrays are calculated using zmu,zmd,zeu.
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

c Local variables:
      integer ncnst
      logical progCwat
c-----------------------------------------------------------------------

c     Zhang convective adjustments.
      if ( arconvccm ) then
         call arconvtran( nstep, lat, pdel, zmu, zmd, zeu
     $,                   zmug, zmdg, zdug, zeug
     $,                   zedg, zdpg, dsubcld, zjtg, zjbg
     $,                   ideep, lengath )
      end if

#ifdef CONVTRAN
c     transport constituents EXCEPT cloud water
      ncnst = pcnst
      if ( progCwat() ) ncnst = pcnst - 1
      call convtran
     $        (q       ,ncnst   ,zmug     ,zmdg    ,
     $         zdug    ,zeug    ,zedg    ,zdpg    ,
     $         dsubcld ,zjtg    ,zjbg    ,ideep   ,
     $         1       ,lengath ,nstep   ,lat     ,
     $         .5*tdt  ,fracis  )

      call qneg3( 'APHYSARz', lat, q )
#endif

c     Archived Hack convection.  (The convective adjustments are done in cmfmca
c     when diagnosing Hack convection.)
      if ( arconvccm ) then
         call hkcmfadj( tdt, rpdel, hketa, hkbeta, q )
         call qneg3( 'APHYSARh', lat, q )
      end if

      return
      end
