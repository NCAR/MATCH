#define LCLJT
#define QPERT
#define RESETTOP
#undef RESETTOP
#undef DEBCONV
      module conv_pjr

      use precision

      real(r8) :: one = 1.
      real(r8) :: zero = 0.

      contains
#include <params.h>
      subroutine buoyan_pjr
     $                 (q       ,t       ,p       ,z       ,pf      ,
     $                  tp      ,qstp    ,tl      ,cape    ,
     $                  pblt    ,lcl     ,lel     ,lfc     ,mx      ,
     $                  nstep   ,
     $                  lat     ,tpert   ,qpert   ,tprime  ,qprime  ,
     $                  capblfc )
C-----------------------------------------------------------------------
C this subroutine works on UNGATHERED ARRAYS
C-----------------------------------------------------------------------
C This is contributed code not fully standardized by the CCM core group.
C
c the documentation has been enhanced to the degree that we are able
c
C Original version:  G. Zhang and collaborators
C Standardized:      Core group staff, 1994 and 195
C Reviewed:          P. Rasch, April 1996
C Revised            P. Rasch and Jon Petch 1997
c-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------------------------------------------------
#include <pmgrid.h>
C--------------------------Commons--------------------------------------
#include <guang_pjr.h>
C-----------------------------------------------------------------------
C
C jul 14/92 - guang jun zhang, m.lazare, n.mcfarlane.  as in
C             previous version buoyan except remove pathalogical
C             cases of "zig-zags" in profiles where lel defined
C             far too high (by use of lelten array, which assumes
C             a maximum of five such crossing points).
C feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.  previous
C             version buoyan.
C
C input arguments
C
      real(r8) q(plond,plev)        ! spec. humidity
      real(r8) t(plond,plev)        ! temperature
      real(r8) p(plond,plev)        ! pressure
      real(r8) z(plond,plev)        ! height
      real(r8) pf(plond,plev+1)     ! pressure at interfaces
      real(r8) pblt(plond)          ! index of pbl depth
      real(r8) tpert(plond)         ! perturbation temperature by pbl processes
      real(r8) qpert(plond)         ! perturbation moisture by pbl processes
c
c output arguments
c
      real(r8) tp(plond,plev)       ! parcel temperature
      real(r8) qstp(plond,plev)     ! saturation mixing ratio of parcel
      real(r8) tl(plond)            ! parcel temperature at lcl
      real(r8) cape(plond)          ! convective aval. pot. energy.
      integer lcl(plond)        ! lifting condensation level
      integer lel(plond)        ! cloud top
c     integer lon(plond)        ! level of onset of deep convection
      integer mx(plond)         ! level of max moist static energy
C
C--------------------------Local Variables------------------------------
C
      real(r8) cape2(plond)          ! available pot. energy of column
      real(r8) capeten(plond,5)     ! provisional value of cape
      real(r8) tv(plond,plev)       ! 
      real(r8) tpv(plond,plev)      ! 
      real(r8) buoy(plond,plev)

      real(r8) a1(plond) 
      real(r8) a2(plond) 
      real(r8) estp(plond) 
      real(r8) pl(plond) 
      real(r8) plexp(plond) 
      real(r8) hmax(plond) 
      real(r8) hmn(plond) 
      real(r8) y(plond)

      logical plge600(plond) 
      integer knt(plond) 
      integer lelten(plond,5)

      real(r8) e

      integer i
      integer k
      integer lat
      integer n
      integer nstep
      integer lnb(plond)


      real(r8) rhd

      real(r8) tprime(plond,plev)
      real(r8) qprime(plond,plev)
      integer lfc(plond)
      real(r8) capblfc(plond)
      real(r8) tparc                ! temperature of the parcel at a given level
      real(r8) qparc                ! spec hum of the parcel at given level
      real(r8) btot(plond)          ! integrated bouyance of parcel at a given level
      real(r8) lel2(plond)          ! second guess on cloud top
      real(r8) off                  ! value of negative bouyance above which the parcel cannot rise
      real(r8) bloc
      real(r8) qs, ha, hb
      real(r8) qs2, esn
      integer numbllcl
      integer numablcl
      integer maxlelten
#include <eslookup.h>
C
C-----------------------------------------------------------------------
      off = -4.                ! value of buoy to switch off at
C-----------------------------------------------------------------------
C

      do n = 1,5
        do i = 1,plon
          lelten(i,n) = plev
          capeten(i,n) = 0.
        end do
      end do
C
      do k = 1,plev
         do i = 1,plon
            buoy(i,k) = 0.0
         end do
      end do
c
      do i = 1,plon
c        lon(i) = plev
        knt(i) = 0
        lel(i) = plev
        mx(i) = plev
        cape(i) = 0.
        hmax(i) = 0.
      end do

c     find the level (mx) of most buoyant parcel in each column, 
c     and the lcl associated with it, and the convective inhibition 
c     of that parcel
      call findcin (lat, t, q, tprime, qprime, mx, 
     $              p, pf, z, lcl, capblfc, lfc, lnb, cape)
#ifdef DEBCONV
      if (lat.eq.latlook) then
         write (6,*) ' at findcin ', lat
         i = ilook
         write (6,*) ' findcin returned mx, lcl, capblfc '//
     $        ' lfc, lnb, cape ',
     $        mx(i), lcl(i), capblfc(i), lfc(i), lnb(i), cape(i)
      endif
#endif
C
c      if (.false.) then
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         write (6,*) ' lat, i, pblt ', lat, i, pblt(i)
         write (6,*) ' buoyan_pjr: t, z, q, qs, rh, h, ha, hpen '
         do k = plev,plev-5,-1
            hmn(i) = cpg*t(i,k) + grav*z(i,k) + rl*q(i,k)
c            e = c1*exp((c2* (t(i,k)-tfreez))/
c     $                ((t(i,k)-tfreez)+c3))
c            qs = eps1*e/ (p(i,k)-e)
c            esn = estblf(t(i,k))
c            qs2 = epsqs*esn/(p(i,k)*100.-esn)
c            ha = cpg*t(i,k) + grav*z(i,k) + rl*(q(i,k)-qs)
c            ha =  hmn(i)+rl*(q(i,k)-qs)
c            hb = cpg*t(i,k-1) + grav*z(i,k-1) + rl*q(i,k-1)
c            write (6,77) k, t(i,k),
c     $           z(i,k), q(i,k),
c     $           qs, q(i,k)/qs, hmn(i), hb,  ha
c   77       format (i5,2f10.3,1p,6e15.5)
c            if (hb .le. hmn(i) .and. q(i,k).ge.qs*0.99) then
c               write (6,*) ' unconditionally unstable profile '
c            endif
         end do
      endif
#endif

      do i = 1,plon
        knt(i) = 0
        lel(i) = plev
        cape(i) = 0.
        hmax(i) = 0.
      end do
C
C
C initialize parcel properties in sub-cloud layer below lcl.
C
      do k = plev,msg + 1,-1
         do i=1,plon
            tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
            if (k.ge.lcl(i) .and. k.le.mx(i)) then
               tparc = t(i,mx(i)) + tprime(i,mx(i))
               qparc = q(i,mx(i)) + qprime(i,mx(i))
               qstp(i,k) = qparc
               tp(i,k)   = tparc*(p(i,k)/p(i,mx(i)))**
     $                     (0.2854*(1.0 - 0.28*qparc))
               tpv(i,k)  = tp(i,k)*(1.+ 1.608*qparc)/
     $                     (1. + qparc)
               buoy(i,k) = tpv(i,k) - tv(i,k)
            else
               tp(i,k) = t(i,k) + tprime(i,k)
               qparc = q(i,k) + qprime(i,k)
               tpv(i,k)  = tp(i,k)*(1.+ 1.608*qparc)/
     $                     (1. + qparc)
               buoy(i,k) = 0.
            end if
        end do
      end do
C     
C     Define parcel properties at lcl (i.e. level immediately above pl).
C
      do i=1,plon
         if (lcl(i).ne.0) then
         k = lcl(i)
         e = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
         tl(i) = 2840./ (3.5*log(t(i,mx(i)))-log(e)-4.805) + 55.
         if (tl(i).lt.t(i,mx(i))) then
            plexp(i) = (1./ (0.2854* (1.-0.28*q(i,mx(i)))))
            pl(i) = p(i,mx(i))* (tl(i)/t(i,mx(i)))**plexp(i)
         else
            tl(i) = t(i,mx(i))
            pl(i) = p(i,mx(i))
         end if
c        tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
         qstp(i,k) = q(i,mx(i)) + qprime(i,mx(i))
C
         tp(i,k) = tl(i)* (p(i,k)/pl(i))**
     $                (0.2854* (1.-0.28*qstp(i,k)))
C
         estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/
     $                ((tp(i,k)-tfreez)+c3))
C
         qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
C
         a1(i) = cpg/rl + qstp(i,k)* (1.+qstp(i,k)/eps1)
     $        *rl*eps1/ (rgas*tp(i,k)**2)
C
         a2(i) = .5* (qstp(i,k)* (1.+2./eps1*qstp(i,k))*
     $        (1.+qstp(i,k)/eps1)*eps1**2*rl*rl/
     $        (rgas**2*tp(i,k)**4)-qstp(i,k)*
     $        (1.+qstp(i,k)/eps1)*2.*eps1*rl/
     $        (rgas*tp(i,k)**3))
C        
         a1(i) = 1./a1(i)
         a2(i) = -a2(i)*a1(i)**3
         y(i)    = q(i,mx(i)) + qprime(i,mx(i)) - qstp(i,k)
         tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
C
c     estp(i)  =exp(a-b/tp(i,k))
C
         estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/
     $        ((tp(i,k)-tfreez)+c3))
C
         qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
C
         tpv(i,k) = (tp(i,k))* (1.+1.608*qstp(i,k))/
     $        (1. + q(i,mx(i)) + qprime(i,mx(i)))
         buoy(i,k) = tpv(i,k) - tv(i,k) 
      endif
      end do
C
C main buoyancy calculation.
C
      do k = plev - 1,msg + 1,-1
        do i=1,plon
          if (k.lt.lcl(i)) then
C
c            tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
            qstp(i,k) = qstp(i,k+1)
            tp(i,k) = tp(i,k+1)* (p(i,k)/p(i,k+1))**
     $                (0.2854* (1.-0.28*qstp(i,k)))
C
c          estp(i) = exp(a-b/tp(i,k))
C
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/
     $                ((tp(i,k)-tfreez)+c3))
C
            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
C
            a1(i) = cpg/rl + qstp(i,k)* (1.+qstp(i,k)/eps1)
     $              *rl*eps1/ (rgas*tp(i,k)**2)
C
            a2(i) = .5* (qstp(i,k)* (1.+2./eps1*qstp(i,k))*
     $              (1.+qstp(i,k)/eps1)*eps1**2*rl*rl/
     $              (rgas**2*tp(i,k)**4)-qstp(i,k)*
     $              (1.+qstp(i,k)/eps1)*2.*eps1*rl/
     $              (rgas*tp(i,k)**3))
C
            a1(i) = 1./a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = qstp(i,k+1) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
c          estp(i)  =exp(a-b/tp(i,k))
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/
     $                ((tp(i,k)-tfreez)+c3))
C
            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
c-jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/
cjt            (1.+q(i,mx(i)))
            tpv(i,k) = tp(i,k)* (1.+1.608*qstp(i,k))/
     $                   (1.+ qstp(i,k))
c     $                   (1.+ q(i,mx(i)) + qprime(i,mx(i))) 
            buoy(i,k) = tpv(i,k) - tv(i,k)
          end if
        end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         do k=plev,1,-1
            if (k.le.mx(i)) then
               write (6,*) 'buoyan_pjr: k, tpv, tv, buoy ', 
     $              k, tpv(i,k), tv(i,k), buoy(i,k)
            endif
         end do
      endif
#endif
C
      maxlelten = 5
      do k = msg + 2,plev
         numbllcl = 0
         do i = 1,plon
            if (k.lt.lcl(i) ) then
               if (buoy(i,k+1).gt.0. .and. buoy(i,k).le.0.) then
                  knt(i) = min(maxlelten,knt(i) + 1)
                  lelten(i,knt(i)) = k
               end if
            else
               numbllcl = numbllcl + 1
            endif
         end do
         if (numbllcl.eq.plon) then
            goto 200
         endif
      end do
      write (6,*) 'buoy_pjr: numbllcl problem ', numbllcl
      stop
  200 continue
C
C    Calculate convective available potential energy (cape).
C
      do n = 1,5
        do k = msg + 1,plev
          do i = 1,plon
            if (k.le.mx(i) .and.k.gt.lelten(i,n)) then
              capeten(i,n) = capeten(i,n) +
     $                       rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
          end do
        end do
      end do
C
C     find maximum cape from all possible tentative capes from
C     one sounding,
C     and use it as the final cape, april 26, 1995
C
      do n = 1,maxlelten
        do i = 1,plon
          if (capeten(i,n).gt.cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
          end if
        end do
      end do

      do i = 1,plon
         capblfc(i) = 0.0
         lfc(i) = plev+2
      end do

c     find the first layer above the lcl with positive buoyancy 
c     as the lfc (level of free convection)
      do k = plev,msg + 1,-1
         do i = 1,plon
            if (k.le.lcl(i)
     $           .and.buoy(i,k).gt.0..and.lfc(i).eq.plev+2) then
               lfc(i) = k
            endif
         end do
      end do

c     now calulate cape between cloud base and the lfc
      do k = plev,msg+1,-1
         do i = 1,plon
            if (k.le.mx(i).and.k.gt.lfc(i)) then
               bloc = rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
c               if (bloc .lt. -100. .or. capblfc(i).lt.-1000.) then
c                  write (6,*) ' strange: i, k, lfc, lcl, bloc, capblfc',
c     $                 i, k, lfc(i), lcl(i), bloc, capblfc(i)
c                  stop
c               endif
               capblfc(i) = capblfc(i) + bloc
            end if
         end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
c         do i = 1,plon
            write (6,*) ' buoyan_pjr: capblfc, mx, lcl, lfc, cape ', 
     $           i, capblfc(i), mx(i), lcl(i), lfc(i), cape(i)
c         end do
      endif
#endif
C
C     put lower bound on cape for diagnostic purposes.
C
      do i = 1,plon
        cape(i) = max(cape(i), zero)
      end do
#ifdef RESETTOP
C
C     limit the cloud top, recalculate cape
C
      do i = 1,plon
         btot(i) = 0.0
         lel2(i) = plev+1
      end do
      do k = plev,1,-1
         do i = 1,plon
            bloc = rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            btot(i) = btot(i) + bloc
            if (k.le.lfc(i).and.btot(i).lt.off
     $           .and.lel2(i).eq.plev+1) then
c              write(*,*) 'lel2 found at',k+1,btot(i)
               lel2(i) = k+1
            endif
         end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         write (6,*) ' lel reset from ', lel(i),' to ', lel2(i)
      endif
#endif
      do i = 1,plon
         if (lel(i).ne.lel2(i)) then
            lel(i) = lel2(i)
         endif
         cape2(i) = 0.0
      end do
      do k = msg + 1,plev
         do i = 1,plon
            if (k.le.mx(i).and.k.gt.lel(i)) then
               cape2(i) = cape2(i) +
     $              rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
         end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         write (6,*) ' new, old cape ', cape2(i), cape(i), 
     $        lel(i), mx(i)
         do k = lel(i), mx(i)
            write (6,*)  'cape by level ', k,
     $              rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
         end do
      endif
#endif
      do i = 1,plon
         cape(i) = max(zero,cape2(i))
      end do

#endif

      return
      end subroutine buoyan_pjr
 
      subroutine cldprp_pjr
     $                 (q       ,t       ,u       ,v       ,p       ,
     $                  z       ,s       ,mu      ,eu      ,du      ,
     $                  md      ,ed      ,sd      ,qd      ,ud      ,
     $                  vd      ,mc      ,qu      ,su      ,zf      ,
     $                  qst     ,hmn     ,hsat    ,alpha   ,shat    ,
     $                  ql      ,totpcp  ,totevp  ,cmeg    ,jb      ,
     $                  lel     ,jt      ,jlcl    ,mx      ,j0      ,
     $                  jd      ,il1g    ,il2g    ,
     $                  nstep   ,lat     ,
     $                  pflx    ,evp     ,cu      ,mu2     ,eu2     ,
     $                  du2     ,md2     ,ed2     ,cmfdqr  ,tpert   ,
     $                  qpert   ,tprime  ,qprime  ,ireset  ,ilookg  )
C------------------------------------------------------------------------------
C This is contributed code not fully standardized by the CCM core group.
C
C this code is very much rougher than virtually anything else in the CCM
C there are debug statements left strewn about and code segments disabled
C these are to facilitate future development. We expect to release a
C cleaner code in a future release
C
c the documentation has been enhanced to the degree that we are able
c
C Original version:  G. Zhang and collaborators
C Standardized:      Core group staff, 1994 and 195
C Reviewed:          P. Rasch, April 1996
C
C**** PLEASE NOTE ****
C
C we are aware of a specific problem in this code 
C (identified by the string ---> PROBLEM ONE)
C during the calculation of the updraft cloud properties,
C rather than adding a perturbation to the updraft temperature of 
C half a degree, (there was an inadvertant addition of cp*0.5) degrees
C or about 500 degrees. (This problem was in the code prior to its 
C contribution to the NCAR effort)

C Fortunately, the erroneous values
C are overwritten later in the code. The problem is quite subtle.
C The erroneous values would persist between cloud base and the lifting 
C condensation level. The addition of the very high perturbation to the updraft
C temperature causes the saturation mixing ratio to be set to zero, 
C and later the lcl to be set to one level above cloud base.
C There are therefore no levels between cloud base and the lcl. Therefore
C all erroneous values are overwritten.

C The only manifestation we are aware of with respect to this problem
C is that the lifting condensation level is constrained to be one level above
C cloud base.

C We discovered the problem after too much had been invested in
C very long integrations (in terms of computer time)
C to allow for a modification and model retuning. It is our expectation that
C this problem will be fixed in the next release of the model.
C
C *********** 
c-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
c-----------------------------------------------------------------------
C nov 20/92 - guang jun zhang,m.lazare. now has deeper (more
C             real(r8)istic) downdrafts.
C jul 14/92 - guang jun zhang,m.lazare. add shallow mixing
C             formulation.
C nov 21/91 - m.lazare. like previous cldprop except minimum "f"
C                       now 0.0004 instead of 0.001 (more
C                       real(r8)istic with more deep).
C may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
C             original version cldprop.
C------------------------------------------------------------------------------
#include <implicit.h>
C------------------------------------------------------------------------------
#include <pmgrid.h>
C------------------------------------------------------------------------------
#include <guang_pjr.h>
C------------------------------------------------------------------------------
C
C Input arguments
C
      real(r8) q(plond,plev)        ! spec. humidity of env
      real(r8) t(plond,plev)        ! temp of env
      real(r8) p(plond,plev)        ! pressure of env
      real(r8) z(plond,plev)        ! height of env
      real(r8) s(plond,plev)        ! normalized dry static energy of env
      real(r8) zf(plond,plevp)      ! height of interfaces
      real(r8) u(plond,plev)        ! zonal velocity of env
      real(r8) v(plond,plev)        ! merid. velocity of env
      real(r8) tpert(plond)         ! perturbation temperature by pbl processes
      real(r8) qpert(plond)         ! perturbation moisture by pbl processes
      real(r8) tprime(plond,plev)   ! level dependent tpert
      real(r8) qprime(plond,plev)   ! level dependent qpert
C
      integer jb(plond)         ! updraft base level
      integer lel(plond)        ! updraft launch level
      integer jt(plond)         ! updraft plume top
      integer jlcl(plond)       ! updraft lifting cond level
      integer mx(plond)         ! updraft base level (same is jb)
      integer j0(plond)         ! level where updraft begins detraining
      integer jd(plond)         ! level of downdraft
c
c output
c
      real(r8) alpha(plond,plev)    !
      real(r8) cmfdqr(plond,plev)   ! rate of production of precip at that layer
      real(r8) du(plond,plev)       ! detrainement rate of updraft
      real(r8) ed(plond,plev)       ! entrainment rate of downdraft
      real(r8) eu(plond,plev)       ! entrainment rate of updraft
      real(r8) hmn(plond,plev)      ! moist stat energy of env
      real(r8) hsat(plond,plev)     ! sat moist stat energy of env
      real(r8) mc(plond,plev)       ! net mass flux
      real(r8) md(plond,plev)       ! downdraft mass flux
      real(r8) mu(plond,plev)       ! updraft mass flux
      real(r8) pflx(plond,plevp)    ! precipitation flux thru layer
      real(r8) qd(plond,plev)       ! spec humidity of downdraft
      real(r8) ql(plond,plev)       ! liq water of updraft
      real(r8) qst(plond,plev)      ! saturation spec humidity of env.
      real(r8) qu(plond,plev)       ! spec hum of updraft
      real(r8) sd(plond,plev)       ! normalized dry stat energy of downdraft
      real(r8) shat(plond,plev)     ! interface values of dry stat energy
      real(r8) su(plond,plev)       ! normalized dry stat energy of updraft
      real(r8) ud(plond,plev)       ! downdraft u
      real(r8) vd(plond,plev)       ! downdraft v
c
c     these version of the mass fluxes conserve mass (used in tracer transport)
c
      real(r8) mu2(plond,plev)      ! updraft mass flux
      real(r8) eu2(plond,plev)      ! updraft entrainment
      real(r8) du2(plond,plev)      ! updraft detrainment
      real(r8) md2(plond,plev)      ! downdraft mass flux
      real(r8) ed2(plond,plev)      ! downdraft entrainment

      integer il1g              !CORE GROUP REMOVE
      integer il2g              !CORE GROUP REMOVE

      integer nstep             ! time step index
      integer lat               ! lat index
C
C Local workspace
C
      real(r8) gamma(plond,plev)  
      real(r8) dz(plond,plev)  
      real(r8) iprm(plond,plev)  
      real(r8) hu(plond,plev)  
      real(r8) hd(plond,plev)  
      real(r8) eps(plond,plev)  
      real(r8) f(plond,plev)  
      real(r8) k1(plond,plev)  
      real(r8) i2(plond,plev)  
      real(r8) ihat(plond,plev)  
      real(r8) i3(plond,plev)  
      real(r8) idag(plond,plev)  
      real(r8) i4(plond,plev)  
      real(r8) qsthat(plond,plev)  
      real(r8) hsthat(plond,plev)  
      real(r8) gamhat(plond,plev)  
      real(r8) cu(plond,plev)  
      real(r8) evp(plond,plev)  
      real(r8) cmeg(plond,plev)  
      real(r8) qds(plond,plev) 
      real(r8) hmin(plond)  
      real(r8) expdif(plond)  
      real(r8) expnum(plond)  
      real(r8) ftemp(plond)  
      real(r8) eps0(plond)  
      real(r8) rmue(plond)  
      real(r8) zuef(plond)  
      real(r8) zdef(plond)  
      real(r8) epsm(plond)  
      real(r8) ratmjb(plond)  
      real(r8) est(plond)  
      real(r8) totpcp(plond)  
      real(r8) totevp(plond)  
      real(r8) alfa(plond) 
      real(r8) beta
      real(r8) c0
      real(r8) ql1
      real(r8) weight
      real(r8) tu
      real(r8) estu
      real(r8) qstu

      real(r8) small
      real(r8) mdt  
      real(r8) cu2

      integer khighest
      integer klowest  
      integer kount 
      integer i,k

      logical doit(plond)
      logical done(plond)
      integer ireset(plond)
      integer ilookg
      real(r8) tmp
C
C------------------------------------------------------------------------------
C
      do i = 1,il2g
         ftemp(i) = 0.
         expnum(i) = 0.
         expdif(i) = 0.
      end do
c
cjr Change from msg+1 to 1 to prevent blowup
c
      do k = 1,plev
        do i = 1,il2g
          dz(i,k) = zf(i,k) - zf(i,k+1)
        end do
      end do

c
c initialize many output and work variables to zero
c
      do k = 1,plev
        do i = 1,il2g
          k1(i,k) = 0.
          i2(i,k) = 0.
          i3(i,k) = 0.
          i4(i,k) = 0.
          mu(i,k) = 0.
          f(i,k) = 0.
          eps(i,k) = 0.
          eu(i,k) = 0.
          du(i,k) = 0.
          ql(i,k) = 0.
          cu(i,k) = 0.
          evp(i,k) = 0.
          cmeg(i,k) = 0.
          qds(i,k) = q(i,k)
          md(i,k) = 0.
          ed(i,k) = 0.
          sd(i,k) = s(i,k)
          qd(i,k) = q(i,k)
          ud(i,k) = u(i,k)
          vd(i,k) = v(i,k)
          mc(i,k) = 0.
          qu(i,k)   = q(i,k) + qprime(i,k)
          su(i,k)   = s(i,k) + tprime(i,k)
c        est(i)=exp(a-b/t(i,k))
          est(i)     = c1*exp((c2*(t(i,k) - tfreez))
     $                   /((t(i,k) - tfreez) + c3))
c++bee
          if ( p(i,k)-est(i) .gt. 0. ) then
             qst(i,k) = eps1*est(i)/ (p(i,k)-est(i))
          else
             qst(i,k) = 1.0
          end if
c--bee
          gamma(i,k) = qst(i,k)*(1. + qst(i,k)/eps1)*eps1*rl/
     $                 (rgas*t(i,k)**2)*rl/cpg
          hmn(i,k) = cpg*t(i,k) + grav*z(i,k) + rl*q(i,k)
          hsat(i,k) = cpg*t(i,k) + grav*z(i,k) + rl*qst(i,k)
c*old          hu(i,k) = hmn(i,k)
          hu(i,k)    = cpg*su(i,k) + rl*qu(i,k)
          hd(i,k) = hmn(i,k)
          mu2(i,k) = 0.
          eu2(i,k) = 0.
          du2(i,k) = 0.
          md2(i,k) = 0.
          ed2(i,k) = 0.
          pflx(i,k) = 0.
          cmfdqr(i,k) = 0.
        end do
      end do
c
c
c interpolate the layer values of qst, hsat and gamma to
c layer interfaces
c
      do i = 1,il2g
        hsthat(i,msg+1) = hsat(i,msg+1)
        qsthat(i,msg+1) = qst(i,msg+1)
        gamhat(i,msg+1) = gamma(i,msg+1)
        totpcp(i) = 0.
        totevp(i) = 0.
        pflx(i,plevp) = 0.
        dz(i,plev) = zf(i,plev) - zf(i,plev+1)
      end do
      do k = msg + 2,plev
        do i = 1,il2g
          if (abs(qst(i,k-1)-qst(i,k)).gt.1.E-6) then
            qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*
     $                    qst(i,k)/ (qst(i,k-1)-qst(i,k))
          else
            qsthat(i,k) = qst(i,k)
          end if
C
          hsthat(i,k) = cpg*shat(i,k) + rl*qsthat(i,k)
C
          if (abs(gamma(i,k-1)-gamma(i,k)).gt.1.E-6) then
            gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*
     $                    gamma(i,k-1)*gamma(i,k)/
     $                    (gamma(i,k-1)-gamma(i,k))
          else
            gamhat(i,k) = gamma(i,k)
          end if
        end do
      end do
C
C     Initialize cloud top to highest plume top.
C
      do i = 1,il2g
        jt(i) = max(lel(i),4)
        jd(i) = plev
        jlcl(i) = lel(i)
        hmin(i) = 1.E6
        ireset(i) = 0
      end do
#ifdef DEBCONV      
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' cldprp: jt is initialized to ', jt(i)
         write (6,*) ' lel is ', lel(i)
      endif
#endif
C
C     Find the level of minimum hsat, where detrainment starts
C
      do k = msg + 1,plev
        do i = 1,il2g
            if (hsat(i,k).le.hmin(i) .and.
     $         (k.ge.jt(i).and.k.le.jb(i))) then
               hmin(i) = hsat(i,k)
               j0(i) = k
            end if
         end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' j0 set 1 ', j0(i)
      endif
#endif
C
      do i = 1,il2g
        j0(i) = min(j0(i),jb(i)-2)
        j0(i) = max(j0(i),jt(i)+2)
c
c Fix from Guang Zhang to address out of bounds array reference
c
        j0(i) = min(j0(i),plev)
      end do
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' j0 set 2 ', j0(i)
      endif
#endif

c no longer need this since perturbations are applied above
C
C Initialize certain arrays inside cloud
C
c      do k = msg + 1,plev
c        do i = 1,il2g
c          if (k.ge.jt(i) .and. k.le.jb(i)) then
c            hu(i,k) = hmn(i,mx(i)) + cpg*0.5
c            su(i,k) = s(i,mx(i)) + 0.5
c          end if
c        end do
c      end do
C
C *********************************************************
C compute taylor series for approximate eps(z) below
C *********************************************************
C
      do k = plev - 1,msg + 1,-1
        do i = 1,il2g
          if (k.lt.jb(i) .and. k.ge.jt(i)) then
           tmp = min(hmn(i,k),hsat(i,k))
            k1(i,k) = k1(i,k+1) 
     $           + (hmn(i,mx(i))-tmp)*dz(i,k)
c     $           + (hmn(i,mx(i))-hmn(i,k))*dz(i,k)
            ihat(i,k) = 0.5* (k1(i,k+1)+k1(i,k))
            i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
            idag(i,k) = 0.5* (i2(i,k+1)+i2(i,k))
            i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
            iprm(i,k) = 0.5* (i3(i,k+1)+i3(i,k))
            i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)
          end if
        end do
      end do
C
C re-initialize hmin array for ensuing calculation.
C
      do i = 1,il2g
        hmin(i) = 1.E6
      end do
C
      do k = msg + 1,plev
        do i = 1,il2g
c          tmp = hmn(i,k)
c          added to deal with supersaturated layers
           tmp = min(hmn(i,k),hsat(i,k))
          if (k.ge.j0(i).and.k.le.jb(i) .and. tmp.le.hmin(i)) then
            hmin(i) = tmp
            expdif(i) = hmn(i,mx(i)) - hmin(i) 
          end if
        end do
      end do
C
C *********************************************************
C     compute approximate eps(z) using above taylor series
C *********************************************************
C
      do k = msg + 2,plev
        do i = 1,il2g
          expnum(i) = 0.
          ftemp(i) = 0.
          if (k.lt.jt(i) .or. k.ge.jb(i)) then
            k1(i,k) = 0.
            expnum(i) = 0.
          else
            expnum(i) = hmn(i,mx(i)) - (hsat(i,k-1)*(zf(i,k)-z(i,k)) +
     $                  hsat(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))
          end if
#ifdef DEBCONV
          if (lat.eq.latlook.and.i.eq.ilookg) then
             write (6,*) ' expnum, k1 ', k, expnum(i), k1(i,k),
     $            expnum(i)*dz(i,k), expdif(i),
     $            jt(i), jb(i)
          endif
#endif
          if ((expdif(i).gt.100..and.expnum(i).gt.0.) .and.
     $        k1(i,k).gt.expnum(i)*dz(i,k)) then
            ftemp(i) = expnum(i)/k1(i,k)
            f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 +
     $                (2.*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2*
     $                ftemp(i)**3 + (-5.*k1(i,k)*i2(i,k)*i3(i,k)+
     $                5.*i2(i,k)**3+k1(i,k)**2*i4(i,k))/
     $                k1(i,k)**3*ftemp(i)**4
            f(i,k) = max(f(i,k),zero)
            f(i,k) = min(f(i,k),0.0002_r8)
          end if
        end do
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' f ', f(i,k), ftemp(i)
      endif
#endif
      end do
      do i = 1,il2g
        if (j0(i).lt.jb(i)) then
          if (f(i,j0(i)).lt.1.E-6 .and. 
     $         f(i,j0(i)+1).gt.f(i,j0(i))) then
             j0(i) = j0(i) + 1
          endif
        end if
      end do
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' j0 reset ', j0(i)
      endif
#endif
C
      do k = msg + 2,plev
        do i = 1,il2g
          if (k.ge.jt(i) .and. k.le.j0(i)) then
            f(i,k) = max(f(i,k),f(i,k-1))
          end if
        end do
      end do
      do i = 1,il2g
        eps0(i) = f(i,j0(i))
        eps(i,jb(i)) = eps0(i)
      end do
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' set eps0 and epsjb ', j0(i), eps0(i),
     $        eps(i,jb(i))
      endif
#endif
C
      do k = plev,msg + 1,-1
         do i = 1,il2g
            if (k.ge.j0(i) .and. k.le.jb(i)) then
               eps(i,k) = f(i,j0(i))
            end if
         end do
      end do
C
      do k = plev,msg + 1,-1
         do i = 1,il2g
            if (k.lt.j0(i) .and. k.ge.jt(i)) eps(i,k) = f(i,k)
         end do
      end do
C
C     specify the updraft mass flux mu, entrainment eu, detrainment du
C     and moist static energy hu.
C     here and below mu, eu,du, md and ed are all normalized by mb
C
      do i = 1,il2g
        if (eps0(i).gt.0.) then
c          mu(i,jb(i)) = 1.
c          eu(i,jb(i)) = eps0(i)/2.
c**pjr NOTE TO CORE GROUP mu and mu2 (and eu and eu2 should have the same vals now)
          mu2(i,jb(i)) = 1.
          eu2(i,jb(i)) = mu2(i,jb(i))/dz(i,jb(i))
          mu(i,jb(i)) = mu2(i,jb(i))
          eu(i,jb(i)) = eu2(i,jb(i))
#ifdef DEBCONV
          if (lat.eq.latlook.and.i.eq.ilookg) then
             write (6,*) ' jb, mujb eps0 ', jb(i), mu(i,jb(i)), eps0(i)
          endif
#endif
        end if
      end do
      do k = plev,msg + 1,-1
        do i = 1,il2g
           if (eps0(i).gt.0. .and. k.ge.jt(i).and.k.lt.jb(i)) then
            zuef(i) = zf(i,k) - zf(i,jb(i))
            rmue(i) = (1./eps0(i))* (exp(eps(i,k+1)*zuef(i))-1.)/zuef(i)
            mu(i,k) = (1./eps0(i))* (exp(eps(i,k)*zuef(i))-1.)/zuef(i)
#ifdef DEBCONV
          if (lat.eq.latlook.and.i.eq.ilookg) then
               write (6,*) ' muk, epsk  ', k, mu(i,k), eps(i,k)
          endif
#endif
            eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
            du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
            mu2(i,k) = mu(i,k)
            eu2(i,k) = eu(i,k)
            du2(i,k) = du(i,k)
          end if
        end do
      end do
C
      khighest = plevp
      klowest = 1
      do i=1,il2g
        khighest = min(khighest,lel(i))
        klowest = max(klowest,jb(i))
      end do
      do k = klowest-1,khighest,-1
cdir$ ivdep
        do i = 1,il2g
          if (k.le.jb(i)-1 .and. k.ge.lel(i) .and. eps0(i).gt.0.) then
#ifdef DEBCONV
             if (lat.eq.latlook.and.i.eq.ilookg) then
                write (6,*) ' small muk, hu, hmn ', k, mu(i,k),
     $               hu(i,k), hmn(i,k), hsat(i,k)
             endif
#endif
            if (mu(i,k).lt.0.01) then
              hu(i,k) = hu(i,jb(i))
              mu(i,k) = 0.
              mu2(i,k) = mu(i,k)
              eu2(i,k) = 0.
              du2(i,k) = mu2(i,k+1)/dz(i,k)
              eu(i,k) = eu2(i,k)
              du(i,k) = du2(i,k)
            else
              hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) +
     $                 dz(i,k)/mu(i,k)* (eu(i,k)*hmn(i,k)-
     $                 du(i,k)*hsat(i,k))
            end if
#ifdef DEBCONV
             if (lat.eq.latlook.and.i.eq.ilookg) then
                write (6,*) ' revised muk, hu, hmn ', k, mu(i,k),
     $               hu(i,k), hmn(i,k)
             endif
#endif
          end if
        end do
      end do
c#ifdef OFFIT
C
C reset cloud top index beginning from two layers above the
C cloud base (i.e. if cloud is only one layer thick, top is not reset
C
c**pjr there are diffs here, and I hope to god they dont matter
      do i=1,il2g
         doit(i) = .true.
         ireset(i) = 0
      end do
      do k=klowest-2,khighest-1,-1
         do i=1,il2g
            if (doit(i) .and. k.le.jb(i)-2 .and. k.ge.lel(i)-1) then
#ifdef DEBCONV
               if (lat.eq.latlook.and.i.eq.ilookg) then
                  write (6,*) ' hub, hstb, hua, hsta, mu doit ', k,
     $                 hu(i,k), hsthat(i,k), hu(i,k+1),
     $                 hsthat(i,k+1), mu(i,k), doit(i)
               endif
#endif
               if (hu(i,k  ).le.hsthat(i,k) .and.
     $              hu(i,k+1).gt.hsthat(i,k+1) .and. 
     $              mu(i,k).ge.0.02) then
                  if (hu(i,k)-hsthat(i,k).lt.-2000.) then
                     jt(i) = k + 1
                     doit(i) = .false.
                     ireset(i) = 1
                  else
                     jt(i) = k
                     ireset(i) = 2
                     if (eps0(i).le.0.) then
                        doit(i) = .false.
                        ireset(i) = 3
                     endif
                  end if
               else if (hu(i,k).gt.hu(i,jb(i)) .or. mu(i,k).lt.0.01
     $             .or. (hu(i,k)-hsthat(i,k).lt.-2000.)) then
                  jt(i) = k + 1
                  doit(i) = .false.
                  ireset(i) = 4
               end if
            end if
#ifdef DEBCONV
            if (lat.eq.latlook.and.i.eq.ilookg) then
               write (6,*) ' jt reset section ', k,
     $              jt(i), ireset(i), doit(i)
            endif
#endif
        end do
      end do
C
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' jt final section ', jt(i), ireset(i)
      endif
#endif
c#endif
C
C
      do k = plev,msg + 1,-1
cdir$ ivdep
        do i = 1,il2g
          if (k.ge.lel(i) .and. k.le.jt(i) .and. 
     $         eps0(i).gt.0.) then
            mu(i,k) = 0.
            eu(i,k) = 0.
            du(i,k) = 0.
            mu2(i,k) = 0.
            eu2(i,k) = 0.
            du2(i,k) = 0.
C     CPGJR/JP               hu(i,k) = hu(i,jb(i)) !commen out?
          end if
          if (k.eq.jt(i) .and. eps0(i).gt.0.) then
            du(i,k) = mu(i,k+1)/dz(i,k)
            du2(i,k) = mu2(i,k+1)/dz(i,k)
            eu2(i,k) = 0.
            mu2(i,k) = 0.
            eu(i,k) = 0.
            mu(i,k) = 0.
          end if
        end do
      end do
C
C specify downdraft properties (no downdrafts if jd.ge.jb).
C scale down downward mass flux profile so that net flux
C (up-down) at cloud base in not negative.
C
      do i = 1,il2g
C
C in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
C
        alfa(i) = dd
        jt(i) = min(jt(i),jb(i)-1)
        jd(i) = max(j0(i),jt(i)+1)
        jd(i) = min(jd(i),jb(i))
        hd(i,jd(i)) = hmn(i,jd(i)-1)
        ud(i,jd(i)) = u(i,jd(i)-1)
        vd(i,jd(i)) = v(i,jd(i)-1)
        if (jd(i).lt.jb(i) .and. eps0(i).gt.0.) then
          epsm(i) = eps0(i)
c          alfa(i)=2.*epsm(i)*( zf(i,jd(i))-zf(i,jb(i)) )/
c     1         (  exp(2.*epsm(i)*( zf(i,jd(i))-
c               zf(i,jb(i)) ))-1.  )
          md(i,jd(i)) = -alfa(i)*epsm(i)/eps0(i)
          md2(i,jd(i)) = md(i,jd(i))
        end if
      end do
      do k = msg + 1,plev
        do i = 1,il2g
          if ((k.gt.jd(i).and.k.le.jb(i)) .and. eps0(i).gt.0.) then
            zdef(i) = zf(i,jd(i)) - zf(i,k)
            md(i,k) = -alfa(i)/ (2.*eps0(i))*
     $               (exp(2.*epsm(i)*zdef(i))-1.)/zdef(i)
            md2(i,k) = md(i,k)
          end if
        end do
      end do
      do k = msg + 1,plev
cdir$ ivdep
        do i = 1,il2g
          if ((k.ge.jt(i).and.k.le.jb(i)) .and. eps0(i).gt.0. .and.
     $         jd(i).lt.jb(i)) then
            ratmjb(i) = min(abs(mu2(i,jb(i))/md2(i,jb(i))),one)
            md2(i,k) = md2(i,k)*ratmjb(i)
c            ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1.)
c            md(i,k) = md(i,k)*ratmjb(i)
               md(i,k) = md2(i,k)
          end if
        end do
      end do

      small = 1.e-20
      do k = msg + 1,plev
        do i = 1,il2g
          if ((k.ge.jt(i).and.k.le.jb(i)) .and. eps0(i).gt.0.) then
            ed2(i,k-1) = (md2(i,k-1)-md2(i,k))/dz(i,k-1)
               ed(i,k-1) = ed2(i,k-1)
            mdt = min(md2(i,k),-small)
            hd(i,k) = (md(i,k-1)*hd(i,k-1) -
     $                 dz(i,k-1)*ed(i,k-1)*hmn(i,k-1))/mdt
          end if
        end do
      end do
C
C calculate updraft and downdraft properties.
C
      do k = msg + 2,plev
        do i = 1,il2g
          if ((k.ge.jd(i).and.k.le.jb(i)) .and. eps0(i).gt.0. .and.
     $        jd(i).lt.jb(i)) then
c         sd(i,k) = shat(i,k)
c    1             +              (hd(i,k)-hsthat(i,k))/
c    2               (cpg    *(1.+gamhat(i,k)))
            qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k))/
     $                 (rl*(1. + gamhat(i,k)))
          end if
        end do
      end do
C
      do i = 1,il2g
         done(i) = .false.
      end do
      kount = 0
      do k = plev,msg + 2,-1
        do i = 1,il2g
          if ((.not.done(i) .and. k.gt.jt(i) .and. k.lt.jb(i)) .and. 
     $         eps0(i).gt.0.) then
            su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) +
     $                dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s(i,k)
            qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) +
     $                dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)-
     $                du(i,k)*qst(i,k))
            tu = su(i,k) - grav/cpg*zf(i,k)
            estu = c1*exp((c2* (tu-tfreez))/ ((tu-tfreez)+c3))
            qstu = eps1*estu/ ((p(i,k)+p(i,k-1))/2.-estu)
            if (qu(i,k).ge.qstu) then
              jlcl(i) = k
              kount = kount + 1
              done(i) = .true.
            end if
          end if
        end do
        if (kount.ge.il2g) goto 690
      end do
 690  continue
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' cldprp: jlcl found to be ', jlcl(i)
         write (6,*) ' jt and jb are currently ', jt(i), jb(i)
      endif
#endif
C
#ifdef LCLJT
C*PJR*XXX    A Change in SCM JP/PJR 
C
      do i = il1g, il2g
         if (jlcl(i).le.jt(i)) then
c            write(*,*) ' ***jlcl =',jlcl(1),'jt =', jt(1),nstep
c            write(*,*) 'Convection should be supressed by phil'
            jt(i)   = jb(i)     !jpnew
            jlcl(i) = jb(i)     !jpnew
            ireset(i) = 5
         endif
      end do
#endif

      do k = msg + 2,plev
        do i = 1,il2g
c          if (k.eq.jb(i) .and. eps0(i).gt.0.) then
c            qu(i,k) = q(i,mx(i))
c            su(i,k) = (hu(i,k)-rl*qu(i,k))/cpg
c          end if
          if ((k.gt.jt(i).and.k.le.jlcl(i)) .and. eps0(i).gt.0.) then
            su(i,k) = shat(i,k) + (hu(i,k)-hsthat(i,k))/
     $               (cpg* (1.+gamhat(i,k)))
            qu(i,k) = qsthat(i,k) + gamhat(i,k)*
     $               (hu(i,k)-hsthat(i,k))/
     $               (rl* (1.+gamhat(i,k)))
          end if
        end do
      end do
C
      do k = plev,msg + 2,-1
        do i = 1,il2g
          if (k.ge.jt(i) .and. k.lt.jb(i) .and. eps0(i).gt.0.) then
            cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/
     $               dz(i,k)- (eu(i,k)-du(i,k))*s(i,k))/
     $               (rl/cpg)
            if (k.eq.jt(i)) cu(i,k) = 0.
c               cu(i,k) = max(zero,cu(i,k))
c               cu2     = max(zero,
c               cu2     = max(-1.e99,
c     $                   +(eu(i,k)*q(i,k) - du(i,k)*qst(i,k))
c     $                   -(mu(i,k)*qu(i,k)-mu(i,k+1)*qu(i,k+1))/dz(i,k)
c     $                   )
c                  
c               if (abs(cu(i,k)-cu2)/(abs(cu(i,k))+abs(cu2)+1.e-50)
c     $              .gt.0.0000001) then
c                  write (6,*) ' inconsistent condensation rates ', 
c     $                 i, k, lat,
c     $                 cu(i,k), cu2, jt(i), jb(i), jlcl(i), lel(i)
c     $                 ,mu(i,k)
c               endif
               cu(i,k) = max(zero,cu(i,k))
          end if
        end do
      end do
C
      beta = 0.
      c0 = 2.E-3
      do k = plev,msg + 2,-1
        do i = 1,il2g
          cmfdqr(i,k) = 0.
C this modification is for test3 run, modified on 6/20/1995
ccc        if(t(i,jt(i) ).gt.tfreez)    c0=0.
ccc        if(t(i,jt(i) ).le.tfreez   )    c0=2.e-3
          if (k.ge.jt(i) .and. k.lt.jb(i) .and. eps0(i).gt.0. .and.
     $         mu(i,k).ge.0.0) then
            if (mu(i,k).gt.0.) then
              ql1 = 1./mu(i,k)* (mu(i,k+1)*ql(i,k+1)-
     $             dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*cu(i,k))
              ql(i,k) = ql1/ (1.+dz(i,k)*c0)
            else
              ql(i,k) = 0.
            end if
            totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*
     $                  (beta*ql(i,k) + (1. - beta)*ql(i,k+1)))
            cmfdqr(i,k) = c0*mu(i,k)*ql(i,k)
          end if
        end do
      end do
C
      do i = 1,il2g
        qd(i,jd(i)) = qds(i,jd(i))
        sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cpg
      end do
C
      do k = msg + 2,plev
        do i = 1,il2g
          if (k.ge.jd(i).and.k.lt.jb(i) .and. eps0(i).gt.0.) then
            qd(i,k+1) = qds(i,k+1)
            evp(i,k) = -ed(i,k)*q(i,k) +
     $                (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
            evp(i,k) = max(evp(i,k),zero)
            mdt = min(md(i,k+1),-small)
            sd(i,k+1) = ((rl/cpg*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) +
     $                    md(i,k)*sd(i,k))/mdt
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
          end if
        end do
      end do
      do i = 1,il2g
c*guang         totevp(i) = totevp(i) + md(i,jd(i))*q(i,jd(i)-1) -
        totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) -
     $              md(i,jb(i))*qd(i,jb(i))
      end do
      if (.true.) then
        do i = 1,il2g
          k = jb(i)
          if (eps0(i).gt.0.) then
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k))/dz(i,k)
            evp(i,k) = max(evp(i,k),zero)
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
          end if
        end do
      endif

      do i = 1,il2g
        totpcp(i) = max(totpcp(i),zero)
        totevp(i) = max(totevp(i),zero)
      end do
C
      weight = 1.0
      do k = msg + 2,plev
        do i = 1,il2g
          if (totevp(i).gt.0. .and. totpcp(i).gt.0.) then
            md2(i,k) = md2(i,k)*min(one,weight*totpcp(i)/
     $                (totevp(i)+weight*totpcp(i)))
            ed2(i,k) = ed2(i,k)*min(one,weight*totpcp(i)/
     $                (totevp(i)+weight*totpcp(i)))
               evp(i,k) = evp(i,k)*min(one,
     $                     weight*totpcp(i)/ (totevp(i)+
     $                     weight*totpcp(i)))
          else
            md2(i,k) = 0.
            ed2(i,k) = 0.
            evp(i,k) = 0.
          end if
            md(i,k) = md2(i,k)
            ed(i,k) = ed2(i,k)
            cmeg(i,k) = cu(i,k) - evp(i,k)
c           cmeg is the cloud water condensed - rain water evaporated
c           cmfdqr  is the cloud water converted to rain - (rain evaporated)
            cmfdqr(i,k) = cmfdqr(i,k)-evp(i,k)
        end do
      end do
      do k = 2,plevp
        do i = 1,il2g
          pflx(i,k) = pflx(i,k-1) + cmfdqr(i,k-1)*dz(i,k-1)
        end do
      end do
      do i = 1,il2g
        if (totevp(i).gt.0. .and. totpcp(i).gt.0.) then
          totevp(i) = totevp(i)*min(one,
     $               weight*totpcp(i)/(totevp(i) + weight*totpcp(i)))
        else
          totevp(i) = 0.
        end if
      end do
C
      do k = msg + 1,plev
        do i = 1,il2g
          mc(i,k) = mu(i,k) + md(i,k)
        end do
      end do
C
      return
      end subroutine cldprp_pjr

      subroutine closure_pjr
     $                  (q       ,t       ,p       ,z       ,s       ,
     $                   tp      ,qs      ,qu      ,su      ,mc      ,
     $                   du      ,mu      ,md      ,qd      ,sd      ,
     $                   alpha   ,qhat    ,shat    ,dp      ,qstp    ,
     $                   zf      ,ql      ,dsubcld ,mb      ,cape    ,
     $                   tl      ,lcl     ,lel     ,jt      ,mx      ,
     $                   il1g    ,il2g    ,
     $                   capelmt ,nstep   ,lat     )
C-----------------------------------------------------------------------
C
C This is contributed code not fully standardized by the CCM core group.
C
C this code is very much rougher than virtually anything else in the CCM
C We expect to release cleaner code in a future release
C
c the documentation has been enhanced to the degree that we are able
c
C Original version:  G. Zhang and collaborators
C Standardized:      Core group staff, 1994 and 195
C Reviewed:          P. Rasch, April 1996
C-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------------------------------------------------
#include <pmgrid.h>
C--------------------------Commons--------------------------------------
#include <guang_pjr.h>
C-----------------------------------------------------------------------
C
C may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
C
C-----------------------------Arguments---------------------------------
C
      real(r8) q(plond,plev)        ! spec humidity
      real(r8) t(plond,plev)        ! temperature
      real(r8) p(plond,plev)        ! pressure (mb)
      real(r8) z(plond,plev)        ! height (m)
      real(r8) s(plond,plev)        ! normalized dry static energy 
      real(r8) tp(plond,plev)       ! parcel temp
      real(r8) qs(plond,plev)       ! sat spec humidity
      real(r8) qu(plond,plev)       ! updraft spec. humidity
      real(r8) su(plond,plev)       ! normalized dry stat energy of updraft
      real(r8) mc(plond,plev)       ! net convective mass flux 
      real(r8) du(plond,plev)       ! detrainment from updraft
      real(r8) mu(plond,plev)       ! mass flux of updraft
      real(r8) md(plond,plev)       ! mass flux of downdraft
      real(r8) qd(plond,plev)       ! spec. humidity of downdraft
      real(r8) sd(plond,plev)       ! dry static energy of downdraft
      real(r8) alpha(plond,plev)
      real(r8) qhat(plond,plev)     ! environment spec humidity at interfaces
      real(r8) shat(plond,plev)     ! env. normalized dry static energy at intrfcs
      real(r8) dp(plond,plev)       ! pressure thickness of layers
      real(r8) qstp(plond,plev)     ! spec humidity of parcel
      real(r8) zf(plond,plev+1)     ! height of interface levels
      real(r8) ql(plond,plev)       ! liquid water mixing ratio

      real(r8) mb(plond)            ! cloud base mass flux
      real(r8) cape(plond)          ! available pot. energy of column
      real(r8) tl(plond)
      real(r8) dsubcld(plond)       ! thickness of subcloud layer

      integer lcl(plond)        ! index of lcl
      integer lel(plond)        ! index of launch leve
      integer jt(plond)         ! top of updraft
      integer mx(plond)         ! base of updraft
C
C--------------------------Local variables------------------------------
C
      real(r8) dtpdt(plond,plev)
      real(r8) dqsdtp(plond,plev)
      real(r8) dtmdt(plond,plev)
      real(r8) dqmdt(plond,plev)
      real(r8) dboydt(plond,plev)
      real(r8) thetavp(plond,plev)
      real(r8) thetavm(plond,plev)

      real(r8) dtbdt(plond),dqbdt(plond),dtldt(plond)
      real(r8) beta
      real(r8) capelmt
      real(r8) dadt
      real(r8) debdt
      real(r8) dltaa
      real(r8) eb

      integer i
      integer il1g
      integer il2g
      integer k
      integer lat
      integer nstep

C-----------------------------------------------------------------------
C change of subcloud layer properties due to convection is
C related to cumulus updrafts and downdrafts.
C mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
C to define betau, betad and f(z).
C note that this implies all time derivatives are in effect
C time derivatives per unit cloud-base mass flux, i.e. they
C have units of 1/mb instead of 1/sec.
C
      do i = il1g,il2g
         mb(i) = 0.
         eb = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
         dtbdt(i) = (1./dsubcld(i))* (mu(i,mx(i))*
     $               (shat(i,mx(i))-su(i,mx(i)))+
     $               md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
         dqbdt(i) = (1./dsubcld(i))* (mu(i,mx(i))*
     $               (qhat(i,mx(i))-qu(i,mx(i)))+
     $               md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
         debdt = eps1*p(i,mx(i))/ (eps1+q(i,mx(i)))**2*dqbdt(i)
         dtldt(i) = -2840.* (3.5/t(i,mx(i))*dtbdt(i)-debdt/eb)/
     $               (3.5*log(t(i,mx(i)))-log(eb)-4.805)**2
      end do
C
C   dtmdt and dqmdt are cumulus heating and drying.
C
      do k = msg + 1,plev
         do i = il1g,il2g
            dtmdt(i,k) = 0.
            dqmdt(i,k) = 0.
         end do
      end do
C
      do k = msg + 1,plev - 1
         do i = il1g,il2g
            if (k.eq.jt(i)) then
               dtmdt(i,k) = (1./dp(i,k))*
     $                       (mu(i,k+1)* (su(i,k+1)-shat(i,k+1)-
     $                       rl/cpg*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-
     $                       shat(i,k+1)))
               dqmdt(i,k) = (1./dp(i,k))*(mu(i,k+1)* (qu(i,k+1)-
     $                       qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*
     $                       (qd(i,k+1)-qhat(i,k+1)))
            end if
         end do
      end do
C
      beta = 0.
      do k = msg + 1,plev - 1
         do i = il1g,il2g
            if (k.gt.jt(i) .and. k.lt.mx(i)) then
               dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+
     $                       mc(i,k+1)* (s(i,k)-shat(i,k+1)))/
     $                       dp(i,k) - rl/cpg*du(i,k)*
     $                       (beta*ql(i,k)+ (1-beta)*ql(i,k+1))
c          dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))
c     1                +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/dp(i,k)
c     2                +du(i,k)*(qs(i,k)-q(i,k))
c     3                +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

               dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+
     $                       cpg/rl* (su(i,k+1)-s(i,k)))-
     $                       mu(i,k)* (qu(i,k)-qhat(i,k)+cpg/rl* 
     $                       (su(i,k)-s(i,k)))+md(i,k+1)* 
     $                       (qd(i,k+1)-qhat(i,k+1)+cpg/rl* 
     $                       (sd(i,k+1)-s(i,k)))-md(i,k)* 
     $                       (qd(i,k)-qhat(i,k)+cpg/rl* 
     $                       (sd(i,k)-s(i,k))))/dp(i,k) +
     $                       du(i,k)* (beta*ql(i,k)+
     $                       (1-beta)*ql(i,k+1))
            end if
         end do
      end do
C
      do k = msg + 1,plev
         do i = il1g,il2g
            if (k.ge.lel(i) .and. k.le.lcl(i)) then
               thetavp(i,k) = tp(i,k)* (1000./p(i,k))** (rgas/cpg)*
     $                         (1.+1.608*qstp(i,k)-q(i,mx(i)))
               thetavm(i,k) = t(i,k)* (1000./p(i,k))** (rgas/cpg)*
     $                         (1.+0.608*q(i,k))
               dqsdtp(i,k) = qstp(i,k)* (1.+qstp(i,k)/eps1)*eps1*rl/
     $                         (rgas*tp(i,k)**2)
C
C dtpdt is the parcel temperature change due to change of
C subcloud layer properties during convection.
C
               dtpdt(i,k) = tp(i,k)/ (1.+
     $                       rl/cpg* (dqsdtp(i,k)-qstp(i,k)/tp(i,k)))*
     $                        (dtbdt(i)/t(i,mx(i))+
     $                       rl/cpg* (dqbdt(i)/tl(i)-q(i,mx(i))/
     $                       tl(i)**2*dtldt(i)))
C
C dboydt is the integrand of cape change.
C
               dboydt(i,k) = ((dtpdt(i,k)/tp(i,k)+1./ 
     $                        (1.+1.608*qstp(i,k)-q(i,mx(i)))* 
     $                        (1.608 * dqsdtp(i,k) * dtpdt(i,k) -
     $                        dqbdt(i))) - (dtmdt(i,k)/t(i,k)+0.608/ 
     $                        (1.+0.608*q(i,k))*dqmdt(i,k)))*grav*
     $                        thetavp(i,k)/thetavm(i,k)
            end if
         end do
      end do
C
      do k = msg + 1,plev
         do i = il1g,il2g
            if (k.gt.lcl(i) .and. k.lt.mx(i)) then
               thetavp(i,k) = tp(i,k)* (1000./p(i,k))** (rgas/cpg)*
     $                         (1.+0.608*q(i,mx(i)))
               thetavm(i,k) = t(i,k)* (1000./p(i,k))** (rgas/cpg)*
     $                         (1.+0.608*q(i,k))
C
C dboydt is the integrand of cape change.
C
               dboydt(i,k) = (dtbdt(i)/t(i,mx(i))+
     $                        0.608/ (1.+0.608*q(i,mx(i)))*dqbdt(i)-
     $                        dtmdt(i,k)/t(i,k)-
     $                        0.608/ (1.+0.608*q(i,k))*dqmdt(i,k))*
     $                        grav*thetavp(i,k)/thetavm(i,k)
            end if
         end do
      end do

C
C buoyant energy change is set to 2/3*excess cape per 3 hours
C
      do i = il1g,il2g
         dadt = 0.
         do k = lel(i),mx(i) - 1
            dadt = dadt + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
         end do
C
         dltaa = -1.* (cape(i)-capelmt)
         if (dadt.ne.0.) mb(i) = max(dltaa/tau/dadt,zero)
      end do
C
      return
      end subroutine closure_pjr
 
      subroutine conv_ccm_pjr
     $                   (t       ,qh      ,pcpc    ,jctop   ,jcbot   ,
     $                    pblh    ,zm      ,geos    ,zi      ,qtg     ,
     $                    ttg     ,pap     ,paph    ,dpp     ,ts      ,
     $                    delt    ,mcon    ,cme     ,nstep   ,lat     ,
     $                    tpert   ,qpert   ,dlf     ,pflx    ,zdu     ,
c++bee
     $                    cmfdqr  ,mup     ,mdn     ,wpert3d ,capblfc ,
     $                    pblt    ,ireset  ,
     $                    mu2     ,md2     ,du2     ,eu2     ,ed2     ,
     $                    dp      ,dsubcld ,jtg     ,maxg    ,ideep   ,
     $                    lengath ,ql      )
c--bee
C
C ##### MAIN DRIVER FOR ZHANG-MCFARLANE CONVECTION SCHEME #####
C
C-----------------------------------------------------------------------
C This is contributed code not fully standardized by the CCM core group.
C All variables have been typed, where most are identified in comments
C The current procedure will be reimplemented in a subsequent version 
C of the CCM where it will include a more straightforward formulation 
C and will make use of the standard CCM nomenclature
C-----------------------------------------------------------------------
#include <implicit.h>
C******************************************************************
#include <pmgrid.h>
C******************************************************************

C same as conv.up except saturation vapor pressure is calculated
C in a different way.
C
C jul 17/92 - guang jun zhang, m.lazare. calls new buoyan, q1q2
C             and moment (several new work fields added for later).
C
C nov 21/91 - m.lazare. like previous conv except calls new
C                       clpdprp.
C feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.
C             previous version conv.
C performs deep convective adjustment based on mass-flux closure
C algorithm.
C
C ************************ index of variables **********************
C  i      => input arrays.
C  i/o    => input/output arrays.
C  w      => work arrays.
C  wg     => work arrays operating only on gathered points.
C  ic     => input data constants.
C  c      => data constants pertaining to subroutine itself.
C
C  wg * alpha    array of vertical differencing used (=1. for upstream).
C  wg * betad    downward mass flux at cloud base.
C  wg * betau    upward   mass flux at cloud base.
C  w  * cape     convective available potential energy.
C  wg * capeg    gathered convective available potential energy.
C  c  * capelmt  threshold value for cape for deep convection.
C  ic  * cpg    specific heat at constant pressure in j/kg-degk.
C  i  * dpp      local sigma half-level thickness (i.e. dshj).
C  ic  * delt     length of model time-step in seconds.
C  wg * dp       layer thickness in mbs (between upper/lower interface).
C  wg * dqdt     mixing ratio tendency at gathered points.
C  wg * dsdt     dry static energy ("temp") tendency at gathered points.
C  wg * dudt     u-wind tendency at gathered points.
C  wg * dvdt     v-wind tendency at gathered points.
C  wg * dsubcld  layer thickness in mbs between lcl and maxi.
C  ic  * grav     acceleration due to gravity in m/sec2.
C  wg * du       detrainment in updraft. specified in mid-layer
C  wg * ed       entrainment in downdraft.
C  wg * eu       entrainment in updraft.
C  wg * hmn      moist static energy.
C  wg * hsat     saturated moist static energy.
C  w  * ideep    holds position of gathered points vs longitude index.
C  ic  * plev     number of model levels.
C  ic  * ilg      lon+2 = size of grid slice.
C  wg * j0       detrainment initiation level index.
C  wg * jd       downdraft   initiation level index.
C  ic  * jlat     gaussian latitude index.
C  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
C  wg * jtg      top  level index of deep cumulus convection.
C  ic  * kount    current model timestep number.
C  w  * lcl      base level index of deep cumulus convection.
C  wg * lclg     gathered values of lcl.
C  w  * lel      index of highest theoretical convective plume.
C  wg * lelg     gathered values of lel.
C  w  * lon      index of onset level for deep convection.
C  wg * long     gathered values of lon.
C  ic  * lev      plev+1.
C  w  * maxi     index of level with largest moist static energy.
C  wg * maxg     gathered values of maxi.
C  wg * mb       cloud base mass flux.
C  wg * mc       net upward (scaled by mb) cloud mass flux.
C  wg * md       downward cloud mass flux (positive up).
C  wg * mu       upward   cloud mass flux (positive up). specified 
C                at interface
C  ic  * msg      number of missing moisture levels at the top of model.
C  c  * kups     number of points undergoing deep convection.
C  w  * p        grid slice of ambient mid-layer pressure in mbs.
C  i  * pblt     row of pbl top indices.
C  i/o * pcp      row of precipitable water in metres.
C  w  * pcpdh    scaled surface pressure.
C  w  * pf       grid slice of ambient interface pressure in mbs.
C  wg * pg       grid slice of gathered values of p.
C  i  * pressg   row of surface pressure in pa.
C  w  * q        grid slice of mixing ratio.
C  wg * qd       grid slice of mixing ratio in downdraft.
C  wg * qdb      row of qd at cloud base.
C  wg * qg       grid slice of gathered values of q.
C  i/o * qh       grid slice of specific humidity.
C  w  * qh0      grid slice of initial specific humidity.
C  wg * qhat     grid slice of upper interface mixing ratio.
C  wg * ql       grid slice of cloud liquid water.
C  wg * qs       grid slice of saturation mixing ratio.
C  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
C  wg * qstpg    grid slice of gathered values of qstp.
C  wg * qug      grid slice of mixing ratio in updraft.
C  ic  * rd      dry air gas constant.
C  wg * rl       latent heat of vaporization.
C  w  * s        grid slice of scaled dry static energy (t+gz/cpg).
C  wg * sd       grid slice of dry static energy in downdraft.
C  wg * sdb      row of sd at cloud base.
C  wg * sg       grid slice of gathered values of s.
C  wg * shat     grid slice of upper interface dry static energy.
C  i  * shbj     grid slice of local bottom interface sigma values.
C  i  * shj      grid slice of local half-level sigma values.
C  i  * shtj     row of local top interfaces of first level.
C  wg * sug      grid slice of dry static energy in updraft.
C  wg * sumde    row of vertically-integrated moist static energy 
C                change.
C  wg * sumdq    row of vertically-integrated scaled mixing ratio 
C                change.
C  wg * sumdt    row of vertically-integrated dry static energy change.
C  wg * sumq     row of vertically-integrated mixing ratio change.
C  i/o * t        grid slice of temperature at mid-layer.
C  o  * jctop    row of top-of-deep-convection indices passed out.
C  o  * jcbot    row of base of cloud indices passed out.
C  w  * tf       grid slice of temperature at interface.
C  wg * tg       grid slice of gathered values of t.
C  w  * tl       row of parcel temperature at lcl.
C  wg * tlg      grid slice of gathered values of tl.
C  w  * tp       grid slice of parcel temperatures.
C  wg * tpg      grid slice of gathered values of tp.
C  i/o * u        grid slice of u-wind (real(r8)).
C  wg * ug       grid slice of gathered values of u.
C  i/o * utg      grid slice of u-wind tendency (real(r8)).
C  i/o * v        grid slice of v-wind (real(r8)).
C  w  * va       work array re-used by called subroutines.
C  wg * vg       grid slice of gathered values of v.
C  i/o * vtg      grid slice of v-wind tendency (real(r8)).
C  i  * w        grid slice of diagnosed large-scale vertical velocity.
C  w  * z        grid slice of ambient mid-layer height in metres.
C  w  * zf       grid slice of ambient interface height in metres.
C  wg * zfg      grid slice of gathered values of zf.
C  wg * zg       grid slice of gathered values of z.
C
C-----------------------------------------------------------------------
C
C multi-level i/o fields:
C
C input/output arguments:
C
      real(r8) t(plond,plev)
c
cjr Added pcnst index to moisture array for consistency
c
!     real(r8) qh(plond,plev,pcnst) 
      real(r8) qh(plond,plev) 
      real(r8) u(plond,plev) 
      real(r8) v(plond,plev) 
c      real(r8) utg(plond,plev) 
c      real(r8) vtg(plond,plev) 
      real(r8) qtg(plond,plev) 
      real(r8) ttg(plond,plev)
C
C       input arguments
      real(r8) pap(plond,plev) 
      real(r8) paph(plond,plev+1) 
      real(r8) dpp(plond,plev) 
      real(r8) zm(plond,plev) 
      real(r8) geos(plond) 
      real(r8) zi(plond,plev+1)
      real(r8) pblh(plond) 
      real(r8) zs(plond) 
      real(r8) tpert(plond) 
      real(r8) qpert(plond)
C
C output arguments
C
      real(r8) pcpck(plond,plev)
      real(r8) mcon(plond,plev) 
      real(r8) dlg(plond,plev)     ! gathrd version of the detraining cld h2o tend
      real(r8) dlf(plond,plev)     ! scattrd version of the detraining cld h2o tend
      real(r8) pflx(plond,plevp)   ! scattered precip flux at each level
      real(r8) pflxg(plond,plevp)  ! gather precip flux at each level
      real(r8) cug(plond,plev)     ! gathered condensation rate 
      real(r8) evpg(plond,plev)    ! gathered evap rate of rain in downdraft
      real(r8) mumax(plond) 
      real(r8) cme(plond,plev)
      real(r8) zdu(plond,plev)
      real(r8) mup(plond,plevp)
      real(r8) mdn(plond,plevp)
      real(r8) wpert3d(plond,plev)
c++bee
      real(r8) cmfdqr(plond,plev)
c--bee
C
C single-level i/o fields:
C
C       input arguments
      real(r8) ts(plond) 
      real(r8) pblt(plond)
C
C       input/output arguments:
C
      real(r8) paprc(plond) 
      real(r8) paprs(plond) 
C
C       output arguments:
C
      real(r8) jctop(plond) 
      real(r8) jcbot(plond)
      real(r8) pcpr(plond) 
      real(r8) pcps(plond) 
      real(r8) pcpc(plond)
C
C-----------------------------------------------------------------------
C
C general work fields (local variables):
C
      real(r8) q(plond,plev) 
      real(r8) p(plond,plev) 
      real(r8) z(plond,plev) 
      real(r8) s(plond,plev) 
      real(r8) qh0(plond,plev) 
      real(r8) tp(plond,plev) 
      real(r8) zf(plond,plev+1) 
      real(r8) pf(plond,plev+1) 
      real(r8) qstp(plond,plev) 

      real(r8) cape(plond) 
      real(r8) tl(plond) 
      real(r8) sumq(plond) 
c      real(r8) sumdt(plond) 
c      real(r8) sumdq(plond) 
c      real(r8) sumde(plond)

      integer lcl(plond) 
      integer lfc(plond)
      integer lel(plond) 
      integer lon(plond) 
      integer maxi(plond) 
      integer ideep(plond) 
      integer index(plond)
      real(r8) precip
C
C gathered work fields:
C
      real(r8) qg(plond,plev) 
      real(r8) tg(plond,plev) 
      real(r8) pg(plond,plev) 
      real(r8) zg(plond,plev) 
      real(r8) sg(plond,plev) 
      real(r8) tpg(plond,plev) 
      real(r8) zfg(plond,plev+1) 
      real(r8) qstpg(plond,plev) 
      real(r8) ug(plond,plev) 
      real(r8) vg(plond,plev) 
      real(r8) cmeg(plond,plev)
      real(r8) cu(plond,plev)
      real(r8) evp(plond,plev)

c++bee
      real(r8) cmfdqrg(plond,plev)
c--bee
      real(r8) capeg(plond) 
      real(r8) tlg(plond)

      integer lclg(plond) 
      integer lelg(plond) 
      integer maxg(plond)
C
C work fields arising from gathered calculations.
C
      real(r8) mu(plond,plev) 
      real(r8) eu(plond,plev) 
      real(r8) dqdt(plond,plev) 
      real(r8) dsdt(plond,plev) 
      real(r8) du(plond,plev) 
      real(r8) md(plond,plev) 
      real(r8) ed(plond,plev) 
      real(r8) alpha(plond,plev) 
      real(r8) sd(plond,plev) 
      real(r8) qd(plond,plev) 
      real(r8) mc(plond,plev) 
      real(r8) qhat(plond,plev) 
      real(r8) qug(plond,plev) 
      real(r8) sug(plond,plev) 
      real(r8) qs(plond,plev) 
      real(r8) shat(plond,plev) 
      real(r8) dp(plond,plev) 
      real(r8) hmn(plond,plev) 
      real(r8) hsat(plond,plev) 
      real(r8) ql(plond,plev) 
c++bee
      real(r8) qlg(plond,plev) 
c--bee
      real(r8) dudt(plond,plev) 
      real(r8) dvdt(plond,plev) 
      real(r8) ud(plond,plev) 
      real(r8) vd(plond,plev)

c      real(r8) deltat(plond,plev) 
c      real(r8) deltaq(plond,plev)

      real(r8) betau(plond) 
      real(r8) betad(plond) 
      real(r8) qdb(plond) 
      real(r8) sdb(plond) 
      real(r8) dsubcld(plond) 
      real(r8) mb(plond) 
      real(r8) totpcp(plond) 
      real(r8) totevp(plond)
      real(r8) mu2(plond,plev) 
      real(r8) eu2(plond,plev) 
      real(r8) du2(plond,plev) 
      real(r8) md2(plond,plev) 
      real(r8) ed2(plond,plev)

      integer jtg(plond) 
      integer jlcl(plond) 
      integer j0(plond) 
      integer jd(plond)

      real(r8) capelmt
      real(r8) delt

      integer i
      integer ii
      integer k
      integer lat
      integer lengath
      integer nstep
      real(r8) psdiss
      real(r8) psevap
      real(r8) psheat
      real(r8) psrain
      real(r8) qdifr
      real(r8) qeff
      real(r8) sdifr
      real(r8) qpmax
      real(r8) tpmax
      real(r8) alpha8
      real(r8) pblhgt
      real(r8) fac1
      real(r8) tprime(plond,plev)
      real(r8) qprime(plond,plev)
      real(r8) tprg(plond,plev)
      real(r8) qprg(plond,plev)
      real(r8) tpertg(plond) 
      real(r8) qpertg(plond)
      real(r8) capblfc(plond)
      real(r8) caped(plond)
      real(r8) wprime
      integer ireset(plond)
      integer irg(plond)
      integer ilookg
C
C--------------------------Commons--------------------------------------
C 
#include <guang_pjr.h>
C#include <comcmf.h>
C
C--------------------------Data statements------------------------------
C
      logical momentm
      data capelmt/70./
      data momentm/.FALSE./
C
C initialize necessary arrays.
C zero out variables not used in ccm
C
      do i = 1,plon
         paprc(i) = 0.
         paprs(i) = 0.
      end do
      psdiss = 0.
      psheat = 0.
      psevap = 0.
      psrain = 0.
C
C initialize convective tendencies
C
      do k = 1,plev
         do i = 1,plon
            dqdt(i,k) = 0.
            dsdt(i,k) = 0.
            dudt(i,k) = 0.
            dvdt(i,k) = 0.
c            deltaq(i,k) = qh(i,k,1)
c            deltat(i,k) = t(i,k)
            mcon(i,k) = 0.
            pcpck(i,k) = 0.
            qtg(i,k) = 0.
            ttg(i,k) = 0.
            pflx(i,k) = 0.
            pflxg(i,k) = 0.
            cme(i,k) = 0.
c++bee
            cmfdqr(i,k) = 0.
            ql(i,k) = 0.
            qlg(i,k) = 0.
c--bee
            zdu(i,k) = 0.
            mup(i,k) = 0.
            mdn(i,k) = 0.
         end do
      end do
      do i = 1,plon
         pcpc(i) = 0.
         pflx(i,plevp) = 0
         pflxg(i,plevp) = 0
      end do
      if (.not.momentm) then
         do k = msg + 1,plev
            do i = 1,plon
               u(i,k) = 0.
               v(i,k) = 0.
            end do
         end do
      end if
C
      do i = 1,plon
         pblt(i) = plev
         pcpr(i) = 0.
         pcps(i) = 0.
         dsubcld(i) = 0.
         sumq(i) = 0.
c         sumdt(i) = 0.
c         sumdq(i) = 0.
         jctop(i) = plev
         jcbot(i) = 1
      end do
C
C calculate local pressure (mbs) and height (m) for both interface
C and mid-layer locations.
C
      do i = 1,plon
         zs(i) = geos(i)*rgrav
         pf(i,plev+1) = paph(i,plev+1)*0.01
         zf(i,plev+1) = zi(i,plev+1) + zs(i)
      end do
      do k = 1,plev
         do i = 1,plon
            p(i,k) = pap(i,k)*0.01
            pf(i,k) = paph(i,k)*0.01
            z(i,k) = zm(i,k) + zs(i)
            zf(i,k) = zi(i,k) + zs(i)
         end do
      end do
C
      do k = plev - 1,msg + 1,-1
         do i = 1,plon
            if (abs(z(i,k)-zs(i)-pblh(i)).lt.
     $          (zf(i,k)-zf(i,k+1))*0.5) pblt(i) = k
         end do
      end do
C
C store incoming specific humidity field for subsequent calculation
C of precipitation (through change in storage).
C convert from specific humidity (bounded by qmin) to mixing ratio.
C define dry static energy (normalized by cpg).
C
      do k = 1,plev
         do i = 1,plon
            qh0(i,k) = qh(i,k)
            qeff = max(qh(i,k),qmin)
            q(i,k) = qeff
            s(i,k) = t(i,k) + (grav/cpg)*z(i,k)
            tp(i,k)=0.0
            shat(i,k) = s(i,k)
            qhat(i,k) = q(i,k)
            dp(i,k) = dpp(i,k)*0.01
            qg(i,k) = q(i,k)
            tg(i,k) = t(i,k)
            pg(i,k) = p(i,k)
            zg(i,k) = z(i,k)
            sg(i,k) = s(i,k)
            tpg(i,k) = tp(i,k)
            zfg(i,k) = zf(i,k)
            qstpg(i,k) = q(i,k)
            ug(i,k) = u(i,k)
            vg(i,k) = v(i,k)
            dlg(i,k) = 0
            dlf(i,k) = 0
            cu(i,k) = 0
            evp(i,k) = 0
         end do
      end do
      do i = 1,plon
         zfg(i,plev+1) = zf(i,plev+1)
         ireset(i) = -1
         capeg(i) = 0.
         lclg(i) = 1
         lelg(i) = plev
         maxg(i) = 1
         tlg(i) = 400.
         dsubcld(i) = 0.
         qdb(i) = 0.
         sdb(i) = 0.
         betau(i) = 0.
         betad(i) = 0.
      end do
C
C     Calculate tprime and qprime as in cmfmc **JP/PJR
C
      tpmax = 1.5
      qpmax = 5.e-4
      qpmax = 1.5e-3
      alpha8 = 2.        ! number of bl depths to propogate perturbations
      alpha8 = 1.        ! number of bl depths to propogate perturbations
C
      do k=1,plev
         do i=1,plon
            pblhgt = max(pblh(i),one)
            fac1   = max(zero,1.0-zm(i,k)*rgrav/(alpha8*pblhgt))
            tprime(i,k) = min(tpert(i),tpmax)*fac1
#ifdef QPERT
            qprime(i,k) = min(qpert(i),qpmax)*fac1
#else
            qprime(i,k) = 0.
#endif
C           write(*,*) tprime(i,k),qprime(i,k),k
         enddo
      enddo
C
      
C
C evaluate covective available potential energy (cape).
C
      call buoyan_pjr
     $           (q       ,t       ,p       ,z       ,pf      ,
     $            tp      ,qstp    ,tl      ,cape    ,
     $            pblt    ,lcl     ,lel     ,lfc     ,maxi    ,
     $            nstep   ,
     $            lat     ,tpert   ,qpert   ,tprime  ,qprime  ,
     $            capblfc )

      if (trigon) then
c++pjr
      do i = 1,plon
         wprime = wpert3d(i,maxi(i))
         caped(i) = cape(i)
         if (max(zero,wprime)**2 + 2*capblfc(i).lt.0.) then
            cape(i) = -.0001
         endif
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         wprime = wpert3d(i,maxi(i))
         write (6,*) 'conv_ccm_pjr: maxi, wprime, capblfc ',
     $        maxi(i), wprime, capblfc(i),
     $        max(zero,wprime)**2 + 2*capblfc(i), wpert3d(i,plev)
         write (6,*) ' caped, cape, lcl, lel, lfc, maxi ',
     $        caped(i), cape(i), lcl(i), lel(i), lfc(i), maxi(i)
      endif
#endif
c--pjr
      endif
C
C determine whether grid points will undergo some deep convection
C (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
C (require cape.gt. 0 and lel<lcl as minimum conditions).
C
         call whenfgt(plon,cape,1,capelmt,ideep,lengath)
         if (lengath.eq.0) return
         
#ifdef DEBCONV
         if (lat.eq.latlook) then
            ilookg = 0
            do i = 1,lengath
               if (ideep(i) .eq. ilook) then
                  ilookg = i
               endif
            end do
            write (6,*) ' ilook ', ilook, 
     $           ' corresponds to gather index', ilookg
         endif
#endif
cc        jyes = 0
cc        jno = plon - 1 + 2
cc        do il = 1,plon
cc           if (cape(il).gt.capelmt) then
cc              jyes = jyes + 1
cc              ideep(jyes) = il
cc           else
cc              jno = jno - 1
cc              ideep(jno) = il
cc           end if
cc        end do
cc        lengath = jyes
cc        if (lengath.eq.0) return
C
C obtain gathered arrays necessary for ensuing calculations.
C
      do k = 1,plev
         do i = 1,lengath
            dp(i,k) = 0.01*dpp(ideep(i),k)
            qg(i,k) = q(ideep(i),k)
            tg(i,k) = t(ideep(i),k)
            pg(i,k) = p(ideep(i),k)
            zg(i,k) = z(ideep(i),k)
            sg(i,k) = s(ideep(i),k)
            tpg(i,k) = tp(ideep(i),k)
            zfg(i,k) = zf(ideep(i),k)
            qstpg(i,k) = qstp(ideep(i),k)
            ug(i,k) = u(ideep(i),k)
            vg(i,k) = v(ideep(i),k)
            tprg(i,k) = tprime(ideep(i),k)
            qprg(i,k) = qprime(ideep(i),k)
         end do
      end do
C
      do i = 1,lengath
         zfg(i,plev+1) = zf(ideep(i),plev+1)
      end do
      do i = 1,lengath
         capeg(i) = cape(ideep(i))
         lclg(i) = lcl(ideep(i))
         lelg(i) = lel(ideep(i))
         maxg(i) = maxi(ideep(i))
         tlg(i) = tl(ideep(i))
         tpertg(i) = tpert(ideep(i))
         qpertg(i) = qpert(ideep(i))
      end do
C
C calculate sub-cloud layer pressure "thickness" for use in
C closure and tendency routines.
C
      do k = msg + 1,plev
         do i = 1,lengath
            if (k.ge.maxg(i)) then
               dsubcld(i) = dsubcld(i) + dp(i,k)
            end if
         end do
      end do
C
C define array of factors (alpha) which defines interfacial
C values, as well as interfacial values for (q,s) used in
C subsequent routines.
C
      do k = msg + 2,plev
         do i = 1,lengath
            alpha(i,k) = 0.5
            sdifr = 0.
            qdifr = 0.
            if (sg(i,k).gt.0. .or. sg(i,k-1).gt. 0.) 
     $          sdifr = abs((sg(i,k)-sg(i,k-1))/
     $          max(sg(i,k-1),sg(i,k)))
            if (qg(i,k).gt.0. .or. qg(i,k-1).gt.0.) 
     $          qdifr = abs((qg(i,k)-qg(i,k-1))/
     $          max(qg(i,k-1),qg(i,k)))
            if (sdifr.gt.1.E-6) then
               shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/
     $                     (sg(i,k-1)-sg(i,k))
            else
               shat(i,k) = 0.5* (sg(i,k)+sg(i,k-1))
            end if
            if (qdifr.gt.1.E-6) then
               qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/
     $                     (qg(i,k-1)-qg(i,k))
            else
               qhat(i,k) = 0.5* (qg(i,k)+qg(i,k-1))
            end if
         end do
      end do
C
C obtain cloud properties.
C
      call cldprp_pjr
     $           (qg      ,tg      ,ug      ,vg      ,pg      ,
     $            zg      ,sg      ,mu      ,eu      ,du      ,
     $            md      ,ed      ,sd      ,qd      ,ud      ,
     $            vd      ,mc      ,qug     ,sug     ,zfg     ,
     $            qs      ,hmn     ,hsat    ,alpha   ,shat    ,
     $            qlg     ,totpcp  ,totevp  ,cmeg    ,maxg    ,
     $            lelg    ,jtg     ,jlcl    ,maxg    ,j0      ,
     $            jd      ,1       ,lengath ,
     $            nstep   ,lat     ,
     $            pflxg   ,evpg    ,cug     ,mu2     ,eu2     ,
     $            du2     ,md2     ,ed2     ,cmfdqrg ,tpertg  ,
     $            qpertg  ,tprg     ,qprg   ,irg     ,ilookg  )
C
C determine cloud base mass flux.
C
      do i = 1,lengath
         qdb(i) = qd(i,maxg(i))
         sdb(i) = sd(i,maxg(i))
         betad(i) = md(i,maxg(i))
         betau(i) = mu(i,maxg(i))
      end do
C
C convert detrainment from units of "1/m" to "1/mb".
C
      do k = msg + 1,plev
         do i = 1,lengath
            du(i,k) = du(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            eu(i,k) = eu(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            ed(i,k) = ed(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            cug(i,k) = cug(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
c++bee, bug fix from /home/pjr/ccm2/omega0.10.1/fspj01/.
            cmeg(i,k) = cmeg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            cmfdqrg(i,k) = cmfdqrg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
c--bee
            evpg(i,k) = evpg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            du2(i,k) = du2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            eu2(i,k) = eu2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            ed2(i,k) = ed2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         end do
      end do

      call closure_pjr
     $            (qg      ,tg      ,pg      ,zg      ,sg      ,
     $             tpg     ,qs      ,qug     ,sug     ,mc      ,
     $             du      ,mu      ,md      ,qd      ,sd      ,
     $             alpha   ,qhat    ,shat    ,dp      ,qstpg   ,
     $             zfg     ,qlg     ,dsubcld ,mb      ,capeg   ,
     $             tlg     ,lclg    ,lelg    ,jtg     ,maxg    ,
     $             1       ,lengath ,
     $             capelmt ,nstep   ,lat     )
C
C limit cloud base mass flux to theoretical upper bound.
C
      do i=1,lengath
        mumax(i) = 0
      end do
      do k=msg + 2,plev
        do i=1,lengath
          mumax(i) = max(mumax(i), mu(i,k)/dp(i,k))
        end do
      end do
      do i=1,lengath
        if (mumax(i).gt.0.) then
          mb(i) = min(mb(i),0.5/(delt*mumax(i)))
        else
          mb(i) = 0.
        endif
      end do
      do k=msg+1,plev
        do i=1,lengath
          mu(i,k) = mu(i,k)*mb(i)
          md(i,k) = md(i,k)*mb(i)
          mc(i,k) = mc(i,k)*mb(i)
          du(i,k) = du(i,k)*mb(i)
          eu(i,k) = eu(i,k)*mb(i)
          ed(i,k) = ed(i,k)*mb(i)
          cmeg(i,k) = cmeg(i,k)*mb(i)
c++bee
          cmfdqrg(i,k) = cmfdqrg(i,k)*mb(i)
c--bee
          cug(i,k) = cug(i,k)*mb(i)
          evpg(i,k) = evpg(i,k)*mb(i)
          pflxg(i,k+1) = pflxg(i,k+1)*mb(i)*100./grav
          mu2(i,k) = mu2(i,k)*mb(i)
          md2(i,k) = md2(i,k)*mb(i)
          du2(i,k) = du2(i,k)*mb(i)
          eu2(i,k) = eu2(i,k)*mb(i)
          ed2(i,k) = ed2(i,k)*mb(i)
        end do
      end do
      do i = 1,lengath
         betau(i) = betau(i)*mb(i)
         betad(i) = betad(i)*mb(i)
C
C totpcp from cldprp has the dimension of kg/kg, here it is 
C converted to kg/(m^2*s), the precipitation rate
C
         totpcp(i) = totpcp(i)*mb(i)*100./grav
         totevp(i) = totevp(i)*mb(i)*100./grav
      end do
C
C compute temperature and moisture changes due to convection.
C
      call q1q2pjr2
     $             (dqdt    ,dsdt    ,qg      ,sg      ,qs      ,
     $              qug     ,sug     ,mc      ,du      ,alpha   ,
     $              qhat    ,shat    ,dp      ,mu      ,md      ,
     $              sd      ,qd      ,qlg     ,dsubcld ,qdb     ,
     $              sdb     ,betau   ,betad   ,mb      ,lclg    ,
     $              jtg     ,maxg    ,2.*delt ,1       ,lengath ,
     $              nstep   ,lat     ,
     $              dlg     ,evpg    ,cug     ,eu      ,ed      )
C
C compute momentum changes due to convection, if desired (i.e
C if logical switch set).
C
c
c      if(momentm)                                                   
c       then
c        call moment(dudt,dvdt,du,alpha,dp,ed,eu,mc,md,mu,
c     1             pg,qd,qu,qhat,sd,su,shat,ud,vd,tg,ug,vg,zg,zfg,
c     2             dsubcld,maxg,jd,jt,
c     3             msg,2.*delt,plev,1,lengath,plond,lat)
c      endif
c
c++bee, move convective transports outside conv_ccm
c      if (pcnst.gt.1) then
c         call convtran_pjr
c     $                (qh      ,mu2     ,md2     ,du2     ,eu2     ,
c     $                 ed2     ,dp      ,dsubcld ,jtg     ,maxg    ,
c     $                 ideep   ,1       ,lengath ,nstep   ,lat     ,
c     $                 delt    )
c      endif
c--bee
C
C gather back temperature and mixing ratio.
C
      do k = msg + 1,plev
         do i = 1,lengath
            psdiss = psdiss + (dudt(i,k)*u(ideep(i),k)+
     $               dvdt(i,k)*v(ideep(i),k))*dpp(ideep(i),k)/grav
C
C q is updated to compute net precip, and then reset to old value.
C the last line is overwritten. so the input basic variables, i.e.
C q, t, u and v are updated to include convective increments. 
C (5/2/95)
C
            q(ideep(i),k) = q(ideep(i),k) + 2.*delt*dqdt(i,k)
            t(ideep(i),k) = t(ideep(i),k) + 2.*delt*dsdt(i,k)
            u(ideep(i),k) = u(ideep(i),k) + 2.*delt*dudt(i,k)
            v(ideep(i),k) = v(ideep(i),k) + 2.*delt*dvdt(i,k)
            cme(ideep(i),k) = cmeg(i,k)
c++bee
            cmfdqr(ideep(i),k) = cmfdqrg(i,k)
            ql(ideep(i),k) = qlg(i,k)
c--bee
            zdu(ideep(i),k) = du2(i,k)
            mup(ideep(i),k) = mu(i,k)
            mdn(ideep(i),k) = md(i,k)
            mcon(ideep(i),k) = mc(i,k)
            qtg(ideep(i),k) = dqdt(i,k)
            ttg(ideep(i),k) = dsdt(i,k)
c            utg(ideep(i),k) = dudt(i,k)
c            vtg(ideep(i),k) = dvdt(i,k)
            dlf(ideep(i),k) = dlg(i,k)
            pflx(ideep(i),k) = pflxg(i,k)
            cu(ideep(i),k) = cug(i,k)
            evp(ideep(i),k) = evpg(i,k)
            ireset(ideep(i)) = irg(i)
         end do
      end do
C
      do i = 1,lengath
         jctop(ideep(i)) = jtg(i)
c++bee
         jcbot(ideep(i)) = maxg(i)
c--bee
         pflx(ideep(i),plevp) = pflxg(i,plevp)
         psevap = psevap + totevp(i)
         psrain = psrain + totpcp(i)
      end do
C
C convert back to specific humidity from mixing ratio.
C take into account any moisture added to ensure positiveness
C of specific humidity at start of routine.
C
      do k = msg + 1,plev
         do i = 1,plon
            qh(i,k) = q(i,k)
            qh(i,k) = qh(i,k) - max((qmin-qh0(i,k)),zero)
         end do
      end do
      do k = plev,msg + 1,-1
         do i = 1,plon
            sumq(i) = sumq(i) - dpp(i,k)* (qh(i,k)-qh0(i,k))
c
c account for the detraining cloud water in the precip 
c
            sumq(i) = sumq(i) - dpp(i,k)*dlf(i,k)*2*delt
            pcpck(i,k) = max(zero,sumq(i))
         end do
      end do
C
C obtain final precipitation rate.
C
      do i = 1,plon
c         llo1 = ts(i) .ge. tfreez
C
C here pcpr and pcps are in units of kg/m^2, ie. precip per
C time step
C
c         pcpr(i) = cvmgt(pcpdh(i)*max(sumq(i),zero),0.,llo1)
c         pcps(i) = cvmgt(0.,pcpdh(i)*max(sumq(i),zero),llo1)
         precip = rgrav*max(sumq(i),zero)
         if (ts(i) .ge. tfreez) then
           pcpr(i) = precip
           pcps(i) = 0.
         else
           pcpr(i) = 0.
           pcps(i) = precip
         end if
      end do

C
C accumulate precipitation, the 1000. is the density of water, so
C paprc and paprs are now in units of meters.
C
      do i = 1,plon
         paprc(i) = paprc(i) + (pcpr(i)+pcps(i))/1000.
         paprs(i) = paprs(i) + pcps(i)/1000.
      end do
C
C convert precipitation to m/s, ie, precip rate.
C
      do i = 1,plon
         pcpr(i) = pcpr(i)/ (2.*delt)/1000.
         pcps(i) = pcps(i)/ (2.*delt)/1000.
         pcpc(i) = pcpr(i) + pcps(i)
         psheat = psheat + (pcps(i)+pcpr(i))*rl
      end do
      do k = msg + 1,plev
         do i = 1,plon
            pcpck(i,k) = rgrav*pcpck(i,k)/ (2.*delt)
         end do
      end do
C
C calculate conservation of quantities.
C

C       if(lat.eq.jlatpr)then
c        do l=msg+1,plev
c        do i=1,lengath
c          sumdq(i) = sumdq(i) + 2.*delt*(rl/cp)*dpp(ideep(i),l)*
c     1                          dqdt(i,l)
c          sumdt(i) = sumdt(i) + 2.*delt*dpp(ideep(i),l)*dsdt(i,l)
c        end do
c        end do
c
c        print *,'sumdq,sumdt,sumde in convection subroutine########'
c        do i=1,lengath
c          sumde(i) = sumdt(i) + sumdq(i)
c          write(6, 901) sumdq(i), sumdt(i),sumde(i), i, ideep(i)
c        end do
cc
c        print *,'sumdq,sumdt,sumde ... all points'
c      do i=1,plon
c         sumdq(i) = 0.0
c         sumdt(i) = 0.0
c      end do
cc
c      do l=msg+1,plev
c      do i=1,plon
c        deltaq(i,l) = qh(i,l) - deltaq(i,l)
c        deltat(i,l) = t (i,l) - deltat(i,l)
c        sumdq(i) = sumdq(i) + (rl/cp)*dpp(i,l)*deltaq(i,l)
c        sumdt(i) = sumdt(i) + dpp(i,l)*deltat(i,l)
c      end do
c      end do
c      do i=1,plon
c        sumde(i) = sumdt(i) + sumdq(i)
c      end do
c        write(6, 902) (i,sumdq(i),sumdt(i),sumde(i),i=1,plon)
c
c  901   format(1x,3e20.12, i10, i10)
c  902   format(1x,i10, 3e20.12)
c
c      endif

      return
      end subroutine conv_ccm_pjr

      subroutine convtran_pjr
     $                   (q     ,mu     ,md      ,du     ,eu     ,
     $                    ed    ,dp     ,dsubcld ,jt     ,mx     ,
     $                    ideep ,il1g   ,il2g    ,nstep  ,lat    ,
     $                    delt    )
C-----------------------------------------------------------------------
C
C Convective transport of trace species
C
C Note that we are still assuming that the tracers are in a moist mixing ratio
C this will change soon
C
C-------------------------Code History----------------------------------
C
C Original version:  P. Rasch, Jan 1996 
C Standardized:      L. Buja,  Feb 1996
C Reviewed:          P. Rasch, Feb 1996      
C 
C-----------------------------------------------------------------------
      implicit none
C--------------------------Commons--------------------------------------
#include <pmgrid.h>
C-----------------------------Arguments---------------------------------
C 
C Input
C
      real(r8) mu(plond,plev)       ! Mass flux up
      real(r8) md(plond,plev)       ! Mass flux down
      real(r8) du(plond,plev)       ! Mass detraining from updraft
      real(r8) eu(plond,plev)       ! Mass entraining from updraft
      real(r8) ed(plond,plev)       ! Mass entraining from downdraft
      real(r8) dp(plond,plev)       ! Delta pressure between interfaces
      real(r8) dsubcld(plond)       ! Delta pressure from cloud base to sfc

      integer jt(plond)         ! Index of cloud top for each column
      integer mx(plond)         ! Index of cloud top for each column
      integer ideep(plond)      ! Gathering array
      integer il1g              ! Gathered min lon indices over which to operate
      integer il2g              ! Gathered max lon indices over which to operate
      integer lat               ! Latitude
      integer nstep             ! Time step index

      real(r8) delt                 ! Time step

C input/output

      real(r8) q(plond,plev,pcnst)  ! Tracer array including moisture

C--------------------------Local Variables------------------------------

      integer i                 ! Work index
      integer k                 ! Work index
      integer kbm               ! Highest altitude index of cloud base
      integer kk                ! Work index
      integer kkp1              ! Work index
      integer km1               ! Work index
      integer kp1               ! Work index
      integer ktm               ! Highest altitude index of cloud top
      integer m                 ! Work index

      real(r8) cabv                 ! Mix ratio of constituent above
      real(r8) cbel                 ! Mix ratio of constituent below
      real(r8) cdifr                ! Normalized diff between cabv and cbel
      real(r8) chat(plond,plev)     ! Mix ratio in env at interfaces
      real(r8) cond(plond,plev)     ! Mix ratio in downdraft at interfaces
      real(r8) const(plond,plev)    ! Gathered tracer array 
      real(r8) conu(plond,plev)     ! Mix ratio in updraft at interfaces
      real(r8) dcondt(plond,plev)   ! Gathered tend array 
      real(r8) small                ! A small number
      real(r8) mbsth                ! Threshold for mass fluxes
      real(r8) mupdudp              ! A work variable
      real(r8) minc                 ! A work variable
      real(r8) maxc                 ! A work variable
      real(r8) qn                   ! A work variable
      real(r8) fluxin               ! A work variable
      real(r8) fluxout              ! A work variable
      real(r8) netflux              ! A work variable

C-----------------------------------------------------------------------
C
      small = 1.e-36
c mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
      mbsth = 1.e-15

c Find the highest level top and bottom levels of convection
      ktm = plev
      kbm = plev
      do i = il1g, il2g
         ktm = min(ktm,jt(i))
         kbm = min(kbm,mx(i))
      end do

c Loop ever each constituent
      do m = 2,pcnst

c Gather up the constituent and set tend to zero
         do k = 1,plev
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
            end do
         end do

c From now on work only with gathered data

c Interpolate environment tracer values to interfaces
         do k = 1,plev
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc.lt.0) then
                  cdifr = 0.
               else
                  cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
               endif

c If the two layers differ significantly use a geometric averaging
c procedure
               if (cdifr.gt.1.E-6) then
                  cabv = max(const(i,km1),maxc*1.e-12)
                  cbel = max(const(i,k),maxc*1.e-12)
                  chat(i,k) = log(cabv/cbel)
     $                          /(cabv-cbel)
     $                          *cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5* (const(i,k)+const(i,km1))
               end if

c Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

c              provisional tends
               dcondt(i,k) = 0.

            end do
         end do

#ifdef CHECKNEG
         do k = 1,plev
            km1 = max(1,k-1)
            do i = il1g, il2g
               if (chat(i,k).lt.0.) then
                  write (6,*) ' negative chat ', i, k, lat, chat(i,k),
     $                 const(i,km1), const(i,k)
                  stop
               endif
            end do
         end do
#endif


c Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = plev 
         do i = il1g,il2g
            mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
            if (mupdudp.gt.mbsth) then
               conu(i,kk) = ( 
     $                       +eu(i,kk)*const(i,kk)*dp(i,kk)
     $                      )/mupdudp
            endif
            if (md(i,k).lt.-mbsth) then
               cond(i,k) =  (  
     $                       -ed(i,km1)*const(i,km1)*dp(i,km1)
     $                      )/md(i,k)
            endif
         end do

c Updraft from bottom to top
         do kk = plev-1,1,-1
            kkp1 = min(plev,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
               if (mupdudp.gt.mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1) 
     $                           +eu(i,kk)*const(i,kk)*dp(i,kk)
     $                         )/mupdudp
               endif
            end do
         end do

c Downdraft from top to bottom
         do k = 3,plev
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k).lt.-mbsth) then
                  cond(i,k) =  (  md(i,km1)*cond(i,km1) 
     $                           -ed(i,km1)*const(i,km1)*dp(i,km1)
     $                         )/md(i,k)
               endif
            end do
         end do


         do k = ktm,plev
            km1 = max(1,k-1)
            kp1 = min(plev,k+1)
            do i = il1g,il2g

c version 1 hard to check for roundoff errors
c               dcondt(i,k) = 
c     $                  +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))
c     $                    -mu(i,k)*   (conu(i,k)-chat(i,k))
c     $                    +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))
c     $                    -md(i,k)*   (cond(i,k)-chat(i,k))
c     $                   )/dp(i,k)

c version 2 har to limit fluxes
c               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k) 
c     $                 -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
c               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)
c     $                 -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

c version 3 limit fluxes outside convection to mass in appropriate layer
c these limiters are probably only safe for positive definite quantitities
c it assumes that mu and md already satify a courant number limit of 1
               fluxin =  mu(i,kp1)*conu(i,kp1) 
     $                 + mu(i,k)*min(chat(i,k),const(i,km1)) 
     $                 -(md(i,k)  *cond(i,k)   
     $                   + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
               fluxout = mu(i,k)*conu(i,k)     
     $                  +mu(i,kp1)*min(chat(i,kp1),const(i,k))
     $                 -(md(i,kp1)*cond(i,kp1) 
     $                   + md(i,k)*min(chat(i,k),const(i,k)))

               netflux = fluxin - fluxout
               if (abs(netflux).lt.max(fluxin,fluxout)*1.e-12) then
                  netflux = 0.
               endif
               dcondt(i,k) = netflux/dp(i,k)
            end do
         end do
C
         do k = kbm,plev             
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k.eq.mx(i)) then

c version 1
c                  dcondt(i,k) = (1./dsubcld(i))*
c     $              (-mu(i,k)*(conu(i,k)-chat(i,k))
c     $               -md(i,k)*(cond(i,k)-chat(i,k))
c     $              )

c version 2
c                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
c                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
c version 3
                  fluxin =  mu(i,k)*min(chat(i,k),const(i,km1))
     $                    - md(i,k)*cond(i,k)
                  fluxout = mu(i,k)*conu(i,k) 
     $                 - md(i,k)*min(chat(i,k),const(i,k))

                  netflux = fluxin - fluxout
                  if (abs(netflux).lt.max(fluxin,fluxout)*1.e-12) then
                     netflux = 0.
                  endif
c                  dcondt(i,k) = netflux/dsubcld(i)
                  dcondt(i,k) = netflux/dp(i,k)
               else if (k.gt.mx(i)) then
c                  dcondt(i,k) = dcondt(i,k-1)
                  dcondt(i,k) = 0.
               end if
            end do
         end do

c Update and scatter data back to full arrays

         do k = 1,plev
            kp1 = min(plev,k+1)
            do i = il1g,il2g
               qn = const(i,k)+dcondt(i,k)*2.*delt
               q(ideep(i),k,m) = qn

            end do
         end do
         
      end do                    ! m = 2,pcnst
C
      return
      end subroutine convtran_pjr

      subroutine q1q2pjr2
     $               (dqdt    ,dsdt    ,q       ,s       ,qs      ,
     $                qu      ,su      ,mc      ,du      ,alpha   ,
     $                qhat    ,shat    ,dp      ,mu      ,md      ,
     $                sd      ,qd      ,ql      ,dsubcld ,qdb     ,
     $                sdb     ,betau   ,betad   ,mb      ,lcl     ,
     $                jt      ,mx      ,dt      ,il1g    ,il2g    ,
     $                nstep   ,lat     ,
     $                dl      ,evp     ,cu      ,eu      ,ed)

      implicit none

#include <pmgrid.h>
#include <guang_pjr.h>

c rewritten by phil rasch dec 19 1995

      real(r8) dqdt(plond,plev),
     $     dsdt(plond,plev)

      real(r8) fact
      integer il1g
      integer il2g
      integer i
      integer k
      integer lat
      integer nstep
      real(r8) dt
      real(r8) q(plond,plev),
     $     s(plond,plev),
     $     qs(plond,plev),
     $     qu(plond,plev),
     $     su(plond,plev),
     $     mc(plond,plev),
     $     du(plond,plev),
     $     alpha(plond,plev),
     $     qhat(plond,plev),
     $     shat(plond,plev),
     $     dp(plond,plev),
     $     mu(plond,plev),
     $     md(plond,plev),
     $     sd(plond,plev),
     $     qd(plond,plev),
     $     ql(plond,plev)
      real(r8)
     $     dl(plond,plev),
     $     evp(plond,plev),
     $     cu(plond,plev),
     $     eu(plond,plev),
     $     ed(plond,plev)
      integer kbm
      integer ktm
      real(r8) emc
      real(r8) dsubcld(plond),
     $     qdb(plond),
     $     sdb(plond),
     $     betau(plond),
     $     betad(plond),
     $     mb(plond)

      integer lcl(plond),
     $     jt(plond),
     $     mx(plond)
C
C work fields:
C
C-------------------------------------------------------------------
      do k = msg + 1,plev
         do i = il1g,il2g
            dsdt(i,k) = 0.
            dqdt(i,k) = 0.
            dl(i,k) = 0.
         end do
      end do
c
c find the highest level top and bottom levels of convection
c
      ktm = plev
      kbm = plev
      do i = il1g, il2g
         ktm = min(ktm,jt(i))
         kbm = min(kbm,mx(i))
      end do
C
      fact = 0.
c      fact = 1.

      do k = ktm,plev-1
         do i = il1g,il2g
c            fact = 1.
c            if (q(i,k).gt.0.8*qs(i,k) .and. k.lt.plev-3) fact = 0.

            emc = +fact*du(i,k)*ql(i,k+1) ! evaporating cloud detraining to env
     $            -cu(i,k)                ! condensation in updraft
     $            +evp(i,k)               ! evaporating rain in downdraft
c            emc = 0

            dsdt(i,k) = -rl/cpg*emc 
     $                  + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1))
     $                     -mu(i,k)*   (su(i,k)-shat(i,k))
     $                     +md(i,k+1)* (sd(i,k+1)-shat(i,k+1))
     $                     -md(i,k)*   (sd(i,k)-shat(i,k))
     $                    )/dp(i,k)
     $     

            dqdt(i,k) = emc
     $                  +(+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1))
     $                    -mu(i,k)*   (qu(i,k)-qhat(i,k))
     $                    +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1))
     $                    -md(i,k)*   (qd(i,k)-qhat(i,k))
     $                   )/dp(i,k)
     $                  

            dl(i,k) = (1-fact)*du(i,k)*ql(i,k+1)

         end do
      end do

C
      do k = kbm,plev             
         do i = il1g,il2g
            if (k.eq.mx(i)) then
               dsdt(i,k) = (1./dsubcld(i))*
     $              (-mu(i,k)* (su(i,k)-shat(i,k))
     $               -md(i,k)* (sd(i,k)-shat(i,k))
     $              )
               dqdt(i,k) = (1./dsubcld(i))*
     $              (-mu(i,k)*(qu(i,k)-qhat(i,k))
     $               -md(i,k)*(qd(i,k)-qhat(i,k))
     $              )
            else if (k.gt.mx(i)) then
               dsdt(i,k) = dsdt(i,k-1)
               dqdt(i,k) = dqdt(i,k-1)
            end if
         end do
      end do
C
      return
      end subroutine q1q2pjr2

      subroutine conv_ini (gravx, cpx, rgasx, hypi, xtrigon)
c
c     a start on a subroutine to initialize the constants in guangs code
c
C-----------------------------------------------------------------------
#include <implicit.h>
C--------------------------Commons--------------------------------------
#include <guang_pjr.h>
#include <pmgrid.h>
C-----------------------------Arguments---------------------------------
      real(r8) gravx                ! input val of gravity
      real(r8) cpx                  ! input val of heat capacity of dry air
      real(r8) rgasx                ! input val of gas constant for dry air
      real(r8) hypi(plevp)          ! reference pressures at interfaces
      logical xtrigon           ! true => triggers turned on

C
      integer limcnv
      integer i
      integer k
c
      cpg = cpx
      grav = gravx
      rgas = rgasx
      rgrav = 1./grav
      trigon = xtrigon
c
      a = 21.656
      b = 5418.
      c1 = 6.112
      c2 = 17.67
      c3 = 243.5
      eps1 = 0.622
      qmin = 1.E-20
      tfreez = 273.16
      rl = 2.5104E6
c
c Limit convection to regions below 40 mb
c
      if (hypi(1) .ge. 4.e3) then
        limcnv = 1
      else
        do k=1,plev
          if (hypi(k).lt.4.e3 .and. hypi(k+1).ge.4.e3) then
            limcnv = k
            goto 10
          end if
        end do
        limcnv = plevp
      end if
   10 continue
c
c Set internal variable "msg" (convection limit) to "limcnv-1"
c
      msg = limcnv - 1
C
C tau=4800. were used in canadian climate center. however, when it
C is used here in echam3 t42, convection is too weak, thus 
C adjusted to 2400. i.e the e-folding time is 1 hour now.
C
      tau = 7200.
      dd = 0.1 ! fraction of updraft mass flux going into downdraft
c      dd = 0.3 ! fraction of updraft mass flux going into downdraft

      write (6,*) ' zhang/mcfarlane convection initialization '
      write (6,*) ' relaxation time tau is set to ', tau
      write (6,*) ' downdraft param dd is set to ', dd
      write (6,*) ' top index for zhang convection is ', msg
      write (6,*) ' convective trigger is set ', trigon
      write (6,*) ' done '

#ifdef DEBCONV
#ifdef LINUX
      latlook = 1
      ilook = 1
#else
      latlook = 1
      ilook = 69
c      ilook = 1
c      latlook = .true.
c      latlook = .false.
#endif
         write (6,*) ' debuging latlook, ilook ', latlook, ilook
#endif
      return
      end subroutine conv_ini

      subroutine findcin (lat, t, q, tprime, qprime, 
     $     mx, pm, pi, zm, lclo, bmax, lfc, lnb, lcape)

c     find the convective inhibition for each parcel
c     define the dry and moist adiabats as lines of constant s and h
c     where s = c_p T + gz, h = s + Lq

c     note that the variations of c_p and R (for hydrostatic equation) with
c     water vapor are ignored

      implicit none
#include <pmgrid.h>
c     
      integer lat                       ! lat index for this slice
      real(r8) t(plond,plev)                ! temperature of env
      real(r8) q(plond, plev)               ! moisture of env
      real(r8) tprime(plond,plev)           ! tpert for updraft plume
      real(r8) qprime(plond,plev)           ! qpert for updraft plume
      real(r8) pm(plond,plev)               ! pressume at layer midpoints
      real(r8) pi(plond,plevp)              ! pressure at layer interfaces
      real(r8) zm(plond,plev)               ! height at layer midpoints
      integer lclo(plond)               ! lcl of parcel leaving layer mx
      integer mx(plond)                 ! level of parcel with max cin
      real(r8) bmax(plond)                  ! buoyancy between mx and lcl
      integer klcl(plond,plev)
      integer klfc(plond,plev)
      integer klnb(plond,plev)
      integer lnb(plond)
      integer lfc(plond)
      real(r8) lcape(plond)
      real(r8) capen

c     local workspace
      real(r8) tubase(plond,plev)
      real(r8) qu(plond,plev)
      real(r8) bloc
      real(r8) bloc2
      real(r8) btot(plond,plev)
      real(r8) dq
      real(r8) ql(plond)
      real(r8) tl(plond)
      real(r8) tvl
      real(r8) tev(plond,plev)
      real(r8) es
      real(r8) buoy
      real(r8) c1n
      real(r8) c2n
      real(r8) e
      real(r8) tlcl(plond), plcl(plond)
      real(r8) plexp, qlcl(plond), esd
      real(r8) dp
      real(r8) bpos(plond)
      real(r8) tparc(plond,plev)
      real(r8) qparc(plond,plev)
      real(r8) hold(plond)
      real(r8) hnew(plond)
      real(r8) tsat(plond)
      real(r8) qs(plond)
      real(r8) tnew
      real(r8) qsold
      real(r8) dtd
      real(r8) dtw
      real(r8) err
      real(r8) capgam
      real(r8) dz
      real(r8) capblfc(plond,plev)
      real(r8) cape(plond,plev)
      real(r8) emax
      real(r8) qold(plond), told(plond)
      real(r8) ppa(plond), esv(plond)
      real(r8) gam(plond)
      integer imax
      real(r8) tlast(plond)
      real(r8) qlast(plond)
      real(r8) gmax


      integer kb
      integer i
      integer itmax
      integer k
      integer lcl(plond,plev)
      integer numlcl, numhuge
      integer kt
      integer iter
      integer ktop
      integer nmax

#include <guang_pjr.h>
#include <eslookup.h>

      do k = plev,1,-1
         do i = 1,plon
            qu(i,k) = q(i,k) + qprime(i,k) ! moisture for parcels
            tubase(i,k) = t(i,k) + tprime(i,k) ! temperature for parcels
c            tev(i,k) = t(i,k)*(1.+.608*q(i,k)) ! virtual temperature for env
            tev(i,k) = t(i,k)
     $               * (1.+1.608*q(i,k))/ (1.+q(i,k))
            lcl(i,k) = 0
            btot(i,k) = -1.e36
            klcl(i,k) = 0
            cape(i,k) = 0
            capblfc(i,k) = 0
            klfc(i,k) = 0
            klnb(i,k) = plev+1
         end do
      end do

      do i = 1,plon
         mx(i) = 1
         lclo(i) = 0
         bmax(i) = -1.e36
         bpos(i) = 0.
         lcape(i) = -1.e36
      end do

      c1n = grav/cpg
      c2n = hlatv/cpg
      capgam = -c1n
      itmax = 5
      ktop = plev+1                     ! highest level to lift a parcel
      kt = plev+1                       ! highest level to start parcels from
      do k = 1,plev
         do i = 1,plon
            if (pm(i,k).gt.100.) then
               ktop = min(k-1,ktop)
            endif
c            if (pm(i,k).gt.850.) then
            if (pi(i,plev)-pm(i,k).lt.150) then
               kt = min(k-1,kt)
            endif
         end do
      end do
#ifdef DEBCONV
      write (6,*) ' lat ', lat
      write (6,*) ' kt highest possible level of cloud base ', kt
      write (6,*) ' ktop highest possible level of cloud top ', ktop
#endif
            
c     first calculate the properties of parcel at base of lifting
      do kb = plev,kt,-1
#ifdef DEBCONV
         if (lat.eq.latlook) then
            write (6,*) ' kb ', kb
            write (6,*) ' k, lnb, lfc, lcl, tpv, tev, tp, t, p'//
     $              ' bouy, bloc2, capblfc, cape '
         endif
#endif
         do i = 1,plon
c            write (6,*) ' override ', i
c            tubase(i,kb) = 190.1707794843005
c            qu(i,kb) = 9.833855485759975E-6
c            pm(i,kb) =  647.8308110287944
c            zm(i,kb) = 3117.237132057446
            ppa(i) = pm(i,kb)*100.
         end do
c         tubase(1,kb) = 300.
         call vqsatd(tubase(1,kb), ppa, esv, qs, gam, plon)
c         qu(1,kb) = qs(1)*1.04
         do i = 1,plon

            plcl(i) = 0.

c           start assuming the parcels are unsaturated
            tparc(i,kb) = tubase(i,kb)
            qparc(i,kb) = qu(i,kb)
            hold(i) = cpg*tparc(i,kb)+grav*zm(i,kb)+hlatv*qparc(i,kb)

c           check whether they are saturated
            tsat(i) = tparc(i,kb)

c           if saturated, then first assume water condenses out and
c           releases all latent heat
            if (qs(i) .le. qparc(i,kb)) then
               plcl(i) = max(plcl(i),pm(i,kb))
               klcl(i,kb) = max(klcl(i,kb),kb)
               gmax = max(gam(i),0.1_r8)   ! limit the size of dtw for cold temp
               dtw = hlatv*(qparc(i,kb)-qs(i))/(cpg*gmax)
               tsat(i) = tparc(i,kb) + dtw
c               es = estblf(tsat(i))*0.01
c               qs(i) = eps1*es/(pm(i,kb)-es)
            endif
         end do
#ifdef DEBCONV
         if (lat.eq.latlook) then
            i = ilook
            write (6,*) ' deb 1 ', tsat(i), qs(i), tparc(i,kb), gam(i)
         endif
#endif

c        now correct it so parcel is just saturated and conserves h
         do iter=1,itmax
            do i = 1,plon
               told(i) = tsat(i)
               qold(i) = qs(i)
            end do
            call vqsatd(tsat, ppa, esv, qs, gam, plon)
            nmax = 0
            do i = 1,plon
               if (klcl(i,kb).ne.0) then
                  hnew(i) = cpg*tsat(i)+grav*zm(i,kb)+hlatv*qs(i)
                  err = hnew(i) - hold(i)
                  tsat(i) = tsat(i) - err/(cpg*(1+gam(i)))
                  es = estblf(tsat(i))*0.01
                  qs(i) = eps1*es/(pm(i,kb)-(1.-eps1)*es)
C
C The following check is to avoid the generation of negative
C values that can occur in the upper stratosphere and mesosphere
C
                  qs(i) = min(one,qs(i))
C
                  if (qs(i) .lt. 0.0) then
                     qs(i) = 0.0
                     es = pm(i,k)
                  end if
                  hnew(i) = cpg*tsat(i)+grav*zm(i,kb)+hlatv*qs(i)
                  err = abs((hnew(i)-hold(i))/(hnew(i)+hold(i)))
                  if (err.gt.1.e-4) then
                     nmax = nmax+1
                  endif
               endif
            end do
#ifdef DEBCONV
            if (lat.eq.latlook) then
               i = ilook
               write (6,*) ' deb 2 ',
     $              iter,hnew(i),hold(i),tsat(i), qs(i)
            endif
#endif
            if (nmax.eq.0) goto 10
         end do

c        now check to make sure it converged
         emax = 0.
         imax = 0
         do i = 1,plon
            if (klcl(i,kb).ne.0) then
               hnew(i) = cpg*tsat(i)+grav*zm(i,kb)+hlatv*qs(i)
               tparc(i,kb) = tsat(i)
               qparc(i,kb) = qs(i)
               err = abs((hnew(i)-hold(i))/(hnew(i)+hold(i)))
               if (err.gt.emax) then
                  emax = err
                  imax = i
               endif
            endif
         end do
#ifdef DEBCONV
         if (lat.eq.latlook) then
            i = ilook
            write (6,*) ' deb 3 ', 
     $           iter, hnew(i), hold(i), tsat(i), qs(i)
         endif
#endif

         if (emax.gt.1.e-4) then
            i = imax
            write (6,*) ' problems initializing moist k base ',
     $           emax, lat, i, klcl(i,kb), hold(i), hnew(i)
            write (6,*) ' told, tnew ', 
     $           told(i), tsat(i)
            write (6,*) ' qold, qnew ', 
     $           qold(i), qs(i)
            write (6,*) ' tubase, qu ', tubase(i,kb), qu(i,kb)
            write (6,*) ' pm, zm ', pm(i,kb), zm(i,kb)
            stop
         endif
   10    continue
c         stop 111
         
c        now check it for lfc
         do i = 1,plon
c            tvl = tparc(i,kb)*(1.+.608*qparc(i,kb))
            tvl = tparc(i,kb)
     $               * (1.+1.608*qparc(i,kb))/ (1.+qparc(i,kb))
            buoy = tvl-tev(i,kb)
            if (klcl(i,kb).gt.0.and.buoy.gt.0) then
               klfc(i,kb) = max(klfc(i,kb),kb)
            endif
            bloc2 = rgas*buoy*(pi(i,kb+1)-pi(i,kb))/pm(i,kb)
c            bloc2 = rgas*buoy*log(pi(i,kb+1)/pi(i,kb))
            capen = cape(i,kb) + bloc2
            if (klfc(i,kb).eq.0) then
               capblfc(i,kb) = capblfc(i,kb) + bloc2
            endif
            if (buoy.ge.0) then
               klnb(i,kb) = min(klnb(i,kb),kb)
            endif
            if (klfc(i,kb).eq.0.or.buoy.ge.0) then
               cape(i,kb) = capen
            endif
#ifdef DEBCONV
            if (lat.eq.latlook.and.i.eq.ilook) then
               write (6,69) kb, klnb(i,kb), klfc(i,kb),
     $               klcl(i,kb), tvl, tev(i,kb), 
     $              tparc(i,kb), t(i,kb), pm(i,kb), 
     $              buoy, bloc2,
     $              capblfc(i,kb), cape(i,kb)
            endif
#endif
         end do

#ifdef DEBCONV
         write (6,*) ' completed cloud base test for level ',kb
#endif
c        move the parcel up
         do k = kb-1,ktop,-1
            do i = 1,plon
               hold(i) = cpg*tparc(i,k+1)
     $              +grav*zm(i,k+1)+hlatv*qparc(i,k+1)

c              first assume it moves up a dry adiabat
               dz = (zm(i,k)-zm(i,k+1))
               dtd = capgam*dz
               td = tparc(i,k+1) + dtd
c              es = exp( -6763.6/td - 4.9283*log(td) + 54.2190 )
c               es = estblf(td)*0.01
c               qs(i) = eps1*es/(pm(i,k)-es)
               tparc(i,k) = tparc(i,k+1) + dtd
               qparc(i,k) = qparc(i,k+1)

               hnew(i) = cpg*tparc(i,k)+grav*zm(i,k)+hlatv*qparc(i,k)
               ppa(i) = pm(i,k)*100.
            end do

c           if it hits saturation, note it, and move up a wet adiabat
c           note this is very approximate because
c           1) qs is assumed linear
c           2) the lcl was assumed occurred precisely at level k+1

            call vqsatd(tparc(1,k), ppa, esv, qs, gam, plon)
            do i = 1,plon
               if (qs(i) .lt. qparc(i,k+1)) then
                  plcl(i) = max(plcl(i),pm(i,k))
                  klcl(i,kb) = max(klcl(i,kb),k)
                  dtw = capgam*dz/(1+gam(i))
                  tsat(i) = tparc(i,k+1) + dtw
                  es = estblf(tsat(i))*0.01
                  qs(i) = eps1*es/(pm(i,k)-(1.-eps1)*es)
C
C The following check is to avoid the generation of negative
C values that can occur in the upper stratosphere and mesosphere
C
                  qs(i) = min(one,qs(i))
C
                  if (qs(i) .lt. 0.0) then
                     qs(i) = 0.0
                     es = pm(i,k)
                  end if
                  tparc(i,k) = tsat(i)
               endif
            end do
            
c           now correct any parcels that have hit their lcl
c           so they conserve h, and are saturated.
c           using a newton method iteration
c           this remove the error in assumptions mentioned above
            do iter = 1,itmax
               call vqsatd(tparc(1,k), ppa, esv, qs, gam, plon)
               nmax = 0
               do i = 1,plon
                  if (klcl(i,kb).ne.0) then
                     hnew(i) = cpg*tsat(i)+grav*zm(i,k)+hlatv*qs(i)
                     err = hnew(i) - hold(i)
                     dtw =  - err/(cpg*(1+gam(i)))
                     tsat(i) = tsat(i) + dtw
                     es = estblf(tsat(i))*0.01
                     qs(i) = eps1*es/(pm(i,k)-(1.-eps1)*es)
C
C The following check is to avoid the generation of negative
C values that can occur in the upper stratosphere and mesosphere
C
                     qs(i) = min(one,qs(i))
C
                     if (qs(i) .lt. 0.0) then
                        qs(i) = 0.0
                        es = pm(i,k)
                     end if
#ifdef DEBCONV
                     if (lat.eq.latlook.and.i.eq.ilook) then
                        write (6,*) ' deb 4 ', k,hnew(i), hold(i), dtw,
     $                       tsat(i), pm(i,k),gam(i)
                     endif
#endif
                     tparc(i,k) = tsat(i)
                     qparc(i,k) = qs(i)
                  end if
                  hnew(i) = cpg*tparc(i,k)+grav*zm(i,k)+hlatv*qparc(i,k)
                  err = abs((hnew(i)-hold(i))/(hnew(i)+hold(i)))
                  if (err.gt.1.e-4) then
                     nmax = nmax + 1
                  endif
               end do
               if (nmax.eq.0) goto 20
            end do

c           now check to see whether it converged.
            emax = 0.
            imax = 0
            do i = 1,plon
               hnew(i) = cpg*tparc(i,k)+grav*zm(i,k)+hlatv*qparc(i,k)
               err = abs((hnew(i)-hold(i))/(hnew(i)+hold(i)))
               if (err.gt.emax) then
                  emax = err
                  imax = i
               endif
            end do
            if (emax.gt.1.e-4) then
               i = imax
               write (6,*) ' problems initializing moist k mid',lat,k,
     $              i, hold(i), hnew(i), tparc(i,k), zm(i,k), gam(i)
               stop
            end if
   20       continue

c           now check it for lfc
            do i = 1,plon
c               tvl = tparc(i,k)*(1.+.608*qparc(i,k))
               tvl = tparc(i,k)
     $               * (1.+1.608*qparc(i,k))/ (1.+qparc(i,k))
               buoy = tvl-tev(i,k)
               if (klcl(i,kb).gt.0.and.buoy.gt.0) then
                  klfc(i,kb) = max(klfc(i,kb),k)
               endif
               bloc2 = rgas*buoy*(pi(i,k+1)-pi(i,k))/pm(i,k)
               capen = cape(i,kb) + bloc2
               if (klfc(i,kb).eq.0) then
                  capblfc(i,kb) = capblfc(i,kb) + bloc2
               endif
               if (buoy.ge.0) then
                  klnb(i,kb) = min(klnb(i,kb),k)
               endif
               if (klfc(i,kb).eq.0.or.buoy.ge.0) then
                  cape(i,kb) = capen
               endif
#ifdef DEBCONV
               if (lat.eq.latlook.and.i.eq.ilook) then
                  write (6,*) ' hnew ', hnew(i)
                  write (6,69) k, klnb(i,kb), klfc(i,kb),
     $                 klcl(i,kb), tvl, tev(i,k),
     $                 tparc(i,k), t(i,k), pm(i,k), 
     $                 buoy, bloc2,
     $                 capblfc(i,kb), cape(i,kb)
   69              format (4i5,10f11.3)
               endif
#endif
            end do
         end do

      end do

c     now choose the cloud base
      do kb = plev,kt,-1
         do i = 1,plon
c            if (capblfc(i,kb).gt.bmax(i)) then
c            if (capblfc(i,kb).ge.bmax(i)
c     $           .and.cape(i,kb).gt.lcape(i)) then
c           choose a layer with the least inhibition
c           till we find a positive capeblfc, then stop
            if (capblfc(i,kb).ge.bmax(i)
     $           .and.bmax(i).lt.0.) then
               mx(i) = kb
               lclo(i) = klcl(i,kb)
               bmax(i) = capblfc(i,kb)
               lfc(i) = klfc(i,kb)
               lnb(i) = klnb(i,kb)
               lcape(i) = cape(i,kb)
            endif
         end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         write (6,*) ' properties of best parcel ', 
     $        'kb, mx(i), lclo(i), bmax(i), lfc(i), lnb(i), lcape(i)',
     $        kb, mx(i), lclo(i), bmax(i), lfc(i), lnb(i), lcape(i)
      endif
#endif

      return
      end subroutine findcin

      end module conv_pjr
