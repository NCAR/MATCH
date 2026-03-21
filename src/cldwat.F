      module cldwat

      use precision

      real(r8) :: one = 1.
      real(r8) :: zero = 0.

      contains

#define SHORTMASK
#undef SHORTMASK

c $Id$
      subroutine inimc( etamid )
c
c initialize the common block variables for the prognostic condensate 
c parameterization
c
c phil rasch april 1997      
c
#include <implicit.h> 
#include <pmgrid.h> 
#include <mc.h>

      real(r8) etamid(plev)
      integer k

#ifdef CRAY
      real(r8) signgam              ! variable required by cray gamma function
      external gamma
#endif
      rhonot = 1.275e-3         ! air density at surface (gm/cm3)
      rhos = .1                 ! assumed snow density (gm/cm3)
      rhow = 1.                 ! water density
      rhoi = 1.                 ! ice density
      esi = 1.0                 ! collection efficient for ice by snow
      esw = 0.1                 ! collection efficient for water by snow
      t0 = 273.16               ! approximate freezing temp
      cldmin = 0.02             ! assumed minimum cloud amount 
      small = 1.e-22            ! a small number compared to unity
      c = 152.93                ! constant for graupel like snow cm**(1-d)/s
      d = 0.25                  ! constant for graupel like snow
      nos = 3.e-2                ! particles snow / cm**4
      pi = 4.*atan(1.0)         
      prhonos = pi*rhos*nos
      thrpd = 3. + d
#ifdef CRAY
      call gamma(3.+d, signgam, gam3pd)
      gam3pd = sign(exp(gam3pd),signgam)
      call gamma(4.+d, signgam, gam4pd)
      gam4pd = sign(exp(gam4pd),signgam)
      write (6,*) ' d, gamma(3+d), gamma(4+d) =', gam3pd, gam4pd
#else
      if (d.eq.0.25) then
         gam3pd = 2.549256966718531 ! only right for d = 0.25
         gam4pd = 8.285085141835282
      else
         write (6,*) ' can only use d ne 0.25 on a cray '
         stop
      endif
#endif
      mcon01 = pi*nos*c*gam3pd/4.
      mcon02 = 1./(c*gam4pd*sqrt(rhonot)/(6*prhonos**(d/4.)))
      mcon03 = -(0.5+d/4.)
      mcon04 = 4./(4.+d)
      mcon05 = (3+d)/(4+d)
      mcon06 = (3+d)/4.
      mcon07 = mcon01*sqrt(rhonot)*mcon02**mcon05/prhonos**mcon06
      mcon08 = -0.5/(4.+d)

c     find the level where the etamids are all greater than 10 hPa
      ktop = 1
      do k = 1,plev
         if (etamid(k).lt.0.01) then
            ktop = k
         endif
      end do

      write (6,*) 'INIMC: prognostic cloud water capped at level ',ktop
     $,           ' whose midpoint is ',etamid(ktop)*1.e5,' pascals.'

      return
      end subroutine inimc

c##############################################################################

      subroutine inimland (oro, clat)

      implicit none

#include <pmgrid.h>      
#include <mc.h>
      real(r8) oro(plond,plat)
      real(r8) clat(plat) ! lat in radians
      real(r8) clon(plond)
      real(r8) cont(plond,plat)
      real(r8) temp(plond,plat)
      real(r8) dmax
      real(r8) arad
      real(r8) dist
      real(r8) sum
      real(r8) cs(plat)
      real(r8) ss(plat)
      real(r8) c1
      real(r8) s1
      real(r8) c2
      real(r8) s2
      real(r8) dx
      real(r8) dy
      real(r8) term

      integer i
      integer j
      integer ii
      integer jj
      integer iplm1
      integer jof
      integer iof
      integer itmp
      integer jmin, jmax

      dmax = 2.e6 ! distance to carry the mask
      arad = 6.37e6
      do i = 1,plon
         clon(i) = 2.*(i-1)*pi/plon
      end do
c
c first isolate the contenents 
c as land points not surrounded by ocean or ice
c
      do j = 1,plat
         cs(j) = cos(clat(j))
         ss(J) = sin(clat(j))
         do i = 1,plon
            cont(i,j) = 0.
            if (nint(oro(i,j)) .eq. 1) then
               cont(i,j) = 1.
            endif
         end do
         temp(1,j) = cont(1,j)
         temp(plon,j) = cont(plon,j)
      end do

      do i = 1,plon
         temp(i,1) = cont(i,1)
         temp(i,plat) = cont(i,plat)
      end do
c
c get rid of one and two point islands
c
      do j = 2,plat-1
         do i = 2,plon-1
            sum =  cont(i,j+1) + cont(i,j-1) 
     $           + cont(i+1,j+1) + cont(i+1,j-1)
     $           + cont(i-1,j+1) + cont(i-1,j-1) 
     $           + cont(i+1,j) + cont(i-1,j)
     $           + cont(i,j)
            if (sum.le.2.) then
               temp(i,j) = 0.
            else
               temp(i,j) = 1.
            endif
         enddo
      end do

      do j = 1,plat
         do i = 1,plon
            cont(i,j) = temp(i,j)
         end do
      end do
#ifdef SHORTMASK
      write (6,*) ' FAST LANDMASK SELECTED '
      do j = 1,plat
         do i = 1,plon
            if (cont(i,j).gt.0.9) then
               landm(i,j) = 1
            else
               landm(i,j) = 0
            endif
         end do
      end do
#else
c
c construct a function which is one over land, 
c zero over ocean points beyond dmax from land
c
      iplm1 = 2*plon - 1
      dy = pi*arad/plat
      jof = dmax/dy + 1
c      write (6,*) ' lat bands to check ', 2*jof+1
      do j = 1,plat
         c1 = cs(j)
         s1 = ss(j)
         dx = 2*pi*arad*cs(j)/plon
         iof = min(int(dmax/dx) + 1, plon)
         do i = 1,plon
            temp(i,j) = 0.
            landm(i,j) = 0.
            jmin = max(1,j-jof)
            jmax = min(plat,j+jof)
            do jj = jmin, jmax
               s2 = ss(jj)
               c2 = cs(jj)
               do itmp = -iof,iof
                  ii = mod(i+itmp+iplm1,plon)+1
                  term = s1*s2 + c1*c2*cos(clon(ii)-clon(i))
                  if (term.gt.0.9999999) term = 1.
                  dist = arad*acos(term)
                  landm(i,j) = max(landm(i,j), 
     $                             (1.-dist/dmax)*cont(ii,jj)
     $                            )
c                  if (dist.lt.dmax .and. cont(ii,jj).eq.1) then
c                     landm(i,j) = max(landm(i,j), 1.-dist/dmax)
c                  endif
               end do
            end do
         end do
      end do
#endif      
c
      return
      end subroutine inimland

c##############################################################################

      subroutine pcond (qtend   ,ttend   ,omega   ,lat    , 
     $                  q       ,cwat    , t      ,p      ,pdel    ,
     $                  oro     ,cldn    ,cldo    ,cme    ,evapr   ,
     $                  prain   ,rmelt   ,deltat  ,pcflx  , 
     $                  fwaut   ,fsaut   ,fracw   ,fsacw  ,fsaci   )
C-----------------------------------------------------------------------
C
C calculate the prognostic condensate amount and tendencies
C     P. Rasch and J.E. Kristjansson, April 1997
C
C-----------------------------------------------------------------------
#include <implicit.h>
#include <pmgrid.h>
#include <mc.h> 
c---------------------------------------------------------------------
C
C Input Arguments
C
      integer lat

      real(r8) cldn(plond,plevp)    ! new value of cloud fraction    (fraction)
      real(r8) cldo(plond,plevp)    ! old value of cloud fraction    (fraction)
      real(r8) cwat(plond,plev)     ! cloud water (kg/kg)
      real(r8) omega(plond,plev)    ! vert pressure vel (Pa/s)
      real(r8) oro(plond)           ! flag for land/ice/ocean
      real(r8) p(plond,plev)        ! pressure          (K)
      real(r8) pcflx(plond,plevp)   ! convective precip level by level (kg/m2/s) 
                                !   (DISABLED)
      real(r8) pdel(plond,plev)     ! pressure thickness (Pa)
      real(r8) q(plond,plev)        ! water vapor (kg/kg)
      real(r8) qtend(plond,plev)    ! mixing ratio tend (kg/kg/s)
      real(r8) t(plond,plev)        ! temperature       (K)
      real(r8) ttend(plond,plev)    ! temp tend         (K/s)
C
C Output Arguments
C
      real(r8) cme(plond,plev)      ! rate of cond-evap within the cloud
      real(r8) deltat               ! time step to advance solution over
      real(r8) evapr(plond,plev)    ! rate of evaporation of falling precip (1/s)
      real(r8) prain(plond,plev)    ! rate of conversion of condensate to precip 
                                !   (1/s)
      real(r8) rmelt(plond,plev)    ! heating rate due to precip phase change 
                                !   (K/s) (DISABLED)
c
c Local workspace
c
      integer i                 ! work variable
      integer iter              ! #iterations for precipitation calculation
      integer k                 ! work variable
      integer l                 ! work variable

      real(r8) cldm(plond)          ! mean cloud fraction over the time step
      real(r8) cldmax(plond)        ! max cloud fraction above
      real(r8) cldnt(plond,plevp)   ! sane new cloud fraction    (temp , fraction)
      real(r8) cldot(plond,plevp)   ! sane old cloud fraction    (temp , fraction)
      real(r8) coef(plond)          ! conversion time scale for condensate to rain
      real(r8) conke                ! rate of evaporation of precipitation:
      real(r8) convm(plond,plev)    ! moistening rate   (kg/kg/s)
      real(r8) cwm(plond)           ! cwat mixing ratio at midpoint
      real(r8) cwn(plond)           ! cwat mixing ratio at end
      real(r8) denom                ! work variable
      real(r8) dqsdp                ! change in sat spec. hum. wrt pressure
      real(r8) dqsdt                ! change in sat spec. hum. wrt temperature
      real(r8) es(plond)            ! sat. vapor pressure
      real(r8) fice(plond)          ! fraction of cwat that is ice
      real(r8) fracw(plond,plev)    ! relative importance of collection of 
                                !   liquid by rain
      real(r8) fsaci(plond,plev)    ! relative importance of collection of ice by 
                                !   snow
      real(r8) fsacw(plond,plev)    ! relative importance of collection of liquid 
                                !   by snow
      real(r8) fsaut(plond,plev)    ! relative importance of ice auto conversion
      real(r8) fwaut(plond,plev)    ! relative importance of warm cloud auto 
                                !   conversion
      real(r8) gamma(plond)         ! d qs / dT
      real(r8) grav
      real(r8) iceab(plond)         ! rate of ice only from above
      real(r8) icwc                 ! in-cloud water content (kg/kg)
      real(r8) mincld               ! a small cloud fraction to avoid / zero
      real(r8) omeps                ! 1 minus epsilon
      real(r8) omsm                 ! a number just less than unity (for rounding)
      real(r8) pcme(plond)          ! provisional condensation minus evaporation 
      real(r8) pcme1                ! work variable
      real(r8) precab(plond)        ! rate of precipitation (kg / (m**2 * s))
      real(r8) prect(plond)         ! rate of precipitation including convection 
                                !   (kg / (m**2 * s))
      real(r8) prprov(plond)        ! provisional value of precip at btm of layer
      real(r8) prtmp                ! work variable
      real(r8) qn(plond,plev)       ! mixing rat at end of time step ignoring 
                                !   condensate
      real(r8) qs(plond)            ! spec. hum. of water vapor
      real(r8) qsn, esn             ! work variable
      real(r8) qsp(plond,plev)      ! sat pt mixing ratio
      real(r8) qtl(plond)           ! tendency which would saturate the grid box 
                                !   in deltat
      real(r8) qtmp, ttmp           ! work variable
      real(r8) relhum(plond)        ! relative humidity
      real(r8) tc                   ! crit temp of transition to ice
      real(r8) tn(plond,plev)       ! temp at end of time step ignoring condensate
      real(r8) tsp(plond,plev)      ! sat pt temperature
      real(r8) pol                  ! work variable
      real(r8) cdt                  ! work variable
c
c------------------------------------------------------------
#include <eslookup.h>
c------------------------------------------------------------
c
      omeps = 1.0 - epsqs
#ifdef PERGRO
      mincld = 1.e-4
      iter = 1   ! number of times to iterate the precipitation calculation
#else
      mincld = 1.e-10
      iter = 2
#endif
      grav = 9.80664
      omsm = 0.99999
c
c Constant for computing rate of evaporation of precipitation:
c
      conke = 1.e-5 
c
c initialize a few single level fields
c
      do i = 1,plon
         precab(i) = 0
         prect(i) = 0
         iceab(i) = 0		! latent heat of precip above
         cldmax(i) = 0.
      end do
c
c initialize some multi-level fields
c
      do k = 1,plev
         do i = 1,plon
            qn(i,k) = q(i,k) + deltat*qtend(i,k)
            tn(i,k) = t(i,k) + deltat*ttend(i,k)
         end do
      end do

c++bee, add initialization for top layers of output variables
      do k = 1, ktop-1
         do i = 1, plon
            cme(i,k) = 0.0
            prain(i,k) = 0.0
            evapr(i,k) = 0.0
            rmelt(i,k) = 0.0
         end do
      end do
c--
c
c find the saturation point for the provisional t and q without condensation
c
      call findsp (qn, tn, p, tsp, qsp)

      do k = ktop,plev
         call vqsatd (t(1,k), p(1,k), es, qs, gamma, plon)
         do i = 1,plon
            relhum(i) = q(i,k)/qs(i)
c
c sane limits on cloud fractions
c
            cldnt(i,k) = max(zero,min(one,cldn(i,k),relhum(i)))
            cldot(i,k) = max(zero,min(one,cldo(i,k),relhum(i)))
c
c the mean cloud fraction over the time step
c
            cldm(i) = max((cldnt(i,k) + cldot(i,k))*0.5,mincld)
c
c the max cloud fraction above this level
c
            cldmax(i) = max(cldmax(i), cldm(i))
         end do
         dqsdp = 0
         do i = 1,plon
c
c fractions of ice at this level
c
            tc = t(i,k) - 273.16
            fice(i) = max(zero,min(-tc*0.05,one))
            dqsdt = cp*gamma(i)/hlatv
c
c the moistening term
c
            convm(i,k) = qtend(i,k) - cldm(i)*dqsdt*ttend(i,k)
     $                              + cldm(i)*dqsdp*omega(i,k)
c     
c amount of water which must be evaporated to change the cloud
c from fraction cldo to cldn keeping same incloud amount
c
cjr this 0.0001 should real(r8)ly be mincld--fix later
            icwc = cwat(i,k)/max(cldot(i,k),0.0001_r8)
#ifdef PERGRO
            pcme1 = min((max(cldnt(i,k),mincld) - 
     $                   max(cldot(i,k),mincld))*icwc/deltat,zero)
#else
            pcme1 = min((cldnt(i,k) - cldot(i,k))*icwc/deltat,zero)
#endif
c
c***** THE ABOVE SHOULD BE CHANGED TO ALLOW EITHER EVAP OR COND TO TAKE PLACE
c
c first guess on q-e, have not yet checked to see that there
c is enough cloud water to support evaporation
c     
            denom = 1. + cldm(i) * hlatv/cp * dqsdt
            pcme(i) = (cldm(i)*convm(i,k) + pcme1) / denom
c
c calculate the cooling due to a phase change of the rainwater
c from above
c
            if (t(i,k).ge.273.16) then
c              rmelt(i,k) =  -hlatf/cp*iceab(i)*grav/pdel(i,k) 
               rmelt(i,k) = 0.
               iceab(i) = 0.
            else
               rmelt(i,k) = 0.
            endif
         end do
c
c put reasonable bounds on the first guess
c
         do i = 1,plon
            qtmp = q(i,k) + ( qtend(i,k) - pcme(i) )*deltat
            ttmp = t(i,k) + deltat*(ttend(i,k)+rmelt(i,k)
     $           +hlatv/cp*pcme(i))
c
c qtl has the tend required to bring to saturation
c
            qtl(i) =  max((qsp(i,k) - qtmp)/deltat,zero)
c
c make sure we condense enough to bring it back to saturation
c this guards (partially) against inconsistent cloud fractions
c if no cloud, then make sure we evaporate all the water
c
            if (qtmp.gt.qsp(i,k)) then
               pcme(i) = pcme(i) + (qtmp-qsp(i,k))/deltat
            else
               if (cldm(i).gt.mincld) then
c                    limit the estimated evaporation of cloud water
c                    to that available 
c                    or that which brings the box to saturation
c                     pcme(i) = max(-cwat(i,k)/deltat,pcme(i),-qtl(i))
               else
                  pcme(i) = - min(cwat(i,k)/deltat,qtl(i))
               endif
            endif
            pcme(i) = max(-cwat(i,k)/deltat,pcme(i),-qtl(i))*omsm
         end do
         do i = 1,plon
c
c a provisional value of the cloud water from cond-evap at midpoint of timestep
c
            cwm(i) = max(cwat(i,k) + pcme(i)*deltat*0.5,zero)
c
c a provisional value at endpoint 
c
            cwn(i) = max(cwat(i,k) + pcme(i)*deltat,zero)
c
c a provisional value of precip
c
            prain(i,k) = 0.0
c
c move provisional cond-evap into final location
c
            cme(i,k) = pcme(i) 
         end do
c
c calculate the formation of precip. Since this is a highly nonlinear
c calculation, we do it iteratively, using values from the midpoint of
c the time step
c
         do l = 1,iter
            do i = 1,plon
               prprov(i) = prect(i) + prain(i,k)*pdel(i,k)/grav
            end do

            call findmcnew (lat, k, prprov, t, p, cwm, cldm, cldmax,
     $           fice, oro, coef, fwaut(1,k), fsaut(1,k),
     $           fracw(1,k), fsacw(1,k), fsaci(1,k))
c
c calculate the precip rate
c
#ifdef PERGRO
            do i = 1,plon
               if (cldm(i).gt.0) then
                  cdt = coef(i)*deltat
c                 first predict the cloud water
                  if (cdt.gt.0.1) then
                     pol = cme(i,k)/coef(i) ! production over loss
                     cwn(i) = max(
     $                             zero,
     $                            (cwat(i,k)-pol)*exp(-cdt)
     $                             + pol 
     $                           )
                  else
                     cwn(i) = max(zero,
     $                            (cwat(i,k) + cme(i,k)*deltat)
     $                            /(1+cdt) 
     $                           )
                  endif
c
c now back out the tendency
c
                  prain(i,k)=max(zero,cme(i,k)-(cwn(i)-cwat(i,k))/deltat)
               else
                  prain(i,k) = 0.0
                  cwn(i) = 0.
               endif
c
c update any remaining  provisional values
c
               cwm(i) = (cwn(i) + cwat(i,k))*0.5
            end do              ! end of do i = 1,plon

#else

            do i = 1,plon
               if (cldm(i).gt.0) then
c
c first predict the cloud water
c
                  if (coef(i).gt.0.) then
                     cwn(i) = max(zero,(cwat(i,k)-cme(i,k)/coef(i))
     $                             *exp(-coef(i)*deltat)
     $                             + cme(i,k)/coef(i))
                  else
                     cwn(i) = max(zero,cwat(i,k) + cme(i,k)*deltat)
                  endif
c
c now back out the tendency
c
                  prain(i,k)=max(zero,cme(i,k)-(cwn(i)-cwat(i,k))/deltat)
               else
                  prain(i,k) = 0.0
                  cwn(i) = 0.
               endif
c
c update any remaining  provisional values
c
               cwm(i) = (cwn(i) + cwat(i,k))*0.5
            end do              ! end of do i = 1,plon
#endif
         end do                 ! end of do l = 1,iter
c
c we now have estimates of the condensation-evaporation terms, and
c the precipitation
c so given decent estimates of pcme, calculate provisional value
c of cloud water for evapr calculation
c
         do i = 1,plon
            qtmp = q(i,k) + (qtend(i,k) - cme(i,k))*deltat
            ttmp = t(i,k) + (ttend(i,k)+rmelt(i,k))*deltat
     $                    + hlatv/cp*deltat*cme(i,k)
            esn = estblf(ttmp)
            qsn = min(epsqs*esn/(p(i,k) - omeps*esn),one)
            qtl(i) = max((qsn - qtmp)/deltat,zero)
            relhum(i) = qtmp/qsn
         end do
c
c evaporation of rain
c
         do i = 1,plon
#ifdef PERGRO
            evapr(i,k) = conke*(1. - max(cldot(i,k),mincld))*
     $         sqrt(precab(i))*(1. - min(relhum(i),1.))
#else
            evapr(i,k) = conke*(1. - cldot(i,k))*sqrt(precab(i))
     $                        *(1. - min(relhum(i),one))
#endif
c
c limit the evaporation to the amount which is entering the box
c or saturates the box
c
            prtmp = precab(i)*grav/pdel(i,k)
            evapr(i,k) = min(evapr(i,k), prtmp, qtl(i))*omsm
         end do
c
c precipitation 
c
         do i = 1,plon
           prtmp = pdel(i,k) / grav *(prain(i,k) - evapr(i,k))
           iceab(i) = iceab(i) + fice(i)*prtmp
           precab(i) = precab(i) + prtmp
           prect(i) = prect(i) + prtmp + pcflx(i,k+1)
           if (abs(precab(i)).lt.1.e-10) then
             precab(i) = 0.
           endif
           if (abs(prect(i)).lt.1.e-10) then
             prect(i) = 0.
           endif
         end do
      end do                    ! level loop (k=1,plev)

      return
      end subroutine pcond

c##############################################################################

      subroutine findmcnew (lat, k, precab, t, p, cwm, cldm, cldmax,
     $                      fice, oro, coef, fwaut, fsaut, 
     $                      fracw, fsacw, fsaci)
c
c calculate the conversion of condensate to precipitate
c
c     written by phil rasch april 1997
c
#include <implicit.h>
#include <pmgrid.h>
#include <mc.h>
c
c input args
c
      integer k                 ! level index 
      integer lat               ! lat index (for debug stuff)

      real(r8) precab(plond)        ! rate of precipitation from above (kg / (m**2 * s))
      real(r8) t(plond,plev)        ! temperature       (K)
      real(r8) p(plond,plev)        ! pressure          (Pa)
      real(r8) cldm(plond)          ! cloud fraction
      real(r8) cldmax(plond)        ! max cloud fraction above this level
      real(r8) cwm(plond)           ! condensate mixing ratio (kg/kg)
      real(r8) fice(plond)          ! fraction of cwat that is ice
      real(r8) oro(plond)           ! oro flag (1 is land)

c output arguments
      real(r8) coef(plond)          ! conversion rate (1/s)
      real(r8) fwaut(plond)         ! relative importance of liquid autoconversion (a diagnostic)
      real(r8) fsaut(plond)         ! relative importance of ice autoconversion (a diagnostic)
      real(r8) fracw(plond)         ! relative  importance of rain accreting liquid (a diagnostic)
      real(r8) fsacw(plond)         ! relative  importance of snow accreting liquid (a diagnostic)
      real(r8) fsaci(plond)         ! relative  importance of snow accreting ice (a diagnostic)

c work variables

      integer i
      integer ii
      integer ind(plond)
      integer nlons

      real(r8) alpha                ! ratio of 3rd moment radius to 2nd
      real(r8) capc                 ! constant for autoconversion
      real(r8) capn                 ! local cloud particles / cm3
      real(r8) capnc                ! cold and oceanic cloud particles / cm3
      real(r8) capnw                ! warm continental cloud particles / cm3
      real(r8) ciaut                ! coefficient of autoconversion of ice (1/s)
      real(r8) ciautb               ! coefficient of autoconversion of ice (1/s)
      real(r8) cldloc(plond)        ! non-zero amount of cloud
      real(r8) cldpr(plond)         ! assumed cloudy volume occupied by rain and cloud
      real(r8) con1                 ! work constant
      real(r8) con2                 ! work constant
      real(r8) convfw               ! constant used for fall velocity calculation
      real(r8) cracw                ! constant used for rain accreting water
      real(r8) critpr               ! critical precip rate collection efficiency changes
      real(r8) csacx                ! constant used for snow accreting liquid or ice
      real(r8) dtice                ! interval for transition from liquid to ice
      real(r8) effc                 ! collection efficiency
      real(r8) icemr(plond)         ! in-cloud ice mixing ratio
      real(r8) icrit                ! threshold for autoconversion of ice
      real(r8) icritc               ! threshold for autoconversion of cold ice
      real(r8) icritw               ! threshold for autoconversion of warm ice
      real(r8) kconst               ! const for terminal velocity (stokes regime)
      real(r8) liqmr(plond)         ! in-cloud liquid water mixing ratio
      real(r8) pracw                ! rate of rain accreting water
      real(r8) prlloc(plond)        ! local rain flux in mm/day
      real(r8) prscgs(plond)        ! local snow amount in cgs units
      real(r8) psaci                ! rate of collection of ice by snow (lin et al 1983)
      real(r8) psacw                ! rate of collection of liquid by snow (lin et al 1983)
      real(r8) psaut                ! rate of autoconversion of ice condensate
      real(r8) ptot                 ! total rate of conversion
      real(r8) pwaut                ! rate of autoconversion of liquid condensate
      real(r8) r3l                  ! volume radius
      real(r8) r3lcrit              ! critical radius at which autoconversion become efficient
      real(r8) rainmr(plond)        ! in-cloud rain mixing ratio
      real(r8) rat1                 ! work constant
      real(r8) rat2                 ! work constant
      real(r8) rdtice               ! recipricol of dtice
      real(r8) rho(plond)           ! density (mks units)
      real(r8) rhocgs               ! density (cgs units)
      real(r8) snowfr               ! fraction of precipate existing as snow
      real(r8) totmr(plond)         ! in-cloud total condensate mixing ratio
      real(r8) vfallw               ! fall speed of precipitate as liquid
      real(r8) wp                   ! weight factor used in calculating pressure dep of autoconversion
      real(r8) wt                   ! fraction of ice

c      real(r8) work1
c      real(r8) csaci
c      real(r8) csacw
c      real(r8) cwaut
c      real(r8) efact
c      real(r8) lamdas
c      real(r8) lcrit
c      real(r8) rcwm
c      real(r8) r3lc2
c      real(r8) snowmr(plond)
c      real(r8) vfalls


c     inline statement functions
      real(r8) heavy, heavym, a1, a2, heavyp, heavymp
      heavy(a1,a2) = max(zero,sign(one,a1-a2))  ! heavyside function
      heavym(a1,a2) = max(0.01_r8,sign(one,a1-a2))  ! modified heavyside function
!
! New heavy functions to perhaps address error growth problems
!
      heavyp(a1,a2) = a1/(a2+a1+1.e-36)
      heavymp(a1,a2) = (a1+0.01*a2)/(a2+a1+1.e-36)

c critical precip rate at which we assume the collector drops can change the
c drop size enough to enhance the auto-conversion process (mm/day)
      critpr = 0.5
      
      convfw = 1.94*2.13*sqrt(rhow*1000.*9.81*2.7e-4)

c liquid microphysics
c      cracw = 6                 ! beheng
      cracw = .884*sqrt(9.81/(rhow*1000.*2.7e-4)) ! tripoli and cotton

c ice microphysics
c      ciautb = 6.e-4
c      ciautb = 1.e-3
      ciautb = 5.e-4
c      icritw = 1.e-5
c      icritw = 5.e-5
      icritw = 4.e-4
c      icritc = 1.5e-6
      icritc = 5.e-6
c      icritc = 7.5e-6

      dtice = 20.
      rdtice = 1./dtice

      capnw = 400.              ! warm continental cloud particles / cm3
      capnc =  80.              ! cold and oceanic cloud particles / cm3
      capnc = 150.              ! cold and oceanic cloud particles / cm3

      kconst = 1.18e6           ! const for terminal velocity

c      effc = 1.                 ! autoconv collection efficiency following boucher 96
c      effc = .55*0.05           ! autoconv collection efficiency following baker 93
      effc = 0.55                ! autoconv collection efficiency following tripoli and cotton

      alpha = 1.1**4
      capc = pi**(-.333)*kconst*effc ! constant for autoconversion
     $       *(0.75)**(1.333)*alpha

c      r3lcrit = 0.e-6            ! 0u  crit radius where liq conversion begins
c      r3lcrit = 6.e-6            ! 6u  crit radius where liq conversion begins
c      r3lcrit = 7.e-6            ! 7u  crit radius where liq conversion begins
c      r3lcrit = 8.e-6            ! 8u crit radius where liq conversion begins
c      r3lcrit = 10.e-6           ! 10u crit radius where liq conversion begins
      r3lcrit = 5.e-6            ! 5u  crit radius where liq conversion begins
c
c find all the points where we need to do the microphysics
c and set the output variables to zero
c
      nlons = 0
      do i = 1,plon
         coef(i) = 0.
         fwaut(i) = 0.
         fsaut(i) = 0.
         fracw(i) = 0.
         fsacw(i) = 0.
         fsaci(i) = 0.
         if (cwm(i).gt.1.e-20) then
            nlons = nlons + 1
            ind(nlons) = i
         endif
      end do

cdir$ ivdep
      do ii = 1,nlons
         i = ind(ii)
c
c the local cloudiness at this level
c
         cldloc(i) = max(cldmin,cldm(i))
c
c a weighted mean between max cloudiness above, and this layer
c
         cldpr(i) = max(cldmin,(cldmax(i)+cldm(i))*0.5)
c
c decompose the suspended condensate into 
c an incloud liquid and ice phase component
c
         totmr(i) = cwm(i)/cldloc(i)
         icemr(i) = totmr(i)*fice(i)
         liqmr(i) = totmr(i)*(1.-fice(i))
c
c density
c
         rho(i) = p(i,k)/(287.*t(i,k))
         rhocgs = rho(i)*1.e-3     ! density in cgs units
c
c decompose the precipitate into a liquid and ice phase 
c
         if (t(i,k).gt.t0) then
            vfallw = convfw/sqrt(rho(i))
            rainmr(i) = precab(i)/(rho(i)*vfallw*cldpr(i))
            snowfr = 0
         else
            snowfr = 1
            rainmr(i) = 0.
         endif
c
c local snow amount in cgs units
c
         prscgs(i) = precab(i)/cldpr(i)*0.1*snowfr
c
c local rain amount in mm/day
c
         prlloc(i) = precab(i)*86400./cldpr(i)
      end do

      con1 = 1./(1.333*pi)**0.333 * 0.01 ! meters
c
c calculate the conversion terms
c
cdir$ ivdep
      do ii = 1,nlons
         i = ind(ii)
         rhocgs = rho(i)*1.e-3     ! density in cgs units
c         
c exponential temperature factor
c
c        efact = exp(0.025*(t(i,k)-t0))
c
c some temperature dependent constants
c
         wt = min(one,max(zero,(t0-t(i,k))*rdtice))
         icrit = icritc*wt + icritw*(1-wt)
c
c linear weight factor in pressure (1 near sfc, 0 at .8 of sfc) 
c
         wp = min(one,max(zero,(p(i,k)-0.8*p(i,plev))/(0.2*p(i,plev))))
c
c near land near sfc raise the number concentration
c
         capn =  landm(i,lat)*(capnw*wp + capnc*(1-wp))
     $          +(1.-landm(i,lat))*capnc
c
c useful terms in following calculations
c
         rat1 = rhocgs/rhow
         rat2 = liqmr(i)/capn
         con2 = (rat1*rat2)**0.333 
c
c volume radius
c
c        r3l = (rhocgs*liqmr(i)/(1.333*pi*capn*rhow))**0.333 * 0.01 ! meters
         r3l = con1*con2 
c
c critical threshold for autoconversion if modified for mixed phase
c clouds to mimic a bergeron findeisen process
c r3lc2 = r3lcrit*(1.-0.5*fice(i)*(1-fice(i)))
c
c autoconversion of liquid
c
c        cwaut = 2.e-4
c        cwaut = 1.e-3
c        lcrit = 2.e-4
c        lcrit = 5.e-4
c        pwaut = max(zero,liqmr(i)-lcrit)*cwaut
c
c pwaut is following tripoli and cotton (and many others)
c we reduce the autoconversion below critpr, because these are regions where
c the drop size distribution is likely to imply much smaller collector drops than
c those relevant for a cloud distribution corresponding to the value of effc = 0.55 
c suggested by cotton (see austin 1995 JAS, baker 1993)

c easy to follow form
c        pwaut = capc*liqmr(i)**2*rhocgs/rhow
c    $           *(liqmr(i)*rhocgs/(rhow*capn))**(.333)
c    $           *heavy(r3l,r3lcrit)
c    $           *max(0.10_r8,min(one,prlloc(i)/critpr))
c somewhat faster form
#ifdef PERGRO
           pwaut = capc*liqmr(i)**2*rat1
     $             *con2
     $             *heavymp(r3l,r3lcrit)
     $             *max(0.10_r8,min(one,prlloc(i)/critpr))
#else
           pwaut = capc*liqmr(i)**2*rat1
     $             *con2
     $             *heavym(r3l,r3lcrit)
     $             *max(0.10_r8,min(one,prlloc(i)/critpr))
#endif
c
c autoconversion of ice
c
c        ciaut = ciautb*efact
         ciaut = ciautb
c        psaut = capc*totmr(i)**2*rhocgs/rhoi
c     $           *(totmr(i)*rhocgs/(rhoi*capn))**(.333)
c
c autoconversion of ice condensate
c
#ifdef PERGRO
           psaut = heavyp(icemr(i),icrit)*icemr(i)*ciaut
#else
           psaut = max(zero,icemr(i)-icrit)*ciaut
#endif
c
c collection of liquid by rain 
c
c        pracw = cracw*rho(i)*liqmr(i)*rainmr(i) !(beheng 1994)
         pracw = cracw*rho(i)*sqrt(rho(i))*liqmr(i)*rainmr(i) !(tripoli and cotton)
c
c the following lines calculate the slope parameter and snow mixing ratio
c from the precip rate using the equations found in lin et al 83
c in the most natural form, but it is expensive, so after some tedious
c algebraic manipulation you can use the cheaper form found below
c            vfalls = c*gam4pd/(6*lamdas**d)*sqrt(rhonot/rhocgs)
c     $               *0.01   ! convert from cm/s to m/s
c            snowmr(i) = snowfr*precab(i)/(rho(i)*vfalls*cldpr(i))
c            snowmr(i) = ( prscgs(i)*mcon02 * (rhocgs**mcon03) )**mcon04
c            lamdas = (prhonos/max(rhocgs*snowmr(i),small))**0.25
c            csacw = mcon01*sqrt(rhonot/rhocgs)/(lamdas**thrpd) 
c
c coefficient for collection by snow independent of phase
c
         csacx = mcon07*rhocgs**mcon08*prscgs(i)**mcon05
c
c collection of liquid by snow (lin et al 1983)
c
         psacw = csacx*liqmr(i)*esw
c
c collection of ice by snow (lin et al 1983)
c
         psaci = csacx*icemr(i)*esi
c
c total conversion of condensate to precipitate 
c
         ptot = pwaut + psaut + pracw + psacw + psaci
c
c the recipricol of cloud water amnt (or zero if no cloud water)
c
c         rcwm =  totmr(i)/(max(totmr(i),small)**2)
c
c turn the tendency back into a loss rate (1/seconds)
c
         if (totmr(i).gt.0.) then
            coef(i) = ptot/totmr(i)
         else
            coef(i) = 0.
         endif

         fwaut(i) = pwaut/max(ptot,small)
         fsaut(i) = psaut/max(ptot,small)
         fracw(i) = pracw/max(ptot,small)
         fsacw(i) = psacw/max(ptot,small)
         fsaci(i) = psaci/max(ptot,small)
 
      end do

      return
      end subroutine findmcnew

c##############################################################################

      subroutine findsp (q, t, p, tsp, qsp)
c
c     find the saturation point for a given t and q
c     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
c     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
c
#include <implicit.h>
#include <pmgrid.h>
#include <mc.h>
c-----------------------------------------------------------------------
c
c     input arguments
c
      real(r8) q(plond,plev)        ! water vapor (kg/kg)
      real(r8) t(plond,plev)        ! temperature (K)
      real(r8) p(plond,plev)        ! pressure    (Pa)
c
c output arguments
c
      real(r8) tsp(plond,plev)      ! saturation temp (K)
      real(r8) qsp(plond,plev)      ! saturation mixing ratio (kg/kg)
c
c local variables
c
      integer i                 ! work variable
      integer k                 ! work variable
      logical lflg              ! work variable
      integer iter              ! work variable
      integer l                 ! work variable

      real(r8) omeps                ! 1 minus epsilon
      real(r8) trinv                ! work variable
      real(r8) es                   ! sat. vapor pressure
      real(r8) desdt                ! change in sat vap pressure wrt temperature
c     real(r8) desdp                ! change in sat vap pressure wrt pressure
      real(r8) dqsdt                ! change in sat spec. hum. wrt temperature
      real(r8) dgdt                 ! work variable
      real(r8) g                    ! work variable
      real(r8) weight(plond)        ! work variable
      real(r8) hlatsb               ! (sublimation)
      real(r8) hlatvp               ! (vaporization)
      real(r8) hltalt(plond,plev)   ! lat. heat. of vap.
      real(r8) tterm                ! work var.
      real(r8) qs                   ! spec. hum. of water vapor
      real(r8) tc                   ! crit temp of transition to ice

      real(r8) t1, q1, dt, dq
      real(r8) dtm, dqm
      real(r8) qvd, a1, tmp
      real(r8) rair
      real(r8) r1b, c1, c2, c3
c-drb
      integer idt,idq
c-drb

#include <eslookup.h>

      omeps = 1.0 - epsqs
      trinv = 1.0/ttrice
      a1 = 7.5*log(10.)
      rair =  287.04
      c3 = rair*a1/cp
      dtm = 0.
      dqm = 0.

c number of times to iterate the calculation
c-drb      iter = 6
      iter = 8

c     dont even bother calculating the saturation point above ktop
      do k = 1,ktop
         do i = 1,plon
            tsp(i,k) = t(i,k)
            qsp(i,k) = 1.0
         end do
      end do

      do k = ktop,plev
c
c first guess on the saturation point
c
         do i = 1,plon

            es = estblf(t(i,k))
c
c Saturation specific humidity
c
            qs = min(epsqs*es/(p(i,k) - omeps*es),one)
c
            if (qs .lt. 0.0) then
               qsp(i,k) = 1.0
               tsp(i,k) = t(i,k)
            else
c
c "generalized" analytic expression for t derivative of es
c  accurate to within 1 percent for 173.16 < t < 373.16
c
c Weighting of hlat accounts for transition from water to ice
c polynomial expression approximates difference between es over
c water and es over ice from 0 to -ttrice (C) (min of ttrice is
c -40): required for accurate estimate of es derivative in transition 
c range from ice to water also accounting for change of hlatv with t 
c above 273.16 where const slope is given by -2369 j/(kg c) = cpv - cw
c
               tc     = t(i,k) - 273.16
               lflg   = (tc.ge.-ttrice .and. tc.lt.0.0)
               weight(i) = min(-tc*trinv,one)
               hlatsb = hlatv + weight(i)*hlatf
               hlatvp = hlatv - 2369.0*tc
               if (t(i,k).lt.273.16) then
                  hltalt(i,k) = hlatsb
               else
                  hltalt(i,k) = hlatvp
               end if

               tmp =  q(i,k) - qs
               c1 = hltalt(i,k)*c3
               c2 = (t(i,k) + 36.)**2
               r1b    = c2/(c2 + c1*qs)
               qvd   = r1b*tmp
               tsp(i,k) = t(i,k) + ((hltalt(i,k)/cp)*qvd)
               es = estblf(tsp(i,k))
               qsp(i,k) = min(epsqs*es/(p(i,k) - omeps*es),one)
            end if
         end do
c
c now iterate on first guess
c
         do l = 1, iter
            dtm = 0
            dqm = 0
            idt = 0
            idq = 0
            do i = 1,plon
               if ( qsp(i,k) .lt. 1.0 ) then
                  es = estblf(tsp(i,k))
c
c Saturation specific humidity
c
                  qs = min(epsqs*es/(p(i,k) - omeps*es),one)
c
c "generalized" analytic expression for t derivative of es
c accurate to within 1 percent for 173.16 < t < 373.16
c
c Weighting of hlat accounts for transition from water to ice
c polynomial expression approximates difference between es over
c water and es over ice from 0 to -ttrice (C) (min of ttrice is
c -40): required for accurate estimate of es derivative in transition 
c range from ice to water also accounting for change of hlatv with t 
c above 273.16 where const slope is given by -2369 j/(kg c) = cpv - cw
c
                  tc     = tsp(i,k) - 273.16
                  lflg   = (tc.ge.-ttrice .and. tc.lt.0.0)
                  weight(i) = min(-tc*trinv,one)
                  hlatsb = hlatv + weight(i)*hlatf
                  hlatvp = hlatv - 2369.0*tc
                  if (tsp(i,k).lt.273.16) then
                     hltalt(i,k) = hlatsb
                  else
                     hltalt(i,k) = hlatvp
                  end if
                  if (lflg) then
                     tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) 
     $                    + tc*(pcf(4) + tc*pcf(5))))
                  else
                     tterm = 0.0
                  end if
c
c The following check is to avoid the generation of negative
c values that can occur in the upper stratosphere and mesosphere
c
                  if (qs .lt. 0.0) then
                     qs = 1.0
                     es = p(i,k)
                     dqsdt = 0.0
                     qs = q(i,k)
                  else
                     desdt = hltalt(i,k)*es
     $                    /(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
                     dqsdt = (epsqs + omeps*qs)
     $                    /(p(i,k) - omeps*es)*desdt
                  endif
                  g = cp*(t(i,k)-tsp(i,k)) + hltalt(i,k)*q(i,k)
     $                 - hltalt(i,k)*qsp(i,k)
                  dgdt = -(cp + hltalt(i,k)*dqsdt)
                  t1 = tsp(i,k) - g/dgdt
                  dt = abs(t1 - tsp(i,k))/t1
                  tsp(i,k) = t1
                  es = estblf(tsp(i,k))
                  q1 = min(epsqs*es/(p(i,k) - omeps*es),one)
                  dq = abs(q1 - qsp(i,k))/max(q1,1.e-12_r8)
                  qsp(i,k) = q1

c                 write (6,*) ' rel chg lev, iter, t, q ', k, l, dt, dq
c-drb                  dtm = max(dtm,dt)
c-drb                  dqm = max(dqm,dq)
                  if (dt.gt.dtm) then
                     idt = i 
                     dtm = dt
                  endif
                  if (dq.gt.dqm) then
                     idq = i 
                     dqm = dq
                  endif
               end if
            end do              ! do i = 1,plon

            if (dtm.lt.1.e-4.and.dqm.lt.1.e-4) then
               go to 10
            endif

         end do                 ! do l = 1,iter
   10    continue



         if (iter.gt.0 .and. (dtm.gt.1.e-4.or.dqm.gt.1.e-4)
c++match/ncep hack
     $       .and. k.gt.2 ) then
c--
            write (6,*) ' findsp not converging ', dtm, dqm, k
c           write (6,*) 'q', q
c-drb
            write (6,*) 'worst t points are: ',k,idt,t(idt,k),tsp(idt,k),dtm
            write (6,*) 'worst q points are: ',k,idq,q(idq,k),qsp(idq,k),dqm
            call abort()
         endif

      end do                    ! level loop (k=1,plev)

      return
      end subroutine findsp

      end module cldwat
