c#include <misc.h>
#include <params.h>
      subroutine initu10 (gravitx, cpairx, cpwvx, rairx, zvirx)

c     initialize the common blocks for getu10
c     author: Phil Rasch


      use precision

      implicit none

#include <pmgrid.h>
      
      real(r8) gravitx
      real(r8) cpairx
      real(r8) cpwvx
      real(r8) rairx
      real(r8) zvirx

#include <u10com.h>

      gravit = gravitx
      cpair = cpairx
      cpvir  = cpwvx/cpairx - 1.
      rair = rairx
      zvir =  zvirx

      write (6,*) ' initu10 complete '

      return
      end

      subroutine getu10 (ioro, pmid, u, v, t, sh, zm,
     $     ts, u10)

c     calculate a 10m wind, given the meteorology in the model atmosphere
c     author: Phil Rasch
c     


      use precision

      implicit none

#include <pmgrid.h>

      integer ioro(plond) ! nint(oro)

      real(r8)
     $  pmid(plond,plev)       ! pressure at layer midpoints
     $, u(plond, plev) ! zonal wind
     $, v(plond, plev) ! meridional wind
     $, t(plond, plev) ! temperature
     $, sh(plond, plev) ! specific humidity
     $, zm(plond,plev) !  height at midpoints
     $, ts(plond)       ! surface temperature
     $, u10(plond)      ! wind speed at 10 meters

      real(r8) thbot(plond)
      integer indx(plond)
      integer i
      integer ii
      real(r8) p0

      p0 = 1.e5

c     find the points over ocean
c     initialize u10 and thbot at all points
      ii = 1
      do i = 1,plon
         if (ioro(i) .eq. 0) then
            indx(ii) = i
            ii = ii + 1
         endif
         u10(i) = 0.
         thbot(i) = t(i,plev)*(p0/pmid(i,plev))**0.286
      end do
      ii = ii-1

      call flxoce (indx, ii, pmid, u(1,plev), v(1,plev),
     $             t(1,plev), sh(1,plev), thbot, zm(1,plev), ts,
     $             u10)

      return
      end

      subroutine flxoce(indx    ,npts    ,pmidm1  ,ubot    ,vbot    ,
     $                  tbot    ,qbot    ,thbot   ,zbot    ,ts      ,
     $                  u10)
c     $                  ltheat  ,shf     ,lhf     ,taux    ,tauy    ,
c     $                  lwup    ,tref    )
C
c     code borrowed from the CSM
c     modified by Phil Rasch to return the 10m neutral wind
c
C Compute ocean to atmosphere surface fluxes of sensible, latent heat
C and stress components:
C
C Assume:
C   1) Neutral 10m drag coeff: 
C         cdn = .0027/U10N + .000142 + .0000764 U10N
C   2) Neutral 10m stanton number: 
C         ctn = .0327 sqrt(cdn), unstable
C         ctn = .0180 sqrt(cdn), stable
C   3) Neutral 10m dalton number:  
C         cen = .0346 sqrt(cdn)
C   4) The saturation humidity of air at T(K): 
C         qsat(T)  (kg/m^3)
C
C Note:
C   1) here, tstar = <WT>/U*, and qstar = <WQ>/U*.
C   2) wind speeds should all be above a minimum speed (umin)
C
C--------------------------Code History---------------------------------
C
C Original:      Bill Large/M.Vertenstein, Sep. 1995
C Standardized:  L. Buja,     Feb 1996
C Reviewed:      B. Briegleb, March 1996
C
C-----------------------------------------------------------------------

      use precision

#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
#include <parpbl.h>
C------------------------------Commons----------------------------------
#include <u10com.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer indx(plond) ! Longitude index array (land)
      integer npts        ! Number of land points
C                         
      real(r8) pmidm1(plond)  ! Bottom level pressure
      real(r8) ubot(plond)    ! Bottom level u wind
      real(r8) vbot(plond)    ! Bottom level v wind
      real(r8) tbot(plond)    ! Bottom level temperature
      real(r8) qbot(plond)    ! Bottom level specific humidity
      real(r8) thbot(plond)   ! Bottom level potential temperature
      real(r8) zbot(plond)    ! Bottom level height above surface
      real(r8) ts(plond)      ! Surface temperature
c     real(r8) ltheat(plond)  ! Latent heat for given srf conditions
C                         
C Output arguments        
C                         
      real(r8) u10(plond)      ! 10 meter wind assuming neutral BL
c      real(r8) shf(plond)     ! Initial sensible heat flux (W/m2)
c      real(r8) lhf(plond)     ! Initial latent heat flux (W/m2)
c      real(r8) taux(plond)    ! X surface stress (N/m2)
c      real(r8) tauy(plond)    ! Y surface stress (N/m2)
c      real(r8) lwup(plond)    ! Longwave up flux at surface (W/m2)
c      real(r8) tref(plond)    ! 2m reference temperature 
C
C---------------------------Local variables-----------------------------
C
      integer i,ii        ! Longitude indices
      real(r8) ssq            ! Surface saturation specific humidity
      real(r8) ustar          ! ustar
      real(r8) tstar          ! tstar
      real(r8) qstar          ! qstar
      real(r8) u10n           ! neutral 10 m wind speed over ocean
      real(r8) vmag           ! Surface wind magnitude
      real(r8) thvbot         ! Bottom lev virtual potential temp
      real(r8) delt           ! potential T difference (K)
      real(r8) delq           ! specific humidity difference (kg/kg)
      real(r8) rdn            ! sqrt of neutral exchange coeff (momentum)
      real(r8) rhn            ! sqrt of neutral exchange coeff (heat)
      real(r8) ren            ! sqrt of neutral exchange coeff (tracers)          
      real(r8) rd             ! sqrt of exchange coeff (momentum)
      real(r8) rh             ! sqrt of exchange coeff (heat)
      real(r8) re             ! sqrt of exchange coeff (tracers)
      real(r8) hol            ! Ref hgt (10m) / monin-obukhov length
      real(r8) xsq            ! Temporary variable
      real(r8) xqq            ! Temporary variable
      real(r8) alz            ! ln(zbot/z10)
      real(r8) cp             ! Specific heat of moist air
      real(r8) tau            ! Reference height stress
      real(r8) psimh          ! Stability funct at ref lev (momentum)
      real(r8) psixh          ! Stability funct at ref lev (heat & tracers) 
      real(r8) stable         ! Stability factor
      real(r8) rbot(plond)    ! Density at bottom model level
      real(r8) bn             ! exchange coef funct for interpolation
      real(r8) bh             ! exchange coef funct for interpolation
      real(r8) fac            ! interpolation factor
      real(r8) ln0            ! log factor for interpolation
      real(r8) ln3            ! log factor for interpolation
      real(r8) ztref          ! reference height for air temperature
      parameter (ztref=2.0)

      real(r8) :: half = 0.5
C
C--------------------------Statement functions--------------------------
C
      real(r8) psimhu         ! Unstable part of psimh
      real(r8) psixhu         ! Unstable part of psixh
      real(r8) qsat           ! Saturation specific humidty of air
      real(r8) cdn            ! Neutral drag coeff at bottom model level
      real(r8) xd             ! Dummy argument
      real(r8) Tk             ! Temperature (K)
      real(r8) Umps           ! Wind velocity (m/sec)
C
      qsat(Tk)   = 640380. / exp(5107.4/Tk)
      cdn(Umps)  = 0.0027 / Umps + .000142 + .0000764 * Umps
      psimhu(xd) = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.)
     $              - 2.*atan(xd) + 1.571
      psixhu(xd) = 2. * log((1. + xd*xd)/2.)
C
C-----------------------------------------------------------------------
C
      do ii=1,npts
         i = indx(ii)
C
C---------------------------------------------------------------
C Set up necessary variables
C---------------------------------------------------------------
C
         rbot(i)= pmidm1(i) / (rair*tbot(i))
         vmag   = max(umin, sqrt(ubot(i)**2+vbot(i)**2))
         thvbot = thbot(i) * (1.0 + zvir*qbot(i))
         ssq    = 0.98 * qsat(ts(i)) / rbot(i)
         delt   = thbot(i) - ts(i)
         delq   = qbot(i) - ssq  
         alz    = log(zbot(i)/zref) 
         cp     = cpair*(1. + cpvir*ssq) 
C
C---------------------------------------------------------------
C First iteration to converge on Z/L and hence the fluxes
C---------------------------------------------------------------
C
C Initial guess for roots of neutral exchange coefficients, 
C assume z/L=0. and u10n is approximated by vmag.
C Stable if (thbot > ts ).
C
         stable = half + sign(half , delt)
         rdn  = sqrt(cdn(vmag))
         rhn  = (1.-stable) * 0.0327 + stable * 0.018 
         ren  = 0.0346 
C
C Initial guess of ustar,tstar and qstar
C
         ustar = rdn*vmag
         tstar = rhn*delt
         qstar = ren*delq
C
C Compute stability and evaluate all stability functions
C Stable if (thbot > ts or hol > 0 )
C
         hol = xkar *gravit *zbot(i)
     $         * (tstar/thvbot + qstar/(1./zvir+qbot(i)))
     $         / ustar**2
         hol = sign( min(abs(hol),10._r8), hol )
         stable = half + sign(half , hol)
         xsq   = max(sqrt(abs(1. - 16.*hol)) , 1._r8)
         xqq   = sqrt(xsq)
         psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
         psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)
C
C Shift 10m neutral wind speed using old rdn coefficient
C
         rd   = rdn / (1.+rdn/xkar*(alz-psimh))
         u10n = vmag * rd / rdn
C
C Update the neutral transfer coefficients at 10m and neutral stability
C
         rdn = sqrt(cdn(u10n))
         ren = 0.0346
         rhn = (1.-stable) * 0.0327 + stable * 0.018 
C
C Shift all coeffs to measurement height and stability
C
         rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
         rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
         re = ren / (1.+ren/xkar*(alz-psixh))
C
C Update ustar, tstar, qstar using updated, shifted coeffs 
C
         ustar = rd * vmag 
         tstar = rh * delt 
         qstar = re * delq 
C
C---------------------------------------------------------------
C Second iteration to converge on Z/L and hence the fluxes
C---------------------------------------------------------------
C
C Recompute stability & evaluate all stability functions  
C Stable if (thbot > ts or hol > 0 )
C 
         hol = xkar * gravit * zbot(i)
     $        * (tstar/thvbot + qstar/(1./zvir+qbot(i)))
     $        / ustar**2
         hol = sign( min(abs(hol),10._r8), hol )
         stable = half + sign(half , hol)
         xsq   = max(sqrt(abs(1. - 16.*hol)) , 1._r8)
         xqq   = sqrt(xsq)
         psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
         psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)
C
C Shift 10m neutral wind speed using old rdn coefficient
C
         rd   = rdn / (1.+rdn/xkar*(alz-psimh))
         u10n = vmag * rd / rdn
c++pjr 
#ifdef V2
c        borrowed from temp interpolation at ztref
c        converted to zref
c        assume wind goes to zero at sfc
c
C---------------------------------------------------------------
C Following Geleyn(1988), interpolate ts to fixed height zref
C---------------------------------------------------------------
C
C Compute function of exchange coefficients. Assume that 
C cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
C 1/sqrt(cn(i))=1/rdn and sqrt(cm(i))/ch(i)=1/rh 
C
         bn = xkar/rdn
         bh = xkar/rh
C
C Interpolation factor for stable and unstable cases
C
         ln0 = log(1.0 + (zref/zbot(i))*(exp(bn) - 1.0))
         ln3 = log(1.0 + (zref/zbot(i))*(exp(bn - bh) - 1.0))
         fac = (ln0 - zref/zbot(i)*(bn - bh))/bh * stable
     $       + (ln0 - ln3)/bh * (1.-stable)
         fac = min(max(fac,0.),1.)
C
C Actual interpolation
C
c
         u10(i) = 0. + vmag*fac
#else
         u10(i) = u10n
#endif
c--pjr
#ifdef TURNEDOFF
C
C Update the neutral transfer coefficients at 10m and neutral stability
C
         rdn = sqrt(cdn(u10n))
         ren = 0.0346
         rhn = (1.-stable) * 0.0327 + stable * 0.018 
C
C Shift all coeffs to measurement height and stability

         rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
         rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
         re = ren / (1.+ren/xkar*(alz-psixh))
C
C---------------------------------------------------------------
C Compute the fluxes
C---------------------------------------------------------------
C
C Update ustar, tstar, qstar using updated, shifted coeffs 
C
         ustar = rd * vmag 
         tstar = rh * delt 
         qstar = re * delq 
C
C Compute surface stress components
C
         tau     =  rbot(i) * ustar * ustar 
         taux(i) = -tau * ubot(i) / vmag 
         tauy(i) = -tau * vbot(i) / vmag 
C
C Compute heat flux components at current surface temperature
C (Define positive latent and sensible heat as upwards into the atm)
C
         shf(i) = -cp * tau * tstar / ustar 
         lhf(i) = -ltheat(i) * tau * qstar / ustar
         lwup(i) = stebol * ts(i)**4 
C
C---------------------------------------------------------------
C Following Geleyn(1988), interpolate ts to fixed height zref
C---------------------------------------------------------------
C
C Compute function of exchange coefficients. Assume that 
C cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
C 1/sqrt(cn(i))=1/rdn and sqrt(cm(i))/ch(i)=1/rh 
C
         bn = xkar/rdn
         bh = xkar/rh
C
C Interpolation factor for stable and unstable cases
C
         ln0 = log(1.0 + (ztref/zbot(i))*(exp(bn) - 1.0))
         ln3 = log(1.0 + (ztref/zbot(i))*(exp(bn - bh) - 1.0))
         fac = (ln0 - ztref/zbot(i)*(bn - bh))/bh * stable
     $       + (ln0 - ln3)/bh * (1.-stable)
         fac = min(max(fac,0.),1.)
C
C Actual interpolation
C
         tref(i) = ts(i) + (tbot(i) - ts(i))*fac
#endif
C
      end do
C
      return
      end

 
