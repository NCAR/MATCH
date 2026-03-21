c $Id$
#include <params.h>

      subroutine esinti(epslon  ,latvap  ,latice  ,rh2o    ,cpair, xip )
C-----------------------------------------------------------------------
C
C Initialize es lookup tables 
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Hack
C Standardized:      L. Buja, June 1992
C Reviewed:          J. Hack, G. Taylor, August 1992
C
C-----------------------------------------------------------------------
      use precision

      implicit none
#if ( defined RS6000 )
      implicit automatic (a-z)
#endif
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) epslon,         ! Ratio of h2o to dry air molecular weights 
     $     latvap,         ! Latent heat of vaporization
     $     latice,         ! Latent heat of fusion
     $     rh2o,           ! Gas constant for water vapor
     $     cpair           ! Specific heat of dry air
c++bee
      logical xip          ! Ice phase (true or false)
c--bee
C
C---------------------------Local workspace-----------------------------
C
      real(r8) tmn,            ! Minimum temperature entry in table
     $     tmx,            ! Maximum temperature entry in table
     $     trice           ! Trans range from es over h2o to es over ice
      logical ip           ! Ice phase (true or false)
      external gestbl
C
C-----------------------------------------------------------------------
C
C Specify control parameters first
C
      tmn   = 173.16
      tmx   = 375.16
      trice =  20.00
c++bee
      ip    = xip
c--bee
C
C Call gestbl to build saturation vapor pressure table.
C
      call gestbl(tmn     ,tmx     ,trice   ,ip      ,epslon  ,
     $            latvap  ,latice  ,rh2o    ,cpair   )
C
      return
      end

c##############################################################################

      subroutine aqsat(t       ,p       ,es      ,qs        ,ii      ,
     $                 ilen    ,kk      ,kstart  ,kend      )
C-----------------------------------------------------------------------
C
C Utility procedure to look up and return saturation vapor pressure from 
C precomputed table, calculate and return saturation specific humidity 
C (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
C This routine is useful for evaluating only a selected region in the 
C vertical.
C 
C---------------------------Code history--------------------------------
C
C Original version: J. J. Hack, February 1990
C Reviewed:         J. J. Hack, August 1992 
C
C-----------------------------------------------------------------------
      use precision

      implicit none
#if ( defined RS6000 )
      implicit automatic (a-z)
#endif
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer ii,            ! I dimension of arrays t, p, es, qs
     $        kk             ! K dimension of arrays t, p, es, qs 
      real(r8) t(ii,kk),         ! Temperature
     $     p(ii,kk)          ! Pressure
      integer ilen,          ! length of vectors in I direction
C                            ! which are assumed to start at 1
     $        kstart,        ! starting location in K direction
     $        kend           ! ending location in K direction
C
C Output arguments
C
      real(r8) es(ii,kk),        ! Saturation vapor pressure
     $     qs(ii,kk)         ! Saturation specific humidity
C
C---------------------------Local workspace-----------------------------
C
      real(r8) omeps             ! 1 - 0.622
      integer i,             ! i index
     $        k              ! k index
C
C-----------------------------------------------------------------------
#include <eslookup.h>
C-----------------------------------------------------------------------
C
      omeps = 1.0 - epsqs
      do k=kstart,kend
         do i=1,ilen
            es(i,k) = estblf(t(i,k))
C
C Saturation specific humidity
C
            qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
C   
C The following check is to avoid the generation of negative values
C that can occur in the upper stratosphere and mesosphere
C
            qs(i,k) = min(1._r8,qs(i,k))
C
            if (qs(i,k) .lt. 0.0) then
               qs(i,k) = 1.0
               es(i,k) = p(i,k)
            end if
         end do
      end do
C
      return
      end

c##############################################################################

      subroutine aqsatd(t       ,p       ,es      ,qs      ,gam     ,
     $                  ii      ,ilen    ,kk      ,kstart  ,kend    )
C-----------------------------------------------------------------------
C
C Utility procedure to look up and return saturation vapor pressure from 
C precomputed table, calculate and return saturation specific humidity 
C (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT) for input 
C arrays of temperature and pressure (dimensioned ii,kk).
C
C---------------------------Code history--------------------------------
C
C Original version: J. J. Hack, February 1990
C Reviewed:         J. J. Hack, August 1992 
C
C-----------------------------------------------------------------------
      use precision

      implicit none

#if ( defined RS6000 )
      implicit automatic (a-z)
#endif
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer ii,           ! I dimension of arrays t, p, es, qs
     $        kk            ! K dimension of arrays t, p, es, qs 
      real(r8) t(ii,kk),        ! Temperature
     $     p(ii,kk)         ! Pressure
      integer ilen,         ! vector length in I direction
     $        kstart,       ! starting location in K direction
     $        kend          ! ending location in K direction
C 
C Output arguments
C
      real(r8) es(ii,kk),       ! Saturation vapor pressure
     $     qs(ii,kk),       ! Saturation specific humidity
     $     gam(ii,kk)       ! (l/cp)*(d(qs)/dt)
C
C---------------------------Local workspace-----------------------------
C
      logical lflg          ! true if in temperature transition region
      integer i,            ! i index for vector calculations
     $        k             ! k index 
      real(r8) omeps,           ! 1. - 0.622
     $     trinv,           ! reciprocal of ttrice (transition range)
     $     tc,              ! temperature (in degrees C)
     $     weight,          ! weight for es transition from water to ice
     $     hltalt,          ! appropriately modified hlat for T derivatives
     $     hlatsb,          ! hlat weighted in transition region
     $     hlatvp,          ! hlat modified for t changes above 273.16
     $     tterm,           ! account for d(es)/dT in transition region
     $     desdt            ! d(es)/dT
C-----------------------------------------------------------------------
#include <eslookup.h>
C-----------------------------------------------------------------------
      omeps = 1.0 - epsqs
      do k=kstart,kend
         do i=1,ilen
            es(i,k) = estblf(t(i,k))
C
C Saturation specific humidity
C
            qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
C
C The following check is to avoid the generation of negative
C values that can occur in the upper stratosphere and mesosphere
C
            qs(i,k) = min(1._r8,qs(i,k))
C
            if (qs(i,k) .lt. 0.0) then
               qs(i,k) = 1.0
               es(i,k) = p(i,k)
            end if
         end do
      end do
C
C "generalized" analytic expression for t derivative of es
C accurate to within 1 percent for 173.16 < t < 373.16
C
      trinv = 0.0
      if ((.not. icephs) .or. (ttrice.eq.0.0)) go to 10
      trinv = 1.0/ttrice
C
      do k=kstart,kend
         do i=1,ilen
C
C Weighting of hlat accounts for transition from water to ice
C polynomial expression approximates difference between es over
C water and es over ice from 0 to -ttrice (C) (min of ttrice is
C -40): required for accurate estimate of es derivative in transition 
C range from ice to water also accounting for change of hlatv with t 
C above 273.16 where constant slope is given by -2369 j/(kg c) =cpv - cw
C
            tc     = t(i,k) - 273.16
            lflg   = (tc.ge.-ttrice .and. tc.lt.0.0)
            weight = min(-tc*trinv,1._r8)
            hlatsb = hlatv + weight*hlatf
            hlatvp = hlatv - 2369.0*tc
            if (t(i,k).lt.273.16) then
               hltalt = hlatsb
            else
               hltalt = hlatvp
            end if
            if (lflg) then
               tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) +
     $                 tc*pcf(5))))
            else
               tterm = 0.0
            end if
            desdt  = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k)) + tterm*trinv
            gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/
     $               (cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
            if(qs(i,k).eq.1.0) gam(i,k) = 0.0
         end do
      end do
C
      go to 20
C
C No icephs or water to ice transition
C
   10 do k=kstart,kend
         do i=1,ilen
C
C Account for change of hlatv with t above 273.16 where
C constant slope is given by -2369 j/(kg c) = cpv - cw
C
            hlatvp = hlatv - 2369.0*(t(i,k)-273.16)
            if (icephs) then
               hlatsb = hlatv + hlatf
            else
               hlatsb = hlatv
            end if
            if (t(i,k).lt.273.16) then
               hltalt = hlatsb
            else
               hltalt = hlatvp
            end if
            desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k))
            gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/
     $                 (cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
            if (qs(i,k) .eq. 1.0) gam(i,k) = 0.0
         end do
      end do
C
   20 return
      end

c##############################################################################

      subroutine gestbl(tmn     ,tmx     ,trice   ,ip      ,epsil   ,
     $                  latvap  ,latice  ,rh2o    ,cpair   )
C-----------------------------------------------------------------------
C
C Builds saturation vapor pressure table for later lookup procedure.
C Uses Goff & Gratch (1946) relationships to generate the table
C according to a set of free parameters defined below.  Auxiliary
C routines are also included for making rapid estimates (well with 1%)
C of both es and d(es)/dt for the particular table configuration.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Hack
C Standardized:      L. Buja, June 1992
C Reviewed:          J. Hack, G. Taylor, August 1992
C
C-----------------------------------------------------------------------
      use precision

      implicit none
#if ( defined RS6000 )
      implicit automatic (a-z)
#endif
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) tmn,          ! Minimum temperature entry in es lookup table
     $     tmx,          ! Maximum temperature entry in es lookup table
     $     epsil,        ! Ratio of h2o to dry air molecular weights
     $     trice,        ! Transition range from es over range to es over ice
     $     latvap,       ! Latent heat of vaporization
     $     latice,       ! Latent heat of fusion
     $     rh2o,         ! Gas constant for water vapor
     $     cpair         ! Specific heat of dry air
C
C---------------------------Local variables-----------------------------
C
      real(r8) t             ! Temperature
      integer n,         ! Increment counter
     $        lentbl,    ! Calculated length of lookup table
     $        itype      ! Ice phase: 0 -> no ice phase
                         !            1 -> ice phase, no transition
                         !           -x -> ice phase, x degree transition
      logical ip         ! Ice phase logical flag
C
C------------------------------Externals--------------------------------
C
      external endrun    ! Abnormal termination
C
C---------------------------Statement function--------------------------
C
#include <eslookup.h>
C-----------------------------------------------------------------------
C
C Set es table parameters
C
      tmin   = tmn       ! Minimum temperature entry in table
      tmax   = tmx       ! Maximum temperature entry in table
      ttrice = trice     ! Trans. range from es over h2o to es over ice
      icephs = ip        ! Ice phase (true or false)
C
C Set physical constants required for es calculation
C
      epsqs  = epsil
      hlatv  = latvap
      hlatf  = latice
      rgasv  = rh2o
      cp     = cpair
C
c-drb      lentbl = ifix(tmax-tmin+2.000001)
      lentbl = int(tmax-tmin+2.000001)
      if (lentbl .gt. plenest) then
         write(6,9000) tmax, tmin, plenest
         call endrun
      end if
C
C Begin building es table.
C Check whether ice phase requested.
C If so, set appropriate transition range for temperature
C
      if (icephs) then
         if(ttrice.ne.0.0) then
            itype = -ttrice
         else
            itype = 1
         end if
      else
         itype = 0
      end if
C
      t = tmin - 1.0
      do n=1,lentbl
         t = t + 1.0
         call gffgch(t,estbl(n),itype)
      end do
C
      do n=lentbl+1,plenest
         estbl(n) = -99999.0
      end do
C
C Table complete -- Set coefficients for polynomial approximation of
C difference between saturation vapor press over water and saturation
C pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
C is valid in the range -40 < t < 0 (degrees C).
C
C                  --- Degree 5 approximation ---
C
      pcf(1) =  5.04469588506e-01
      pcf(2) = -5.47288442819e+00
      pcf(3) = -3.67471858735e-01
      pcf(4) = -8.95963532403e-03
      pcf(5) = -7.78053686625e-05
C
C                  --- Degree 6 approximation ---
C
C-----pcf(1) =  7.63285250063e-02
C-----pcf(2) = -5.86048427932e+00
C-----pcf(3) = -4.38660831780e-01
C-----pcf(4) = -1.37898276415e-02
C-----pcf(5) = -2.14444472424e-04
C-----pcf(6) = -1.36639103771e-06
C
      write(6,*)' ***** SATURATION VAPOR PRESSURE TABLE COMPLETED *****'
      return
C
 9000 format('GESTBL: FATAL ERROR *********************************',/,
     $     ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH',
     $     ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/,
     $     ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
C 
      end

c##############################################################################

      subroutine vqsatd(t       ,p       ,es      ,qs      ,gam      ,
     $                  len     )
C-----------------------------------------------------------------------
C
C Utility procedure to look up and return saturation vapor pressure from 
C precomputed table, calculate and return saturation specific humidity 
C (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
C function as qsatd, but operates on vectors of temperature and pressure
C
C----------------------------Code History-------------------------------
C
C Original version:  J. Hack
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, August 1992
C
C-----------------------------------------------------------------------
      use precision

      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     len       ! vector length
      real(r8) t(len),   ! temperature
     $     p(len)    ! pressure
C 
C Output arguments
C
      real(r8) es(len),  ! saturation vapor pressure
     $     qs(len),  ! saturation specific humidity
     $     gam(len)  ! (l/cp)*(d(qs)/dt)
C
C---------------------------Local workspace-----------------------------
C
      logical lflg   ! true if in temperature transition region
      integer
     $     i         ! index for vector calculations
      real(r8) omeps,    ! 1. - 0.622
     $     trinv,    ! reciprocal of ttrice (transition range)
     $     tc,       ! temperature (in degrees C)
     $     weight,   ! weight for es transition from water to ice
     $     hltalt,   ! appropriately modified hlat for T derivatives  
     $     hlatsb,   ! hlat weighted in transition region
     $     hlatvp,   ! hlat modified for t changes above 273.16
     $     tterm,    ! account for d(es)/dT in transition region
     $     desdt     ! d(es)/dT
C-----------------------------------------------------------------------
#include <eslookup.h>
C-----------------------------------------------------------------------
      omeps = 1.0 - epsqs
      do i=1,len
         es(i) = estblf(t(i))
C
C Saturation specific humidity
C
         qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
C
C The following check is to avoid the generation of negative
C values that can occur in the upper stratosphere and mesosphere
C
         qs(i) = min(1._r8,qs(i))
C
         if (qs(i) .lt. 0.0) then
            qs(i) = 1.0
            es(i) = p(i)
         end if
      end do
C
C "generalized" analytic expression for t derivative of es
C accurate to within 1 percent for 173.16 < t < 373.16
C
      trinv = 0.0
      if ((.not. icephs) .or. (ttrice.eq.0.0)) go to 10
      trinv = 1.0/ttrice
      do i=1,len
C
C Weighting of hlat accounts for transition from water to ice
C polynomial expression approximates difference between es over
C water and es over ice from 0 to -ttrice (C) (min of ttrice is
C -40): required for accurate estimate of es derivative in transition 
C range from ice to water also accounting for change of hlatv with t 
C above 273.16 where const slope is given by -2369 j/(kg c) = cpv - cw
C
         tc     = t(i) - 273.16
         lflg   = (tc.ge.-ttrice .and. tc.lt.0.0)
         weight = min(-tc*trinv,1._r8)
         hlatsb = hlatv + weight*hlatf
         hlatvp = hlatv - 2369.0*tc
         if (t(i).lt.273.16) then
            hltalt = hlatsb
         else
            hltalt = hlatvp
         end if
         if (lflg) then
            tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) +
     $              tc*pcf(5))))
         else
            tterm = 0.0
         end if
         desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
         gam(i) = hltalt*qs(i)*p(i)*desdt/
     $           (cp*es(i)*(p(i) - omeps*es(i)))
         if(qs(i).eq.1.0) gam(i) = 0.0
      end do
      return
C
C No icephs or water to ice transition
C
   10 do i=1,len
C
C Account for change of hlatv with t above 273.16 where
C constant slope is given by -2369 j/(kg c) = cpv - cw
C
         hlatvp = hlatv - 2369.0*(t(i)-273.16)
         if (icephs) then
            hlatsb = hlatv + hlatf
         else
            hlatsb = hlatv
         end if
         if (t(i).lt.273.16) then
            hltalt = hlatsb
         else
            hltalt = hlatvp
         end if
         desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
         gam(i) = hltalt*qs(i)*p(i)*desdt/
     $            (cp*es(i)*(p(i) - omeps*es(i)))
         if (qs(i) .eq. 1.0) gam(i) = 0.0
      end do
C
      return
      end

c##############################################################################

      subroutine gffgch(t       ,es      ,itype   )
C-----------------------------------------------------------------------
C
C Computes saturation vapor pressure over water and/or over ice using
C Goff & Gratch (1946) relationships.  T (temperature), and itype are
C input parameters, while es (saturation vapor pressure) is an output
C parameter.  The input parameter itype serves two purposes: a value of
C zero indicates that saturation vapor pressures over water are to be
C returned (regardless of temperature), while a value of one indicates
C that saturation vapor pressures over ice should be returned when t is
C less than 273.16 degrees k.  If itype is negative, its absolute value
C is interpreted to define a temperature transition region below 273.16
C degrees k in which the returned saturation vapor pressure is a
C weighted average of the respective ice and water value.  That is, in
C the temperature range 0 => -itype degrees c, the saturation vapor
C pressures are assumed to be a weighted average of the vapor pressure
C over supercooled water and ice (all water at 0 c; all ice at -itype
C c).  Maximum transition range => 40 c
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Hack
C Standardized:      L. Buja, June 1992
C Reviewed:          J. Hack, G. Taylor, August 1992
C
C-----------------------------------------------------------------------
      use precision

      implicit none
#if ( defined RS6000 )
      implicit automatic (a-z)
#endif
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) t          ! Temperature
      integer itype   ! Flag for ice phase and associated transition
C
C Output arguments
C
      real(r8) es         ! Saturation vapor pressure
C
C---------------------------Local variables-----------------------------
C
      real(r8) e1,        ! Intermediate scratch variable for es over water
     $     e2,        ! Intermediate scratch variable for es over water
     $     eswtr,     ! Saturation vapor pressure over water
     $     f,         ! Intermediate scratch variable for es over water
     $     f1,        ! Intermediate scratch variable for es over water
     $     f2,        ! Intermediate scratch variable for es over water
     $     f3,        ! Intermediate scratch variable for es over water
     $     f4,        ! Intermediate scratch variable for es over water
     $     f5,        ! Intermediate scratch variable for es over water
     $     ps,        ! Reference pressure (mb)
     $     t0,        ! Reference temperature (freezing point of water)
     $     term1,     ! Intermediate scratch variable for es over ice
     $     term2,     ! Intermediate scratch variable for es over ice
     $     term3,     ! Intermediate scratch variable for es over ice
     $     tr,        ! Transition range for es over water to es over ice
     $     ts,        ! Reference temperature (boiling point of water)
     $     weight     ! Intermediate scratch variable for es transition
      integer itypo   ! Intermediate scratch variable for holding itype
C
C------------------------------Externals--------------------------------
C
      external endrun ! Abnormal termination
C
C-----------------------------------------------------------------------
C      
C Check on whether there is to be a transition region for es
C
      if (itype.lt.0) then
         tr    = abs(float(itype))
         itypo = itype
         itype = 1
      else
         tr    = 0.0
         itypo = itype
      end if
      if (tr .gt. 40.0) then
         write(6,900) tr
        call endrun
      end if
C
      if(t .lt. (273.16 - tr) .and. itype.eq.1) go to 10
C
C Water
C
      ps = 1013.246
      ts = 373.16
      e1 = 11.344*(1.0 - t/ts)
      e2 = -3.49149*(ts/t - 1.0)
      f1 = -7.90298*(ts/t - 1.0)
      f2 = 5.02808*log10(ts/t)
      f3 = -1.3816*(10.0**e1 - 1.0)/10000000.0
      f4 = 8.1328*(10.0**e2 - 1.0)/1000.0
      f5 = log10(ps)
      f  = f1 + f2 + f3 + f4 + f5
      es = (10.0**f)*100.0
      eswtr = es
C
      if(t.ge.273.16 .or. itype.eq.0) go to 20
C
C Ice
C
   10 continue
      t0    = 273.16
      term1 = 2.01889049/(t0/t)
      term2 = 3.56654*log(t0/t)
      term3 = 20.947031*(t0/t)
      es    = 575.185606e10*exp(-(term1 + term2 + term3))
C
      if (t.lt.(273.16 - tr)) go to 20
C
C Weighted transition between water and ice
C
      weight = min((273.16 - t)/tr,1._r8)
      es = weight*es + (1.0 - weight)*eswtr
C
   20 continue
      itype = itypo
      return
C
  900 format('GFFGCH: FATAL ERROR ******************************',/,
     $       'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR',
     $       ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF',
     $       ' 40.0 DEGREES C',/, ' TR = ',f7.2)
C
      end
 

c##############################################################################

      subroutine limsh( ps, t, q )

c Limit maximum specific humidity to the saturation value.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  ps(plond,plat)       ! surface pressure
     $, t(plond,plev,plat)   ! temperature at model levels

c Input/Output arguments:
      real(r8)
     $  q(plond,plev,plat)  ! specific humidity

c Local variables:
      integer i, j, k
      real(r8)
     $  pint(plond,plevp)    ! Pressure at model interfaces
     $, pmid(plond,plev)     ! Pressure at model levels
     $, pdel(plond,plev)     ! Layer thickness (pint(k+1) - pint(k))
     $, es(plond,plev)       ! saturation vapor pressure
     $, qs(plond,plev)       ! saturation specific humidity
c-----------------------------------------------------------------------

      do j = 1, plat
         call plevs0( ps(1,j), pint, pmid, pdel )
         call aqsat( t(1,1,j), pmid, es, qs, plond, plon, plev, 1, plev)
         do k = 1, plev
            do i = 1, plon
               q(i,k,j) = min( q(i,k,j), qs(i,k) )
            end do
         end do
      end do

      return
      end

c##############################################################################

      subroutine relhum( ps, t, q, rh )

c Calculate relative humidity for a latitude slice.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  ps(plond)       ! surface pressure
     $, t(plond,plev)   ! temperature at model levels
     $, q(plond,plev)   ! specific humidity

c Output argument:
      real(r8)
     $  rh(plond,plev)  ! specific humidity

c Local variables:
      integer i, j, k
      real(r8)
     $  pint(plond,plevp)    ! Pressure at model interfaces
     $, pmid(plond,plev)     ! Pressure at model levels
     $, pdel(plond,plev)     ! Layer thickness (pint(k+1) - pint(k))
     $, es(plond,plev)       ! saturation vapor pressure
     $, qs(plond,plev)       ! saturation specific humidity
c-----------------------------------------------------------------------

      call plevs0( ps, pint, pmid, pdel )
      call aqsat( t, pmid, es, qs, plond, plon, plev, 1, plev)
      do k = 1, plev
         do i = 1, plon
            rh(i,k) = q(i,k)/qs(i,k)
         end do
      end do

      return
      end

c##############################################################################

