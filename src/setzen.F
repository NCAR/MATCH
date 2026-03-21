c $Id$
      subroutine setzen(clat, pi, calday, plond, plon, coszrs)

c setzen computes coszrs by the same algorithm as used in radinp (CCM2).


      use precision

      implicit none
c-----------------------------------------------------------------------

      real(r8) clat                 ! current latitude
      real(r8) pi                   ! 3.14159...
      real(r8) calday               ! time at Greenwich
      integer plond
      integer plon
      real(r8) coszrs(plond)
C---------------------------local variables-----------------------------
C
      integer i
c
      real(r8) delta                ! 
      real(r8) theta                ! 
      real(r8) cosd                 ! 
      real(r8) cosc                 ! 
      real(r8) sind                 ! 
      real(r8) sinc                 ! 
      real(r8) loctim               ! 
      real(r8) dayspy               ! 
      parameter( dayspy=365. )
C-----------------------------------------------------------------------
C
C
C     !  Set cos(zenith angle) for diurnal or daily avg insolation.
C     !  Determine values for frac and coszrs.
C     !  Replace 'radzen' for diurnal cycle.
C
C     !  compute three solar related quantities necessary to compute inc
C     !  solar radiation in 'radcsw':
C
C     !  (5) coszrs = longitude, latitude, and time dependent cosine sol
C     !  zenith angle: positive for sun above horizon, and negative for
C     !  sun below horizon.
C
C     !  (6) frac   = fraction of solar day (useful only when no diurnal
C     !  cycle is assumed; set to 1 when diurnal cycle is assumed).
C
C     !  the computation of these three quantities depends only
C     !  on the day of the year, so that their values repeat
C     !  from one year to the next during model run.
C
C     !  compute solar distance factor and cosine solar zenith angle usi
C     !  day value where a round day (such as 213.0) refers to 0z at
C     !  greenwich longitude.
C
C     !  use formulas from paltridge and platt 1976  p. 57, p. 62,63.
c
c      pi = 4.*atan(1.0)  ! get pi from input arg. so that we use the same
c                         ! value as used in other parts of code.
      theta = 2.*pi*calday/dayspy
C
C     !  Solar declination in radians:
C
      delta = .006918 - .399912*cos(theta) + .070257*sin(theta) -
     $        .006758*cos(2.*theta) + .000907*sin(2.*theta) -
     $        .002697*cos(3.*theta) + .001480*sin(3.*theta)
C
C     !  Compute local time based on longitude and day of year.
C
      sinc = sin(clat)
      sind = sin(delta)
      cosc = cos(clat)
      cosd = cos(delta)
c
      do 10 i=1,plon
         loctim = calday + (real(i-1)/real(plon))
C
C     !  Compute local cosine solar zenith angle.
C
           coszrs(i) = sinc*sind - cosc*cosd*cos(2.*pi*loctim)
   10 continue
C
      return
      end
