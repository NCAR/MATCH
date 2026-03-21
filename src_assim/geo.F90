MODULE geo

use floats

PUBLIC :: EarthDist, InsideLons

CONTAINS

    LOGICAL FUNCTION InsideLons( Lon1, Lon2, Lon )
!+
!
!-
      use precision
  IMPLICIT NONE
!
! INPUTS:
      REAL(r8), INTENT(IN) ::  Lon1, Lon2, Lon
!
! INTERNALS

!
      IF ( Lon1 .GT. Lon2 ) THEN
         IF ( ((Lon .GE. 0)    .AND. (Lon .LE. Lon2)) .OR. & 
              ((Lon .GE. Lon1) .AND. (Lon .LE. 360.0) ) ) THEN
           InsideLons = .True.
         ELSE
           InsideLons = .False.
         ENDIF  
      ELSE
         IF ( (Lon .GE. Lon1) .AND. (Lon .LE. Lon2) ) THEN
            InsideLons = .True. 
         ELSE
           InsideLons = .False.
         ENDIF
      ENDIF   
!
    RETURN
    END FUNCTION InsideLons


      
    FUNCTION EarthDist( Lat1,Lon1, Lat2,Lon2 )
!+
! Computes distance in km between two points on
! the Earth surface
!-
        use precision
      IMPLICIT NONE
!
! INPUTS:
      REAL(r8) Lat1, Lat2, Lon1, Lon2
!
! INTERNALS
      REAL(r8) R, Pi
      PARAMETER (R = 6359., Pi = 3.14159045)
      REAL(r8) Lat1_, Lat2_, Lon1_, Lon2_, dLon, dLat, F, Dist, X

      REAL(r8) EarthDist

!
      F = Pi/180.0
      Lat1_ = Lat1 * F
      Lat2_ = Lat2 * F
      Lon1_ = Lon1 * F
      Lon2_ = Lon2 * F
      dLon =  Lon1_ - Lon2_
      dLat =  Lat1_ - Lat2_

      IF ((ABS(dLat) .LE. TINY(0.0)) .AND. &
          (ABS(dLon) .LE. TINY(0.0))) THEN
         EarthDist = 0.0
         RETURN
      ENDIF
   
      X = SIN(Lat1_)*SIN(Lat2_) + COS(Lat1_) * COS(Lat2_) * COS(dLon)

      IF (X .GT.  1.0) THEN
!         PRINT *, 'X, Lat1,Lat2,Lon1,Lon2: ', X, Lat1,Lat2,Lon1,Lon2
         X =  1.0
      ENDIF   
      IF (X .LT. -1.0) X = -1.0

      Dist = R * ACOS( X )

      EarthDist = Dist
!
      RETURN
      END FUNCTION EarthDist



      
      function EarthDist1  ( y1g, x1g, y2g, x2g )
!
!***   compute distance of two points on the globe.
!*
!***   this formula is correct, also for points close together.
!*
!***   x1g, y1g : longitude and latitude of first point  (degrees)
!***   x2g, y2g : longitude and latitude of second point (degrees)
!***   ddg      : distance between these points          (degrees)
!*
      use precision
      implicit  none

!*
      REAL(r8)  x1g, x2g, y1g, y2g, ddg
      REAL(r8)  x1, x2, y1, y2, dd
      REAL(r8)  pi, cc, dx, dy, dx2, dy2, dd2, rearth

      REAL(r8) EarthDist1
!*
      data  pi  / 3.141593 /
      data  rearth  / 6378.5 /
!*
      cc = pi / 180.0
!*
      x1 = x1g * cc
      x2 = x2g * cc
      y1 = y1g * cc
      y2 = y2g * cc
!*
      dy = 0.5 * ( y2 - y1 )
      dy = sin(dy)
      dy2 = dy * dy
!*
      dx = 0.5 * ( x2 - x1 )
      dx = sin(dx)
      dx2 = dx * dx * cos(y1) * cos(y2)
!*
      dd2 = dx2 + dy2
      dd = sqrt(dd2)
      dd = 2.0 * asin(dd)
!*
      ddg = dd / cc
      ddg = (ddg * 2 * pi * rearth) / 360  
!*
      EarthDist1 = ddg
!*
      return
!*
      end FUNCTION  EarthDist1

END MODULE geo
