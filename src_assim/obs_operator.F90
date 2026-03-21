MODULE obs_operator

use mathlib
use avg_kernel
use f2v
use geo
use precision

IMPLICIT NONE

PRIVATE 
PUBLIC ::  Obs_Operator_XY, Obs_Operator_Z


CONTAINS  

 SUBROUTINE Obs_Operator_XY ( LatObs, LonObs, LatGrd, LonGrd, H, nH )
!+
!
! Computes elements of the observational operator -- matrix (in general
! NLat*NLon x NObs) that interpolates from locations of the model
! grid points to the locations of observations. Since we are on a plane
! only three model grid points are needed to define a plane and interpolate
! locally. Therefore H will only have 3 non-zero elements in each raw.
! Thus H is defined as (3,NObs). nH(3,NObs) contains the indices of the
! non-zero elements. 
!
! H acts on a vector of values defined at the model grid points (the
! vector length is NLat*NLon) and produces a vector of values interpolated
! to the observations locations (the vector length is NObs).
!
! We assume that the vector of observations is arranged in the same
! order as that in LatObs, LonObs. 
!
! We also assume that the model values are arranged in a vector as series
! of swaps along meridians. That is latitudes change the fastest. 
! Since observations follow the satellite track and are more or less 
! aligned along meridians, such arrangment should
! results in sparse matrices containing compact "blocks" of non-zero
! elements. Unless, of course, I totally screwed up picturing it all.  
!
! This arrangment is controlled by the external function 
!                      IndexF2V(NLat,NLon, iLat,iLon) 
! which returns the position of an element (iLat,iLon) in the vector
! of model values. The inverse of this is a subroutine 
!                      IndexV2F(NLat,NLon,k,iLat,iLon)
! that computes iLat,iLon given the position k in the vector.  
!
!-
      IMPLICIT NONE
! INPUTS:
      real(r8), INTENT(IN) :: LatObs(:), LonObs(:)
      real(r8), INTENT(IN) :: LatGrd(:), LonGrd(:)
!
! OUTPUTS:
      real(r8), INTENT(OUT) :: H(3, SIZE(LatObs) )
      INTEGER, INTENT(OUT) :: nH(3, SIZE(LatObs) )
!
! INTERNALS: 
      INTEGER :: NObs, NLat, NLon
!      INTEGER :: IndexF2V
      INTEGER :: iObs, i,j, ii1,ii2, jj1,jj2
      INTEGER :: k1, k2, k3, iLat1, iLat2, iLat3, jLon1, jLon2, jLon3
      real(r8) :: Lat1, Lat2, Lat3, Lon1, Lon2, Lon3
      real(r8) :: dR1, dR2, dR3, dR4, alpha, beta, gamma!, EarthDist
      LOGICAL :: LonGap
!
      IF ( SIZE(LatObs) .NE. SIZE(LonObs) ) THEN
         PRINT *, 'Obs_Operator_XY% Incompatible LatObs & LonObs.'
         STOP
      ENDIF
   
      NObs = SIZE(LatObs)
      NLat = SIZE(LatGrd)
      NLon = SIZE(LonGrd)


      H(:,:) = 0.0
      nH(:,:) = 0

!
! First we find which grid box contains a partucular observation:
      DO iObs=1,NObs
         ii1 = 0
         ii2 = 0
         jj1 = 0
         jj2 = 0

         DO i=1,NLat-1
            IF ( (LatObs(iObs) .GE. LatGrd(i)) .AND. &
                (LatObs(iObs) .LT. LatGrd(i+1)) ) THEN
              ii1=i
              ii2=i+1
            ENDIF            
         ENDDO  

         IF ( LatObs(iObs) .LE. LatGrd(1) ) THEN
            PRINT *, 'Obs_Operator_XY% Observation outside model domain: ', &
             LatObs(iObs), LatGrd(1)
            STOP
            ii1 = 1
            ii2 = 1
         ENDIF

         IF ( LatObs(iObs) .GE. LatGrd(NLat) ) THEN
            PRINT *, 'Obs_Operator_XY% Observation outside model domain: ', &
             LatObs(iObs), LatGrd(NLat)
            STOP
            ii1 = NLat
            ii2 = NLat
         ENDIF

         IF ( LatObs(iObs) .EQ. LatGrd(1) ) THEN
            ii1 = 1
            ii2 = 2
         ENDIF

         IF ( LatObs(iObs) .EQ. LatGrd(NLat) ) THEN
            ii1 = NLat-1
            ii2 = NLat
         ENDIF

         DO j=1,NLon-1
            IF ( (LonObs(iObs) .GE. LonGrd(j)) .AND. &
                (LonObs(iObs) .LT. LonGrd(j+1)) ) THEN
              jj1=j
              jj2=j+1
            ENDIF            
         ENDDO       
         LonGap = .False.
         IF ( LonObs(iObs) .LT. LonGrd(1) ) THEN
            jj1 = NLon
            jj2 = 1
            LonGap = .True.
         ENDIF
          IF ( LonObs(iObs) .GE. LonGrd(NLon) ) THEN
            jj1 = NLon
            jj2 = 1
            LonGap = .True.
         ENDIF    

         !
         ! Now we have four model grid points around the observation
         ! so we should pick three nearest ones and triangulate.
         dR1 = EarthDist( LatObs(iObs), LonObs(iObs), LatGrd(ii1), LonGrd(jj1) )
         dR2 = EarthDist( LatObs(iObs), LonObs(iObs), LatGrd(ii1), LonGrd(jj2) )
         dR3 = EarthDist( LatObs(iObs), LonObs(iObs), LatGrd(ii2), LonGrd(jj1) )
         dR4 = EarthDist( LatObs(iObs), LonObs(iObs), LatGrd(ii2), LonGrd(jj2) )

         !
         ! Figure which point is the most distant one and through it away.
            iLat1 = 0
            jLon1 = 0
            iLat2 = 0
            jLon2 = 0
            iLat3 = 0
            jLon3 = 0 

         Lon1 = 0
         Lon2 = 0
         Lon3 = 0

         IF ( dR1 .GE. MAXVAL( (/dR2, dR3, dR4/) ) ) THEN
            iLat1 = ii1
            jLon1 = jj2
            iLat2 = ii2
            jLon2 = jj1
            iLat3 = ii2
            jLon3 = jj2
            IF (LonGap) THEN 
               Lon1 = 360.0 + LonGrd(1)
              jLon1 = 1
               Lon2 = LonGrd(NLon)
              jLon2 = NLon 
               Lon3 = 360.0 + LonGrd(1)
              jLon3 = 1
             ELSE
               Lon1 = LonGrd(jLon1)
               Lon2 = LonGrd(jLon2)
               Lon3 = LonGrd(jLon3)
            ENDIF    
         ENDIF   

         IF ( dR2 .GE. MAXVAL( (/dR1, dR3, dR4/) ) ) THEN
            iLat1 = ii1
            jLon1 = jj1
            iLat2 = ii2
            jLon2 = jj1
            iLat3 = ii2
            jLon3 = jj2 
            IF (LonGap) THEN 
               Lon1 = LonGrd(NLon)
              jLon1 = NLon 
               Lon2 = LonGrd(NLon)
              jLon2 = NLon 
               Lon3 = 360.0 + LonGrd(1)
              jLon3 = 1 
             ELSE
               Lon1 = LonGrd(jLon1)
               Lon2 = LonGrd(jLon2)
               Lon3 = LonGrd(jLon3)
            ENDIF    
         ENDIF   

         IF ( dR3 .GE. MAXVAL( (/dR1, dR2, dR4/) ) ) THEN
            iLat1 = ii1
            jLon1 = jj1
            iLat2 = ii1
            jLon2 = jj2
            iLat3 = ii2
            jLon3 = jj2 
            IF (LonGap) THEN 
               Lon1 = LonGrd(NLon)
              jLon1 = NLon 
               Lon2 = 360.0 + LonGrd(1)
              jLon2 = 1 
               Lon3 = 360.0 + LonGrd(1)
              jLon3 = 1 
             ELSE
               Lon1 = LonGrd(jLon1)
               Lon2 = LonGrd(jLon2)
               Lon3 = LonGrd(jLon3)
            ENDIF    
         ENDIF 

         IF ( dR4 .GE. MAXVAL( (/dR1, dR2, dR3/) ) ) THEN
            iLat1 = ii1
            jLon1 = jj1
            iLat2 = ii1
            jLon2 = jj2
            iLat3 = ii2
            jLon3 = jj1 
            IF (LonGap) THEN 
               Lon1 = LonGrd(NLon)
              jLon1 = NLon 
               Lon2 = 360.0 + LonGrd(1)
              jLon2 = 1 
               Lon3 = LonGrd(NLon)
              jLon3 = NLon 
             ELSE
               Lon1 = LonGrd(jLon1)
               Lon2 = LonGrd(jLon2)
               Lon3 = LonGrd(jLon3)
            ENDIF    
         ENDIF 

         IF ( iLat1*iLat2*iLat3*jLon1*jLon2*jLon3 .EQ. 0) THEN
            PRINT *, 'Obs_Operator_XY% Something is wrong with my logic!'
            PRINT *, 'Obs_Operator_XY% ',iLat1,iLat2,iLat3,jLon1,jLon2,jLon3
            PRINT *, 'Obs_Operator_XY% ',ii1,ii2,jj1,jj2
            STOP
         ENDIF

         Lat1 = LatGrd(iLat1)
         Lat2 = LatGrd(iLat2)
         Lat3 = LatGrd(iLat3)

         !
         ! Now do the triangulation and return the interpolation
         ! coefficients alpha, beta and, gamma. The interpolated
         ! value is f(LatObs,LonObs) 
         !      = alpha * f(Lat1,Lon1) 
         !      +  beta * f(Lat2,Lon2) 
         !      + gamma * f(Lat3,Lon3)
         !
         IF ( LonObs(iObs) .LT. LonGrd(1) ) THEN
           CALL Triangulate( LatObs(iObs), LonObs(iObs)+360.0_r8, &
               Lat1,Lon1, Lat2,Lon2, Lat3,Lon3, alpha, beta, gamma )
          ELSE
           CALL Triangulate( LatObs(iObs), LonObs(iObs), &
               Lat1,Lon1, Lat2,Lon2, Lat3,Lon3, alpha, beta, gamma )    
         ENDIF  


         ! 
         ! Compute positions of the nearest points 1,2,3 in the
         ! vector of the model values.
         k1 = IndexF2V( NLat,NLon, iLat1,jLon1 ) 
         k2 = IndexF2V( NLat,NLon, iLat2,jLon2 ) 
         k3 = IndexF2V( NLat,NLon, iLat3,jLon3 )


         nH(1,iObs) = k1 
         nH(2,iObs) = k2
         nH(3,iObs) = k3
         H (1,iObs) = alpha
         H (2,iObs) = beta
         H (3,iObs) = gamma 
 
!!$         !
!!$         ! We need to arrange them in the increasing order (in k)
!!$         ! because this order is assumed in the sparse matrix
!!$         ! routines that will deal with H later. Pretty messy piece.
!!$         ! Try something more elegant if you can.
!!$         IF (k1 .LE. MIN(k2,k3)) THEN
!!$            nH(1,iObs) = k1
!!$            H (1,iObs) = alpha
!!$            IF (k2 .LE. k3) THEN
!!$               nH(2,iObs) = k2
!!$               nH(3,iObs) = k3
!!$               H (2,iObs) = beta
!!$               H (3,iObs) = gamma
!!$             ELSE
!!$               nH(2,iObs) = k3
!!$               nH(3,iObs) = k2
!!$               H (2,iObs) = gamma
!!$               H (3,iObs) = beta
!!$            ENDIF   
!!$         ENDIF  
!!$ 
!!$         IF (k2 .LE. MIN(k1,k3)) THEN
!!$            nH(1,iObs) = k2
!!$            H (1,iObs) = beta
!!$            IF (k1 .LE. k3) THEN
!!$               nH(2,iObs) = k1
!!$               nH(3,iObs) = k3
!!$               H (2,iObs) = alpha
!!$               H (3,iObs) = gamma
!!$             ELSE
!!$               nH(2,iObs) = k3
!!$               nH(3,iObs) = k1
!!$               H (2,iObs) = gamma
!!$               H (3,iObs) = alpha
!!$            ENDIF   
!!$         ENDIF   
!!$
!!$         IF (k3 .LE. MIN(k1,k2)) THEN
!!$            nH(1,iObs) = k3
!!$             H(1,iObs) = gamma
!!$            IF (k1 .LE. k2) THEN
!!$               nH(2,iObs) = k1
!!$               nH(3,iObs) = k2
!!$               H (2,iObs) = alpha
!!$               H (3,iObs) = beta
!!$             ELSE
!!$               nH(2,iObs) = k2
!!$               nH(3,iObs) = k1
!!$               H (2,iObs) = beta
!!$               H (3,iObs) = alpha
!!$            ENDIF   
!!$         ENDIF  
!!$
!!$         !
!!$         ! And make sure we did not miss anything in this piece above:
!!$         IF ( PRODUCT(nH(:,iObs)) .EQ. 0 ) THEN
!!$           PRINT *, 'Obs_Operator_XY% Failed to set nH - check the code'
!!$           PRINT *, 'Obs_Operator_XY% iObs, nH(:,iObs): ', iObs, nH(:,iObs)
!!$           PRINT *, 'Obs_Operator_XY% k1,k2,k3:', k1,k2,k3
!!$           PRINT *, 'Obs_Operator_XY% iLat1,jLon1:', iLat1,jLon1
!!$           PRINT *, 'Obs_Operator_XY% iLat2,jLon2:', iLat2,jLon2
!!$           PRINT *, 'Obs_Operator_XY% iLat3,jLon3:', iLat3,jLon3
!!$           STOP
!!$         ENDIF  
!!$         IF (.NOT. ((nH(1,iObs) .LE. nH(2,iObs)) .AND. 
!!$     &              (nH(2,iObs) .LE. nH(3,iObs))) )  THEN
!!$           PRINT *, 'Obs_Operator_XY% Failed to sort H, nH - check the code'
!!$           PRINT *, 'Obs_Operator_XY% iObs, nH(:,iObs): ', iObs, nH(:,iObs)
!!$           PRINT *, 'Obs_Operator_XY% k1,k2,k3:', k1,k2,k3
!!$           PRINT *, 'Obs_Operator_XY% iLat1,jLon1:', iLat1,jLon1
!!$           PRINT *, 'Obs_Operator_XY% iLat2,jLon2:', iLat2,jLon2
!!$           PRINT *, 'Obs_Operator_XY% iLat3,jLon3:', iLat3,jLon3
!!$           STOP
!!$         ENDIF           
      ENDDO
!
      RETURN
      END SUBROUTINE Obs_Operator_XY



      SUBROUTINE Triangulate( x,y, x1,y1, x2,y2, x3,y3, alpha, beta, gamma )
!+
!
!-
      IMPLICIT NONE
! INPUTS:
      real(r8), INTENT(IN) :: x,y, x1,y1, x2,y2, x3,y3
!
! OUTPUTS:
      real(r8), INTENT(OUT) :: alpha, beta, gamma
!
! INTERNALS:
      real(r8) x01,x02,y01,y02,x12,x13,y12,y13,x23,y23
      real(r8) :: xMin,xMax, yMin,yMax
      real(r8) :: discriminant
      real(r8) :: d21, d32, d31, sqdiam ! squared distances and set diameter

      xMin = MIN( x1, MIN(x2,x3) )
      xMax = MAX( x1, MAX(x2,x3) )
      yMin = MIN( y1, MIN(y2,y3) )
      yMax = MAX( y1, MAX(y2,y3) )

      IF (x .LT. xMin) THEN
         PRINT *, 'Triangulate% x < MIN(x1,x2,x3) !'
         PRINT *, 'Triangulate% x, x1, x2, x3: ', x, x1, x2, x3
      ENDIF
      IF (x .GT. xMax) THEN
         PRINT *, 'Triangulate% x > MAX(x1,x2,x3) !'
         PRINT *, 'Triangulate% x, x1, x2, x3: ', x, x1, x2,x3
      ENDIF
   
      IF (y .LT. yMin) THEN
         PRINT *, 'Triangulate% y < MIN(y1,y2,y3) !'
         PRINT *, 'Triangulate% y, y1, y2, y3: ', y, y1, y2, y3
      ENDIF
      IF (y .GT. yMax) THEN
         PRINT *, 'Triangulate% y > MAX(y1,y2,y3) !'
         PRINT *, 'Triangulate% y, y1, y2, y3: ', y, y1, y2, y3
      ENDIF  


      x01 = x  - x1
      x02 = x  - x2
      x12 = x1 - x2
      x13 = x1 - x3
      x23 = x2 - x3
      y01 = y  - y1
      y02 = y  - y2
      y12 = y1 - y2
      y13 = y1 - y3
      y23 = y2 - y3

      discriminant = ( y12*x23 - y23*x12 )
      d21 = x12*x12 + y12*y12
      d31 = x13*x13 + y13*y13
      d32 = x23*x23 + y23*y23
      sqdiam = MAX(d21, MAX(d31,d32))

      IF ( abs(discriminant / sqdiam ) .LT. epsilon(1._r8) ) THEN
         PRINT *, 'Triangulate% Points are in a line to machine precision.'
         PRINT *, '1: ', x1,y1, ' 2: ', x2,y2, ' 3: ', x3,y3 
         STOP
      ENDIF   

      alpha = (x23*y02 - y23*x02)/discriminant
      beta  = (x01*y13 - y01*x13)/discriminant
      gamma = (x12*y01 - y12*x01)/discriminant

!      if (( alpha.lt.0 ) .or. ( beta.lt.0 ).or. ( gamma.lt.0 ))  then
      if ( alpha*beta*gamma.lt.0 )  then
         write(6,*)'WARNING: triangulate LT0 ',alpha,beta,gamma
      end if


      RETURN
      END SUBROUTINE Triangulate

 SUBROUTINE Obs_Operator_Z( zObs, z, dz, H )
!+
!
!-
   IMPLICIT NONE
!
! INPUTS:
   real(r8), INTENT(IN) :: zObs(:), z(:), dz(:)
!
! OUTPUTS:
   real(r8), INTENT(OUT) :: H(SIZE(z), SIZE(zObs))
! INTERNALS:
!
   INTEGER :: i,j
   real(r8) :: A(N_Kernel_Grd,SIZE(zObs)), z_A(N_Kernel_Grd)
   real(r8) :: H_(SIZE(z), SIZE(zObs))
!

!
! First get the tabulated averaging kernels on grid z_A:
   !PRINT *,'      Obs_Operator_Z% Calling GET_AvgKernel . . .'
   CALL Get_AvgKernel( A, zObs, z_A )
!
! Now interpolate the averaging kernels to grid z:
   !PRINT *,'      Obs_Operator_Z% Interpolating averaging kernels . . .'
   DO i=1,SIZE(zObs)
        H(:,i) = LInterp1( A(:,i), z_A, z )
   ENDDO
!
! Normalize, i.e. make sure that the integral of the vertical observational operator is 1
! PRINT *,'      Obs_Operator_Z% Normalize interpolated averaging kernels . . .'
!
   DO i=1,SIZE(zObs)
      IF (SUM(H(:,i)) .EQ. 0) THEN
         PRINT *, '      Obs_Operator_Z% No non-zero elements in averaging kernel!'
         PRINT *, '      Obs_Operator_Z% ', i, SIZE(zObs), H(:,i)
         PRINT *, '      Obs_Operator_Z% z:', z
         PRINT *, '      Obs_Operator_Z% z_A:', z_A
         PRINT *, '      Obs_Operator_Z% A(:,i):', A(:,i)
         PRINT *, '      Obs_Operator_Z% STOP.'
         call abort()
!
! WDC -- no normalization necessary
!
!      ELSE
!         H(:,i) = H(:,i)/SUM(dz)
!        H(:,i) = H(:,i)/SUM(H(:,i))
      ENDIF
   ENDDO
!
!  h = 0.
!  h(15,1) = 1.
!
 RETURN
 END SUBROUTINE Obs_Operator_Z



END MODULE obs_operator
