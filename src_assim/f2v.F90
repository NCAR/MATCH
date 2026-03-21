MODULE f2v

use floats

CONTAINS

      SUBROUTINE Field2Vector ( NLat, NLon, F, V )
!+
!
!- 
!
        use precision
      IMPLICIT NONE
! INPUTS:
      INTEGER, INTENT(IN) :: NLat,NLon
      REAL(r8), INTENT(IN) :: F(NLat,NLon)
! OUTPUT:
      REAL(r8), INTENT(OUT) :: V(NLon*NLat)
!
! INTERNALS:
!      INTEGER :: IndexF2V
      INTEGER :: iLat,jLon,k
!
      DO iLat=1,NLat
         DO jLon=1,NLon
            k = IndexF2V( NLat, NLon, iLat, jLon ) 
            V(k) = F(iLat,jLon)
         ENDDO   
      ENDDO   
!
      RETURN
      END  SUBROUTINE Field2Vector 




      SUBROUTINE Vector2Field ( NLat, NLon, V, F )
!+
!
!-
!
        use precision
      IMPLICIT NONE
! INPUTS:
      INTEGER, INTENT(IN) :: NLat,NLon
      REAL(r8), INTENT(IN) :: V(NLat*NLon)
! OUTPUT:
      REAL(r8), INTENT(OUT) :: F(NLat,NLon)
!
! INTERNALS:
!      INTEGER :: IndexF2V
      INTEGER :: iLat,jLon,k
!
      DO iLat=1,NLat
         DO jLon=1,NLon
            k = IndexF2V( NLat, NLon, iLat, jLon ) 
            F(iLat,jLon) = V(k)
         ENDDO   
      ENDDO   
!
      RETURN
      END SUBROUTINE Vector2Field 




      INTEGER FUNCTION IndexF2V( NLat, NLon, iLat, iLon )
!+
!
!
!-
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NLat, NLon, iLat, iLon
! 
      IndexF2V = iLat + (iLon-1)*NLat
!
      RETURN
      END FUNCTION IndexF2V




      SUBROUTINE IndexV2F( NLat, NLon, k, iLat, iLon )
!+
!
!
!-
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NLat, NLon, k
      INTEGER, INTENT(OUT) :: iLat, iLon
!
!      INTEGER :: IndexF2V
!
      iLon = INT((k-0.001)/NLat) + 1
      iLat = k - (iLon-1)*NLat
      IF (iLat .EQ. 0) THEN
         PRINT *, 'IndexV2F% NLat,NLon, k,iLat,iLon: ', NLat,NLon, k, iLat,iLon
         STOP
      ENDIF
      IF (iLon .EQ. 0) THEN
         PRINT *, 'IndexV2F% NLat,NLon, k,iLat,iLon: ', NLat,NLon, k, iLat,iLon
         STOP
      ENDIF   
      IF (IndexF2V(NLat,NLon,iLat,iLon) .NE. k) THEN
         PRINT *, 'IndexV2F% Re-write me! ', iLat, iLon, k, &
         IndexF2V(NLat,NLon,iLat,iLon)
         PRINT *, 'IndexV2F% STOP.'
         STOP
      ENDIF   
!
      RETURN
      END SUBROUTINE IndexV2F


END MODULE f2v
