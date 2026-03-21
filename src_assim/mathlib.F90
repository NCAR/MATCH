
MODULE mathlib

USE floats
use precision

IMPLICIT NONE

PUBLIC :: LInterp1, Interpol, Invert, MatSolve

CONTAINS



FUNCTION LInterp1( y, x, x1, NoData )
!+
! PURPOSE:
!      1-D Linear interpolation and extrapolation.
!
! Input:
!	y:	 The input vector 
!	x:	 The absicissae values for y
!	x1:	 The absicissae values for the result (y1)
!	NoData:  Missing data value for y
!
! Output:
!	Returns vector of interpolated values. No extrapolation
!	is performed if NoData is set - NoData value is assigned to
!	the interpolated data outside of the x range.
!-
  use precision
 IMPLICIT NONE
!
 REAL(r8), INTENT(IN) :: y(:), x(:), x1(:)
 REAL(r8), OPTIONAL, INTENT(IN) :: NoData
 REAL(r8), DIMENSION(SIZE(x1)) :: LInterp1
!
 INTEGER :: i, k, N, Counter(SIZE(x)) 
 REAL(r8), DIMENSION(SIZE(x1)) :: y1
 REAL(r8), ALLOCATABLE, DIMENSION(:) :: xx, yy
 REAL(r8) :: xMin, xMax
!
  IF (SIZE(x) .NE. SIZE(y)) THEN
    PRINT *, 'LInterp1% Incompatible x & y'
    STOP
  ENDIF

  IF (PRESENT(NoData)) THEN 

     Counter = 0

     WHERE (y .NE. NoData)
         Counter = 1
     ENDWHERE 

     N = SUM( Counter )

    IF (N .GE. 2) THEN
         k = 1 
         DO  i=1, SIZE(x)
            IF ( y(i) .NE. NoData ) THEN
             xx(k) = x(i)
             yy(k) = y(i)
             k = k+1
            ENDIF
        ENDDO
     ELSE 
        PRINT *, 'LInterp1% Only one valid point.'
        y1 = NoData
        LInterp1 = y1
        RETURN
     ENDIF

     y1=INTERPOL(yy,xx,x1)

     xMax = MAXVAL(xx)
     xMin = MINVAL(xx)

     WHERE (x1 < xMin)
         y1 = NoData
     ENDWHERE

     WHERE (x1 > xMax)
         y1 = NoData
     ENDWHERE
  ELSE 
     y1=INTERPOL(y,x,x1)
  ENDIF

 LInterp1 = y1

!
RETURN
END FUNCTION LINTERP1




FUNCTION INTERPOL( V, X, U )
!+
! NAME:
!	INTERPOL
!
! PURPOSE:
!	Linearly interpolate and extrapolate vectors with a regular or irregular grid.
!
! CATEGORY:
!	E1 - Interpolation
!
! CALLING SEQUENCE:
!	Result = INTERPOL(V, X, U)	
!
! INPUTS:
!	V:	The input vector can be any type except string.
!
!	X:	The absicissae values for V.  This vecotr must have same # of
!		elements as V.  The values MUST be monotonically ascending 
!		or descending.
!
!	U:	The absicissae values for the result.  The result will have 
!		the same number of elements as U.  U does not need to be 
!		monotonic.
!	
!
! OUTPUTS:
!	INTERPOL returns a floating-point vector of N points determined
!	by linearly interpolating the input vector.
!
! PROCEDURE:
!	Result(i) = V(x) + (x - FIX(x)) * (V(x+1) - V(x))
!
!	where        x = U(i).
!		m = number of points of input vector.
!
! MODIFICATION HISTORY:
!                 IDL routine translated to F90 by Boris Khattatov, 9/25/96
!-
!
  use precision
 IMPLICIT NONE
!
 REAL(r8), INTENT(IN) :: v(:), x(:), u(:)
 REAL(r8), DIMENSION( SIZE(u) ) :: Interpol
!
 REAL(r8), DIMENSION( SIZE(u) ) :: r
 REAL(r8) :: s1, d
 INTEGER :: M, N, ix, i
!
 M = SIZE(v)	!# of input pnts

 IF (SIZE(x) .NE. M) THEN  
    PRINT *, 'INTERPOL %  V and X must have same # of elements'
    STOP
 ENDIF

 N= SIZE(u)	!# of output points

 r = V(1)	!floating, dbl or cmplx result

 IF (x(2) - x(1) .GE. 0) THEN 
    s1 = 1.0 
 ELSE 
    s1 = -1.0 
 ENDIF !Incr or Decr X

 ix = 1			!current point
 DO i=1,N 	                   !point loop
    d = s1 * (u(i)-x(ix))	!difference
    IF (ABS(d) .LT. TINY(d)) THEN 
       r(i)=v(ix) 
    ELSE   !at point
       IF (d > 0) THEN 
          DO WHILE ( (s1*(u(i)-x(ix+1)) > 0) .AND. (ix < m-1) )
             ix=ix+1
          ENDDO
       ELSE 
          DO WHILE ( (s1*(u(i)-x(ix)) < 0) .AND.  (ix > 1) ) 
             ix=ix-1
          ENDDO
       ENDIF
       r(i) = v(ix) + (u(i)-x(ix))*(v(ix+1)-v(ix))/(x(ix+1)-x(ix))
    ENDIF
 ENDDO

 INTERPOL = r
 !
RETURN
END FUNCTION INTERPOL



FUNCTION Invert( A )
!+
!
!-
  use precision
 IMPLICIT NONE
!
  REAL(r8), INTENT(IN) :: A(:,:)
  REAL(r8), DIMENSION(SIZE(A(:,1)), SIZE(A(1,:)) ) :: Invert
!
  REAL(r8), DIMENSION(SIZE(A(:,1)), SIZE(A(1,:)) ) :: AInv, B, BInv
  REAL(r8) :: d
  INTEGER :: i,j, N, Indx(SIZE(A(1,:))), Flag
!
      N = SIZE(A(:,1))

      DO i=1,N
       DO j=1,N
        BInv(i,j) = Zero
        B(i,j) = Zero
       ENDDO
       BInv(i,i) = One
      ENDDO 
 
      DO i=1,N 
       DO j=1,N
         B(i,j)=A(i,j)
       ENDDO 
      ENDDO
 
      CALL LUDCMP( B, N, N, Indx, d, Flag )
      DO i=1,N
         CALL LUBKSB( B, N, N, Indx, BInv(1,i) )
      ENDDO

      DO i=1,N 
       DO j=1,N
         AInv(i,j)=BInv(i,j)
       ENDDO 
      ENDDO
 
      Invert = AInv

RETURN
END FUNCTION Invert



SUBROUTINE MatSolve( A, b, x)
!+
! Solves A*x=b
!-
  use precision
  IMPLICIT NONE
!
! INPUTS:
  REAL(r8), INTENT(IN) :: A(:,:), b(:)
!
! OUTPUTS:
  REAL(r8), INTENT(OUT) :: x(SIZE(b))
!
! INTERNALS:
  INTEGER:: N
  REAL(r8) :: d
  INTEGER :: Indx(SIZE(b)), Flag
  REAL(r8) :: AA(SIZE(b),SIZE(b)), bb(SIZE(b))
!
  real(r8) :: p(size(b))
!
  N = SIZE(b)
!
  AA = A
  bb = b
!
  CALL LUDCMP( AA, N, N, Indx, d, Flag ) 
  CALL LUBKSB( AA, N, N, Indx, bb )
  x = bb
!
! use Cholesky decomposition instead of LU-decomposition
!
! call choldc(aa,n,n,p)
! call cholsl(aa,n,n,p,bb,x)
!
! use conjugate gradient method
!
! call conjugate_gradient(n,aa,b,x)
!
  RETURN
END SUBROUTINE MatSolve



SUBROUTINE lubksb(a,n,np,indx,b)
!+
!
!-
  use precision
  IMPLICIT NONE
!
 INTEGER n,np,indx(n)
 REAL(r8) a(np,np),b(n)
!
 INTEGER i,ii,j,ll
 REAL(r8) sum_
!
      ii=0
      DO 12 i=1,n
        ll=indx(i)
        sum_=b(ll)
        b(ll)=b(i)
        IF (ii  .NE.  0) THEN
          DO  j=ii,i-1
            sum_=sum_-a(i,j)*b(j)
          ENDDO 
        ELSE IF (sum_  .NE.  0.) THEN
          ii=i
        ENDIF
        b(i)=sum_
12    CONTINUE
      DO i=n,1,-1
        sum_=b(i)
         DO  j=i+1,n
            sum_=sum_-a(i,j)*b(j)
         ENDDO
        b(i)=sum_/a(i,i)
      ENDDO
!
RETURN
END SUBROUTINE lubksb



SUBROUTINE ludcmp(a,n,np,indx,d, Flag)
!+
!
!-
  use precision
  IMPLICIT NONE
!
 INTEGER :: n,np,indx(n)
 REAL(r8) :: d,a(np,np)
 INTEGER, INTENT(OUT) :: Flag
!
 INTEGER i,imax,j,k
 REAL(r8) aamax,dum,sum_,vv(NP)
!
     Flag = 0

      d=1.
      DO 12 i=1,n
        aamax=0.
        DO  j=1,n
          IF (ABS(a(i,j)) .GT. aamax) aamax=ABS(a(i,j))
        ENDDO
        IF (aamax .EQ. 0.) THEN
          PRINT *,  'Singular matrix in ludcmp'
          PRINT *, A
          Flag = -1
          RETURN
        ENDIF    
        vv(i)=1./aamax
12    CONTINUE
      DO  19 j=1,n
        DO  i=1,j-1
          sum_=a(i,j)
          DO  k=1,i-1
            sum_=sum_-a(i,k)*a(k,j)
          ENDDO
          a(i,j)=sum_
        ENDDO
        aamax=0.
        DO 16 i=j,n
          sum_=a(i,j)
          DO k=1,j-1
            sum_=sum_-a(i,k)*a(k,j)
          ENDDO
          a(i,j)=sum_
          dum=vv(i)*ABS(sum_)
          IF (dum.GE.aamax) THEN
            imax=i
            aamax=dum
          ENDIF
16      CONTINUE
        IF (j  .NE.  imax) THEN
          DO k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          ENDDO 
          d=-d
          vv(imax)=vv(j)
        ENDIF
        indx(j)=imax
        IF (a(j,j) .EQ. 0.) a(j,j)=Small
        IF (j  .NE.  n) THEN
          dum=1./a(j,j)
          DO  i=j+1,n
            a(i,j)=a(i,j)*dum
          ENDDO
        ENDIF
19    CONTINUE
!
RETURN
END SUBROUTINE ludcmp
!
! the following subroutines for the Cholesky decompostion
! come from the Numerical Recipes (chapter 2.9)
!
      SUBROUTINE choldc(a,n,np,p)
!
      use precision
      implicit none
!
      INTEGER n,np
      REAL(r8) a(np,np),p(n)
      INTEGER i,j,k
      REAL(r8) wk_sum
!
      do 13 i=1,n
        do 12 j=i,n
          wk_sum=a(i,j)
          do 11 k=i-1,1,-1
            wk_sum=wk_sum-a(i,k)*a(j,k)
11        continue
          if(i.eq.j)then
            if(wk_sum.le.0.) then
              print *,'choldc failed ',wk_sum
              do k=i,1,-1
                print *,i,k,a(i,k)
              end do
              write(110) n
              write(110) a
              call abort()
            endif
            p(i)=sqrt(wk_sum)
          else
            a(j,i)=wk_sum/p(i)
          endif
12      continue
13    continue
      return
      END subroutine choldc
!
      SUBROUTINE cholsl(a,n,np,p,b,x)
!
      use precision
      implicit none
!
      INTEGER n,np
      REAL(r8) a(np,np),b(n),p(n),x(n)
      INTEGER i,k
      REAL(r8) wk_sum
      do 12 i=1,n
        wk_sum=b(i)
        do 11 k=i-1,1,-1
          wk_sum=wk_sum-a(i,k)*x(k)
11      continue
        x(i)=wk_sum/p(i)
12    continue
      do 14 i=n,1,-1
        wk_sum=x(i)
        do 13 k=i+1,n
          wk_sum=wk_sum-a(k,i)*x(k)
13      continue
        x(i)=wk_sum/p(i)
14    continue
      return
      END subroutine cholsl
!
      subroutine conjugate_gradient(n,a,b,x)
!
! based on algorithm 10.2.7 from Golub and van Loon
!
      use precision
      implicit none
!
! arguments
!
      integer n
      real(r8) a(n,n),b(n),x(n)
!
! local
!
      integer i,j,niter,niter_max
      real(r8) r(n),p(n),w(n),sum
      real(r8) alpha,beta,epsilon,rho,norm_b
      allocatable :: rho(:)
!
! parameters
!
      niter_max = 100
      epsilon   = 1.e-8

! set x to 0, r to b
!
      do i=1,n
        x(i) = 0.
        r(i) = b(i)
      end do
!
! set counter to 0
!
      niter = 0
!
! iterate
!
 100  continue
!
! compute 2-norm of r
!
      rho(niter) = 0.d0
      do i=1,n
        rho(niter) = rho(niter) + (r(i)*r(i))
      end do
!
! check convergence
!
! ! rho(0) = ||b||^2
!
      if(sqrt(rho(niter)).lt.(epsilon*sqrt(rho(0)))) then
        print    *,'Conjugate-gradient has converged after ',niter &
                  ,' iterations'
        return
      endif
!
      niter = niter + 1
!
! test number of iterations
!
      if(niter.gt.niter_max) then
        print *,'reached maximum number of iterations in CG, stopping'
        call abort()
      endif
!
      if(niter.eq.1) then
        do i=1,n
          p(i) = r(i)
        end do
      else
        beta = rho(niter-1)/rho(niter-2)
        do i=1,n
          p(i) = r(i) + beta * p(i)
        end do
      endif
!
      do i=1,n
        w(i) = 0.
        do j=1,n
          w(i) = w(i) + a(i,j) * p(j)
        end do
      end do
      sum = 0.
      do i=1,n
        sum = sum + p(i) * w(i)
      end do
      alpha = rho(niter-1)/sum
      do i=1,n
        x(i) = x(i) + alpha * p(i)
        r(i) = r(i) - alpha * w(i)
      end do
!
      goto 100
!
      end subroutine conjugate_gradient
!
END MODULE mathlib
