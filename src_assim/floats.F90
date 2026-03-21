MODULE floats

!
! Floating point precision parameters
!
  use precision

 INTEGER, PARAMETER :: low = 4
 INTEGER, PARAMETER :: dbl = 8
 INTEGER, PARAMETER :: ext = 16

 INTEGER, PARAMETER :: std = dbl


 REAL(r8), PARAMETER :: Zero   = 0.0_r8
 REAL(r8), PARAMETER :: Half   = 0.5_r8 
 REAL(r8), PARAMETER :: One    = 1.0_r8
 REAL(r8), PARAMETER :: Two    = 2.0_r8 
 REAL(r8), PARAMETER :: Three  = 3.0_r8
 
 REAL(r8), PARAMETER :: Small = TINY(Zero)

END MODULE floats
