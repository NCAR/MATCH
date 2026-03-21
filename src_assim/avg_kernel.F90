MODULE avg_kernel

  use precision
  use floats

!
! WDC -- the averaging kernel should be unity at all levels, and the
!        operator used on the model is linearly interpolated from this kernel
!        So just define the kernel as a constant (1) at 2 levels at TOA
!        and surface
!
 INTEGER, PARAMETER :: N_Kernel_Grd  = 2   ! number of gridpoints in the kernel definition
 INTEGER, PARAMETER :: N_AOD_Layers =  1   ! number of layers in the vertical profile(s)

 REAL(r8), PARAMETER, DIMENSION(N_Kernel_Grd) :: z_Kernel_Grd = (/1100.0,0.0/)

 PUBLIC :: N_Kernel_Grd, Get_AvgKernel

CONTAINS


 SUBROUTINE Get_AvgKernel( A, zRet, zGrd )
!
   use precision
   IMPLICIT NONE
!
! INPUTS: 
!
   REAL(r8), INTENT(IN) :: zRet(:)
!
! OUTPUTS:
!
   REAL(r8), INTENT(OUT) :: A(N_Kernel_Grd,N_AOD_Layers), zGrd(N_Kernel_Grd)
!
! INTERNALS:
!
   INTEGER :: NRet, i
   REAL(r8) :: A1(N_Kernel_Grd)
!
   zGrd = z_Kernel_Grd 
   NRet = N_AOD_Layers
!
! Define the averaging kernels for each layer:
!
   A( 1,1) = 1.000 ! 1100 mb
   A( 2,1) = 1.000 !    0 mb
!
 RETURN
 END SUBROUTINE Get_AvgKernel

END MODULE avg_kernel
