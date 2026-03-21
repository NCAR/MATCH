module precision
!-------------------------------------------------------------------------------
!
! Purpose:
!	Define the precision to use for floating point and integer operations
!	throughout the model.
!
! Author: CCM Core group
!
!-------------------------------------------------------------------------------
  integer, parameter :: r8 = selected_real_kind(RPREC)
  integer, parameter :: i8 = selected_int_kind(13)
  integer, parameter :: DBLKIND = selected_real_kind( p = 12 )   
end module precision
