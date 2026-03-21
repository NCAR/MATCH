MODULE data_structures

use floats
use precision

IMPLICIT NONE

INTEGER, PARAMETER :: Profile_N_Lev = 1

TYPE ProfileType
   REAL(r8) :: Lon
   REAL(r8) :: Lat
   integer :: Date
   integer :: Secs
   REAL(r8), DIMENSION(Profile_N_Lev) :: V
   REAL(r8), DIMENSION(Profile_N_Lev) :: sigmaV
END TYPE ProfileType


TYPE Orbits
   INTEGER :: N
   REAL(r8) :: NoData
   REAL(r8), DIMENSION(Profile_N_Lev) :: Profile_Grid
   TYPE (ProfileType), DIMENSION(:), POINTER :: Profile 
END TYPE Orbits  


END MODULE data_structures
