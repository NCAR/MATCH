c $Id$
c Restart info file data.

      integer
     $  lunrst0    ! unit number for restart info file
     $, lunrst1    ! unit number for restart data file
     $, rstrrt     ! retention period (days) for restart file
     $, rstfrq     ! output frequency (steps) for restart files
     $, nestep     ! number of ending timestep

      common / rest_I / 
     $  lunrst0, lunrst1, rstrrt, rstfrq, nestep

      logical
     $  async       ! t => dispose to remote device asynchronously
     $, rmrst       ! t => remove local copy after disposing
     $, diconvccm   ! t => diagnose CCM convection
     $, dicldphys   ! t => diagnose cloud physics
     $, advqth      ! t => advect Q and THETA (potential temperature)

      common / rest_L / 
     $  async, rmrst, diconvccm, dicldphys, advqth

      character
     $  lpthrst0*128   ! local pathname for restart info file
     $, rpthrst0*128   ! remote pathname for restart info file
     $, lpthrst1*128   ! local pathname for restart data
     $, rpthrst1*128   ! remote pathname for restart data

      common / rest_C / 
     $  lpthrst0,  rpthrst0, lpthrst1, rpthrst1
