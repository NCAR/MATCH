c     $Header$ -*-fortran-*-

c     Purpose: Set parameters needed by netCDF routines

c     Dependencies: dsttibds.h MUST have access to params.h to work
#ifndef PARAMS_H
#include <params.h> /* Preprocessor tokens */ 
#endif /* not PARAMS_H */ 

c     Usage: 
c#include <dsttibds.h> /* Time-invariant boundary data sets */

      integer 
     $     pcnst,
     $     plat,
     $     plev,
     $     plon,
     $     plond
      parameter(pcnst=PCNST,
     $     plat=PLAT,
     $     plev=PLEV,
     $     plon=PLON,
     $     plond=PLON+3)

