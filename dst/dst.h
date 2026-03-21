c     $Header$ -*-fortran-*-

c     Purpose: dst.h sets and describes all tokens used in dust parameterization

c     Usage: 
c     #include <dst.h> /* Dust preprocessor tokens */ 

#ifndef DST_H
#define DST_H

#ifdef DST

c     dst.h MUST have access to params.h to work
c     CCM files all include #include params.h, but MATCH files do not
c     The following three lines ensure access to params.h for MATCH files
#ifndef PARAMS_H
#include <params.h>
#endif /* not PARAMS_H */ 

c     Host model:
c     Tokens must be set to enable the correct output routines
c     CCM activates the 5-argument CCM outfld() call
c     When CCM is not defined, the 7-argument MATCH outfld() call is used
c     BXM activates the extensive 1-D netCDF output calls
c     Normally CCM, and BXM are set in the Makefile rather than here
c#define CCM
c#define BXM

c     Set radiative mode:
c     DST_RAD activates passing of dust to radiation routines
c     If DST_RAD is not defined then the code should be bit-for-bit with CCM
c     In this case the dust transport will occur, but optical depth is the only radiative property archived
#define DST_RAD
#ifdef DST_RAD
c     When DST_RAD is defined, either DST_FDB or DST_FRC must also be defined
c     DST_FDB implements radiative feedback mode
c     DST_FDB calls the radiation routines only once
c     The DST_FDB dynamics simulation will differ from CCM
c     The other alternative when DST_RAD is defined is to use forcing mode, DST_FRC
c     DST_FRC calls the radiation routines with and without dust
c     The difference is the diagnostic radiative forcing due to dust
c     Archiving all the DST_FRC fields involves a lot of bookkeeping
c     The DST_FRC dynamics simulation should be bit-for-bit with CCM
#define DST_FRC
c#define DST_FDB
#endif /* not DST_RAD */ 

c     Set land surface mode:
c     DST_LSM activates passing of LSM data into CCM
c     Currently, DST_LSM is required by the mobilization scheme
#define DST_LSM

c     Set chemistry mode:
c     DST_CHM activates input and processing of chemistry arrays
#define DST_CHM

c     Set debugging modes:
c     DST_DBG activates non-vectorized sanity checks and verbose output:
c     1. Tendency procedures (e.g., dst_dps_wet()) will check for unreasonable mixing ratios after each adjustment
c     2. Surface soil properties at the debug gridpoint are sometimes printed
c     3. qneg3() is called after each source/sink routine
c     Production runs are either on CRAY or SGI nowadays so default is no DST_DBG on these machines
c#if (!defined CRAY) && (!defined SGI)
c #define DST_DBG
c#endif /* not CRAY or SGI */
#define DST_DBG

c     DST_MSS_BDG activates global mass budget output into netCDF file
c     The diagnostic mass budget requires running global averages on many quantities each timestep
c     This can be very expensive, but makes checking for mass balance much easier
#if (!defined CRAY) && (!defined SGI)
!#define DST_MSS_BDG
#endif /* not CRAY or SGI */

c     Automated check for erroneous token combinations.
c     No editing should need to be done beneath this line.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Dynamic source scheme requires interactive LSM data
#if (defined DST_DYN) && (!defined DST_LSM)
      "ERROR dst.h defines DST_DYN but not DST_LSM"
#endif /* not (DST_DYN && DST_TMS) */ 

#ifdef DST_RAD
c     Exactly one of DST_FDB and DST_FRC must be defined when DST_RAD is defined 
#if (defined DST_FDB) && (defined DST_FRC)
      "ERROR dst.h defines both DST_FDB and DST_FRC"
#endif /* not (DST_FDB && DST_FRC) */ 
#if (!defined DST_FDB) && (!defined DST_FRC)
      "ERROR dst.h defines neither DST_FDB nor DST_FRC"
#endif /* (DST_FDB && DST_FRC) */ 
#endif /* not DST_RAD */ 

#endif /* not DST */

#endif /* not DST_H */

#define MODIS
