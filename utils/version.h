/*$Id: version.h,v 1.2 1998/04/07 22:32:53 eaton Exp $*/

/* Define version of program for examining a CCM history tape.  The
   floating point representations indicated by CRAY_FLOAT, IEEE_SP, and IEEE_DP
   refer to how numbers are stored internally on the machine running the
   program and NOT to how numbers are represented on the history tape.
   There are conversion routines built into the program so that it can
   examine any format of history tape regardless of how numbers are
   represented internally.
*/

#ifndef VERSION_H
#define VERSION_H

#if defined(CRAY_FLOAT) || defined(IEEE_SP)
   typedef int Hint;
   typedef float Hfloat;
#elif defined(IEEE_DP)
   typedef int Hint;
   typedef double Hfloat;
#endif

#endif
