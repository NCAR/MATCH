c     $Header$ -*-fortran-*- 

c     Purpose: Systematic nomenclature for debugging/verbosity levels

c     Usage:
c#include <dbg.com> /* Debugging constants */ 

      common / dbg /
     $     dbg_lvl,             ! Initialized in main
     $     prg_nm               ! Initialized in cmd_ln_sng()
      integer dbg_lvl
      character prg_nm*80

c     Debugging types are as follows:
      integer dbg_nbr
      parameter(dbg_nbr=9)      ! Number of different debugging levels

      integer dbg_off
      integer dbg_fl
      integer dbg_scl
      integer dbg_crr
      integer dbg_sbr
      integer dbg_vec
      integer dbg_io
      integer dbg_vrb
      integer dbg_old

      parameter(dbg_off=0)      ! Production mode. Debugging is turned off.
      parameter(dbg_fl=1)       ! Filenames
      parameter(dbg_scl=2)      ! Scalars
      parameter(dbg_crr=3)      ! Current task
      parameter(dbg_sbr=4)      ! Subroutine names on entry and exit
      parameter(dbg_io=5)       ! Subroutine I/O
      parameter(dbg_vec=6)      ! Entire vectors
      parameter(dbg_vrb=7)      ! Everything
      parameter(dbg_old=8)      ! Old debugging blocks not used anymore

