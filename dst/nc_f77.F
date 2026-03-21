c     $Header$

c     Purpose: General purpose netCDF routines

c     Usage: 
c     cd $HOME/f; make -W nc_f77.F OPTS=D libcsz_f77; cd -
c     cd $HOME/f; make -W nc_f77.F libcsz_f77; cd -
c     cd $HOME/f; make OPTS=D libcsz_f77; cd -

      subroutine nc_err_exit(rcd,msg)
c     Purpose: Print netCDF error message, routine name, and exit
#include <dbg.com> /* Debugging constants */
#include <netcdf.inc>
c     Input
      integer rcd               ! Return code
      character*(*) msg         ! Diagnostic message
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
c     Main Code
      write (6,'(2a)') prg_nm(1:f_strlen(prg_nm)),nf_strerror(rcd)
      write (6,'(3a)') prg_nm(1:f_strlen(prg_nm)),': ERROR ',msg(1:f_strlen(msg))
      stop 
      end                       ! end nc_err_exit()

      integer function nf_open_csz(fl_nm,mode,nc_id)
c     Purpose: Close netCDF file
c     Prototype:
c     integer nf_open_csz       ! Open existing netCDF file (libcsz_f77)
      implicit none
c     Commons
#include <dbg.com> /* Debugging constants */
#include <netcdf.inc>
c     Input
      character*(*) fl_nm       ! Filename
      integer mode              ! netCDF read/write mode
      integer nc_id             ! netCDF file handle
c     Local
      integer rcd               ! [rcd] Return success code
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
c     Main code
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf_open(fl_nm,mode,nc_id)
      nf_open_csz=rcd
      if (rcd /= NF_NOERR) then
         write (6,'(3a)') prg_nm(1:f_strlen(prg_nm)),': ERROR opening ',fl_nm(1:f_strlen(fl_nm))
         write (6,'(2a,i3,2a)') prg_nm(1:f_strlen(prg_nm)),': rcd = ',rcd,', ',nf_strerror(rcd)
         stop
      endif                     ! endif err
      return
      end                       ! end nf_open_csz()

      integer function nf_close_csz(nc_id,fl_nm,vrb_sng)
c     Purpose: Close netCDF file
c     Prototype:
c     integer nf_close_csz       ! Close netCDF file (libcsz_f77)
      implicit none
c     Commons
#include <dbg.com> /* Debugging constants */
#include <netcdf.inc>
c     Input
      character*(*) fl_nm       ! Filename
      character*(*) vrb_sng     ! Verb string
      integer nc_id             ! netCDF file handle
c     Local
      integer rcd               ! [rcd] Return success code
c     External
      integer f_strlen            ! f_strlen() is defined in libcsz
c     Main code
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf_close(nc_id)
      nf_close_csz=rcd
      if (rcd /= NF_NOERR) then
         write (6,'(3a)') prg_nm(1:f_strlen(prg_nm)),': ERROR closing ',fl_nm(1:f_strlen(fl_nm))
         write (6,'(2a,i3,2a)') prg_nm(1:f_strlen(prg_nm)),': rcd = ',rcd,', ',nf_strerror(rcd)
         stop
      endif                     ! endif err
      write (6,'(a,1x,a)') vrb_sng,fl_nm(1:f_strlen(fl_nm))
      return
      end                       ! end nf_close_csz()

