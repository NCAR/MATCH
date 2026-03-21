c     $Header$

c     Purpose: Library of general purpose utility routines

c     Usage: 
c     cd $HOME/f; make -W utl.F OPTS=D libcsz_f77; cd -
c     cd $HOME/f; make -W utl.F libcsz_f77; cd -
c     cd $HOME/f; make OPTS=D libcsz_f77; cd -

      logical function mnt_chk(grd,grd_nbr)
c     Purpose: Return .true. if array is monotonic, .false. if array is not
c     Prototype:
c     logical mnt_ncr_chk       ! Returns true for monotonic arrays (libcsz_f77)

      use precision

      implicit none
c     Commons
#include <dbg.com> /* Debugging constants */
c     Input
      integer grd_nbr
      real(r8) grd(grd_nbr)
c     Local
      integer idx
      logical dcr
      logical ncr

c     Main Code
      ncr=.true.
      if (grd_nbr > 1) then 
         if (grd(2)-grd(1) < 0.0) ncr=.false.
      endif                     ! endif
      dcr=.not.ncr

      if (ncr) then
         do idx=2,grd_nbr
            if (grd(idx)-grd(idx-1) < 0.0) goto 100
         end do                 ! end loop over characters
      endif                     ! endif
      if (dcr) then
         do idx=2,grd_nbr
            if (grd(idx)-grd(idx-1) > 0.0) goto 100
         end do                 ! end loop over characters
      endif                     ! endif

 100  continue
      if (idx == grd_nbr+1) then
         mnt_chk=.true.
      else
         mnt_chk=.false.
         if (dbg_lvl > dbg_io) then
c     Print non-monotonic points
            write (6,'(a)') 'mnt_chk() reports non-monotonic array:'
            if (ncr) write (6,'(a)') 'Array should be increasing'
            if (dcr) write (6,'(a)') 'Array should be decreasing'
            write (6,'(1p,2(a,i4,a,e10.3),0p)') 
     $           'grd(',idx-1,') = ',grd(idx-1),
     $           ', grd(',idx,') = ',grd(idx)
         endif                  ! endif dbg
      endif                     ! endif

      return
      end                       ! end mnt_chk()

      logical function mnt_ncr_chk(grd,grd_nbr)
c     Purpose: Return .true. if array monotonically increases, .false. if not
c     Prototype:
c     logical mnt_ncr_chk       ! True for monotonic increasing arrays (libcsz_f77)

      use precision

      implicit none
c     Commons
c     Input
      integer grd_nbr
      real(r8) grd(grd_nbr)
c     Local
      logical mnt               ! Monotonicity flag
c     External
      logical mnt_chk
c     Main Code
      mnt=mnt_chk(grd,grd_nbr)
      mnt_ncr_chk=.false.
      if (mnt) then 
         if (grd(2)-grd(1) > 0.0) mnt_ncr_chk=.true.
      endif                     ! endif
      return
      end                       ! end mnt_ncr_chk()

      logical function mnt_dcr_chk(grd,grd_nbr)
c     Purpose: Return .true. if array monotonically decreases, .false. if not
c     Prototype:
c     logical mnt_dcr_chk       ! True for monotonic decreasing arrays (libcsz_f77)

      use precision

      implicit none
c     Commons
c     Input
      integer grd_nbr
      real(r8) grd(grd_nbr)
c     Local
      logical mnt               ! Monotonicity flag
c     External
      logical mnt_chk
c     Main Code
      mnt=mnt_chk(grd,grd_nbr)
      mnt_dcr_chk=.false.
      if (mnt) then 
         if (grd(1)-grd(2) > 0.0) mnt_dcr_chk=.true.
      endif                     ! endif
      return
      end                       ! end mnt_dcr_chk()

      subroutine date_time_get(lcl_date_time)
c     Purpose: Return date and time information

      use precision

      implicit none
c     Parameters
c     Commons
c     Input
c     Output 
      character*26::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
c     Local workspace
      character*8 date
      character*10 time
      character*5 zone
      integer time_val(8)
c Main code

      call date_and_time(date,time,zone,time_val) ! Fortran90 DATE_AND_TIME intrinsic
      write (lcl_date_time,'(
     $     i4.4,a1,i2.2,a1,i2.2,a1,
     $     i2.2,a1,i2.2,a1,i2.2,a1,i3.3
     $     )') 
     $     time_val(1),'/',time_val(2),'/',time_val(3),' ',
     $     time_val(5),':',time_val(6),':',time_val(7),'.',time_val(8)

      return
      end                       ! end date_time_get()
