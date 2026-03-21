c $Id$
      subroutine inirest( rstflg, iunrst0, iunrst1, nrt, xrstfrq
     $,                   xnestep, rrootd, case, lrst0, rrst0
     $,                   lrst1, rrst1, xasync, xrmrst
     $,                   xdiconvccm, xdicldphys, xadvqth
     $,                   itim, date, datesec, day, daysec
     $,                   ntshst, lpthdyn, rpthdyn, lpthhst, rpthhst
     $)

c Initialize /rest0/.  Acquire and read restart file.  Reset namelist
c parameters.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <restart.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  rstflg     ! restart flag
     $, iunrst0    ! unit number for restart info file
     $, iunrst1    ! unit number for restart data file
     $, nrt        ! retention period for restart file
     $, xrstfrq    ! output frequency for restart files
     $, xnestep    ! ending timestep number

      character
     $  rrootd*16  ! Default root directory for files written to remote storage.
     $, case*16    ! case name from the environment variable CASE
     $, lrst0*128  ! local pathname for restart info file
     $, rrst0*128  ! remote pathname for restart info file
     $, lrst1*128  ! local pathname for restart data
     $, rrst1*128  ! remote pathname for restart data

      logical
     $  xasync       ! t => dispose to remote device asynchronously
     $, xrmrst       ! t => remove restart data files after disposing
c     The following flags are passed to the restart module so that it will
c     know to call the corresponding read and write restart file functions
c     for the parameterization module.
     $, xdiconvccm   ! t => diagnose CCM convection
     $, xdicldphys   ! t => cloud fraction calcs
     $, xadvqth      ! t => advect Q and THETA (potential temperature)

c Output arguments:

      integer
     $  itim         ! time index of restart.  (not used on an initial run)
     $, date         ! date of restart data
     $, datesec      ! seconds relative to rstdate
     $, day          ! day of restart data
     $, daysec       ! seconds relative to rstday
     $, ntshst       ! number of time samples on history file

      character
     $  lpthdyn*168  ! local path for dynamics file containing (date,sec).
     $, rpthdyn*168  ! remote path for dynamics file containing (date,sec).
     $, lpthhst*168  ! local path for history file at (date,sec).
     $, rpthhst*168  ! remote path for history file at (date,sec).

c Local variables:

      integer
     $  istat      ! return status

      logical
     $  ldum    ! dummy

c Externals:

      integer attach, glc

      external
     $  attach  ! open file on local disk or acquire from remote storage
     $, glc     ! position of last non-blank or non-null character
c-----------------------------------------------------------------------

c     Set default filenames if values not supplied via namelist
      if ( lrst0(1:1) .eq. ' ' ) lrst0 = 'rstrt'
      if ( rrst0(1:1) .eq. ' ' ) rrst0 = '/' // 
     $   rrootd(:glc(rrootd)) // '/match/' // case(:glc(case)) //
     $   '/rest/rstrt'
      if ( lrst1(1:1) .eq. ' ' ) lrst1 = 'r0000'
      if ( rrst1(1:1) .eq. ' ' ) rrst1 = '/' // 
     $   rrootd(:glc(rrootd)) // '/match/' // case(:glc(case)) //
     $   '/rest/r0000'

      if ( rstflg .eq. 1 ) then

c        Get restart info file.

         istat = attach( lrst0, rrst0, iunrst0, .false., ldum )
         if ( istat .ne. 0 ) then
            write(*,*)'inirest: attach returned ', istat
            stop
         end if

c        Read restart info file.

         read( iunrst0, '(6i10)' ) date, datesec, day,daysec,itim,ntshst
         read( iunrst0, '(a168)' ) lpthdyn
         read( iunrst0, '(a168)' ) rpthdyn
         read( iunrst0, '(a168)' ) lpthhst
         read( iunrst0, '(a168)' ) rpthhst
         read( iunrst0, '(a128)' ) lrst1
         read( iunrst0, '(a128)' ) rrst1

      end if

c     Set variables in /rest0/.
      lunrst0  = iunrst0
      lunrst1  = iunrst1
      rstrrt   = nrt
      rstfrq   = xrstfrq
      nestep   = xnestep
      async    = xasync
      rmrst    = xrmrst
      lpthrst0 = lrst0
      rpthrst0 = rrst0
      lpthrst1 = lrst1
      rpthrst1 = rrst1
      diconvccm= xdiconvccm
      dicldphys= xdicldphys
      advqth   = xadvqth

      if ( rstflg .eq. 0 ) then
c        Open restart info file.
         open( lunrst0, file=lpthrst0(:glc(lpthrst0)), status='unknown')
      end if

      return
      end

c##############################################################################

      subroutine wrrst( nstep, closehst, as, x, psdry )

c Write and dispose restart info and data files if it is time.

      use advect, only: wadvrst
      use dyninp, only: qdyn
      use histout, only: qhstc, qhsti, disphst, appndhst, whstrst
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <restart.com>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  nstep        ! Model timestep number.
      logical
     $  closehst     ! true => close history file (flush output buffers)
      real(r8)
     $  as(plond,plev,pcnst,plat) ! advected species
     $, x(plond,plev,padv,platd)  ! misc advected fields
     $, psdry(plond,plat)         ! surface pressure of dry atm

c Local variables:
      integer
     $  istat   ! return status
     $, ncdate  ! current date in yymmdd format
     $, ncsec   ! seconds relative to ncdate
     $, ncday   ! current day
     $, ncdsec  ! seconds relative to ncday
     $, ntsdyn  ! number of most recently read dynamics time sample
     $, ntshst  ! number of time samples written on current history file
     $, mfilt   ! maximum number of time samples on a history file
      character
     $  lpthdyn*168 ! local path for dynamics file
     $, rpthdyn*168 ! remote path for dynamics file
     $, lpthhst*168 ! local path for history file
     $, rpthhst*168 ! remote path for history file
     $,  rpthrst0lbl*168  ! remote path for rstrt with label
     $, com*80      ! comment field
     $, lndir
     $, msvol
     $, wpasswd
      logical
     $  isNC     ! .t. => filename has '.nc' extension
     $, remove   ! .true. => remove local file after disposing to remote device
     $, link     ! .true. => link local file to lndir

c Externals:
      integer glc, incstr
      logical
     $  adj_winds       ! adjust winds to enforce mass conservation
c-----------------------------------------------------------------------

c     Update restart info file.
      rewind( lunrst0 )

c     Current date/time.
      call getcurdate( ncdate, ncsec )
      call getcurday( ncday, ncdsec )

c     Get name of dynamics file for current time.  If the most recently read
c     time sample was the first on the file then must decrement file number
c     to get correct starting file.
      call qdyn( lpthdyn, rpthdyn, ntsdyn )
      if ( ntsdyn .eq. 1 ) then
c        Is dynamics from a netCDF file?
         isNC = lpthdyn(glc(lpthdyn)-2:glc(lpthdyn)) == '.nc'
         if ( isNC ) then
            istat = incstr( lpthdyn(:glc(lpthdyn)-3), -1 )
         else
            istat = incstr( lpthdyn, -1 )
         end if
         if ( istat .ne. 0 ) then
            write(*,*)'wrrst: incstr returned ', istat
            write(*,*)'       while trying to decrement ', lpthdyn
            stop
         end if
         if ( glc( rpthdyn ) .gt. 0 ) then
            if ( isNC ) then
               istat = incstr( rpthdyn(:glc(rpthdyn)-3), -1 )
            else
               istat = incstr( rpthdyn, -1 )
            end if
            if ( istat .ne. 0 ) then
               write(*,*)'wrrst: incstr returned ', istat
               write(*,*)'       while trying to decrement ', rpthdyn
               stop
            end if
         end if
      end if

c     Get history file paths and number of time samples currently written.
      istat = qhstc( 'LPATH', lpthhst )
      if ( istat .ne. 0 ) then
         write(*,*)'wrrst: qhstc returned error for LPTHHST'
         stop
      end if
      istat = qhstc( 'RPATH', rpthhst )
      if ( istat .ne. 0 ) then
         write(*,*)'wrrst: qhstc returned error for RPTHHST'
         stop
      end if
      istat = qhsti( 'NTSHST', ntshst )
      if ( istat .ne. 0 ) then
         write(*,*)'wrrst: qhsti returned error for NTSHST'
         stop
      end if
      istat = qhsti( 'MXTSHST', mfilt )
      if ( istat .ne. 0 ) then
         write(*,*)'wrrst: qhsti returned error for MXTSHST'
         stop
      end if

c     If this history file has not been closed, then close it but don't
c     remove it (noremove=.true.).  The history file must be re-opened and
c     positioned at the end.
      if ( .not. closehst ) then
         call disphst( .false., .true. )
         call appndhst()
      end if

c     Increment the restart data filenames.
      istat = incstr( lpthrst1, 1 )
      if ( istat .ne. 0 ) then
         write(*,*)'wrrst: incstr returned ', istat
         write(*,*)'       while trying to increment ', lpthrst1
         stop
      end if
      istat = incstr( rpthrst1, 1 )
      if ( istat .ne. 0 ) then
         write(*,*)'wrrst: incstr returned ', istat
         write(*,*)'       while trying to increment ', rpthrst1
         stop
      end if

c     Write restart info file.

      write(*,*)'wrrst: writing restart file at date = ',ncdate
     $,         '  seconds = ', ncsec

      write( lunrst0, '(6i10)' ) ncdate, ncsec,ncday,ncdsec,nstep,ntshst
      write( lunrst0, '(a168)' ) lpthdyn
      write( lunrst0, '(a168)' ) rpthdyn
      write( lunrst0, '(a168)' ) lpthhst
      write( lunrst0, '(a168)' ) rpthhst
      write( lunrst0, '(a128)' ) lpthrst1
      write( lunrst0, '(a128)' ) rpthrst1

c     Dispose restart info file.

      remove = .false.
      link = .false.
      lndir = ' '
      msvol = ' '
      wpasswd = ' '
      com = ' '

c     Label the remote restart information file with the name of the restart data file
      rpthrst0lbl = trim(rpthrst0) // '_' // rpthrst1(glc(rpthrst1)-4:) 
      write(6,*)'remote rstrt filename ', rpthrst0lbl

      call dispose( lunrst0, lpthrst0, rpthrst0lbl, rstrrt
     $,             async, com, wpasswd, remove )

c     Re-open restart info file (it is closed in dispose).

      open( lunrst0, file=lpthrst0(:glc(lpthrst0)), status='old' )

c     Open new restart data file.

      open( lunrst1, file=lpthrst1(:glc(lpthrst1)), status='unknown'
     $,     form='unformatted' )

c     Write restart data.

c     Advection algorithm.
      call wadvrst( lunrst1 )

c     History buffer and header info.
      call whstrst( lunrst1 )

c     Advected species.
      write( lunrst1 ) as

c     Dry mass of atmosphere
      call wdryrst( lunrst1 )

c     Miscellaneous advected fields.
      if ( advqth ) write( lunrst1 ) x

c     Surface pressure.
      if ( adj_winds() ) write( lunrst1 ) psdry

c     Dispose restart data file.
      remove = rmrst
      if ( nstep .eq. nestep ) remove = .false.
      call datecom( ncdate, ncsec, com )
      call dispose( lunrst1, lpthrst1, rpthrst1, rstrrt
     $,             async, com, wpasswd, remove )

      return
      end

c##############################################################################

      subroutine rdrst( as, x, psdry )

c Read restart data file.

      use advect, only: radvrst
      use histout, only: rhstrst
      use precision


      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <restart.com>
c-----------------------------------------------------------------------

c Output argument:
      real(r8)
     $  as(plond,plev,pcnst,plat) ! advected species
     $, x(plond,plev,padv,platd)  ! miscellaneous advected fields
     $, psdry(plond,plat)         ! surface pressure of dry atm

c Local variables:
      integer
     $  istat   ! return status
      character
     $  cmd*256     !  Command string
      logical
     $  cosbrst  ! .true. => file is COS blocked.  Not used since can only
c                ! do a restart on machine that is binary compatible with
c                ! machine that wrote the restart file.

c Externals:
      integer attach, glc, fsystem
      logical
     $  adj_winds       ! adjust winds to enforce mass conservation
c-----------------------------------------------------------------------

c     Get restart data file.

      istat = attach( lpthrst1, rpthrst1, lunrst1, .true., cosbrst )
      if ( istat .ne. 0 ) then
         write(*,*)'rdrst: attach returned ', istat
         stop
      end if

c     Read restart data.

c     Advection algorithm.
      call radvrst( lunrst1 )

c     History buffer and header info.
      call rhstrst( lunrst1 )

c     Advected species.
      read( lunrst1 ) as

c     Dry mass of atmosphere
      call rdryrst( lunrst1 )

c     Miscellaneous advected fields.
      if ( advqth ) read( lunrst1 ) x

c     Surface pressure.
      if ( adj_winds() ) read( lunrst1 ) psdry


      return
      end

c##############################################################################

      integer function qrsti( name, ival )

c Query restart module for integer data.

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <restart.com>
c-----------------------------------------------------------------------

c Input argument:

      character
     $  name*(*)  ! variable name of requested integer data

c Output argument:

      integer
     $  ival      ! requested value

c Return value:

c  0 => success
c -1 => variable name not found
c-----------------------------------------------------------------------

      qrsti = 0

      if ( name .eq. 'RSTFRQ' ) then
         ival = rstfrq
         return
      end if

      qrsti = -1
      return
      end

c##############################################################################

      subroutine prtrst( )

c Print contents of /rest0/.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <restart.com>
c-----------------------------------------------------------------------

      write(*,*)'/rest0/ contains:'
      write(*,*)'lunrst0=',lunrst0
      write(*,*)'lunrst1=',lunrst1
      write(*,*)'rstrrt=',rstrrt
      write(*,*)'rstfrq=',rstfrq
      write(*,*)'nestep=',nestep
      write(*,*)'async=',async
      write(*,*)'rmrst=',rmrst
      write(*,*)'lpthrst0=',lpthrst0
      write(*,*)'rpthrst0=',rpthrst0
      write(*,*)'lpthrst1=',lpthrst1
      write(*,*)'rpthrst1=',rpthrst1

      return
      end
