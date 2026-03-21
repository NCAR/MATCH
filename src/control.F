c $Id$
      subroutine inictl( itimrst, rstdate, rstdatesec, rstday
     $,                  rstdaysec )

c Initialize control variables.
c N.B.1.  The time variables for the initial conditions are initialized in 
c         inirun() and potentially modified by the namelist input.
c N.B.2.  Date calculations are performed which assumes that the calendar
c         module has been initialized.

      use histout, only: qhsti
      use calendar
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------

c Input argument:
      integer
     $  itimrst      ! time index of restart data
     $, rstdate      ! date of restart data
     $, rstdatesec   ! seconds relative to rstdate
     $, rstday       ! day of restart data
     $, rstdaysec    ! seconds relative to rstday

c Local variables:
      integer ret
c-----------------------------------------------------------------------

c     Time index at beginning of either initial or restart run.
      itim0 = 0
      if ( rstflg .eq. 1 ) itim0 = itimrst

c     Init current timestep to be the one that ended at time index itim0.
      curstep = itim0

c     Current date/time info.
      if ( rstflg .eq. 0 ) then
         curdate = icdate
         curdatesec = icdatesec
         curday = icday
         curdaysec = icdaysec
      else
         curdate = rstdate
         curdatesec = rstdatesec
         curday = rstday
         curdaysec = rstdaysec
      end if

c     Set base date.
      call addsec2dat( -(icday*86400+icdaysec), icdate, icdatesec
     $,                basedate, basedatesec )

      return
      end

c#############################################################################

      subroutine ctlstep( nstep )

c Set control variables for this timestep.

      use histout, only: qhsti
      use calendar
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  nstep       ! current timestep

c Local variables:
      integer
     $  itim        ! time index at the end of timestep nstep
     $, istat       ! return status
     $, nhtfrq      ! frequency of history time sample writes in timesteps
     $, ntshst      ! number of time samples currently in history file
     $, mxtshst     ! maximum number of time samples in history file
     $, rstfrq      ! frequency of restart file writes in timesteps
     $, dd          ! month day of curdate
      logical
     $  shutdown    ! true => shut down run this timestep

c Externals:
      integer qrsti
c------------------------------------------------------------------------------

c     Update current date/day variables.
      call setcurtime( nstep )

      itim = nstep

      begstep = .false.
      if ( nstep .eq. 1 ) begstep = .true.

      laststep = .false.
      if ( nstep .eq. nestep ) laststep = .true.

c     endofrun is set true either because all
c     requested timesteps have been completed, or because the job needs to
c     shut down for external reasons (e.g. operator could set a switch
c     requesting the job to shut down, or the job might check to see if enough
c     cpu time remains to allow another timestep to complete)

      shutdown = .false.
      endofrun = .false.
      if ( laststep .or. shutdown ) endofrun = .true.

c     Is it time to write a time sample to the history file?
      istat = qhsti( 'NHTFRQ', nhtfrq )
      if ( istat .ne. 0 ) then
         write(*,*)'ctlstep: qhsti returned error for NHTFRQ'
         stop
      end if
      wrhstts = .false.
      if ( nhtfrq .gt. 0 ) then
         if ( mod( itim, nhtfrq ).eq.0 ) wrhstts = .true.
      else
         dd = mod( curdate, 100 )
         if ( dd .eq. 1 .and. curdatesec .eq. 0 ) wrhstts = .true.
      end if

c     Is it time to close the history file?  Only set this flag when a
c     time sample is to be written.

      istat = qhsti( 'NTSHST', ntshst )
      if ( istat .ne. 0 ) then
         write(*,*)'ctlstep: qhsti returned error for NTSHST'
         stop
      end if
      istat = qhsti( 'MXTSHST', mxtshst )
      if ( istat .ne. 0 ) then
         write(*,*)'ctlstep: qhsti returned error for MXTSHST'
         stop
      end if

      fullhst = .false.
      if ( wrhstts .and. ntshst .eq. mxtshst-1 ) fullhst = .true.
      closehst = .false.
      if ( fullhst  .or.  laststep ) closehst = .true.

c     Is it time to write a restart file?
      istat = qrsti( 'RSTFRQ', rstfrq )
      if ( istat .ne. 0 ) then
         write(*,*)'ctlstep: qrsti returned error for RSTFRQ'
         stop
      end if
      wrestart = .false.
      if ( rstfrq .ne. 0 ) then
         if ( mod( itim, rstfrq ).eq.0 .or. endofrun ) wrestart = .true.
      else
         if ( wrhstts .or. endofrun ) wrestart = .true.
      end if

      return
      end

c#############################################################################

      subroutine setcurtime( nstep )

c Set the date/time variables valid at the mid-point and the end of the
c input timestep.  Then update the current timestep number.

      use calendar
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  nstep       ! current timestep

c Local varibles:
      integer
     $  totsec      ! Elapsed time to end point of current timestep
     $, mpsec       ! Elapsed time to mid-point of current timestep

c Externals:
!      integer
!     $  newdate
c-----------------------------------------------------------------------

      totsec = (nstep-curstep) * delt
      mpsec = totsec - delt/2

c     Update mid-point date
      mpdate = newdate( curdate, (curdatesec + mpsec)/86400 )
      mpdatesec = mod( curdatesec + mpsec, 86400 )

c     Update current day
      curday = curday + (curdaysec + totsec)/86400
      curdaysec = mod( curdaysec + totsec, 86400 )

c     Update current date
      curdate = newdate( curdate, (curdatesec + totsec)/86400 )
      curdatesec = mod( curdatesec + totsec, 86400 )

c     Update current timestep number.
      curstep = nstep

      return
      end

c#######################################################################

      subroutine getmpdate( date, datesec )

c Return date at mid-point of current timestep.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------

c Output arguments:
      integer
     $  date     ! date at mid-point of current timestep
     $, datesec  ! seconds relative to date
c-----------------------------------------------------------------------
      date = mpdate
      datesec = mpdatesec
      return
      end

c#######################################################################

      subroutine getcurdate( cdate, cdatesec )

c Return date at end point of current timestep.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------

c Output arguments:
      integer
     $  cdate     ! date at end of current timestep
     $, cdatesec  ! seconds relative to cdate
c-----------------------------------------------------------------------
      cdate = curdate
      cdatesec = curdatesec
      return
      end

c#######################################################################

      subroutine getcurday( cday, cdaysec )

c Return current time.

      use calendar
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------

c Output arguments:
      integer
     $  cday      ! day at end of current timestep
     $, cdaysec   ! seconds relative to curday
c-----------------------------------------------------------------------
      cday = curday
      cdaysec = curdaysec
      return
      end

c#######################################################################

      integer function qctli( name, ival )

c Query control module for integer data.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <control.com>
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

      qctli = 0

      if ( name .eq. 'CURDAY' ) then
         ival = curday
         return
      else if ( name .eq. 'CURDAYSEC' ) then
         ival = curdaysec
         return
      else if ( name .eq. 'CURDATE' ) then
         ival = curdate
         return
      else if ( name .eq. 'CURDATESEC' ) then
         ival = curdatesec
         return
      end if

      qctli = -1
      return
      end

c#######################################################################

      logical function progCwat()
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------
      progCwat = dicldphys
      return
      end

c#######################################################################

      logical function dryDep()

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------
      dryDep = dodrydep
      return
      end

c#######################################################################

      logical function altConv()
      use precision
      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------
      altConv = useAltConv
      return
      end

c#######################################################################

      logical function outputTS()
      use precision
      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------
      outputTS = wrhstts
      return
      end

c#######################################################################

      logical function adj_winds()
      use precision
      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------
      adj_winds = wind_ps_adj
      return
      end

c#######################################################################

      logical function calc_etadot()
      use precision
      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------
      calc_etadot = calcedot
      return
      end


c#######################################################################

      logical function AODasm()
      use precision
      implicit none
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------
      AODasm = doAODasm
      return
      end

c#######################################################################
