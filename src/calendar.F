      module calendar

      use precision

      contains
c $Id$
      subroutine inicalendar( xtype )

c Initialize calendar module common block.

      implicit none
c-----------------------------------------------------------------------
#include <calendar.com>
c-----------------------------------------------------------------------

c Input argument:
      character*16
     $  xtype  ! calendar type.  Currently '365' or 'gregorian'
c-----------------------------------------------------------------------
      type = xtype

      if ( type(1:3) .ne. '365' .and. type(1:9) .ne. 'gregorian' ) then
         write(*,*)'inicalendar:FATAL ERROR: unrecognized calendar type'
         stop
      end if

      return
      end subroutine inicalendar

c##############################################################################

      real(DBLKIND) function diffdat( dat1, sec1, dat2, sec2 )

c Compute the difference: (dat2,sec2) - (dat1,sec1)  in days.

      implicit none
c-----------------------------------------------------------------------
#include <calendar.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  dat1   ! date in yyyymmdd format
     $, sec1   ! seconds relative to dat1
     $, dat2   ! date in yyyymmdd format
     $, sec2   ! seconds relative to dat2

c Return value:
c (dat2,sec2) - (dat1,sec1)  in days.

c Externals:
!      real(DBLKIND) diffdat365, diffdatgrg
!      external diffdat365, diffdatgrg
c-----------------------------------------------------------------------

      if ( type(1:3) .eq. '365' ) then
         diffdat = diffdat365( dat1, sec1, dat2, sec2 )
      else if ( type(1:9) .eq. 'gregorian' ) then
         diffdat = diffdatgrg( dat1, sec1, dat2, sec2 )
      else
         write(*,*)'diffdat: FATAL ERROR: calendar not initialized'
         stop
      end if

      return
      end function diffdat

c##############################################################################

      real(DBLKIND) function diffdat365( dat1, sec1, dat2, sec2 )

c Compute the difference: (dat2,sec2) - (dat1,sec1)  in days.
c N.B. Assume 1 year = 365 days.

      implicit none

c Input arguments:

      integer
     $  dat1   ! date in yyyymmdd format
     $, sec1   ! seconds relative to dat1
     $, dat2   ! date in yyyymmdd format
     $, sec2   ! seconds relative to dat2

c Return value:

c (dat2,sec2) - (dat1,sec1)  in days.

c Local variables:

      integer lastd, lasts, firstd, firsts
     $, yr2, mo2, dy2, doy2, yr1, mo1, dy1, doy1, ndays
      real(DBLKIND) sign, days

c Externals:

!     integer doy

!      external
!     $  doy     ! day of year (between 1 and 365)
c-----------------------------------------------------------------------

c     Dates equal?

      if ( dat1 .eq. dat2 .and. sec1 .eq. sec2 ) then
         diffdat365 = 0.
         return
      end if

c     Which date is later?

      if ( dat2 .gt. dat1 ) then
         sign = 1.
         lastd = dat2
         lasts = sec2
         firstd = dat1
         firsts = sec1
      else if ( dat2 .lt. dat1 ) then
         sign = -1.
         lastd = dat1
         lasts = sec1
         firstd = dat2
         firsts = sec2
      else
         if ( sec2 .gt. sec1 ) then
            sign = 1.
            lastd = dat2
            lasts = sec2
            firstd = dat1
            firsts = sec1
         else
            sign = -1.
            lastd = dat1
            lasts = sec1
            firstd = dat2
            firsts = sec2
         end if
      end if

c     Compute number of days between lastd and firstd

      yr2 = lastd / 10000
      mo2 = mod( lastd, 10000 ) / 100
      dy2 = mod( lastd, 100 )
      doy2 = doy( mo2, dy2 )

      yr1 = firstd / 10000
      mo1 = mod( firstd, 10000 ) / 100
      dy1 = mod( firstd, 100 )
      doy1 = doy( mo1, dy1 )

      ndays = 365*(yr2-yr1) + doy2 - doy1

c     Adjust for remaining seconds

      days = float( ndays ) + float( lasts - firsts )/86400.

c     Adjust sign

      diffdat365 = sign * days

      return
      end function diffdat365

c##############################################################################

      real(DBLKIND) function diffdatgrg( dat1, sec1, dat2, sec2 )

c Compute the difference: (dat2,sec2) - (dat1,sec1)  in days.
c N.B. Assume Gregorian calendar.

      implicit none

c Input arguments:

      integer
     $  dat1   ! date in yyyymmdd format
     $, sec1   ! seconds relative to dat1
     $, dat2   ! date in yyyymmdd format
     $, sec2   ! seconds relative to dat2

c Return value:

c (dat2,sec2) - (dat1,sec1)  in days.

c Local variables:

      integer lastd, lasts, firstd, firsts
     $, ndays
      real(DBLKIND) sign, days

c Externals:

!     integer greg2jday

!      external
!     $  greg2jday
c-----------------------------------------------------------------------

c     Dates equal?

      if ( dat1 .eq. dat2 .and. sec1 .eq. sec2 ) then
         diffdatgrg = 0.
         return
      end if

c     Which date is later?

      if ( dat2 .gt. dat1 ) then
         sign = 1.
         lastd = dat2
         lasts = sec2
         firstd = dat1
         firsts = sec1
      else if ( dat2 .lt. dat1 ) then
         sign = -1.
         lastd = dat1
         lasts = sec1
         firstd = dat2
         firsts = sec2
      else
         if ( sec2 .gt. sec1 ) then
            sign = 1.
            lastd = dat2
            lasts = sec2
            firstd = dat1
            firsts = sec1
         else
            sign = -1.
            lastd = dat1
            lasts = sec1
            firstd = dat2
            firsts = sec2
         end if
      end if

c     Compute number of days between lastd and firstd
      ndays = greg2jday( lastd ) - greg2jday( firstd )

c     Adjust for remaining seconds
      days = float( ndays ) + float( lasts - firsts )/86400.

c     Adjust sign
      diffdatgrg = sign * days

      return
      end function diffdatgrg

c##############################################################################

      subroutine addsec2dat( sec, dat1, sec1, dat2, sec2 )

c Add seconds to date and compute new date.
c (dat1,sec1) + sec = (dat2,sec2).

      implicit none

c Input arguments:
      integer
     $  sec    ! seconds to add to (dat1,sec1) (can be negative)
     $, dat1   ! date in yyyymmdd format
     $, sec1   ! seconds relative to dat1

c Output arguments:
      integer
     $  dat2   ! date in yyyymmdd format
     $, sec2   ! seconds relative to dat2

c Externals:
!      integer newdate
!      external
!     $  newdate
c-----------------------------------------------------------------------

      if ( sec .eq. 0 ) then
         dat2 = dat1
         sec2 = sec1
         return
      end if

      sec2 = mod( sec1+sec, 86400 )
      dat2 = newdate( dat1, (sec1+sec)/86400 )

      if ( sec2 .lt. 0 ) then
         sec2 = sec2 + 86400
         dat2 = newdate( dat2, -1 )
      end if

      return
      end subroutine addsec2dat

c##############################################################################

      integer function doy( mon, cday )

c Compute day of year ignoring leap years.  Returns values in the range
c [1,365].

      integer mon, cday

      integer jdcon(12)
      data jdcon/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
      save jdcon

      doy = jdcon(mon) + cday

      return
      end function doy

c##############################################################################

      integer function newdate( date, iday )

c Find the date a specified number of days (possibly negative) from the
c given date.

      implicit none
c-----------------------------------------------------------------------
#include <calendar.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $   date  ! starting date in yyyymmdd format (can be negative)
     $,  iday  ! number of days to increment date (can be negative)

c Externals:
!      integer newdate365, newdategrg
!      external
!     $  newdate365, newdategrg
c-----------------------------------------------------------------------

      if ( type(1:3) .eq. '365' ) then
         newdate = newdate365( date, iday )
      else if ( type(1:9) .eq. 'gregorian' ) then
         newdate = newdategrg( date, iday )
      else
         write(*,*)'newdat: FATAL ERROR: calendar not initialized'
         stop
      end if

      return
      end function newdate

c##############################################################################

      integer function newdate365( date, iday )

c Find the date a specified number of days (possibly negative) from the
c given date.
c N.B. Assume 1 year = 365 days.

      implicit none

c Input arguments:

      integer
     $   date  ! starting date in yyyymmdd format (can be negative)
     $,  iday  ! number of days to increment date (can be negative)

c Local variables:

      integer
     $  idate  ! intial value of date (may change sign)
     $, year   ! year of date
     $, month  ! month of date
     $, imonth ! save the initial value of month
     $, day    ! day of date
     $, ndays  ! running count of days to increment date by
     $, i      ! index
     $, cmonth ! current month as we increment date
     $, nd     ! days to end of month

      integer mdays(24)
      data mdays / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
     $             31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
      save mdays

c Externals:

!      integer nd2endm
!      external
!     $  nd2endm  ! number of days to end of month
c-----------------------------------------------------------------------
      
      if ( iday .eq. 0 ) then
         newdate365 = date
         return
      end if

      idate = date
      year = idate / 10000
      if ( idate .lt. 0 ) idate = -idate
      month = mod( idate, 10000 ) / 100
      imonth = month
      day = mod( idate, 100 )
      
      if ( iday  .gt. 0 ) then
         ndays = iday
         year = year + ndays/365
         ndays = mod( ndays, 365 )
         do i = imonth, imonth+12
            cmonth = i  
            if ( i .gt. 12 ) cmonth = i - 12
            nd = nd2endm( cmonth, day )
            if ( ndays .gt. nd ) then
               month = month + 1
               if ( month .gt. 12 ) then
                  month = 1
                  year = year + 1
               end if
               day = 1
               ndays = ndays - nd - 1
               if ( ndays .eq. 0 ) go to 10
            else
               day = day + ndays
               go to 10
            end if
         end do
   10    continue

      else if ( iday .lt. 0 ) then
         ndays = -iday
         year = year - ndays/365
         ndays = mod( ndays, 365 )
         imonth = month
         do i = imonth+12, imonth, -1
            if ( ndays .ge. day ) then
               month = month - 1
               if ( month .lt. 1 ) then
                  month = 12
                  year = year - 1
               end if
               ndays = ndays - day
               day = mdays( month )
               if ( ndays .eq. 0 ) go to 20
            else
               day = day - ndays
               go to 20
            end if
         end do
   20    continue

      end if

      if ( year .ge. 0 ) then
         newdate365 = year*10000 + month*100 + day
      else
         newdate365 = -year*10000 + month*100 + day
         newdate365 = -newdate365
      end if

      return
      end function newdate365

c##############################################################################

      integer function newdategrg( date, iday )

c Find the date a specified number of days (possibly negative) from the
c given date.
c N.B. Assume Gregorian calendar

      implicit none

c Input arguments:

      integer
     $   date  ! starting date in yyyymmdd format (can be negative)
     $,  iday  ! number of days to increment date (can be negative)

c Local variables:
      integer
     $  curjday, newjday

c Externals:
!      integer greg2jday, jday2greg
!      external
!     $  greg2jday, jday2greg
c-----------------------------------------------------------------------
      
      if ( iday .eq. 0 ) then
         newdategrg = date
         return
      end if

      newdategrg = jday2greg( greg2jday( date ) + iday )

      return
      end function newdategrg

c##############################################################################

      integer function nd2endm( mon, day )

c Returns the number of days to the end of the month.  This number added to
c the input arguement day gives the last day of month mon.

      implicit none

      integer mon, day, mdays(12)

      data mdays / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
      save mdays

      if ( mon .lt. 1  .or.  mon .gt. 12 ) then
         write(*,*)'nd2endm: illegal input... mon= ', mon
         stop
      end if

      nd2endm = mdays(mon) - day

      if ( day .lt. 0 ) then
         write(*,*)'nd2endm: illegal input... mon= ',mon,' day= ',day
         stop
      end if

      return
      end function nd2endm

c##############################################################################

      real(r8) function caldayr( idate, isec )

c Calendar day with fractional part.  Returns values in the range
c [1., 366.)

      implicit none
      integer idate, isec, mon, day

      mon = mod( idate, 10000 ) / 100
      day = mod( idate, 100 )

      caldayr = doy( mon, day ) + float( isec )/86400.

      return
      end function caldayr

c##############################################################################

c Currently the routines greg2jday and jday2greg are working in sense of
c j = greg2jday( jday2greg( j ) ) down to jday 1684595 which corresponds
c approx. to the date -1000228.

      integer function greg2jday( date )

c Return Julian day number given Gregorian date.

c Algorithm from Hatcher,D.A., "Simple Formulae for Julian Day Numbers
c and Calendar Dates," Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.

      implicit none
c-----------------------------------------------------------------------
      integer date
      integer yy, mm, dd
      integer ap, mp
      integer y, d, n, g
c-----------------------------------------------------------------------

c     extract year, month, and day from date
      yy = date / 10000
      mm = mod( abs(date), 10000 ) / 100
      dd = mod( abs(date), 100 )

c     modify year and month numbers
      ap = yy - (12-mm)/10
      mp = mod( mm-3, 12 )
      if ( mp .lt. 0 ) mp = mp + 12

c     Julian day
      y = int( 365.25*( ap + 4712 ) )
      d = int( 30.6*mp + .5 )
      n = y + d + dd  + 59
      g = int( .75*int( ap/100 + 49 ) ) - 38
      greg2jday = n - g

      return
      end function greg2jday

c##############################################################################

      integer function jday2greg( day )

c Return Gregorian date given Julian day number.

c Algorithm from Hatcher,D.A., "Simple Formulae for Julian Day Numbers
c and Calendar Dates," Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.

      implicit none
c-----------------------------------------------------------------------
      integer day
      integer g, n
      integer yy, dp, mm, dd
c-----------------------------------------------------------------------

      g = int( .75*int( (day-4479.5)/36524.25 ) + .5 ) - 37
      n = day + g

      yy = int( n/365.25 ) - 4712
      dp = int( mod( n-59.25, 365.25 ) )
      mm = mod( int( (dp+.5)/30.6 )+2, 12 ) + 1
      dd = int( mod( dp+.5, 30.6 ) ) + 1

      jday2greg = sign( abs(yy)*10000 + mm*100 +dd, yy )

      return
      end function jday2greg

c##############################################################################

      end module calendar
