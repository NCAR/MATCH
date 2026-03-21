c $Id$
      integer function toklen( cs )

c Return token length, i.e., starting from begining of string return length
c of first token as delimited by either a blank or a null byte.

      implicit none

c Input argument:

      character*(*) cs       !  Input character string

c Local variables:

      integer m, n
c-----------------------------------------------------------------------

      m = 0
      n = 1

   10 continue

      if ( cs(n:n) .eq. ' ' .or. cs(n:n) .eq. char(0) ) then
         toklen = m
         return
      else
         m = m + 1
         n = m + 1
         go to 10
      end if

      return
      end                                                   ! TOKLEN
c==============================================================================

      integer function lastsl( cs )

c Position of last slash (/) in the first input token.

      implicit none

c Input argument:

      character*(*) cs       !  Input character string

c Return values:
c     > 0  => position of last slash
c     = 0  => slash not found in token

c Local variables:

      integer n

c Externals:

      integer toklen
      external toklen
c-----------------------------------------------------------------------

      n = toklen( cs )
      if ( n .eq. 0 ) then
         lastsl = 0
         return
      end if

   10 continue

      if ( cs(n:n) .eq. '/' ) then
         lastsl = n
         return
      else
         n = n - 1
         if ( n .gt. 0 ) then
            go to 10
         else
            lastsl = 0
            return
         end if
      end if

      return
      end                                               ! LASTSL
c==============================================================================

      integer function glc( cs )

c Position of last significant character in string.  Here significant means
c non-blank or non-null.

      implicit none

c Input argument:

      character*(*) cs       !  Input character string

c Return values:
c     > 0  => position of last significant character
c     = 0  => no significant characters in string

c Local variables:

      integer i, n
c-----------------------------------------------------------------------

      n = len( cs )
      if ( n .eq. 0 ) then
         glc = 0
         return
      end if

   10 continue

      if ( cs(n:n) .ne. ' ' .and. cs(n:n) .ne. char(0) ) then
         glc = n
         return
      else
         n = n - 1
         if ( n .gt. 0 ) then
            go to 10
         else
            glc = 0
            return
         end if
      end if

      return
      end                                                 ! GLC
c==============================================================================

      character*20 function i2char(kint,kilen)

c Convert integer kint to character type, left-justify and return number
C of characters in i2char in kilen

      implicit none

c Input arguments:

      integer
     $  kint   ! integer to be converted

c Output arguments:

      integer
     $  kilen  ! length of the ascii representation of kint

c Return value:
c The ascii representation of kint

c Local variables:

      character
     $  chtem*20
      integer
     $  i

c Externals:

      integer glc
      external
     $  glc     ! position of last significant character
c-----------------------------------------------------------------------

      write( chtem, '(i20)' ) kint
      i=20
10    continue
      if ( chtem(i:i) .ne. ' ' ) then
         i = i - 1
         if ( i .lt. 1 ) then
            i2char = chtem
            go to 20
         endif
         go to 10
      else
         i2char = chtem(i+1:20)
      endif
20    continue
      kilen = glc( i2char )
      return
      end                                                  ! I2CHAR
c==============================================================================

      character*(*) function upcase(lstring)

c Convert character string lstring to upper case

      implicit none
      character*(*) lstring
      integer glc, i
      external glc
c-----------------------------------------------------------------------
      upcase = ' '
      do i = 1, glc(upcase)
         upcase(i:i) = ' '
      end do
      do i = 1, glc(lstring)
         upcase(i:i) = lstring(i:i)
         if (ichar(lstring(i:i)).ge.97.and.ichar(lstring(i:i)).le.122)
     $       upcase(i:i) = char(ichar(lstring(i:i)) - 32)
      end do
      return
      end                                                   ! UPCASE
c==============================================================================


      integer function incstr( s, inc )


c Increment a string whose ending characters are digits.  The incremented
c integer must be in the range [0 - (10**n)-1] where n is the number of
c trailing digits.

      implicit none

c Input variable:

      integer
     $  inc    ! value to increment string (may be negative)

c Input/Output variable:

      character
     $  s*(*)  ! string with trailing digits

c Return values:

c  0 success
c -1 error: no trailing digits in string
c -2 error: incremented integer is out of range

c Local variables:

      integer
     $  i      ! index
     $, lstr   ! number of significant characters in string
     $, lnd    ! position of last non-digit
     $, ndigit ! number of trailing digits
     $, ival   ! integer value of trailing digits
     $, pow    ! power of ten
     $, digit  ! integer value of a single digit

c Externals:

      integer glc, lastnd, incstrdate
      logical isnum,temp

      external
     $  glc     ! position of last non-blank or non-null character
     $, lastnd  ! position of last non-digit in string
     $, isnum   ! true if tail of string is all integer (to treat as a date)
     $, incstrdate ! called instead to increment date if isnum is true

c-----------------------------------------------------------------------

      if ( isnum(s) ) then
c-drb         write(6,*)'incstr: ',s,'incremeneted as a date'
         incstr = incstrdate(s,inc)
         return
      end if

c-drb      write(6,*)'incstr: file name incremented as a non-date integer'
      lstr = glc( s )
      lnd = lastnd( s )
      ndigit = lstr - lnd

      if ( ndigit .eq. 0 ) then
         incstr = -1
         return
      end if

c     Calculate integer corresponding to trailing digits.

      ival = 0
      pow = 0
      do i = lstr, lnd+1, -1
         digit = ichar(s(i:i)) - ichar('0')
         ival = ival + digit * 10**pow
         pow = pow + 1
      end do

c     Increment the integer.

      ival = ival + inc
      if ( ival .lt. 0 .or. ival .gt. 10**ndigit-1 ) then
         incstr = -2
         return
      end if

c     Overwrite trailing digits.

      pow = ndigit
      do i = lnd+1, lstr
         digit = mod( ival, 10**pow ) / 10**(pow-1)
         s(i:i) = char( ichar('0') + digit )
         pow = pow - 1
      end do

      incstr = 0

      return
      end                                                   ! INCSTR
c==============================================================================

      integer function incstrdate( s, inc )

c incrememt the date at the end of a string.

      use calendar

      implicit none

c Input variable:

      integer
     $  inc    ! value to increment string (may be negative)
               ! corresponds to number of DAYS


c Input/Output variable:

      character
     $  s*(*)  ! string with trailing digits

c Return values:

c  0 success
c -1 error: no trailing digits in string
c -2 error: incremented integer is out of range ??

c Local variables:

      integer
     $  i      ! index
     $, lstr   ! number of significant characters in string
     $, lnd    ! position of last non-digit
     $, ndigit ! number of trailing digits  = 8 here
     $, ival   ! integer value of trailing digits
     $, pow    ! power of ten
     $, digit  ! integer value of a single digit

c Externals:
      integer glc, lastnd
      logical isnum,temp

      external
     $  glc     ! position of last non-blank or non-null character
     $, lastnd  ! position of last non-digit in string
     $, isnum   
c-----------------------------------------------------------------------



      lstr = glc( s )
      lnd = lastnd( s )
      ndigit = lstr - lnd

      if ( ndigit .eq. 0 ) then
         incstrdate = -1
         return
      end if

c     Calculate integer corresponding to trailing digits.
      ival = 0
      pow = 0
      do i = lstr, lnd+1, -1
         digit = ichar(s(i:i)) - ichar('0')
         ival = ival + digit * 10**pow
         pow = pow + 1
      end do

c     Increment the integer.

      ival = newdate( ival, inc )
      if ( ival .lt. 0 .or. ival .gt. 10**ndigit-1 ) then
         incstrdate = -2
         return
      end if

c     Overwrite trailing digits.

      pow = ndigit
      do i = lnd+1, lstr
         digit = mod( ival, 10**pow ) / 10**(pow-1)
         s(i:i) = char( ichar('0') + digit )
         pow = pow - 1
      end do

      incstrdate = 0



      return
      end


c==============================================================================

      logical function isnum(str)
c     returns true if the end of the character string is interpretable as an integer
c     otherwise returns false

      implicit none

      character*(*) str
      character*(32) str2
      integer iwork,glc 

c Externals
      integer lastsl

      str2 = str(lastsl(str)+1:glc(str))

      read(unit=str2,fmt='(i20)',err=10) iwork      
   55 format(i20)
      isnum = .true.
      return
   10 continue
      isnum = .false.
      return
      end

c==============================================================================

      integer function lastnd( cs )

c Position of last non-digit in the first input token.

      implicit none

c Input argument:

      character*(*) cs       !  Input character string

c Return values:
c     > 0  => position of last non-digit
c     = 0  => token is all digits (or empty)

c Local variables:

      integer n, nn, digit

c Externals:

      integer glc

      external
     $  glc     ! position of last non-blank or non-null character
c-----------------------------------------------------------------------

      n = glc( cs )
      if ( n .eq. 0 ) then     ! empty string
         lastnd = 0
         return
      end if

      do nn = n, 1, -1

         digit = ichar(cs(nn:nn)) - ichar('0')

         if ( digit .lt. 0 .or. digit .gt. 9 ) then
            lastnd = nn
            return
         end if
      end do

      lastnd = 0    ! all characters are digits

      return
      end

c###############################################################################

      character*10 function idate2char( idate )

c Convert integer date in YYYYMMDD to a character array in the form
c 'YYYY-MM-DD'

      implicit none

c Input arguments:

      integer
     $  idate

c Local variables:

      integer
     $  yy, mm, dd
      character
     $  dash
c-----------------------------------------------------------------------

      dash = '-'

c     extract year, month, and day from date
      yy = idate / 10000
      mm = mod( abs(idate), 10000 ) / 100
      dd = mod( abs(idate), 100 )

      write( idate2char, '(i4.4,a1,i2.2,a1,i2.2)' ) yy,dash,mm,dash,dd

      return
      end

c###############################################################################

      character*8 function isec2char( isec )

c Convert integer seconds to a character array in the form 'HH:MM:SS'.

      implicit none

c Input arguments:

      integer
     $  isec

c Local variables:

      integer
     $  hh, mm, ss
      character
     $  colon
c-----------------------------------------------------------------------

      colon = ':'

c     extract hour, minutes, seconds
      hh = isec / 3600
      mm = mod( isec, 3600 ) / 60
      ss = mod( isec, 60 )

      write( isec2char, '(i2.2,a1,i2.2,a1,i2.2)' ) hh,colon,mm,colon,ss

      return
      end

c###############################################################################

