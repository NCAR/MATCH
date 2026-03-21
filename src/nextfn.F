c $Id$
      integer function nextfn( fn )

c Build next default history tape or regen file name in sequence.  Default
c formats are: cddd, cdddd, or ccdddd, where the c's are non digit characters
c and the d's are digits.

      implicit none

c Input/Output variable:

      character
     $  fn*(*)  ! filename in default format.  On return, the number part
c               ! of the name is incremented by 1

c Return values:

c  0 success
c -1 filename not in default form
c -2 number associated with filename is too large to increment

c Local variables:

      integer
     $  lname  ! number of significant characters in filename
     $, fnno   ! number associated with filename

      character
     $  cdum   ! hold 1st one or two characters of filename
     $, ctmp*5 ! character string for number associated with filename

c Externals:

      integer glc

      external
     $  glc     ! position of last non-blank or non-null character
c-----------------------------------------------------------------------

      lname = glc( fn )

      if ( lname .lt. 4  .or. lname .gt. 6 ) then
         nextfn = -1
         return
      end if

      if ( lname .eq. 4 ) then
         read( fn, '(a1,i3)' ) cdum, fnno
      elseif ( lname .eq. 5 ) then
         read( fn, '(a1,i4)' ) cdum, fnno
      elseif ( lname .eq. 6 ) then
         read( fn, '(a2,i4)' ) cdum, fnno
      endif

      fnno = fnno + 1

      if ( lname .eq. 4 .and. fnno .gt. 999 ) then
         nextfn = -2
         return
      elseif ( lname .ge. 5 .and. fnno .gt. 9999 ) then
         nextfn = -2
         return
      endif

      fnno = 10000 + fnno
      write( ctmp, '(i5)' ) fnno

      if ( lname .eq. 4 ) then
         fn = fn(1:1) // ctmp(3:5)
      elseif ( lname .eq. 5 ) then
         fn = fn(1:1) // ctmp(2:5)
      elseif ( lname .eq. 6 ) then
         fn = fn(1:2) // ctmp(2:5)
      endif

      nextfn = 0

      return
      end
