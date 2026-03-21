c $Id$
      subroutine datecom( ncdate, ncsec, com )

c Build string containing day and date interval info.

      use precision

      implicit none

c Input arguments:

      integer
     $  ncdate  ! current date
     $, ncsec   ! seconds relative to ncdate
      

c Output variables:

      character
     $  com*(*)  ! comment field for mass storage system

c Local variables:

      integer
     $  i
     $, lenc   ! accumulate length of com string as it is built
     $, kblen  ! length of strings from i2char
     $, khour
      real(r8)
     $  hour
      character
     $  nctemp*20  ! string from i2char
     $, chtem*80   ! string for temporary copies of com

c Externals:

      character i2char*20
      external
     $  i2char  ! convert integer to left justified ascii format
c-----------------------------------------------------------------------

      lenc = 0

      hour = float(ncsec)/3600.
      khour = int(hour)
      nctemp = i2char( khour, kblen )
      com = 'DATE: ' // nctemp(1:kblen) // '.'
      lenc = lenc + kblen + 7
      khour = int((hour - float(khour))*1000. + 1000.)
      nctemp = i2char( khour, kblen )
      chtem = com(1:lenc) // nctemp(2:4) // 'Z '
      lenc = lenc + 5
      nctemp = i2char( ncdate, kblen )
      com = chtem(1:lenc) // nctemp(1:kblen)
      lenc = lenc + kblen

      if ( lenc .lt. len(com) ) then
         do i = lenc+1, len(com)
            com(i:i) = ' '
         enddo
      endif

      return
      end

c##############################################################################

      subroutine ddicom( ndcurf, nscurf, ndcur, nscur, ncdatf, ncsecf
     $,                  ncdate, ncsec, com )

c Build string containing day and date interval info.

      use precision

      implicit none

c Input arguments:

      integer
     $  ndcurf  ! first day
     $, nscurf  ! seconds relative to ndcurf
     $, ndcur   ! current day
     $, nscur   ! seconds relative to ndcur
     $, ncdatf  ! first date
     $, ncsecf  ! seconds relative to ncdatf
     $, ncdate  ! current date
     $, ncsec   ! seconds relative to ncdate
      

c Output variables:

      character
     $  com*(*)  ! comment field for mass storage system

c Local variables:

      integer
     $  i
     $, lenc   ! accumulate length of com string as it is built
     $, kblen  ! length of strings from i2char
     $, kdayno ! 3 digit fractional part of day converted to an integer
     $, khour
      real(r8)
     $  hour
      character
     $  nctemp*20  ! string from i2char
     $, chtem*80   ! string for temporary copies of com

c Externals:

      character i2char*20
      external
     $  i2char  ! convert integer to left justified ascii format
c-----------------------------------------------------------------------

      do i = 1, len( com )
         com(i:i) = ' '
      end do
      lenc = 0

c     Number of first day in file.

      nctemp = i2char( ndcurf, kblen )
      com = 'DAYS: ' // nctemp(1:kblen) // '.'
      lenc = lenc + kblen + 7
      kdayno = int((float(nscurf)/86400.)*1000.) + 1000.
      nctemp = i2char( kdayno, kblen )
      chtem = com(1:lenc) // nctemp(2:4) // '-'
      lenc = lenc + 4

c     Number of last day in file.

      nctemp = i2char( ndcur, kblen )
      com = chtem(1:lenc) // nctemp(1:kblen) // '.'
      lenc = lenc + kblen + 1
      kdayno = int((float(nscur)/86400.)*1000.) + 1000.
      nctemp = i2char( kdayno, kblen )
      chtem = com(1:lenc) // nctemp(2:4) // ' DATES: '
      lenc = lenc + 11

c     First date on file.

      hour = float(ncsecf)/3600.
      khour = int(hour)
      nctemp = i2char( khour, kblen )
      com = chtem(1:lenc) // nctemp(1:kblen) // '.'
      lenc = lenc + kblen + 1
      khour = int((hour - float(khour))*1000. + 1000.)
      nctemp = i2char( khour, kblen )
      chtem = com(1:lenc) // nctemp(2:4) // 'Z '
      lenc = lenc + 5
      nctemp = i2char( ncdatf, kblen )
      com = chtem(1:lenc) // nctemp(1:kblen) // ' - '
      lenc = lenc + kblen + 3

c     Last date on file.

      hour = float(ncsec)/3600.
      khour = int(hour)
      nctemp = i2char( khour, kblen )
      chtem = com(1:lenc) // nctemp(1:kblen) // '.'
      lenc = lenc + kblen + 1
      khour = int((hour - float(khour))*1000. + 1000.)
      nctemp = i2char( khour, kblen )
      com = chtem(1:lenc) // nctemp(2:4) // 'Z '
      lenc = lenc + 5
      nctemp = i2char( ncdate, kblen )
      chtem = com(1:lenc) // nctemp(1:kblen)
      lenc = lenc + kblen

      com = chtem(1:lenc)

      return
      end
