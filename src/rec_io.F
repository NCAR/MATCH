      module rec_io

      use precision

      contains
c $Id$
      subroutine rdarr( iu, iln, buf )

c Read an array from a sequential file.

      implicit none

c Input arguments:

      integer
     $  iu    ! Fortran logical unit
     $, iln   

      real(r8)
     $  buf(iln)
c-----------------------------------------------------------------------

      read( iu ) buf

      return
      end subroutine rdarr

c#############################################################################

      subroutine wrarr( iu, iln, buf )

c Write an array to a file.

      implicit none

c Input arguments:

      integer
     $  iu    ! Fortran logical unit
     $, iln   

      real(r8)
     $  buf(iln)
c-----------------------------------------------------------------------

      write(iu) buf

      return
      end subroutine wrarr

c#############################################################################

      subroutine wrintrec( iu, iln, buf )

c Write an array to a file.

      implicit none

c Input arguments:

      integer
     $  iu    ! Fortran logical unit
     $, iln   

      integer
     $  buf(iln)
c-----------------------------------------------------------------------

      write(iu) buf

      return
      end subroutine wrintrec

c#############################################################################

      subroutine wrch8rec( iu, iln, buf )

c Write an array to a file.

      implicit none

c Input arguments:

      integer
     $  iu    ! Fortran logical unit
     $, iln   

      character
     $  buf(iln)*8
c-----------------------------------------------------------------------

      write(iu) buf

      return
      end subroutine wrch8rec

c#############################################################################

      subroutine readda( lun, irec, lenbuf, buf )
      integer lun, irec, lenbuf
      real(r8) buf(lenbuf)
      read( lun, rec=irec ) buf
      return
      end subroutine readda

c#############################################################################

      subroutine writeda( lun, irec, lenbuf, buf )
      integer lun, irec, lenbuf
      real(r8) buf(lenbuf)
      write( lun, rec=irec ) buf
      return
      end subroutine writeda

c#############################################################################

      end module rec_io

      
