#ifdef AIX
#define SYSCALL fsystem
#else
#define SYSCALL system
#endif
c $Id$
      integer function getfile( lpath, rpath )

c Check to see if specified local file exists.  If not, acquire it from
c remote storage.
c 
c This routine calls locally written NCAR specific code for
c acquiring files from the NCAR Mass Storage System (msrcp).
      
      implicit none

c Input arguments.
      character(len=*), intent(in) ::
     $  lpath  ! pathname for local file
     $, rpath  ! pathname for remote file

c Return values:
c   0 => successful return
c  -2 => error reading remote file

c Local variables:
      integer ::
     $  slash      ! position of last slash in a pathname
     $, istat      ! return value from system call
      logical ::
     $  there      ! test for existence of file
      character(len=128) ::
     $  ctmp       ! tmp string for commands and file names

c Externals:
      integer :: lastsl, SYSCALL
c-----------------------------------------------------------------------

c     Check for the existence of the local file.
      inquire( file=trim( lpath ), exist=there )

      if ( there ) then
         write(*,*)'getfile: using local disk copy of ', trim( lpath )
      else

c        Create the local directory if necessary.
         slash = lastsl( lpath )
         if ( slash .gt. 1 ) then
            ctmp = 'mkdir -p ' // lpath(:slash)
            istat = SYSCALL( ctmp )
c           Ignore error return... directory may already exist.
         end if

c        Acquire remote file from NCAR mass store.
c        This code will need modification at non-NCAR sites.
         ctmp  = 'msrcp mss:' // trim( rpath ) // ' ' // trim( lpath )
         write(*,*)'getfile: SYSCALL issuing command - '
         write(*,*) trim( ctmp )
         istat = SYSCALL( ctmp )
         if ( istat .ne. 0 ) then
            write(*,*)'getfile: SYSCALL: can''t read file from mss'
            getfile = -2
            call abort()
            return
         end if
      end if

      getfile = 0

      return
      end

c#######################################################################

      subroutine dispose( lun, lpath, rpath, irt, async, comment, wpwd
     $,                   remove )

c Close history file, dispose to mass store if irt.gt.0, and remove if requested.
c This routine builds a command string which is issued to the system() command.
c The rdsp command is a perl script for managing writing files to the NCAR MSS.
c *** N.B. *** This routine assumes that the rdsp command is in the PATH of the
c              current environment.

      implicit none

c Input arguments.
      integer, intent(in) ::
     $  lun           ! Logical unit number
     $, irt           ! Mass Store retention time
      character(len=*), intent(in) ::
     $  lpath         ! Local filename
     $, rpath         ! Mass Store full pathname
     $, comment       ! Mass Store comment field
     $, wpwd          ! write password
      logical, intent(in) ::
     $  async         ! true => do not wait for successful file transfer
     $, remove        ! true => remove local file after copy to MSS

c Local variables:
      integer :: ier
      character(len=4) :: krt
      character(len=256) :: cmd

c Externals:
      integer SYSCALL  ! execute shell command
c-----------------------------------------------------------------------

c     Close the output file to flush the buffers.
      close( lun )

c     Dispose dataset to Mass Store only if irt.gt.0
      if ( irt .ne. 0 ) then
         write( krt, '(i4)' ) irt
         cmd = 'rdsp -p ' // krt
         if ( remove ) then
            cmd = trim( cmd ) // ' -d'
         end if
         if ( wpwd(:1) .ne. ' ' ) then
            cmd = trim( cmd ) // ' -w ' // trim( wpwd )
         end if
         if ( comment(:1) .ne. ' ' ) then
            cmd = trim( cmd ) // ' -c "' // trim( comment ) // '"'
         end if
         cmd = trim( cmd ) // ' ' // trim( lpath )
         cmd = trim( cmd ) // ' ' // trim( rpath )
         if ( async ) then
            cmd = trim( cmd ) // ' & '
         end if
         write(*,*)'dispose: SYSCALL issuing command - '
         write(*,*) trim( cmd )
         ier = SYSCALL( trim( cmd ) )
      end if

      return
      end
c#######################################################################

      integer function get_ftn_unit()

c Return next available Fortran I/O unit number.  Return -1 if unit not
c available.

      implicit none

c Local variables

      integer :: lun
      logical :: connected
c-----------------------------------------------------------------------

      do lun = 1, 99
         inquire( unit=lun, opened=connected )
         if ( .not. connected ) exit
      end do
      if ( lun == 100 ) then
         get_ftn_unit = -1
      else
         get_ftn_unit = lun
      end if

      end function get_ftn_unit

c#######################################################################

