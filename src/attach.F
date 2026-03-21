c $Id$
      integer function attach( lpath, rpath, lun, binary, cosb )

c Check to see if specified local file exists.  If not, acquire it from
c remote storage.  Check to see if binary file is COS blocked.  Then open
c as unformatted file (if binary) or formatted file (if not binary).

c This routine contains locally written NCAR specific code for:
c 1) Acquiring files that don't reside on the local disk from the NCAR Mass
c    Storage System (msread).
c 2) Reading CRAY COS blocked datasets on NCAR's IBM RS6000 cluster nodes
c    (crayread).
      
      implicit none

c Input arguments:

      character
     $  lpath*(*)  ! pathname for local file
     $, rpath*(*)  ! pathname for remote file

      logical
     $  binary     ! .true. => file is binary (unformatted),
c                  !           else assume ascii (formatted).

c Input/Output argument:

      integer
     $     lun  ! On a Cray, this is the unit number (input) to open file on.
c               ! On a non-Cray machine, if the file is not COS blocked, then
c               ! it is the unit number (input) to open file on.  If the file
c               ! is COS blocked, then it is the return value from crayopen
c               ! (output).

c Output argument:

      logical
     $  cosb    ! true => file is COS blocked

c Return values:

c   0 => successful return
c  -2 => error reading remote file
c  -3 => error determining if file is COS blocked
c  -4 => error opening COS blocked file on a non-Cray machine

c Local variables:

      integer
     $  slash      ! position of last slash in a pathname
     $, istat      ! return value from ishell
     $, i
     $, cid

      logical
     $  there      ! test for existence of file
      
      character
     $  ctmp*128   ! tmp string for commands and file names

c Externals:

      integer lastsl, glc, fsystem, iscosb

      external
     $  lastsl     ! returns position of last '/' in a string.
     $, glc        ! position of last non-blank or non-null character
     $, fsystem    ! execute a shell command
     $, iscosb     ! determine if file is COS blocked.

c Library functions:

      integer
     $  crayopen   ! open a COS blocked file on non-Cray machine
      external crayopen
c-----------------------------------------------------------------------

c     Check for the existence of the local file.

      inquire( file=lpath(:glc(lpath)), exist=there )

      if ( there ) then
         write(*,*)'attach: using local disk copy of '
     $,            lpath(:glc(lpath))

      else

c        If no access to mss, abort if file not on local filesystem
c         write(*,*)'attach: ',lpath(:glc(lpath))
c     $,            ' not on local disk... stopping'
c         stop

c        Create the local directory if necessary.

         slash = lastsl( lpath )
         if ( slash .gt. 1 ) then
            do i = 1, len(ctmp)
               ctmp(i:i) = ' '
            end do
            ctmp = 'mkdir -p ' // lpath(:slash)
            istat = fsystem( ctmp )
c           Ignore error return... directory may already exist.
         end if

c        Acquire remote file from NCAR mass store.
c        This code will need modification at non-NCAR sites.

         do i = 1, len(ctmp)
            ctmp(i:i) = ' '
         end do
         ctmp  = 'msread ' // lpath(:glc(lpath)) // ' '
     $           // rpath(:glc(rpath)) // char(0)
         write(*,*)'attach: fsystem issuing command - '
         write(*,*) ctmp(:glc(ctmp))
         istat = fsystem( ctmp )
         if ( istat .ne. 0 ) then
            write(*,*)'attach - fsystem: can''t read file from mss'
            attach = -2
            return
         end if
      end if

      if ( binary ) then

c        Check for COS blocked file.

         do i = 1, len(ctmp)
            ctmp(i:i) = ' '
         end do
         ctmp = lpath(:glc(lpath)) // char(0)
         istat = iscosb( ctmp )
         if ( istat .lt. 0 ) then
            write(*,*)'attach: error -- iscosb returned ',istat
            attach = -3
            return
         end if

         cosb = .true.
         if ( istat .eq. 0 ) cosb = .false.

c        Open binary file.

#ifdef CRAY
         if ( cosb ) then
            write(*,*)'attach: opening COS blocked file'
         else
            write(*,*)'attach: opening ieee file'
            call asnunit( lun, '-F f77 -N ieee', istat )
         end if
         open( lun, file=lpath(:glc(lpath)),
     $         status='OLD', form='UNFORMATTED' )
#else
         if ( .not. cosb ) then
            write(*,*)'attach: opening ieee file'
            open( lun, file=lpath(:glc(lpath)),
     $            status='OLD', form='UNFORMATTED' )
         else
            write(*,*)'attach: opening COS blocked file'
            cid = o'644'
            lun = crayopen( lpath, 0, cid )
            if ( lun .lt. 0 ) then
               write(*,*)'attach: crayopen failed'
               attach = -4
               return
            end if
         end if
#endif

      else

c        Open ASCII file.

         write(*,*)'attach: opening ascii file'
         open( lun, file=lpath(:glc(lpath)), status='OLD' )

      end if

      attach = 0

      return
      end

c##############################################################################

