c $Id$
      subroutine initrac( lpth, rpth, icdate, icsec, lun, remove
     $,                   xreadNAS, xreadQIC, as, qth )

c Initialize /tracinp/, and tracer arrays.   The tracer fields must come
c directly from a time sample that exists on the input file (though not
c necessarily the first time sample).

      use ccmhist
      use calendar

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <tracinp.com>
c-----------------------------------------------------------------------

c Input arguments:
      character
     $  lpth*(*)   ! local path for tracer initial conditions
     $, rpth*(*)   ! remote path for tracer initial conditions
      integer
     $  icdate     ! date of initial conditions in yymmdd format
     $, icsec      ! seconds relative to date for initial conditions
     $, lun        ! unit number for input
      logical
     $  remove      ! t => remove file after initialization complete
     $, xreadNAS  ! t => chemistry scheme enabled
     $, xreadQIC    ! t => read initial conditions for the advected Q field
c                   !      rather than initializing from meteorology.

c Output arguments:
      real(r8)
     $  as(plond,plev,pcnst,platd) ! advected species
     $, qth(plond,plev,padv,platd) ! Q and THETA

c Local variables:
      character
     $  ctmp*128    ! temp string for system command
      integer
     $  istat   ! return status
     $, idum
     $, i, j
     $, ncraywds
     $, mpsig   ! position in real header of first word of vertical coord list
     $, mpflds  ! position in integer header of field info list
     $, mpcfld  ! position in character header of field info list
      real(r8)
     $  hyam(plev)
     $, hybm(plev)
     $, hyai(plev+1)
     $, hybi(plev+1)
!     real(DBLKIND), dimension(:), allocatable :: dat
      real(r8), dimension(:), allocatable :: dat

c Externals:
      integer glc, lastsl, attach, system
!     real(DBLKIND) diffdat
      integer
     $  crayread  ! read a COS blocked file on a non-Cray machine
     $, crayclose ! close a COS blocked file on a non-Cray machine
c-----------------------------------------------------------------------
c     Set variables in /tracinp/.

c     If local path not specified, take filename from remote path and put
c     it in subdirectory TRA
      if ( lpth(1:1) .eq. ' ' ) then
         lpath = 'TRA/' // rpth(lastsl(rpth)+1:glc(rpth))
      else
         lpath = lpth
      end if
      rpath    = rpth
      luntra   = lun
      readNAS  = xreadNAS
      readQIC  = xreadQIC

c     Get tracer initial conditions file.
      istat = attach( lpath, rpath, luntra, .true., cosbtra )
      if ( istat .ne. 0 ) then
         write(*,*)'initrac: attach returned ', istat
         stop
      end if

c     Read fixed part of first header record to get size and position info.
      call inirdhd( luntra, cosbtra, lenhdi, lenhdc, lenhdr, maxsiz
     $,             nlon, nlonw, nlat, nlev, idum, idum, idum, nfld
     $,             tracdate, tracsec, idum, mpflds, mpcfld, ndens )

c     Check that data resolution agrees with user specified parameters.
      if ( nlon .ne. plon .or. nlat .ne. plat .or. nlev .ne. plev ) then
         write(*,*)'initrac: nlon, nlat, nlev= ', nlon, nlat, nlev
         stop
      end if

c     Make sure first time sample is not later than initial time of run.
      if ( diffdat(  icdate, icsec, tracdate, tracsec ) .gt. 0. ) then
         write(*,*)'initrac: first time sample= ', tracdate, tracsec
         stop
      end if

c     Make sure fixed length arrays are large enough.
      if ( nfld .gt. mxfldsIn ) then
         write(*,*)'initrac: MXFLDSIN must be increased to ', nfld
         stop
      end if

c     Allocate memory for data records.
      allocate( dat(maxsiz) )

   10 continue  ! jump to here after skipping data records

      write (6,*) 'initrac: searching for icdate, icsec ', icdate, icsec
c     Read all 3 header records.

!      write(*,*)'initrac:L115 passing to rdhd: ',luntra, cosbtra, lenhdi
!     $,             lenhdc, lenhdr, nfld
!     $,             tracdate, tracsec, mflds, mcflds
!     $,             nlev
      istat = rdhd( luntra, cosbtra, lenhdi, lenhdc, lenhdr, nfld
     $,             tracdate, tracsec, mflds, mcflds
     $,             nlev, hyam, hybm, hyai, hybi )
      if ( istat .ne. 0 ) then
         write(*,*)'initrac: rdhd returned ', istat
         stop
      end if

      if ( tracdate < 19000000 ) tracdate = tracdate + 19000000

      write (6,*) 'initrac: found tracdate, tracsec ', tracdate, tracsec

c     Read initial conditions if correct time found.  If not, skip data
c     records and read next time sample header.

      if ( diffdat(  icdate, icsec, tracdate, tracsec ) .eq. 0. ) then

         call gettrac( as, qth )

      else  ! skip the data records in this time sample
         write (6,*) ' initrac: skipping data records '
         do j = 1, nlat
#ifdef CRAY
            read( luntra ) dat
#else
            if ( cosbtra ) then
#if (REALSIZE == 4)
               write (6,*) ' crayread disabled for 32 bit precision '
               stop
#else
               ncraywds = crayread( luntra, dat, maxsiz, 0 )
               if ( ncraywds .le. 0 ) then
                  write(*,*)'initrac: crayread returned ', ncraywds
                  stop
               end if
#endif

            else
               read( luntra ) dat
            end if
#endif  
         end do
         go to 10
      end if

      deallocate( dat )

c     Close file.

#ifdef CRAY
      close( luntra )
#else
      if ( cosbtra ) then
         istat = crayclose( luntra )
      else
         close( luntra )
      end if
#endif

c     Remove if requested.

      if ( remove ) then
         do i = 1, len(ctmp)
            ctmp(i:i) = ' '
         end do
         ctmp = 'rm -f ' // lpath(:glc(lpath))
         write(*,*)'initrac: system issuing command - '
         write(*,*) ctmp(:glc(ctmp))
         istat = system( ctmp )
      endif

      return
      end

c##############################################################################

      subroutine gettrac( as, qth )

c Read constituent data from history tape.
c *** N.B. *** This code will break if it is built for single precision reals 
c              and the data is cray data because the variable drecin(maxsiz)
c              won't be large enough to hold the full data record which will
c              be 8*maxsize bytes.

      use ccmhist

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <tracinp.com>
c-----------------------------------------------------------------------
#include <tracnam.h>
c-----------------------------------------------------------------------

c Output arguments:
      real(r8) as(plond,plev,pcnst,platd) ! advected species
      real(r8) qth(plond,plev,padv,platd) ! Q, THETA

c Local variables:
      integer
     $  j     ! do loop index over latitudes
     $, jcen  ! index into extended arrays corresponding to "lat"
     $, lat   ! latitude index of current data record
     $, m     ! constituent index

      real(r8) fldsin(maxsiz)  ! storage for a data record from the hist tape
c-----------------------------------------------------------------------

      do j = 1, plat

c        Read the data record for the current latitude.
         call rddat( luntra, cosbtra, fldsin, maxsiz, ndens )

         lat = int( fldsin(1) + .01 )     ! type 43 history tape
         jcen   = lat + nxpt + jintmx

c        Extract appropriate slice and transfer to model arrays
         do m = 1, pcnst
            call mkslic( tracnam(m), as(i1,1,m,jcen), nlev, mflds
     $,                  mcflds, nfld, fldsin, maxsiz, nlon, nlonw
     $,                  plond, nlev )
         end do
         if ( readQIC ) then
            call mkslic( 'Q', qth(i1,1,1,jcen), nlev, mflds
     $,                  mcflds, nfld, fldsin, maxsiz, nlon, nlonw
     $,                  plond, nlev )
         end if

         if ( readNAS ) then

c           Set nas slice pointer to current latitude slice in incore array.
c            call nasptr( lat, pnasls )
c
c            do m = 1, pnats
c               call mkslic( natsnam(m), nasls(1,1,m), nlev, mflds
c     $,                     mcflds, nfld, fldsin, maxsiz, nlon, nlonw
c     $,                     plond, nlev )
c            end do

         end if

      end do

      return
      end
