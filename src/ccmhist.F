      module ccmhist

      use precision

      contains
c $Id$
#include <params.h>

      subroutine inirdhd( lun, cosb, lenhdi, lenhdc, lenhdr, maxsiz
     $,                   nlon, nlonw, nlat, nlev, ntrm, ntrn, ntrk
     $,                   nfld, ncdate, ncsec, mpsig, mpflds, mpcfld
     $,                   ndens )

c Read part of the integer header record of a CCM2 history tape.  Extract
c information about sizes and pointers to data in the header records, then
c rewind the file.

      implicit none

c Input arguments:

      integer
     $  lun   ! Fortran logical unit number.

      logical
     $  cosb  ! .true. => file is COS blocked

c Output arguments:

      integer
     $  lenhdi  ! length of integer header record (number of integers)
     $, lenhdc  ! length of character header record (number of character*8)
     $, lenhdr  ! length of real header record (number of reals)
     $, maxsiz  ! length of integer header record (number of reals)
     $, nlon    ! number of longitudes
     $, nlonw   ! number of longitudes written on history tape
     $, nlat    ! number of latitudes
     $, nlev    ! number of model layers
     $, ntrm    ! spectral truncation parameter M
     $, ntrn    ! spectral truncation parameter N
     $, ntrk    ! spectral truncation parameter K
     $, nfld    ! number of fields in data records
     $, ncdate  ! current date (in yyyymmdd format)
     $, ncsec   ! seconds relative to ncdate
     $, mpsig   ! position in real header of first word of vertical coord list
     $, mpflds  ! position in integer header of field info list
     $, mpcfld  ! position in character header of field info list
     $, ndens   ! packing density for data records

c Local variables:

      integer
     $  nint      ! number of integers to read from record
      parameter( nint=40 )

      integer
     $  ihd(nint) ! storage for nint integers from header record
     $, ncraywds  ! number of cray words read if successful, -1 otherwise
     $, iret      ! return value
     $, i         ! index

c Functions:

      integer
     $  crayread  ! read COS blocked records on non-Cray machine
     $, crayrew   ! rewind a COS blocked file on a non-Cray machine
c-----------------------------------------------------------------------

#ifdef CRAY
      read( lun ) (ihd(i),i=1,nint)
      rewind( lun )
#else
      if ( cosb ) then
         ncraywds = crayread( lun, ihd, nint, 3 )
         if ( ncraywds .le. 0 ) then
            write(*,*)'inirdhd: crayread returned ', ncraywds
            stop
         end if
         iret = crayrew( lun )
         if ( iret .ne. 0 ) then
            write(*,*)'inirdhd: crayrew returned ', iret
            stop
         end if
      else
         read( lun ) (ihd(i),i=1,nint)
         rewind( lun )
      end if
#endif

      lenhdi = ihd(1)
      lenhdc = ihd(31)
      lenhdr = ihd(32)
      maxsiz = ihd(6)
      nlon   = ihd(9)
      nlonw  = ihd(10)
      nlat   = ihd(11)
      nlev   = ihd(12)
      ntrm   = ihd(13)
      ntrn   = ihd(14)
      ntrk   = ihd(15)
      nfld   = ihd(16)
      ncdate = ihd(26)
      ncsec  = ihd(27)
      mpsig  = ihd(33)
      mpflds = ihd(36)
      mpcfld = ihd(37)
      ndens  = ihd(40)

#ifdef DEBUG
      write(*,*)'inirdhd: lenhdi,lenhdc,lenhdr,maxsiz,nlon,nlonw,nlat'
      write(*,*)'   nlev,ntrm,ntrn,ntrk,nfld,ncdate,ncsec,mpsig,mpflds'
      write(*,*)'   mpcfld,ndens:'
      write(*,10) lenhdi,lenhdc,lenhdr,maxsiz,nlon,nlonw,nlat
     $,           nlev,ntrm,ntrn,ntrk,nfld,ncdate,ncsec,mpsig,mpflds
     $,           mpcfld,ndens
   10 format(2(10i8/))
#endif

      return
      end subroutine  inirdhd

c##############################################################################

      integer function rdhd( lun, cosb, lenhdi, lenhdc, lenhdr, nfld
     $,                      ncdate, ncsec, mflds, mcflds
     $,                      nlev, hyam, hybm, hyai, hybi )

c Read header records of a CCM2 history tape.  Extract the current date.
c Return values:
c  0 => success
c -1 => eof encountered reading integer header
c  1 => eof encountered reading 2nd or 3rd header, or non recoverable i/o error

      implicit none

c Input arguments:
      integer
     $  lun     ! Fortran logical unit number.
     $, lenhdi  ! length of integer header record (number of integers)
     $, lenhdc  ! length of character header record (number of character*8)
     $, lenhdr  ! length of real header record (number of reals)
     $, nfld    ! number of fields
     $, nlev    ! number of layers
      logical
     $  cosb    ! .true. => file is COS blocked

c Output arguments:
      integer
     $  ncdate            ! current date (in yyyymmdd format)
     $, ncsec             ! seconds relative to ncdate
     $, mflds(3,nfld)     ! field info from integer header
      character(len=8) ::
     $  mcflds(2,nfld)    ! field info from character header
      real(r8)
     $  hyam(nlev)        ! hybrid A coeffs at midpoints
     $, hybm(nlev)        ! hybrid B coeffs at midpoints
     $, hyai(nlev+1)      ! hybrid A coeffs at interfaces
     $, hybi(nlev+1)      ! hybrid B coeffs at interfaces

c Local variables:
      integer
     $  ncraywds      ! number of cray words read if successful, -1 otherwise
     $, ios           ! i/o status
     $, mpflds        ! start index for field info in integer header
     $, mpcfld        ! start index for field info in character header
     $, mpsig         ! start index for hybrid coords in real header
     $, a1            ! start index for hybrid A coeffs in real header
     $, b1            ! start index for hybrid B coeffs in real header
     $, hdi(lenhdi)   ! integer header
      character(len=8) ::
     $  hdc(lenhdc)   ! character header
      real(r8)
     $  hdr(lenhdr)   ! real header

c Externals:
      integer
     $  crayread  ! read COS blocked records on non-Cray machine
c-----------------------------------------------------------------------

c     Read integer header record.

#ifdef CRAY
      read( lun, iostat=ios ) hdi
      if ( ios .lt. 0 ) then
         rdhd = -1
         return
      else if ( ios .gt. 0 ) then
         rdhd = 1
         return
      end if
#else
      if ( cosb ) then
         ncraywds = crayread( lun, hdi, lenhdi, 3 )
         if ( ncraywds .le. 0 ) then
            rdhd = -1
            return
         end if
      else
         read( lun, iostat=ios ) hdi
         if ( ios .lt. 0 ) then
            rdhd = -1
            return
         else if ( ios .gt. 0 ) then
            rdhd = ios
            return
         end if
      end if
#endif

      ncdate = hdi(26)
      ncsec  = hdi(27)
      mpsig  = hdi(33)
      mpflds = hdi(36)
      mpcfld = hdi(37)
      mflds = reshape( hdi(mpflds:lenhdi), (/3,nfld/) )

c     Read character header record.

#ifdef CRAY
      read( lun, iostat=ios ) hdc
      if ( ios .ne. 0 ) then
         rdhd = 1
         return
      end if
#else
      if ( cosb ) then
         ncraywds = crayread( lun, hdc, lenhdc, 0 )
         if ( ncraywds .le. 0 ) then
            rdhd = 1
            return
         end if
      else
         read( lun, iostat=ios ) hdc
         if ( ios .ne. 0 ) then
            rdhd = 1
            return
         end if
      end if
#endif

      mcflds = reshape( hdc(mpcfld:lenhdc), (/2,nfld/) )

c     Read real header record.

#ifdef CRAY
      read( lun, iostat=ios ) hdr
      if ( ios .ne. 0 ) then
         rdhd = 1
         return
      end if
#else
      if ( cosb ) then
         ncraywds = crayread( lun, hdr, lenhdr, 2 )
         if ( ncraywds .le. 0 ) then
            rdhd = 1
            return
         end if
      else
         read( lun, iostat=ios ) hdr
         if ( ios .ne. 0 ) then
            rdhd = 1
            return
         end if
      end if
#endif

c     Unscramble hybrid coefficients.
      a1 = mpsig + 2*nlev+1
      b1 = a1  + 2*nlev+1
      hyai = hdr(a1  :a1+2*nlev:2)
      hyam = hdr(a1+1:a1+2*nlev:2)
      hybi = hdr(b1  :b1+2*nlev:2)
      hybm = hdr(b1+1:b1+2*nlev:2)

      rdhd = 0

      return
      end function rdhd

c##############################################################################

      subroutine rddat( kunit, cosblk, pa, klen, ndens )

c Read data from unit kunit to array pa.  If the file is COS blocked,
c and the data is packed, then it is read without conversion to ieee
c and the first two words of the record are converted to ieee separately.
c The packed data should be passed to the unpacker just as it was written.
c If the COS blocked record is unpacked, then it is converted to ieee.

c Input arguments:

      integer
     $  kunit  ! Fortran unit or Cray file specifier for data file
     $, klen   ! number of reals to read
     $, ndens  ! packing flag, 1 => unpacked, 2 or 4 => packed.

      logical
     $  cosblk    ! true => file is COS blocked
                  ! false => file is ieee

c Output arguments:

      real(r8)
     $  pa(klen)  ! contains data read from kunit

c Local variables:

      integer ncraywds
      real(r8) data(2)

c Functions:

      integer crayread
c-----------------------------------------------------------------------

#ifdef CRAY
      read( kunit ) pa
#else
      if ( cosblk ) then
         if ( ndens .gt. 1 ) then
            ncraywds = crayread( kunit, pa, klen, 0 )
            if ( ncraywds .le. 0 ) then
               write(*,*)'rddat: crayread returned ', ncraywds
               stop
            end if
            call ctodpf( pa, data, 2 )
            pa(1) = data(1)
            pa(2) = data(2)
         else
            ncraywds = crayread( kunit, pa, klen, 2 )
            if ( ncraywds .le. 0 ) then
               write(*,*)'rddat: crayread returned ', ncraywds
               stop
            end if
         end if
      else
         read( kunit ) pa
      end if
#endif

      return
      end subroutine rddat

c##############################################################################

      subroutine mkslic (field, slice, nlevf,
     $     mflds, mcflds, nfldh, hrec,lhrec,
     $     nlon, nlonh, nlond, nlev)
c+
c Mkslic extracts the data for a specified field from an input history
c tape record and returns it as a vertical slice.
c
c The input data record may be packed.
c
c Byron A Boville, October 16, 1990, adated from earlier code by Joe Wag
c-
      implicit none

c+
c Input arguments
c-
      integer
     $     lhrec,               ! length of input record
     $     nfldh,               ! number of fields in input record
     $     nlon,                ! number of longitudes (distinct)
     $     nlond,               ! number of longitudes (dimension)
     $     nlonh,               ! number of longitudes (input record)
     $     nlev,                ! number of levels (multilevel fields)
     $     nlevf,               ! number of levels in current field
     $     mflds(3,nfldh)       ! integer field list from input header

      character
     $  mcflds(2,nfldh)*8       ! character field list (header)
     $, field*(*)               ! desired field

      real(r8) slice(nlond,nlevf),  ! unpacked data (vertical slice)
     $     hrec(lhrec)          ! packed data slice read
c+
c Local storage
c-
      integer
     $     if, i, k, l,         ! counters
     $     ksiz,                ! packed size of data per level
     $     ndens,               ! packing of current field
     $     nlevfi,              ! number of levels in input field
     $     start                ! starting point of current field
      real(r8)
     $  amn, amx                ! dummy args for unpkaf
c+
c Action
c
c Loop through fields on input record until the required field is found.
c The input pointer cannot be used because we allowing for packed data.
c-
      start = 3
      do 100 if = 1, nfldh
c+
c Determine the packing density, and number of words per level for curre
c field.
c-
         ndens = mflds(3,if)
         if (ndens .gt. 1) then
            ksiz = 3 + (nlonh - 1) / ndens
         else
            ksiz = nlonh
         endif
c+
c Determine the number of levels for the current field
c-
         if (mod(mflds(1,if),10).lt.1) then
            nlevfi = 1
         else
            nlevfi = nlev
         endif
c+
c Is this the required field?
c-
         if (field .eq. mcflds(1,if)) then
c+
c This is the required field, check the number of levels.
c-
            if (nlevf .ne. nlevfi) then
               print *, 'Level mismatch for field: ', field
               if (nlevf .eq. 1) then
                  print *, 'The array is single level but the'
     $                 //' data is multilevel'
               else
                  print *, 'The array is multilevel but the'
     $                 //' data is single level'
               end if
               stop 'badlevs'
            end if
c+
c Unpack and load the data.
c-
            l = start
            do 20 k=1,nlevfi
               if (ndens.gt.1) then
                  call unpkaf( slice(1,k), nlon, hrec(l), ndens, amn
     $,                        amx )
               else
                  do 50 i = 1, nlon
                     slice(i,k) = hrec(l+i-1)
   50             continue
               endif
               l = l + ksiz
   20       continue
            return
         else
c+
c Not the required field, skip it.
c-
            start = start + ksiz * nlevfi

         end if

  100 continue
c+
c The required field was not found in the input record.
c-
      print *, 'The required input field ', field,
     $     ' was not on the initial data file'
      print *, 'The available fields are:'
      print *, (mcflds(1,if), if = 1, nfldh)
      stop 'nofield'

      end subroutine mkslic

      end module ccmhist
