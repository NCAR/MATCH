module sat_obs

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Provide satellite observations of aerosol optical depth to the
   ! Khattatov/Lamarque assimilation code.
   !
   ! Author: B. Eaton
   !-----------------------------------------------------------------------

   use pmgrid
   use precision
   use calendar

   implicit none


   save
   private
   public :: &
      sat_obs_ini,   &! initialization
      satDataExists, &! returns true when there is satellite data to assimilate
      getNSSAOD       ! reads satellite AOD data from netCDF files,  subtracts
                      ! sea salt AOD from total.

   ! Local

   integer, parameter ::&
      mxaodfl = 32768,              &! Maximum number of aerosol OD files
      mxfllen = 300                  ! Maximum length of OD file names (in char)

   integer ::&
      naodfl = 0,                   &! Number of aerosol OD files to examine
      aodate(mxaodfl),              &! Dates for AOD data (format=YYYYMMDD)
      aosec(mxaodfl),               &! Seconds rel. to aodate for AOD data
      iaosrt(mxaodfl),              &! Indices to AOD files in ascending chronological order
      iaod                           ! Index into iaosrt for 1st AOD data after current
                                     !    model time
   character(len=mxfllen) ::&
      aofl(mxaodfl)                  ! File names of aerosol files

   real(DBLKIND) ::&
      time(mxaodfl)                  ! aodate + aosec/86400.  (used for sorting)

!##############################################################################
contains
!##############################################################################

   subroutine sat_obs_ini( read_sat_obs, satdir, ncdate, ncsec )

      ! Initialize sat_obs module, and identify AOD data 
      ! closest to but later than the current time slice
      ! Setting read_sat_obs=.false. effectively turns off the assimilation
      ! since naodfl=0 causes satDataExists to return .false.

     use precision
     implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      logical, intent(in) ::&
         read_sat_obs     ! t ==> read the satellite obs files in satdir

      integer, intent(in) ::&
         ncdate,         &! current date
         ncsec            ! current time in seconds relative to ncdate

      character(len=*), intent(in) ::&
         satdir           ! name of assimilation directory

      ! Local variables

      character(len=13) :: routine = 'sat_obs_ini: '

      integer ::&
         i,        &! longitude index
         j,        &! latitude index
         k,        &! level index
         l,        &! index into AOD sorted indices
         laod       ! AOD time slice

      integer ::&
         ierror,   &! error return from directory-handling routines
         ilen,     &! length of file name
         EEND,     &! return code signifying end of directory
         jdirent,  &! structure handle for directory entry
         idirid     ! unique directory identifier

      character(len=mxfllen) ::&
         asfile     ! name of file in assimilation directory

      integer ::&
         ncid,     &! NetCDF handle for assimilation file
         ncstatus, &! Return code from NetCDF call
      
         syear,    &! starting year of satellite data
         sday,     &! starting day of satellite data
         shour,    &! starting hour of satellite data
         smin,     &! starting minute of satellite data
         ehour,    &! ending hour of satellite data
         emin,     &! ending minute of satellite data
      
         syearid,  &! NetCDF handle for starting year of satellite data
         sdayid,   &! NetCDF handle for starting day of satellite data
         shourid,  &! NetCDF handle for starting hour of satellite data
         sminid,   &! NetCDF handle for starting minute of satellite data
         ehourid,  &! NetCDF handle for ending hour of satellite data
         eminid,   &! NetCDF handle for ending minute of satellite data
      
         sdate,    &! Starting date corresponding to (syear,sday)
         ssec,     &! Seconds rel. to sdate from (shour,smin)
         esec,     &! Seconds rel. to sdate from (ehour,emin)
         mdate,    &! Date corresponding to sdate + midpoint(ssec,esec)
         msec       ! Secconds rel. to mdate

      real(r8) ::&
         diff

      character(len=256) :: cmd
      integer :: lun

      ! externals
      integer :: get_ftn_unit
      integer :: fsystem       ! execute shell command
!      integer :: newdate
!      real(DBLKIND) :: diffdat
      !-----------------------------------------------------------------------------

      ! Determine list of file names for assimilation 

      ! Issue system command to run the perl script list_nc_files.
      ! cmd = '/home/bluesky/fillmore/match/run/list_nc_files ' // trim( satdir ) // ' sat_obs_files'
      cmd = '/homedir/dfillmor/MATCH/build/list_nc_files ' // trim( satdir ) // ' sat_obs_files'
      ! cmd = 'list_nc_files ' // trim( satdir ) // ' sat_obs_files'
      write(*,*) routine//'fsystem issuing command - '
      write(*,*) trim( cmd )
      ierror = fsystem( trim( cmd ) )
      write(*,*) 'done with system command list_nc_files'

      ! Get next availble Fortran I/O unit number.
      lun = get_ftn_unit()
      if ( lun < 0 ) then
         write(*,*) routine//'get_ftn_unit returned error'
         call endrun()
      end if

      open( unit=lun, file='sat_obs_files', iostat=ierror )
      if (ierror /= 0) then
         write(*,*) routine//'error opening file sat_obs_files'
         call endrun
      endif

      read( unit=lun, fmt=*, iostat=ierror ) naodfl
      if (ierror /= 0) then
         write(*,*) routine//'error reading naodfl from file sat_obs_files'
         call endrun
      endif
      if ( naodfl > mxaodfl ) then
         write(*,*) routine//'must increase size of mxaodfl to at least ', naodfl
         call endrun
      endif

      do l = 1, naodfl
         read( unit=lun, fmt='(a)', iostat=ierror ) aofl(l)
         if (ierror /= 0) then
            write(*,*) routine//'error reading file sat_obs_files'
            call endrun
         endif
      end do
         
      close( lun )

      ! Determine dates and times of assimilation files
      write(*,*) routine//'reading in sat_obs files to determine dates and times'
      write(*,*) 'number of files',naodfl
      if (naodfl .gt. 0) then 
         do l = 1, naodfl
            write(*,*) routine//'reading in file ',l,TRIM(aofl(l))
            ! Load satellite time information

            ncstatus = nf_open(aofl(l), NF_NOWRITE, ncid)
            call handle_ncerr( ncstatus, routine//'error opening '//TRIM(aofl(l)) )

            ncstatus = nf_inq_varid(ncid, 'syear', syearid)
            call handle_ncerr( ncstatus, routine//'error inquiring about syear' )
            ncstatus = nf_get_var_int(ncid, syearid, syear)
            call handle_ncerr( ncstatus, routine//'error reading syear' )

            ncstatus = nf_inq_varid(ncid, 'sday', sdayid)
            call handle_ncerr( ncstatus, routine//'error inquiring about sday' )
            ncstatus = nf_get_var_int(ncid, sdayid, sday)
            call handle_ncerr( ncstatus, routine//'error reading sday' )

            ncstatus = nf_inq_varid(ncid, 'shour', shourid)
            call handle_ncerr( ncstatus, routine//'error inquiring about shour' )
            ncstatus = nf_get_var_int(ncid, shourid, shour)
            call handle_ncerr( ncstatus, routine//'error reading shour' )

            ncstatus = nf_inq_varid(ncid, 'smin', sminid)
            call handle_ncerr( ncstatus, routine//'error inquiring about smin' )
            ncstatus = nf_get_var_int(ncid, sminid, smin)
            call handle_ncerr( ncstatus, routine//'error reading smin' )

            ncstatus = nf_inq_varid(ncid, 'ehour', ehourid)
            call handle_ncerr( ncstatus, routine//'error inquiring about ehour' )
            ncstatus = nf_get_var_int(ncid, ehourid, ehour)
            call handle_ncerr( ncstatus, routine//'error reading shour' )

            ncstatus = nf_inq_varid(ncid, 'emin', eminid)
            call handle_ncerr( ncstatus, routine//'error inquiring about emin' )
            ncstatus = nf_get_var_int(ncid, eminid, emin)
            call handle_ncerr( ncstatus, routine//'error reading smin' )

            ncstatus = nf_close(ncid)
            call handle_ncerr( ncstatus, routine//'error closing '//TRIM(aofl(l)) )

            ! Assign midpoint between starting and ending times to assim time.

!            write(6,* )l,' sat_obs file ',trim(aofl(l))            
            if (syear < 50 ) then  
               ! Assume 0-50 = 2000-2050, 51-99 = 1951-1999
               syear = syear + 2000;
               sdate = newdate(syear*10000 + 101, sday - 1)   !files use 2 digit year
!               write(6,*)' year < 50, modified to ',syear, 'date ',sdate
            else if (syear < 100 ) then 
               !came in with 2 digit year, must be 51-99 = 1951-1999
               syear = syear + 1900;
               sdate = newdate(syear*10000 + 101, sday - 1)   !files use 2 digit year
!               write(6,*)' year < 100, modified to ',syear, 'date ',sdate
            else
               sdate = newdate(syear*10000 + 101, sday - 1)   !files use 4 digit year
!               write(6,*)' year okay. Not modified. ',syear, 'date ',sdate
            endif

            ssec = shour * 3600 + smin * 60
            esec = ehour * 3600 + emin * 60

            if (esec .lt. ssec) then
               esec = esec + 86400
            endif

            call addsec2dat((esec - ssec) / 2, sdate, ssec, mdate, msec)
            aodate(l) = mdate
            aosec(l) = msec
            time(l) = real(mdate,r8) + real(msec,r8)/real(86400,r8)

!            write(6,999) mdate,msec,time(l)
999         format('date sec  time',i10,i10,f20.8)

         end do
      endif


      ! Determine time slice closest to current time slice

      iaod = -1

      ! First sort the files by date and time

      if (naodfl .eq. 0) then 
         return
      else if (naodfl .eq. 1) then 
         iaosrt(1) = 1
      else if (naodfl .gt. 1) then
         do l = 1, naodfl
            iaosrt(l) = l
         end do
         call sort_aod_times()
!         write(6,*)'DRB sort: ',time(iaosrt(1))
      endif
      
      l = 1
      do while (iaod .eq. -1 .and. l .le. naodfl)
         laod = iaosrt(l)
         diff = diffdat( ncdate, ncsec, aodate(laod), aosec(laod) )
!         write(6,*)'drb: ncdate, ncsec, aodate(laod), aosec(laod),diff'
!         write(6,*)'drb: ',ncdate, ncsec, aodate(laod), aosec(laod),diff
         if (diff .gt. 0) then
            iaod = l
         end if
         l = l + 1
      end do
     
   end subroutine sat_obs_ini

!#####################################################################################

   subroutine sort_aod_times()

      ! sort the times of the aod files using selection sort

      integer :: i, j, jj, min, t
      !--------------------------------------------------------------------------

      do i = 1, naodfl-1
         min = i
         do j = i+1, naodfl
            if ( time(iaosrt(j)) < time(iaosrt(min)) ) min = j
         end do
         t = iaosrt(min)
         iaosrt(min) = iaosrt(i)
         iaosrt(i) = t
      end do

   end subroutine sort_aod_times

!#####################################################################################

   logical function satDataExists( ncdate, ncsec )

      ! Return .true. if BOTH 1. and 2. are true.
      !      (1. There are any files (left) to assimilate
      !      (2. If there are, is the current MATCH timestep .ge. closest
      !          satellite overpass (to within idt).

      use precision

      implicit none

      integer, intent(in) ::&
         ncdate, ncsec   ! current model date and sec

      ! Local
!      real(DBLKIND) :: diffdat
      integer laod
!-----------------------------------------------------------------------

      satDataExists = .false.

      if (naodfl .eq.  0     .or.   &
          iaod   .lt.  1     .or.   &
          iaod   .gt. naodfl) then
         write(*,*)'AOD_ASSIM_3D: no suitable sat file. n files ',naodfl,'index of next file',iaod
         return
      endif

      laod = iaosrt(iaod)
      if (diffdat( ncdate,ncsec,aodate(laod),aosec(laod) ) .gt. 0) then
         write(6,*)'AOD_ASSIM_3D: no file for this time step',ncdate,ncsec
         write(6,*)'            : nearest file at', aodate(laod),aosec(laod),trim( aofl(laod) )

!         write(*,*)'drb:sat_obs: returning, diffdat',diffdat( ncdate,ncsec,aodate(laod),aosec(laod) )
!         write(*,*)'drb:sat_obs ncdate,ncsec,aodate(laod),aosec(laod)',ncdate,ncsec,aodate(laod),aosec(laod)

         return
      end if

      write(*,*)"AOD_ASSIM_3D: MATCH date/time = ", ncdate, " / ", ncsec
      write(*,*)"              AOD date/time = ", aodate(laod), " / ", aosec(laod)
      write(*,*)"              AOD file      = ", trim( aofl(laod) )

      satDataExists = .true.

   end function satDataExists

!#######################################################################

   subroutine getNSSAOD( obsdate, obssec, data_lat, data_lon, nss_aod )

      ! Read satellite AOD data from netCDF files.  Subtract the sea salt optical
      ! depth.

     use precision
     use netcdf

     implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------
#include <ssltcom.h>
!-----------------------------------------------------------------------

     integer, intent(out) ::&
          obsdate, obssec     ! date and sec of obs
     real(r8), intent(out) ::&
         data_lat(plat),    &! Latitudes of gridded sat. obs
         data_lon(plon),    &
         nss_aod(plon,plat)

      ! Local

      character(len=11) :: routine = 'getNSSAOD: '

      integer ::&
         i, j,        &
         laod,        &! Index to current AOD file
         ncstatus,    &! NetCDF return code
         ncid,        &! NetCDF file handle
         lat_varid,   &! Latitude variable ID
         lon_varid,   &! Longitude variable ID
         aod_varid     ! AOD variable ID

      real(r8) ::&
!         tsslt2d(plond,plat), &! sea salt optical depth
         data_aod(plon,plat), &!
         sslt_od
      !-----------------------------------------------------------------------

      write (*,*) 'TEST - getNSSAOD enter'

      ! Read data

      laod = iaosrt(iaod)
      iaod = iaod + 1

      obsdate = aodate(laod)
      obssec = aosec(laod)

      ncstatus = nf_open(aofl(laod), NF_NOWRITE, ncid)
      call handle_ncerr( ncstatus, routine//'error opening '//TRIM(aofl(laod)) )
      
      ncstatus = nf_inq_varid(ncid, 'lat_ctr', lat_varid)
      call handle_ncerr( ncstatus, routine//'error inquiring about lat_ctr' )

      ncstatus = nf90_get_var(ncid, lat_varid, data_lat)

      call handle_ncerr( ncstatus, routine//'error reading lat_ctr' )

      ncstatus = nf_inq_varid(ncid, 'lon_ctr', lon_varid)
      call handle_ncerr( ncstatus, routine//'error inquiring about lon_ctr' )

      ncstatus = nf90_get_var(ncid, lon_varid, data_lon)

      call handle_ncerr( ncstatus, routine//'error reading lon_ctr' )

      ncstatus = nf_inq_varid(ncid, 'tau_avg', aod_varid)
      call handle_ncerr( ncstatus, routine//'error inquiring about tau_avg' )

      ncstatus = nf90_get_var(ncid, aod_varid, data_aod)

      call handle_ncerr( ncstatus, routine//'error reading tau_avg' )

      ncstatus = nf_close(ncid)
      call handle_ncerr( ncstatus, routine//'error closing '//TRIM(aofl(laod)) )

! tsslt2d stored in ssltcom.h, calculated by optDepth in physlic

      do j = 1, plat
         do i = 1, plon
            if ( data_aod(i,j) .ge. 0.0 ) then
               sslt_od = max( 0.0_r8, tsslt2d(i,j) )
               nss_aod(i,j) = max( 0.0_r8, data_aod(i,j) - sslt_od )
            else
               nss_aod(i,j) = data_aod(i,j)
            end if
         end do
      end do

      write (*,*) 'TEST - getNSSAOD exit'

   end subroutine getNSSAOD

!#######################################################################

end module sat_obs
