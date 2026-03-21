      module dsttvbds

      use precision

      contains

c     $Header$

c     Purpose: dstvbds.F contains subroutines for initializing, reading,
c     interpolating, time-varying boundary data sets into the dust model. 
c     These datasets are currently assumed be annual cycle datasets 
c     containing twelve monthly averages.
c     These routines depend heavily on netCDF I/O

      subroutine dst_tvbds_ini(fl_in,time_doy) ! I

      use netcdf

c     Purpose: Initialize annual cycle of land surface characteristics
c     NB: This routine opens fl_in and leaves it open for the duration of the run
c     Remember to use dst_tvbds_close() to close it
c     dst_tvbds_ini() is called by CCM:control/initext(), MATCH:src/inirun()
      implicit none
#include <netcdf.inc>
#include <dsttvbds.h> /* Time-varying boundary data sets */
c     Parameters
c     Input
      character fl_in*(*)       ! netCDF input file
      real(r8) time_doy             ! Day of year
c     Output (common block dstbnd)
c     Input/Output
c     Local
      integer lat_dim_id        ! [enm] Dimension ID for lat
      integer lat_nbr_in        ! Dimension size
      integer lon_dim_id        ! [enm] Dimension ID for lon
      integer lon_nbr_in        ! Dimension size
      integer rcd               ! [rcd] Return success code
      integer time_dim_id       ! [enm] Dimension ID for time
      integer time_id           ! [enm] Variable ID
      integer lai_dst_id        ! [enm] Variable ID
      integer time_nbr_in       ! Dimension size
      
c     Externals
!     integer idx_glb_get       ! idx_glb_get() is defined in dstnc.F
      external nc_err_exit      ! nc_err_exit() is defined in dstnc.F

c     Read in land surface types from high resolution grid
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf90_open(fl_in,NCNOWRIT,nc_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini()')
      write(6,'(a29,1x,a)') 'dst: Opened netCDF input file',fl_in

c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'lat',lat_dim_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini()')
      rcd=rcd+nf90_inq_dimid(nc_id,'lon',lon_dim_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini()')
      rcd=rcd+nf90_inq_dimid(nc_id,'time',time_dim_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini()')
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,lat_dim_id,len=lat_nbr_in)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini()')
      if (lat_nbr /= lat_nbr_in) stop 'lat_nbr /= lat_nbr_in'
      rcd=rcd+nf90_inquire_dimension(nc_id,lon_dim_id,len=lon_nbr_in)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini()')
      if (lon_nbr /= lon_nbr_in) stop 'lon_nbr /= lon_nbr_in'
      rcd=rcd+nf90_inquire_dimension(nc_id,time_dim_id,len=time_nbr_in)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini()')
      write (6,*) ' time_nbr ', time_nbr
      write (6,*) ' time_nbr_in ', time_nbr_in
      if (time_nbr /= time_nbr_in) stop 'time_nbr /= time_nbr_in'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'time',time_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini():varid time')
      rcd=rcd+nf90_inq_varid(nc_id,'lai_lsm',lai_dst_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini():varid lai_lsm')

c     Get time coordinate
      rcd=rcd+nf90_get_var(nc_id,time_id,time)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini():get_var time')
c     Get indices of bounding time slices
      idx_glb_dsk=idx_glb_get(time,time_nbr,time_doy)
      idx_lub_dsk=mod(idx_glb_dsk,time_nbr)+1

c     Initialize hyperslab indices for disk file
      srt_lon_lat_time(1)=1
      srt_lon_lat_time(2)=1
      srt_lon_lat_time(3)=1
      cnt_lon_lat_time(1)=lon_nbr
      cnt_lon_lat_time(2)=lat_nbr
      cnt_lon_lat_time(3)=1

c     Initialize hyperslab indices for variables in RAM
      idx_glb_ram=1
      idx_lub_ram=2

c     Get glb record
      srt_lon_lat_time(3)=idx_glb_dsk
      rcd=rcd+nf90_get_var(nc_id,lai_dst_id,lai_dst_bnd(:,:,idx_glb_ram),srt_lon_lat_time,cnt_lon_lat_time)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini():get_vara lai_dst')

c     Get lub record
      srt_lon_lat_time(3)=idx_lub_dsk
      rcd=rcd+nf90_get_var(nc_id,lai_dst_id,lai_dst_bnd(:,:,idx_lub_ram),srt_lon_lat_time,cnt_lon_lat_time)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ini():get_vara lai_dst')

      write (6,'(3(a,f9.4))') 'dst: dst_tvbds_ini(): Day of year = ',time_doy,
     $     ': read data for days',time(idx_glb_dsk),' and ',time(idx_lub_dsk)

      return
      end subroutine dst_tvbds_ini                       ! end dst_tvbds_ini()

      subroutine dst_tvbds_ntp(time_doy) ! I
c     Purpose: Interpolate annual cycle data to current time value
c     Input new data as necessary
c     dst_tvbds_ntp() is called by CCM:dynamics/advnce(), MATCH:src/main()
      use netcdf
      implicit none
#include <netcdf.inc>
#include <dsttvbds.h> /* Time-varying boundary data sets */
c     Parameters
c     Input
      real(r8) time_doy             ! [day] Day of year
c     Output (common block dstbnd)
c     Input/Output
c     Local
      integer lai_dst_id        ! [enm] Variable ID
      integer dat_nbr           ! Array size passed to interpolation routine
      integer idx_glb_dsk_old   ! Index of most recent past time slice in input file, old
      integer idx_lub_dsk_old   ! Index of nearest future time slice in input file, old
      integer rcd               ! [rcd] Return success code
      real(r8) crd_max              ! Maximum coordinate value for interpolation routine
      real(r8) crd_min              ! Minimum coordinate value for interpolation routine
      real(r8) crd_ntp              ! Arrival coordinate for interpolation routine
      real(r8) time_dlt             ! [day] Time interval between glb and lub data
      real(r8) time_dlt_1st         ! [day] Time interval between glb and current time
c     Externals
!     integer idx_glb_get       ! idx_glb_get() is defined in dstnc.F
      external nc_err_exit      ! nc_err_exit() is defined in dstnc.F
c     Main code

c     Save current indices of bounding time slices
      idx_lub_dsk_old=idx_lub_dsk
      idx_glb_dsk_old=idx_glb_dsk
c     Get new indices of bounding time slices
      idx_glb_dsk=idx_glb_get(time,time_nbr,time_doy)
      idx_lub_dsk=mod(idx_glb_dsk,time_nbr)+1

c     Is model time still bounded by dataset times?
      if (idx_glb_dsk /= idx_glb_dsk_old.or.idx_lub_dsk /= idx_lub_dsk_old) then

c     Get variable IDs (defensive programming, the dataset might have changed...)
         rcd=NF_NOERR           ! NF_NOERR == 0
         rcd=rcd+nf90_inq_varid(nc_id,'lai_lsm',lai_dst_id)
         if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ntp():varid lai_lsm')

         if (idx_lub_dsk /= idx_lub_dsk_old) then
c     Read in new "future" data, replace old "past" data

            idx_glb_ram=idx_lub_ram ! Swap pointers to "past" and "future" slices
            idx_lub_ram=mod(idx_glb_ram,2)+1 ! "future" index is whatever past is not

            srt_lon_lat_time(3)=idx_lub_dsk
            rcd=rcd+nf90_get_var(nc_id,lai_dst_id,lai_dst_bnd(:,:,idx_lub_ram),srt_lon_lat_time,cnt_lon_lat_time)
            if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ntp():get_vara lai_dst')
            write (6,'(2(a,f9.4))') 'dst: dst_tvbds_ntp(): Day of year = ',time_doy,': read data for day ',time(idx_lub_dsk)
            
c     Update glb
c     NB: This is defensive programming---this code should never need to be executed
            if (idx_glb_dsk /= idx_lub_dsk_old) then
               srt_lon_lat_time(3)=idx_glb_dsk
               rcd=rcd+nf90_get_var(nc_id,lai_dst_id,lai_dst_bnd(:,:,idx_glb_ram),srt_lon_lat_time,cnt_lon_lat_time)
               if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tvbds_ntp():get_vara lai_dst')
               write (6,'(2(a,f9.4))') 'dst: dst_tvbds_ntp(): WARNING Day of year = ',time_doy,
     $              ': read data for day ',time(idx_glb_dsk)
            endif               ! endif reading new glb
         endif                  ! endif reading new lub
      endif                     ! endif updating either glb or lub

c     Renormalize coordinate distance to unity to facilitate interpolation 
      if (time(idx_lub_dsk) < time(idx_glb_dsk)) then ! i.e., Dec 15 -- Jan 15
         time_dlt=365.0-time(idx_glb_dsk)+time(idx_lub_dsk)
         if (time_doy <= time(idx_lub_dsk)) then ! i.e., Jan 1 -- Jan 15
            time_dlt_1st=365.0-time(idx_glb_dsk)+time_doy
         else
            time_dlt_1st=time_doy-time(idx_glb_dsk) ! i.e., Dec 15 -- Dec 31
         end if                 ! endif
      else                      ! i.e., Jan 15 -- Dec 15
         time_dlt=time(idx_lub_dsk)-time(idx_glb_dsk)
         time_dlt_1st=time_doy-time(idx_glb_dsk)
      end if                    ! endif

c     Interpolate to current time
      dat_nbr=lat_nbr*lon_nbr
      crd_min=0.0
      crd_max=1.0
      crd_ntp=time_dlt_1st/time_dlt

      call ntp_arr(dat_nbr,lai_dst_bnd(1,1,idx_glb_ram),lai_dst_bnd(1,1,idx_lub_ram), ! I
     $     crd_min,crd_max,crd_ntp, ! I
     $     lai_dst)             ! O
    
      return
      end subroutine dst_tvbds_ntp                      ! end dst_tvbds_ntp()

      subroutine dst_tvbds_close()
c     Purpose: Close the annual cycle data file
c     dst_tvbds_close() is called by CCM:dynamics/eul/stepon(), MATCH:main()
      use netcdf
      implicit none
#include <netcdf.inc>
#include <dsttvbds.h> /* Time-varying boundary data sets */
c     Local
      integer rcd               ! [rcd] Return success code
c     Main Code
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf90_close(nc_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'Unable to close tvbds netCDF file')
      write (6,'(a25,1x,a)') 'dst: Closed tvbds netCDF file'
      return
      end subroutine dst_tvbds_close                      ! end dst_tvbds_close()

      subroutine dst_tvbds_get(lon_nbr_out,lat_idx, ! I
     $     lai_dst_out) ! O
c     Purpose: Return specifed latitude slice of surface boundary data
c     dst_tvbds_get() is called by dst_mbl()
c     NB: dst_tvbds_get() converts arrays, if necessary, from the netCDF
c     dataset resolution (lat_nbr = plat X lon_nbr = plon) into latitude 
c     slices dimensioned according to the model resolution (i.e., set 
c     lon_nbr_out = plond)
      implicit none
#include <dsttvbds.h> /* Time-varying boundary data sets */
c     Parameters
c     Input
      integer lat_idx           ! Latitude index
      integer lon_nbr_out       ! Dimension size (normally plond)
c     Output
      real(r8) lai_dst_out(lon_nbr_out) ! Seasonal cycle
c     Local
      integer lon_idx           ! [idx] Counting index
c     Externals
c     Main code

      if (lon_nbr_out /= lon_nbr+3) then
         write (6,'(a)') 'dst: ERROR dst_tvbds_get() reports lon_nbr_out /= lon_nbr+3'
         stop
      endif                     ! endif err

      do lon_idx=1,lon_nbr
         lai_dst_out(lon_idx)=lai_dst(lon_idx,lat_idx)
      end do                    ! end loop over lon

      return
      end subroutine dst_tvbds_get                   ! end dst_tvbds_get()

      subroutine ntp_arr(dat_nbr,dat_min,dat_max, ! I
     $     crd_min,crd_max,crd_ntp, ! I
     $     dat_out)             ! O
c     Purpose: Interpolate the input arrays to the given coordinate
c     ntp_arr() is called by dst_tvbds_ntp()
c     Input
      integer dat_nbr           ! Dimension size
      real(r8) crd_max              ! Maximum coordinate value for interpolation routine
      real(r8) crd_min              ! Minimum coordinate value for interpolation routine
      real(r8) crd_ntp              ! Arrival coordinate for interpolation routine
      real(r8) dat_max(dat_nbr)     ! Array values at crd_max
      real(r8) dat_min(dat_nbr)     ! Array values at crd_min
c     Output
      real(r8) dat_out(dat_nbr)     ! Array values at crd_ntp
c     Local
      integer idx               ! [idx] Counting index
      real(r8) crd_dlt              ! Interval between min and max coordinates
      real(r8) crd_dlt_1st          ! Interval between min and arrival coordinates
      real(r8) crd_fct              ! Factor crd_dlt_1st/crd_dlt

      if (crd_max <= crd_min) then
         write (6,'(a)') 'dst: ERROR bnt_ntp() reports crd_max <= crd_min'
         stop
      endif                     ! endif err

      crd_dlt=crd_max-crd_min
      crd_dlt_1st=crd_ntp-crd_min
      crd_fct=crd_dlt_1st/crd_dlt

      do idx=1,dat_nbr
         dat_out(idx)=dat_min(idx)+crd_fct*(dat_max(idx)-dat_min(idx))
      end do                    ! end loop over idx

      return
      end subroutine ntp_arr                      ! end ntp_arr()

      integer function idx_glb_get(dat,dat_nbr,val)
c     Purpose: Return the index of the largest data point (in dat) which is less than the specified value (val), i.e., find the index of the greatest lower bound of the specified value.
c     The data array is assumed to be cyclic so the routine should return dat_nbr when val < dat(1)
c     idx_glb_get() is called by dst_tvbds_ini() and dst_tvbds_ntp()
c     Input
      integer dat_nbr           ! Dimension size
      real(r8) dat(dat_nbr)         ! Monotonically increasing input array
      real(r8) val                  ! Value to bound
c     Local
      integer dat_idx           ! [idx] Counting index

      if (val < dat(1).or.val >= dat(dat_nbr)) then
         idx_glb_get=dat_nbr
         return
      endif                     ! endif
      do dat_idx=2,dat_nbr
         if (val < dat(dat_idx)) then
            idx_glb_get=dat_idx-1
            return
         endif                  ! endif
      end do                    ! end loop over dat
      write (6,'(a)') 'dst: ERROR idx_glb_get() unable to find greatest lower bound'
      stop
      end function idx_glb_get                      ! end idx_glb_get()

      end module dsttvbds
