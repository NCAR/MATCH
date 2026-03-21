c     $Header$

c     Purpose: dsttibds.F contains subroutines for reading time-invariant boundary
c     data sets into the dust model. 
c     These routines depend heavily on netCDF I/O

      subroutine dst_tibds_ini(
     $     fl_in)               ! I

      use netcdf
      use precision

c     Initialize time invariant land surface characteristics
c     dst_tibds_ini() is called by CCM:control/initext(), MATCH:src/inirun()
      implicit none
#include <dsttibds.h> /* Time-invariant boundary data sets */
#include <dstsfc.h> /* 2-D surface fields on PLON x PLAT grid */
#include <netcdf.inc>
#ifdef SCCM
c     The SCCM header file latlon.h contains two floats, datalat and datalon, 
c     which specify the lat and lon of the SCCM profile.
#include <latlon.h>
#endif /* not SCCM */
c     Parameters
      integer dim_nbr_max
      integer lat_nbr_max
      integer lon_nbr_max
      parameter(dim_nbr_max=3,  ! lon, lat, time
     $     lat_nbr_max=plat,    ! NB: This must change to enable SCCM to work
     $     lon_nbr_max=plon)
c     Input
      character fl_in*(*)       ! netCDF input file
c     Output
c     Input/Output
c     Local
      integer cnt_lon_lat(dim_nbr_max) ! Hyperslab size
      integer dst_idx           ! [idx] Counting index
      integer idx               ! [idx] Counting index
      integer lat_dim_id        ! [enm] Dimension ID for lat
      integer lat_id            ! [enm] Variable ID
      integer lat_idx           ! [idx] Counting index
      integer lat_nbr           ! Dimension size
      integer lnd_frc_dry_id    ! [enm] Variable ID
      integer lon_dim_id        ! [enm] Dimension ID for lon
      integer lon_id            ! [enm] Variable ID
      integer lon_idx           ! [idx] Counting index
      integer lon_nbr           ! Dimension size
      integer mbl_bsn_fct_id    ! [enm] Variable ID
      integer mss_frc_CaCO3_id  ! [enm] Variable ID
      integer mss_frc_cly_id    ! [enm] Variable ID
      integer mss_frc_snd_id    ! [enm] Variable ID
      integer nc_id             ! File handle
      integer rcd               ! [rcd] Return success code
      integer sfc_typ_id        ! [enm] Variable ID
      integer srt(dim_nbr_max)  ! Starting offsets
      logical flg_err           ! [flg] Bad input data
      real(r8) :: one = 1.0
      real(r8) :: zero = 0.0
      
c     Entire lat and lon arrays are retrieved from disk
      real(r8) lat(lat_nbr_max)     ! Entire lat and lon arrays are retrieved from disk
      real(r8) lon(lon_nbr_max)
      
c++pjr
      real(r8) fmin                 ! min fraction we can attribute to roundoff problems
      real(r8) fmax                 ! max fraction we can attribute to roundoff problems
c--pjr

c     Externals
      external nc_err_exit
      
#ifdef SCCM
      integer lat_idx_sccm      ! SCCM nearest latitude index into boundary data
      integer lon_idx_sccm      ! SCCM nearest longitude index into boundary data
      real(r8) dst_crr              ! Current least distance from gridpoint to target
      real(r8) dst_new              ! Distance of current gridpoint to target
#endif /* not SCCM */ 
c++mahow natalie mahowalds dust modifications
      real(r8) source(plon,plat), area(plon,plat)
      integer i,j
      common /source/ source
c--mahow


c     Initialize default values
      do idx=1,dim_nbr_max
         srt(idx)=1
      end do                    ! end loop over dim
      cnt_lon_lat(1)=plon
      cnt_lon_lat(2)=plat

c++mahow read in additional source field
      open(72,file='obs.source')
      write(*,*) 'reading obs.source'
      do j=1,plat
      read(72,*) (source(i,j),i=1,plon)
c      write(*,*) (source(i,j),i=1,plon)
      enddo
      write(*,*) 'finished reading source'
      close(72)
      open(72,file='ginouxbasin')
      do j=1,plat
      read(72,*) (area(i,j),i=1,plon)
      enddo
      close(72)
      write(*,*) 'finished reading ginoux basin'
c--mahow
c     Read in land surface types from high resolution grid
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf90_open(fl_in,NCNOWRIT,nc_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini()')
      write(6,'(a29,1x,a)') 'dst: Opened netCDF input file',fl_in




c     Get dimension IDs
      rcd=rcd+nf90_inq_dimid(nc_id,'lat',lat_dim_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini(): dimid lat')
      rcd=rcd+nf90_inq_dimid(nc_id,'lon',lon_dim_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini(): dimid lon')
c     Get dimension sizes
      rcd=rcd+nf90_inquire_dimension(nc_id,lat_dim_id,len=lat_nbr)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini(): dimlen lat')
      if (lat_nbr /= lat_nbr_max.and.lat_nbr /= 1) stop 'lat_nbr /= lat_nbr_max.and.lat_nbr /= 1'
      rcd=rcd+nf90_inquire_dimension(nc_id,lon_dim_id,len=lon_nbr)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini(): dimlen lon')
      if (lon_nbr /= lon_nbr_max.and.lon_nbr /= 1) stop 'lon_nbr /= lon_nbr_max.and.lon_nbr /= 1'
c     Get variable IDs
      rcd=rcd+nf90_inq_varid(nc_id,'lon',lon_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():varid lon')
      rcd=rcd+nf90_inq_varid(nc_id,'lat',lat_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():varid lat')
      rcd=rcd+nf90_inq_varid(nc_id,'lnd_frc_dry',lnd_frc_dry_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():varid lnd_frc_dry')
      rcd=rcd+nf90_inq_varid(nc_id,'mbl_bsn_fct',mbl_bsn_fct_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():varid mbl_bsn_fct')
      rcd=rcd+nf90_inq_varid(nc_id,'mss_frc_CaCO3',mss_frc_CaCO3_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():varid mss_frc_CaCO3')
      rcd=rcd+nf90_inq_varid(nc_id,'mss_frc_snd',mss_frc_snd_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():varid mss_frc_snd')
      rcd=rcd+nf90_inq_varid(nc_id,'mss_frc_cly',mss_frc_cly_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():varid mss_frc_cly')
      rcd=rcd+nf90_inq_varid(nc_id,'sfc_typ',sfc_typ_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():varid sfc_typ')
c     Get data
      rcd=rcd+nf90_get_var(nc_id,lat_id,lat)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():get_var lat')
      rcd=rcd+nf90_get_var(nc_id,lon_id,lon)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():get_var lon')
c     In SCCM, find nearest lat and lon indices into boundary data
#ifdef SCCM
      lat_idx_sccm=1
      dst_crr=abs(datalat-lat(lat_idx_sccm))
      do lat_idx=1,lat_nbr
         dst_new=abs(datalat-lat(lat_idx))
         if(dst_new < dst_crr) then 
            lat_idx_sccm=lat_idx
            dst_crr=dst_new
         endif
      end do                    ! end loop over lat

      lon_idx_sccm=1
      dst_crr=abs(datalon-lon(lon_idx_sccm))
      do lon_idx=1,lon_nbr
         dst_new=abs(datalon-lon(lon_idx))
         if(dst_new < dst_crr) then 
            lon_idx_sccm=lon_idx
            dst_crr=dst_new
         endif
      end do                    ! end loop over lon

c     Overwrite hyperslab indices with coordinates of SCCM column
      srt(1)=lon_idx_sccm
      srt(2)=1at_idx_sccm
#endif /* not SCCM */ 

      rcd=rcd+nf90_get_var(nc_id,lnd_frc_dry_id,lnd_frc_dry,srt,cnt_lon_lat)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():get_vara lnd_frc_dry')
      rcd=rcd+nf90_get_var(nc_id,mbl_bsn_fct_id,mbl_bsn_fct,srt,cnt_lon_lat)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():get_vara mbl_bsn_fct')
      rcd=rcd+nf90_get_var(nc_id,mss_frc_CaCO3_id,mss_frc_CaCO3,srt,cnt_lon_lat)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():get_vara mss_frc_CaCO3')
      rcd=rcd+nf90_get_var(nc_id,mss_frc_snd_id,mss_frc_snd,srt,cnt_lon_lat)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():get_vara mss_frc_snd')
      rcd=rcd+nf90_get_var(nc_id,mss_frc_cly_id,mss_frc_cly,srt,cnt_lon_lat)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():get_vara mss_frc_cly')
      rcd=rcd+nf90_get_var(nc_id,sfc_typ_id,sfc_typ,srt,cnt_lon_lat)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'dst_tibds_ini():get_vara sfc_typ')

      rcd=rcd+nf90_close(nc_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'Unable to close file '//fl_in)
      write (6,'(a25,1x,a)') 'dst: Ingested netCDF file',fl_in

c     Sanity check (can save lots of wasted time if boundary dataset is corrupt)
      flg_err=.false.           ! [flg] Bad input data
      fmin = -1.e-10                    ! minimum fraction we can attribute to roundoff
      fmax = 1. - fmin                  ! maximum fraction we can attribute to roundoff
      do lat_idx=1,lat_nbr
         do lon_idx=1,lon_nbr
            if (sfc_typ(lon_idx,lat_idx) < 0.or.sfc_typ(lon_idx,lat_idx) > 28) then
               write (6,*) ' dst: sanity sfc_typ ',lon_idx,lat_idx, sfc_typ(lon_idx,lat_idx)
            endif
            if (lnd_frc_dry(lon_idx,lat_idx) < fmin.or.lnd_frc_dry(lon_idx,lat_idx) > fmax) then
               write (6,*) ' dst: sanity lnd_frc_dry ',lon_idx,lat_idx, lnd_frc_dry(lon_idx,lat_idx)
               stop
            endif
            if (mbl_bsn_fct(lon_idx,lat_idx) < fmin.or.mbl_bsn_fct(lon_idx,lat_idx) > fmax) then
               write (6,*) ' dst: sanity mbl_bsn_fct ',lon_idx,lat_idx, mbl_bsn_fct(lon_idx,lat_idx)
               stop
            endif
            if (mss_frc_CaCO3(lon_idx,lat_idx) < fmin.or.mss_frc_CaCO3(lon_idx,lat_idx) > fmax) then
               write (6,*) ' dst: sanity mss_frc_CaCO3 ',lon_idx,lat_idx, mss_frc_CaCO3(lon_idx,lat_idx)
               stop
            endif
            if (mss_frc_cly(lon_idx,lat_idx) < fmin.or.mss_frc_cly(lon_idx,lat_idx) > fmax) then
               write (6,*) ' dst: sanity mss_frc_cly ',lon_idx,lat_idx, mss_frc_cly(lon_idx,lat_idx)
               stop
            endif
            if (mss_frc_snd(lon_idx,lat_idx) < fmin.or.mss_frc_snd(lon_idx,lat_idx) > fmax) then
               write (6,*) ' dst: sanity mss_frc_snd ',lon_idx,lat_idx, mss_frc_snd(lon_idx,lat_idx)
               stop
            endif
            sfc_typ(lon_idx,lat_idx) = max(0,min(28,sfc_typ(lon_idx,lat_idx)))
            lnd_frc_dry(lon_idx,lat_idx) = max(zero,min(one,lnd_frc_dry(lon_idx,lat_idx)))
            mbl_bsn_fct(lon_idx,lat_idx) = max(zero,min(one,mbl_bsn_fct(lon_idx,lat_idx)))
            mss_frc_CaCO3(lon_idx,lat_idx) = max(zero,min(one,mss_frc_CaCO3(lon_idx,lat_idx)))
            mss_frc_cly(lon_idx,lat_idx) = max(zero,min(one,mss_frc_cly(lon_idx,lat_idx)))
            mss_frc_snd(lon_idx,lat_idx) = max(zero,min(one,mss_frc_snd(lon_idx,lat_idx)))
         end do                 ! end loop over lon
      end do                    ! end loop over lat
      do j=1,plat
      do i=1,plon
         mbl_bsn_fct(i,j)=area(i,j)
	source(i,j)=area(i,j)
      enddo
      enddo
      write(*,*) ' over-rode mbl_bsn_fct and source with ginoux area '
      return
      end                       ! end dst_tibds_ini()
