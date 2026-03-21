c     $Header$ -*-fortran-*-

c     Purpose: dst_tvbds_cmn stores boundary data needed by dust parameterization

c     dsttvbds.h MUST have access to params.h.
#ifndef PARAMS_H
#include <params.h> /* Preprocessor tokens */ 
#endif /* not PARAMS_H */ 

c     Usage:
c#include <dsttvbds.h> /* Time-varying boundary data sets */

c     Time-varying boundary data sets
      integer lat_nbr
      integer lon_nbr
      integer time_nbr
      parameter(lat_nbr=PLAT,
     $     lon_nbr=PLON,
     $     time_nbr=12) 

c     NB: fl_nm is currently not copied into the common block
c     Pending addition of sng.F library in Dust model, it will be
      common / dst_tvbds_cmn_c /  ! dst_tvbds_cmn is initialized in dst_tvbds_ini()
     $     fl_nm                ! netCDF file containing TVBDS data
      character*80 fl_nm        ! netCDF file containing TVBDS data

      common / dst_tvbds_cmn_i /  ! dst_tvbds_cmn is initialized in dst_tvbds_ini()
     $     cnt_lon_lat_time,
     $     idx_glb_dsk,
     $     idx_glb_ram,
     $     idx_lub_dsk,
     $     idx_lub_ram,
     $     nc_id,
     $     srt_lon_lat_time
      integer cnt_lon_lat_time(3) ! Hyperslab size
      integer idx_glb_dsk       ! Index of most recent past time slice in input file
      integer idx_glb_ram       ! Index of most recent past time slice in memory
      integer idx_lub_dsk       ! Index of nearest future time slice in input file
      integer idx_lub_ram       ! Index of nearest future time slice in memory
      integer nc_id             ! File ID
      integer srt_lon_lat_time(3) ! Starting offsets

      common / dst_tvbds_cmn_r /  ! dst_tvbds_cmn is initialized in dst_tvbds_ntp()
     $     lai_dst,
     $     lai_dst_bnd,
     $     time
      real(r8) lai_dst(lon_nbr,lat_nbr) ! [m2 m-2] LAI, interpolated
      real(r8) lai_dst_bnd(lon_nbr,lat_nbr,2) ! [m2 m-2] LAI, bounding data
      real(r8) time(time_nbr)       ! Time coordinate (day of year)
      


