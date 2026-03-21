      module dstbdg

      use precision

      
      contains

c     $Header$
      
c     Purpose: Mass budget diagnostic subroutines for the dust model
c     Most routines are based on Brian Eaton's massbgt.F
      
#include <params.h>
#include <dst.h> /* Dust preprocessor tokens */
      
      subroutine bdg_cmn_ini(grv_sfc,lat_wgt_tmp,tm_dlt_tmp)
c     Initialize parameter-dependent fields in dust common block
c     bdg_cmn_ini() is called by CCM:physics/inti(), MATCH:inirun()
c     NB: MATCH reads delt as integer, CCM reads dtime as real
c     Be sure calling routine passes a real to bdg_cmn_ini()
      use netcdf
      implicit none
#include <dbg.com> /* Debugging constants */
#include <pmgrid.h>
#include <dstbdg.h> /* Mass budget common block */
#include <netcdf.inc>
c     Input
      real(r8),intent(in)::grv_sfc  ! [m s-2] Gravity at surface of Earth
      real(r8),intent(in)::lat_wgt_tmp(plat) ! [frc] Latitude weights (currently must sum to 2.0)
      real(r8),intent(in)::tm_dlt_tmp ! [s] Timestep duration
c     Local
      integer var_idx           ! [idx] Counting index
      integer lat_idx           ! [idx] Counting index
      integer rcd               ! [enm] Return success code
      integer time_dim_id       ! [enm] Dimension ID for time
      integer junk
      real(r8) lat_wgt_ttl          ! [frc] Total of latitude weights (currently must sum to 2.0)
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
c     Main code
c     Initialize defaults
      uninitialized=.true.      ! [flg] True until first bdg_update() call of run
c     Initializing the time index to 1 seems to leave record 0 blank in dst_mss_bdg.nc
      time_idx=1                ! [idx] Current record index
      rcd=NF_NOERR              ! [enm] NF_NOERR == 0
      fl_out='dst_mss_bdg.nc'   ! [sng] Output file
c     Set common block variables from inputs
      grv_sfc_rcp=1.0/grv_sfc   ! [s2 m-1] Reciprocal of gravity 
      tm_dlt=tm_dlt_tmp         ! [s] Timestep duration
      grv_sfc_plon_rcp=1.0/(plon*grv_sfc) ! [s2 m-1] 1.0/(grv_sfc*plon)
      do lat_idx=1,plat
c     Latitude weights should be (at least) double precision in production runs
         lat_wgt(lat_idx)=lat_wgt_tmp(lat_idx) ! [frc] Latitude weights (currently must sum to 2.0)
      end do                    ! end loop over lat
      call vec_set(mss_lat,plat*var_nbr,0.0_r8) ! [kg m-2]
c     Set up identifiers for mass calculations
      var_lst( 1)='dst_sf_dry'
      var_lst( 2)='dst_sf_mbl'
      var_lst( 3)='dst_sf_grv'
      var_lst( 4)='dst_sf_trb'
      var_lst( 5)='dst_sf_pcp'
      var_lst( 6)='dst_ss_dry'
      var_lst( 7)='dst_ss_pcp'
      var_lst( 8)='dst_ss_evp'
      var_lst( 9)='dst_ss_mbl'
c     Which column mass path diagnostics are output is model-specific:
c     MATCH currently outputs: mpc_dst_mdl,mpc_dst_a_fxr,mpc_dst_a_trn,mpc_dst_a_trn_qng,mpc_dst_a_dry,mpc_dst_a_mbl,mpc_dst_a_wet
c     CCM currently outputs: mpc_dst_mdl,mpc_dst_a_fxr,mpc_dst_a_flt,mpc_dst_a_tphys
c     BXM currently outputs: mpc_dst_mdl,mpc_dst_a_fxr,mpc_dst_a_flt,mpc_dst_a_tphys,mpc_dst_a_dry,mpc_dst_a_mbl,mpc_dst_a_wet
      var_lst(10)='mpc_dst_mdl' ! [kg m-2] Column mass path of dust at end of timestep 
      var_lst(11)='mpc_dst_a_fxr' ! [kg m-2] Column mass path of dust after mass fixer
      var_lst(12)='mpc_dst_a_flt' ! [kg m-2] Column mass path of dust after time filter
      var_lst(13)='mpc_dst_b_tphys' ! [kg m-2] Column mass path of dust after tphys()
      var_lst(14)='mpc_dst_a_trn' ! [kg m-2] Column mass path of dust after transport
      var_lst(15)='mpc_dst_a_trn_qng' ! [kg m-2] Column mass path of dust after qneg check for transport
      var_lst(16)='mpc_dst_a_dry' ! [kg m-2] Column mass path of dust after dry deposition
      var_lst(17)='mpc_dst_a_mbl' ! [kg m-2] Column mass path of dust after mobilization
      var_lst(18)='mpc_dst_a_wet' ! [kg m-2] Column mass path of dust after wet deposition
      var_lst(19)='mpc_H2O'     ! [kg m-2] Column mass path of H2O vapor
      
c     Sanity checks
      if (var_nbr /= 19) stop 'dst: var_nbr /= 19 in bdg_cmn_ini()'
      if (tm_dlt <= 0.0 .or. tm_dlt > 3600.0) then
         write(6,'(a,f9.3,a2)') 'dst: bdg_cmn_ini() reports unreasonable timestep duration tm_dlt = ',tm_dlt,' s'
         stop
      endif                     ! endif err
      lat_wgt_ttl=0.0           ! [frc] Total of latitude weights (currently must sum to 2.0)
      do lat_idx=1,plat
         lat_wgt_ttl=lat_wgt_ttl+lat_wgt(lat_idx) ! [frc] Total of latitude weights (currently must sum to 2.0)
      end do                    ! end loop over lat
      if (abs(2.0-lat_wgt_ttl)/2.0 > 2*epsilon(1.0)) then
         write (6,'(a,f9.7,a)') 'dst: bdg_cmn_ini() reports lat_wgt_ttl = ',lat_wgt_ttl,' != 2.0'
         write (6,*) ' smallest real is ', tiny(1.0), epsilon(1.0)
         stop
      endif                     ! endif err
      
c     Create netCDF output file
      call strnul(fl_out)
      rcd=rcd+nf90_create(fl_out,NF_CLOBBER,nc_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'Unable to open file '//fl_out)
      write(6,'(a,a28,1x,a)')  prg_nm(1:f_strlen(prg_nm)),': Created netCDF output file',fl_out(1:f_strlen(fl_out))
      
c     Define dimensions
      rcd=rcd+nf90_def_dim(nc_id,"time",NF_UNLIMITED,time_dim_id)
c     Define variables
      rcd=rcd+nf90_def_var(nc_id,'nstep',NF90_INT,time_dim_id,var_id(1))
      rcd=rcd+nf90_def_var(nc_id,'date',NF_INT,time_dim_id,var_id(2))
      rcd=rcd+nf90_def_var(nc_id,'datesec',NF_INT,time_dim_id,var_id(3))
      do var_idx=1,var_nbr
         rcd=rcd+nf90_def_var(nc_id,var_lst(var_idx),NF_DOUBLE,time_dim_id,var_id(3+var_idx))
      end do                    ! end loop over var
c     Derived budget fields defined in bdg_update()
      rcd=rcd+nf90_def_var(nc_id,'dst_dry_dlt',NF_DOUBLE,time_dim_id,dst_dry_dlt_id)
      rcd=rcd+nf90_def_var(nc_id,'dst_mbl_dlt',NF_DOUBLE,time_dim_id,dst_mbl_dlt_id)
      rcd=rcd+nf90_def_var(nc_id,'dst_pcp_dlt',NF_DOUBLE,time_dim_id,dst_pcp_dlt_id)
      rcd=rcd+nf90_def_var(nc_id,'dst_sf_dps',NF_DOUBLE,time_dim_id,dst_sf_dps_id)
      rcd=rcd+nf90_def_var(nc_id,'dst_sf_net',NF_DOUBLE,time_dim_id,dst_sf_net_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_mss_dps',NF_DOUBLE,time_dim_id,flx_mss_dps_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_mss_dry',NF_DOUBLE,time_dim_id,flx_mss_dry_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_mss_grv',NF_DOUBLE,time_dim_id,flx_mss_grv_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_mss_mbl',NF_DOUBLE,time_dim_id,flx_mss_mbl_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_mss_pcp',NF_DOUBLE,time_dim_id,flx_mss_pcp_id)
      rcd=rcd+nf90_def_var(nc_id,'flx_mss_trb',NF_DOUBLE,time_dim_id,flx_mss_trb_id)
      rcd=rcd+nf90_def_var(nc_id,'mpc_dst_dgn',NF_DOUBLE,time_dim_id,mpc_dst_dgn_id)
      rcd=rcd+nf90_def_var(nc_id,'mpc_dst_dlt',NF_DOUBLE,time_dim_id,mpc_dst_dlt_id)
      rcd=rcd+nf90_def_var(nc_id,'mpc_dst_mdl_avg',NF_DOUBLE,time_dim_id,mpc_dst_mdl_avg_id)
      rcd=rcd+nf90_def_var(nc_id,'mpc_err',NF_DOUBLE,time_dim_id,mpc_err_id)
      rcd=rcd+nf90_def_var(nc_id,'mpc_err_frc',NF_DOUBLE,time_dim_id,mpc_err_frc_id)
      rcd=rcd+nf90_def_var(nc_id,'mpc_err_tnd_frc',NF_DOUBLE,time_dim_id,mpc_err_tnd_frc_id)
      rcd=rcd+nf90_def_var(nc_id,'tau_dps',NF_DOUBLE,time_dim_id,tau_dps_id)
      rcd=rcd+nf90_def_var(nc_id,'tau_dry',NF_DOUBLE,time_dim_id,tau_dry_id)
      rcd=rcd+nf90_def_var(nc_id,'tau_grv',NF_DOUBLE,time_dim_id,tau_grv_id)
      rcd=rcd+nf90_def_var(nc_id,'tau_mbl',NF_DOUBLE,time_dim_id,tau_mbl_id)
      rcd=rcd+nf90_def_var(nc_id,'tau_pcp',NF_DOUBLE,time_dim_id,tau_pcp_id)
      rcd=rcd+nf90_def_var(nc_id,'tau_trb',NF_DOUBLE,time_dim_id,tau_trb_id)
      rcd=rcd+nf90_def_var(nc_id,'tm_ttl',NF_DOUBLE,time_dim_id,tm_ttl_id)
      
c     Add descriptions
c     rcd=rcd+nf90_put_att(nc_id,var_id(3+var_idx),'units','kilogram meter-2')
c     End definitions
      rcd=rcd+nf90_enddef(nc_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'bdg_cmn_ini(): Unable to end define mode')
      return
      end subroutine bdg_cmn_ini                      ! end bdg_cmn_ini()
      
      subroutine bdg_gam_dry(var_nm,lat_idx,prs_dlt,q_H2O_vpr_mmrm,mmrd)
c     Purpose: Calculate latitude contribution to global average mass of tracer field
c     Input units of tracer are dry mass mixing ratio (kg tracer / kg dry air)
c     Input units of specific humidity moist mass mixing ratio (kg H2O vapor/ kg moist air)
c     Output masses are in [kg m-2]
c     bdg_gam_dry() is called by host model once for each latitude slice for each 3D mmrd diagnostic field each timestep
      implicit none
#include <pmgrid.h>
#include <dstbdg.h> /* Mass budget common block */
c     Input
      character,intent(in)::var_nm*(*) ! [sng] Variable name
      integer,intent(in)::lat_idx ! [idx] Latitude index
      real(r8),intent(in)::prs_dlt(plond,plev) ! [Pa] Pressure thickness
      real(r8),intent(in)::q_H2O_vpr_mmrm(plond,plev) ! [kg kg-1] Specific humidity (kg H2O vapor/kg moist air)
      real(r8),intent(in)::mmrd(plond,plev) ! [kg kg-1] Dry mass mixing ratio
c     Local
      integer lon_idx           ! [idx] Counting index
      integer lev_idx           ! [idx] Counting index
      integer var_idx           ! [idx] Counting index
c     Externals
!      integer bdg_idx_get       ! [idx] Index into field list
!      external bdg_idx_get      ! [idx] Index into field list
c     Main Code
      var_idx=bdg_idx_get(var_nm) ! [idx] Counting index
      mss_lat(lat_idx,var_idx)=0.0 ! [kg m-2] Mass path of tracer (modulo pre-factors)
      do lev_idx=1,plev
         do lon_idx=1,plon
            mss_lat(lat_idx,var_idx)=mss_lat(lat_idx,var_idx)+
     $           mmrd(lon_idx,lev_idx)*(1.0-q_H2O_vpr_mmrm(lon_idx,lev_idx))*prs_dlt(lon_idx,lev_idx) ! [kg m-2] Mass path of tracer (modulo pre-factors)
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      mss_lat(lat_idx,var_idx)=mss_lat(lat_idx,var_idx)*lat_wgt(lat_idx) ! [kg m-2] Mass path of tracer (modulo pre-factors)
      return
      end subroutine bdg_gam_dry                      ! end bdg_gam_dry()
      
      subroutine bdg_gam_wet(var_nm,lat_idx,prs_dlt,mmrm)
c     Purpose: Calculate latitude contribution to global average mass of tracer field
c     Input units of tracer are moist mass mixing ratio (kg tracer / kg moist air)
c     Output masses are in [kg m-2]
c     bdg_gam_wet() is called by host model once for each latitude slice for each 3D mmrm diagnostic field each timestep
      implicit none
#include <pmgrid.h>
#include <dstbdg.h> /* Mass budget common block */
c     Input
      character,intent(in)::var_nm*(*) ! [sng] Variable name
      integer,intent(in)::lat_idx ! [idx] Latitude index
      real(r8),intent(in)::prs_dlt(plond,plev) ! [Pa] Pressure thickness
      real(r8),intent(in)::mmrm(plond,plev) ! [kg kg-1] Moist mass mixing ratio
c     Local
      integer lon_idx           ! [idx] Counting index
      integer lev_idx           ! [idx] Counting index
      integer var_idx           ! [idx] Counting index
c     Externals
!      integer bdg_idx_get       ! [idx] Index into field list
!      external bdg_idx_get      ! [idx] Index into field list
c     Main code
      var_idx=bdg_idx_get(var_nm) ! [idx] Counting index
      mss_lat(lat_idx,var_idx)=0.0 ! [kg m-2] Mass path of tracer (modulo pre-factors)
      do lev_idx=1,plev
         do lon_idx=1,plon
            mss_lat(lat_idx,var_idx)=mss_lat(lat_idx,var_idx)+mmrm(lon_idx,lev_idx)*prs_dlt(lon_idx,lev_idx) ! [kg m-2] Mass path of tracer (modulo pre-factors)
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      mss_lat(lat_idx,var_idx)=mss_lat(lat_idx,var_idx)*lat_wgt(lat_idx) ! [kg m-2] Mass path of tracer (modulo pre-factors)
      return
      end subroutine bdg_gam_wet                      ! end bdg_gam_wet()
      
      subroutine bdg_gam_wet_2D(var_nm,lat_idx,prs_dlt,mmrm)
c     Purpose: Calculate latitude contribution to global average mass of a surface tracer field
c     Input concentration units are moist mixing ratio (kg tracer / kg moist air)
c     Output masses are in [kg m-2]
c     bdg_gam_wet() is called by host model once for each latitude slice for each 2D mmrm diagnostic field each timestep
      implicit none
#include <pmgrid.h>
#include <dstbdg.h> /* Mass budget common block */
c     Input
      character,intent(in)::var_nm*(*) ! [sng] Variable name
      integer,intent(in)::lat_idx ! [idx] Latitude index
      real(r8),intent(in)::prs_dlt(plond) ! [Pa] Pressure thickness
      real(r8),intent(in)::mmrm(plond) ! [kg kg-1] Moist mass mixing ratio
c     Local
      integer lon_idx           ! [idx] Counting index
      integer var_idx           ! [idx] Counting index
c     Externals
!      integer bdg_idx_get       ! [idx] Index into field list
!     external bdg_idx_get      ! [idx] Index into field list
c     Main code
      var_idx=bdg_idx_get(var_nm) ! [idx] Counting index
      mss_lat(lat_idx,var_idx)=0.0 ! [kg m-2] Mass path of tracer (modulo pre-factors)
      do lon_idx=1,plon
         mss_lat(lat_idx,var_idx)=mss_lat(lat_idx,var_idx)+mmrm(lon_idx)*prs_dlt(lon_idx) ! [kg m-2] Mass path of tracer (modulo pre-factors)
      end do                    ! end loop over lon
      mss_lat(lat_idx,var_idx)=mss_lat(lat_idx,var_idx)*lat_wgt(lat_idx) ! [kg m-2] Mass path of tracer (modulo pre-factors)
      return
      end subroutine bdg_gam_wet_2D                      ! end bdg_gam_wet_2D()
      
      subroutine bdg_aal(var_nm,lat_idx,var_dat)
c     Purpose: Store in memory the latitude contribution to area average of field
c     bdg_aal() takes zonal average of latitude slice of a 2D field and stores it
c     bdg_gam() is called by host model once for each latititude slice for each 2D diagnostic field each timestep
      implicit none
#include <pmgrid.h>
#include <dstbdg.h> /* Mass budget common block */
c     Input
      character,intent(in)::var_nm*(*) ! [sng] Variable name
      integer,intent(in)::lat_idx ! [idx] Latitude index
      real(r8),intent(in)::var_dat(plond) ! Field to average
c     Local
      integer lon_idx           ! [idx] Counting index
      integer var_idx           ! [idx] Variable index
c     Externals
!      integer bdg_idx_get       ! [idx] Index into field list
!      external bdg_idx_get      ! [idx] Index into field list
c     Main Code
      var_idx=bdg_idx_get(var_nm) ! [idx] Counting index
      mss_lat(lat_idx,var_idx)=0.0 ! [kg m-2] Mass path of tracer (modulo pre-factors)
      do lon_idx=1,plon
         mss_lat(lat_idx,var_idx)=mss_lat(lat_idx,var_idx)+var_dat(lon_idx) ! [kg m-2] Mass path of tracer (modulo pre-factors)
      end do                    ! end loop over lon
      mss_lat(lat_idx,var_idx)=mss_lat(lat_idx,var_idx)*lat_wgt(lat_idx) ! [kg m-2] Mass path of tracer (modulo pre-factors)
      return
      end subroutine bdg_aal                      ! end bdg_aal()
      
      real(r8) function bdg_gam(var_nm)
c     Purpose: Return global average mass (kg m-2) of a tracer field
c     bdg_gam() is called by bdg_update() once for each 3D field each timestep
c     Each latitude's contribution to the field must be in memory prior to this call
c     Routine usually used to convert 3D MMR field to global average MPC field
      implicit none
#include <pmgrid.h>
#include <dstbdg.h> /* Mass budget common block */
c     Input
      character,intent(in)::var_nm*(*) ! [sng] Variable name
c     Local
      integer lat_idx           ! [idx] Counting index
      integer var_idx           ! [idx] Counting index
c     Externals
!      integer bdg_idx_get       ! [idx] Index into field list
!      external bdg_idx_get      ! [idx] Index into field list
c     Main Code
      var_idx=bdg_idx_get(var_nm) ! [idx] Counting index
      bdg_gam=0.0
      do lat_idx=1,plat
         bdg_gam=bdg_gam+mss_lat(lat_idx,var_idx) ! [kg m-2]
      end do                    ! end loop over lat
c     0.5 factor assumes lat_wgt sums to 2.0
      bdg_gam=bdg_gam*grv_sfc_plon_rcp*0.5 ! [kg m-2]
      return
      end function bdg_gam                      ! end bdg_gam()
      
      real(r8) function bdg_aa(var_nm)
c     Purpose: Return area average of a field
c     Convert 2D MPC to global average MPC
c     bdg_aa() is called by bdg_update() once for each field each timestep
c     Each latitude's contribution to the field must be in memory prior to this call
      implicit none
#include <pmgrid.h>
#include <dstbdg.h> /* Mass budget common block */
c     Input
      character,intent(in)::var_nm*(*) ! [sng] Variable name
c     Local
      integer var_idx           ! [idx] Variable index
      integer lat_idx
c     Externals
!      integer bdg_idx_get       ! [idx] Index into field list
!      external bdg_idx_get      ! [idx] Index into field list
c     Main code
      var_idx=bdg_idx_get(var_nm) ! [idx] Counting index
      bdg_aa=0.0                ! [kg m-2]
      do lat_idx=1,plat
         bdg_aa=bdg_aa+mss_lat(lat_idx,var_idx) ! [kg m-2]
      end do                    ! end loop over lat
c     0.5 factor assumes lat_wgt sums to 2.0
      bdg_aa=bdg_aa*0.5/plon    ! [kg m-2]
      return
      end function bdg_aa                      ! end bdg_aa()
      
      integer function bdg_idx_get(var_nm)
c     Purpose: Return index of var_nm string in var_lst
      implicit none
#include <pmgrid.h>
#include <dstbdg.h> /* Mass budget common block */
c     Input
      character,intent(in)::var_nm*(*) ! [sng] Variable name
c     Local
      integer var_idx           ! [idx] Counting index
c     Main code
      bdg_idx_get=0
      do var_idx=1,var_nbr
         if (var_nm == var_lst(var_idx)) then
            bdg_idx_get=var_idx
            return
         end if                 ! endif
      end do                    ! end loop over var
      write(6,'(2a)') 'ERROR bdg_idx_get(): Variable name not found:',var_nm
      stop
      end function bdg_idx_get                      ! end bdg_idx_get()
      
      subroutine bdg_update(
     $     nstep,               ! I [idx] Timestep of host model
     $     date,                ! I [YYMMDD] Current date in YYMMDD format
     $     datesec)             ! I [s] Seconds past current date at 0Z
c     Purpose: Write global average budget terms to netCDF file
c     bdg_update() is called by CCM:dynamics/eul/stepon(), MATCH:main()
c     bdg_update() should be called each timestep
c     A few diagnostic fields in bdg_update() have hysterisis and are incorrect after restarts
      use netcdf
      use precision
      implicit none
#include <pmgrid.h>
#include <dstbdg.h> /* Mass budget common block */
#include <netcdf.inc>
c     Parameters
      real(r8),parameter::mss_val_dbl=1.0e36 ! [frc] Missing value in real(r8)
c     Input variables
      integer,intent(in)::nstep ! I [idx] Timestep of host model
      integer,intent(in)::date  ! I [YYMMDD] Current date in yymmdd format
      integer,intent(in)::datesec ! I [s] Seconds past current date at 0Z
c     Local
      integer hdr_nbr           ! [nbr] Number of header fields
      integer rcd               ! [rcd] Return success code
      real(r8) wgt_pst  ! [frc] Weight of past in running average (fragile)
      real(r8) wgt_crr  ! [frc] Weight of current timestep in running average (fragile)
c     Fields stored in mass budget arrays
      real(r8) mpc_dst_mdl ! [kg m-2] bdg_gam('DSTQ')
      real(r8) dst_sf_dry ! [kg m-2 s-1] bdg_aa('DSTSFDRY')
      real(r8) dst_sf_grv ! [kg m-2 s-1] bdg_aa('DSTSFGRV')
      real(r8) dst_sf_mbl ! [kg m-2 s-1] bdg_aa('DSTSFMBL')
      real(r8) dst_sf_pcp ! [kg m-2 s-1] bdg_aa('DSTSFPCP')
      real(r8) dst_sf_trb ! [kg m-2 s-1] bdg_aa('DSTSFTRB')
      real(r8) dst_ss_evp ! [kg m-2] bdg_gam('DSTSSEVP')
      real(r8) dst_ss_dry ! [kg m-2] bdg_gam('DSTSSDRY')
      real(r8) dst_ss_mbl ! [kg m-2] bdg_gam('DSTSSMBL')
      real(r8) dst_ss_pcp ! [kg m-2] bdg_gam('DSTSSPCP')
c     Derived mass budget fields
      real(r8) dst_dry_dlt ! [kg m-2] dst_ss_dry-dst_sf_dry
      real(r8) dst_mbl_dlt ! [kg m-2] dst_ss_mbl-dst_sf_mbl
      real(r8) dst_pcp_dlt ! [kg m-2] dst_ss_pcp-dst_ss_evp-dst_sf_pcp
      real(r8) dst_sf_dps ! [kg m-2 s-1] dst_sf_dry+dst_sf_pcp (Positive definite) (fragile)
      real(r8) dst_sf_net ! [kg m-2 s-1] dst_sf_mbl-dst_sf_dps (Positive into atmosphere) (fragile)
      real(r8) mpc_dst_dgn ! [kg m-2] Time integral of dst_sf_net (fragile)
      real(r8) mpc_dst_dlt ! [kg m-2] Approximate error defined as mpc_dst_dgn-mpc_dst_mdl (fragile)
      real(r8) mpc_dst_mdl_avg ! [kg m-2] Running average column mass path (fragile)
      real(r8) flx_mss_dps ! [kg m-2 s-1] Running average deposition (fragile)
      real(r8) flx_mss_dry ! [kg m-2 s-1] Running average dry deposition (fragile)
      real(r8) flx_mss_grv ! [kg m-2 s-1] Running average gravitational deposition (fragile)
      real(r8) flx_mss_mbl ! [kg m-2 s-1] Running average mobilization (fragile)
      real(r8) flx_mss_pcp ! [kg m-2 s-1] Running average wet deposition (fragile)
      real(r8) flx_mss_trb ! [kg m-2 s-1] Running average turbulent deposition (fragile)
      real(r8) mpc_err  ! [kg m-2] Mass imbalance
      real(r8) mpc_err_frc ! [frc] Mass imbalance as fraction of mean column burden
      real(r8) mpc_err_tnd_frc ! [frc] Mean imbalance tendency as fraction of mean source tendency

      real(r8) tau_dps  ! [s] Timescale for deposition (fragile)
      real(r8) tau_dry  ! [s] Timescale for dry deposition (fragile)
      real(r8) tau_grv  ! [s] Timescale for gravitational deposition (fragile)
      real(r8) tau_mbl  ! [s] Timescale for mobilization (fragile)
      real(r8) tau_pcp  ! [s] Timescale for wet deposition (fragile)
      real(r8) tau_trb  ! [s] Timescale for turbulent deposition (fragile)
      real(r8) tm_ttl   ! [s] Elapsed time of simulation (fragile)
c     Statics
      save dst_sf_net           ! [kg m-2 s-1] dst_sf_mbl-dst_sf_dps (Positive into atmosphere)
      save mpc_dst_dgn          ! [kg m-2] Time integral of dst_sf_net
      save mpc_dst_mdl_avg      ! [kg m-2] Running average column mass path (fragile)
      save flx_mss_dps          ! [kg m-2 s-1] Running average deposition (fragile)
      save flx_mss_dry          ! [kg m-2 s-1] Running average dry deposition (fragile)
      save flx_mss_grv          ! [kg m-2 s-1] Running average gravitational deposition (fragile)
      save flx_mss_mbl          ! [kg m-2 s-1] Running average mobilization (fragile)
      save flx_mss_pcp          ! [kg m-2 s-1] Running average wet deposition (fragile)
      save flx_mss_trb          ! [kg m-2 s-1] Running average turbulent deposition (fragile)
      save tm_ttl               ! [s] Elapsed time of simulation (fragile)
c     Externals
!      integer bdg_idx_get       ! [idx] Index into field list
!      real(r8) bdg_aa               ! Convert 2D MPC to global average MPC
!      real(r8) bdg_gam              ! Convert 3D MMR to global average MPC
!      external bdg_aa           ! Convert 2D MPC to global average MPC
!      external bdg_idx_get      ! [idx] Index into field list
!      external bdg_gam          ! Convert 3D MMR to global average MPC
c     Main Code
      
c     Initialize defaults
      hdr_nbr=3                 ! [nbr] Number of header fields
      rcd=NF_NOERR              ! NF_NOERR == 0
      
c     Convert surface fluxes to global average surface fluxes
      dst_sf_mbl=bdg_aa('dst_sf_mbl') ! [kg m-2 s-1]
      dst_sf_grv=bdg_aa('dst_sf_grv') ! [kg m-2 s-1]
      dst_sf_dry=bdg_aa('dst_sf_dry') ! [kg m-2 s-1]
      dst_sf_trb=bdg_aa('dst_sf_trb') ! [kg m-2 s-1]
      dst_sf_pcp=bdg_aa('dst_sf_pcp') ! [kg m-2 s-1]
      mpc_dst_mdl=bdg_gam('mpc_dst_mdl') ! [kg m-2]
      
c     Convert source/sink tendencies to global average surface fluxes
      dst_ss_mbl=bdg_gam('dst_ss_mbl') ! [kg m-2 s-1]
      dst_ss_dry=bdg_gam('dst_ss_dry') ! [kg m-2 s-1]
      dst_ss_pcp=bdg_gam('dst_ss_pcp') ! [kg m-2 s-1]
      dst_ss_evp=bdg_gam('dst_ss_evp') ! [kg m-2 s-1]
      
c     Difference between diagnosed and modeled surface flux tendencies
      dst_mbl_dlt=dst_ss_mbl-dst_sf_mbl ! [kg m-2 s-1]
      dst_dry_dlt=dst_ss_dry-dst_sf_dry ! [kg m-2 s-1]
      dst_pcp_dlt=dst_ss_pcp-dst_ss_evp-dst_sf_pcp ! [kg m-2 s-1]
      
c     Because of timesplitting, position of next two lines is somewhat model-dependent
c     There seems to be no elegant way to code this to work satisfactorily on both MATCH and CCM
c     This is largely due to hairy CCM timesplitting
c     For CCM these lines should appear after mpc_dst_dlt computation so that this timestep's dst_sf_dps and dst_sf_net contribute to next timestep's mpc_dst_dgn and mpc_dst_dlt
c     For MATCH (and box model) these lines should appear before mpc_dst_dlt computation so that this timestep's dst_sf_dps and dst_sf_net contribute to this timestep's mpc_dst_dgn and mpc_dst_dlt
      dst_sf_dps=dst_sf_dry+dst_sf_pcp ! [kg m-2 s-1] dst_sf_dry+dst_sf_pcp (Positive definite)
      dst_sf_net=dst_sf_mbl-dst_sf_dps ! [kg m-2 s-1] dst_sf_mbl-dst_sf_dps (Positive into atmosphere)
      
c     In CCM, mpc_dst_mdl lags dst_mbl() and dst_sf_net by 1 timestep
c     In MATCH, mpc_dst_mdl is in sync with dst_mbl() and dst_sf_net
c     Using static local storage for mpc_dst_dgn is OK only for initial runs
c     Variables with unpredicatble behavior on restart are marked "fragile"
c     On restart runs, the diagnostics which rely on integration will break
c     because they are not stored in the restart file and so are re-initialized
c     to 0.0 in bdg_cmn_ini().
c     Uninitialized flag is .true. only on first call to bdg_update()
      if (uninitialized) then
c     Initialize running sums
c     Uncomment following line for CCM (and move dst_sf_net definition after this block)
c     mpc_dst_dgn=0.0        ! [kg m-2] Time integral of dst_sf_net
         mpc_dst_dgn=tm_dlt*dst_sf_net ! [kg m-2] Time integral of dst_sf_net (fragile)
         mpc_dst_mdl_avg=0.0    ! [kg m-2] Running average column mass path (fragile)
         tm_ttl=0.0             ! [s] Elapsed time of simulation (fragile)
         uninitialized=.false.  ! [flg] True until first bdg_update() call of run
      else
         mpc_dst_dgn=mpc_dst_dgn+tm_dlt*dst_sf_net ! [kg m-2] Time integral of dst_sf_net (fragile)
      endif                     ! endif
c     Difference between diagnosed and modeled atmospheric burden
c     Positive when integral of dst_sf_net exceeds mpc_dst_mdl
      mpc_dst_dlt=mpc_dst_dgn-mpc_dst_mdl ! [kg m-2] Approximate error defined as mpc_dst_dgn-mpc_dst_mdl (fragile)
      
c     Running averages
c     Integration method employed here may accumulate significant roundoff error 
      tm_ttl=tm_ttl+tm_dlt      ! [s] Elapsed time of simulation (fragile)
      wgt_pst=(time_idx-1.0)/time_idx ! [frc] Weight of past in running average (fragile)
      wgt_crr=1.0/time_idx      ! [frc] Weight of current timestep in running average (fragile)
      mpc_dst_mdl_avg=wgt_pst*mpc_dst_mdl_avg+wgt_crr*mpc_dst_mdl ! [kg m-2] Running average column mass path (fragile)
      flx_mss_dps=wgt_pst*flx_mss_dps+wgt_crr*dst_sf_dps ! [kg m-2 s-1] Running average deposition (fragile)
      flx_mss_dry=wgt_pst*flx_mss_dry+wgt_crr*dst_sf_dry ! [kg m-2 s-1] Running average dry deposition (fragile)
      flx_mss_grv=wgt_pst*flx_mss_grv+wgt_crr*dst_sf_grv ! [kg m-2 s-1] Running average gravitational deposition (fragile)
      flx_mss_mbl=wgt_pst*flx_mss_mbl+wgt_crr*dst_sf_mbl ! [kg m-2 s-1] Running average mobilization (fragile)
      flx_mss_pcp=wgt_pst*flx_mss_pcp+wgt_crr*dst_sf_pcp ! [kg m-2 s-1] Running average wet deposition (fragile)
      flx_mss_trb=wgt_pst*flx_mss_trb+wgt_crr*dst_sf_trb ! [kg m-2 s-1] Running average turbulent deposition (fragile)
      
c     Mass conservation
c     Mass imbalance is difference between time integrated mean net surface flux and current atmospheric burden
c     Mass imbalance is deviation from perfect mass conservation
      mpc_err=tm_ttl*(flx_mss_mbl-flx_mss_dps)-mpc_dst_mdl ! [kg m-2] Mass imbalance
      mpc_err_frc=mpc_err/mpc_dst_mdl_avg ! [frc] Mass imbalance as fraction of mean column burden
      mpc_err_tnd_frc=mpc_err/(flx_mss_mbl*tm_ttl) ! [frc] Mean imbalance tendency as fraction of mean source tendency

c     Diagnostic timescales
      tau_dps=mss_val_dbl       ! [s] Timescale for deposition (fragile)
      tau_dry=mss_val_dbl       ! [s] Timescale for dry deposition (fragile)
      tau_grv=mss_val_dbl       ! [s] Timescale for gravitational deposition (fragile)
      tau_mbl=mss_val_dbl       ! [s] Timescale for mobilization (fragile)
      tau_pcp=mss_val_dbl       ! [s] Timescale for wet deposition (fragile)
      tau_trb=mss_val_dbl       ! [s] Timescale for turbulent deposition (fragile)
      if (flx_mss_dps /= 0.0) tau_dps=mpc_dst_mdl_avg/flx_mss_dps ! [s] Timescale for deposition (fragile)
      if (flx_mss_dry /= 0.0) tau_dry=mpc_dst_mdl_avg/flx_mss_dry ! [s] Timescale for dry deposition (fragile)
      if (flx_mss_grv /= 0.0) tau_grv=mpc_dst_mdl_avg/flx_mss_grv ! [s] Timescale for gravitational deposition (fragile)
      if (flx_mss_mbl /= 0.0) tau_mbl=mpc_dst_mdl_avg/flx_mss_mbl ! [s] Timescale for mobilization (fragile)
      if (flx_mss_pcp /= 0.0) tau_pcp=mpc_dst_mdl_avg/flx_mss_pcp ! [s] Timescale for wet deposition (fragile)
      if (flx_mss_trb /= 0.0) tau_trb=mpc_dst_mdl_avg/flx_mss_trb ! [s] Timescale for turbulent deposition (fragile)
      
      rcd=rcd+nf90_put_var(nc_id,var_id(1),nstep)
      rcd=rcd+nf90_put_var(nc_id,var_id(2),date)
      rcd=rcd+nf90_put_var(nc_id,var_id(3),datesec)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_sf_mbl')),dst_sf_mbl)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_sf_dry')),dst_sf_dry)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_sf_grv')),dst_sf_grv)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_sf_trb')),dst_sf_trb)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_sf_pcp')),dst_sf_pcp)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_ss_mbl')),dst_ss_mbl)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_ss_dry')),dst_ss_dry)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_ss_pcp')),dst_ss_pcp)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_ss_evp')),dst_ss_evp)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_mdl')),mpc_dst_mdl)
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_fxr')),bdg_gam('mpc_dst_a_fxr'))
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_flt')),bdg_gam('mpc_dst_a_flt'))
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_b_tphys')),bdg_gam('mpc_dst_b_tphys'))
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_trn')),bdg_gam('mpc_dst_a_trn'))
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_trn_qng')),bdg_gam('mpc_dst_a_trn_qng'))
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_dry')),bdg_gam('mpc_dst_a_dry'))
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_mbl')),bdg_gam('mpc_dst_a_mbl'))
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_wet')),bdg_gam('mpc_dst_a_wet'))
      rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_H2O')),bdg_gam('mpc_H2O'))
      
c     Output derived diagnostic constants
      rcd=rcd+nf90_put_var(nc_id,dst_dry_dlt_id,dst_dry_dlt)
      rcd=rcd+nf90_put_var(nc_id,dst_mbl_dlt_id,dst_mbl_dlt)
      rcd=rcd+nf90_put_var(nc_id,dst_pcp_dlt_id,dst_pcp_dlt)
      rcd=rcd+nf90_put_var(nc_id,dst_sf_dps_id,dst_sf_dps)
      rcd=rcd+nf90_put_var(nc_id,dst_sf_net_id,dst_sf_net)
      rcd=rcd+nf90_put_var(nc_id,flx_mss_dps_id,flx_mss_dps)
      rcd=rcd+nf90_put_var(nc_id,flx_mss_dry_id,flx_mss_dry)
      rcd=rcd+nf90_put_var(nc_id,flx_mss_grv_id,flx_mss_grv)
      rcd=rcd+nf90_put_var(nc_id,flx_mss_mbl_id,flx_mss_mbl)
      rcd=rcd+nf90_put_var(nc_id,flx_mss_pcp_id,flx_mss_pcp)
      rcd=rcd+nf90_put_var(nc_id,flx_mss_trb_id,flx_mss_trb)
      rcd=rcd+nf90_put_var(nc_id,mpc_dst_dgn_id,mpc_dst_dgn)
      rcd=rcd+nf90_put_var(nc_id,mpc_dst_dlt_id,mpc_dst_dlt)
      rcd=rcd+nf90_put_var(nc_id,mpc_dst_mdl_avg_id,mpc_dst_mdl_avg)
      rcd=rcd+nf90_put_var(nc_id,mpc_err_id,mpc_err)
      rcd=rcd+nf90_put_var(nc_id,mpc_err_frc_id,mpc_err_frc)
      rcd=rcd+nf90_put_var(nc_id,mpc_err_tnd_frc_id,mpc_err_tnd_frc)
      rcd=rcd+nf90_put_var(nc_id,tau_dps_id,tau_dps)
      rcd=rcd+nf90_put_var(nc_id,tau_dry_id,tau_dry)
      rcd=rcd+nf90_put_var(nc_id,tau_grv_id,tau_grv)
      rcd=rcd+nf90_put_var(nc_id,tau_mbl_id,tau_mbl)
      rcd=rcd+nf90_put_var(nc_id,tau_pcp_id,tau_pcp)
      rcd=rcd+nf90_put_var(nc_id,tau_trb_id,tau_trb)
      rcd=rcd+nf90_put_var(nc_id,tm_ttl_id,tm_ttl)
      
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'bdg_update()')
      
c     Update record index
      time_idx=time_idx+1       ! [idx] Current record index
      
      return
      end subroutine bdg_update                      ! end bdg_update()
      
      subroutine bdg_close()
c     Purpose: Close mass budget netCDF file
c     bdg_close() is called by CCM:stepon(), MATCH:main()
      use netcdf
      implicit none
#include <dbg.com> /* Debugging constants */
#include <pmgrid.h>
#include <dstbdg.h> /* Mass budget common block */
#include <netcdf.inc>
c     Local
      integer rcd               ! [rcd] Return success code
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
c     Main code
c     Initialize defaults
      rcd=NF_NOERR              ! [enm] NF_NOERR == 0
      rcd=rcd+nf90_close(nc_id)
      if (rcd /= NF_NOERR) call nc_err_exit(rcd,'Unable to close file '//fl_out)
      call strnul(fl_out)
      write (6,'(a,a20,1x,a)') prg_nm(1:f_strlen(prg_nm)),': Closed netCDF file',fl_out(1:f_strlen(fl_out))
      return
      end subroutine bdg_close                      ! end bdg_close()
      
#ifndef BXM /* Routine requires qneg3(), only available in global models */
      subroutine qneg3_dgn(
     $     sbr_nm,              ! I [sng] Name of calling routine
     $     lat_idx,             ! I [idx] Latitude index
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     q_mmr,               ! I/O [kg kg-1] Constituent array
     $     tm_adj)              ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
c     Purpose: Diagnostic plugin replacement for qneg3() routine
c     Computes and archives tendencies due to qneg3()
c     Dependencies: Requires qneg3()
      implicit none
#include <pmgrid.h>
c     Input
      character*8,intent(in)::sbr_nm ! I [sng] Name of calling routine
      integer,intent(in)::lat_idx ! I [idx] Latitude index
      real(r8),intent(in)::tm_adj   ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
      real(r8),intent(in)::prs_dlt(plond,plev) ! I [Pa] Pressure thickness
c     Input/Output
      real(r8),intent(inout)::q_mmr(plond,plev,pcnst) ! I [kg kg-1] Constituent array
c     Local
      real(r8) qneg3_tnd(plond,plev,pcnst) ! [kg kg-1 s-1] Tendency due to qneg3() limiter
c     Main code
c     Same array holds initial value q_mmr and qneg3() tendency to save memory
      qneg3_tnd=q_mmr           ! [kg kg-1] Constituent array before qneg3()
      call qneg3(sbr_nm,lat_idx,q_mmr)
c     Compute tendency due to qneg3()
c     By definition this should be a postive tendency because mass after qneg3()
c     is always greater than mass before
      qneg3_tnd=(q_mmr-qneg3_tnd)/tm_adj ! [kg kg-1 s-1] Tendency due to qneg3() limiter
c     Archive tendency to global mass diagnostics
#ifdef DST_MSS_BDG
      call bdg_gam_wet('dst_ss_qneg3',lat_idx,prs_dlt,qneg3_tnd)
#endif /* not DST_MSS_BDG */
      end subroutine qneg3_dgn                       ! end qneg3_dgn()
#endif /* BXM */

      end module dstbdg
