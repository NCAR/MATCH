c     $Header$ -*-fortran-*-
      
c     Purpose: Common blocks in dstbdg store fields needed by mass diagnostic routines for dust parameterization
      
c     These variables are initialized with routine dst_bdg_cmn_ini() 
c     dst_bdg_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
      
c     Usage: 
c     #include <dstbdg.h> /* Mass budget common block */
      
c     dstbdg.h must have access to pmgrid.h to work
      
      integer,parameter::var_nbr=19 ! [nbr] Number of variables in mass budget file
      
      real(r8) grv_sfc_rcp          ! [s2 m-1] Reciprocal of gravity
      real(r8) grv_sfc_plon_rcp     ! [s2 m-1] 1.0/(grv_sfc*plon)
      real(r8) lat_wgt(plat)        ! [frc] Latitude weights (currently must sum to 2.0)
      real(r8) mss_lat(plat,var_nbr) ! [kg m-2] Mass path of tracer (modulo pre-factors)
      common /bdg_mss_0/ grv_sfc_rcp, grv_sfc_plon_rcp, lat_wgt, mss_lat
      
      integer time_idx          ! [idx] Current record index
      integer nc_id             ! File handle
      integer var_id(3+var_nbr) ! [enm] Variable IDs
      common /bdg_mss_1/ time_idx, nc_id, var_id
      
      character var_lst(var_nbr)*64 ! [sng] Variable names
      character fl_out*80       ! [sng] Output file
      common /bdg_mss_2/ var_lst, fl_out
      
c     tm_dlt must be carried in a common block since it is not passed each timestep
      real(r8) tm_dlt               ! [s] Timestep duration
      logical uninitialized     ! [flg] True until first bdg_update() call of run
      common /bdg_mss_3/        ! bdg_mss_3 holds variables for derived fields
     $     tm_dlt,              ! [s] Timestep duration
     $     uninitialized        ! [flg] True until first bdg_update() call of run
      
c     Variable IDs for derived (non-enumerated) variables are carried separately from var_id array
      integer dst_dry_dlt_id    ! [enm] Variable ID
      integer dst_mbl_dlt_id    ! [enm] Variable ID
      integer dst_pcp_dlt_id    ! [enm] Variable ID
      integer dst_sf_dps_id     ! [enm] Variable ID
      integer dst_sf_net_id     ! [enm] Variable ID
      integer flx_mss_dps_id    ! [enm] Variable ID
      integer flx_mss_dry_id    ! [enm] Variable ID
      integer flx_mss_grv_id    ! [enm] Variable ID
      integer flx_mss_mbl_id    ! [enm] Variable ID
      integer flx_mss_pcp_id    ! [enm] Variable ID
      integer flx_mss_trb_id    ! [enm] Variable ID
      integer mpc_dst_dgn_id    ! [enm] Variable ID
      integer mpc_dst_dlt_id    ! [enm] Variable ID
      integer mpc_dst_mdl_avg_id ! [enm] Variable ID
      integer mpc_err_id        ! [enm] Variable ID
      integer mpc_err_frc_id    ! [enm] Variable ID
      integer mpc_err_tnd_frc_id ! [enm] Variable ID
      integer tau_dps_id        ! [enm] Variable ID
      integer tau_dry_id        ! [enm] Variable ID
      integer tau_grv_id        ! [enm] Variable ID
      integer tau_mbl_id        ! [enm] Variable ID
      integer tau_pcp_id        ! [enm] Variable ID
      integer tau_trb_id        ! [enm] Variable ID
      integer tm_ttl_id         ! [enm] Variable ID
      common /bdg_mss_4/        ! bdg_mss_4 holds variable IDs for derived fields
     $     dst_dry_dlt_id,      ! [enm] Variable ID
     $     dst_mbl_dlt_id,      ! [enm] Variable ID
     $     dst_pcp_dlt_id,      ! [enm] Variable ID
     $     dst_sf_dps_id,       ! [enm] Variable ID
     $     dst_sf_net_id,       ! [enm] Variable ID
     $     flx_mss_dps_id,      ! [enm] Variable ID
     $     flx_mss_dry_id,      ! [enm] Variable ID
     $     flx_mss_grv_id,      ! [enm] Variable ID
     $     flx_mss_mbl_id,      ! [enm] Variable ID
     $     flx_mss_pcp_id,      ! [enm] Variable ID
     $     flx_mss_trb_id,      ! [enm] Variable ID
     $     mpc_dst_dgn_id,      ! [enm] Variable ID
     $     mpc_dst_dlt_id,      ! [enm] Variable ID
     $     mpc_dst_mdl_avg_id,  ! [enm] Variable ID
     $     mpc_err_id,          ! [enm] Variable ID
     $     mpc_err_frc_id,      ! [enm] Variable ID
     $     mpc_err_tnd_frc_id,  ! [enm] Variable ID
     $     tau_dps_id,          ! [enm] Variable ID
     $     tau_dry_id,          ! [enm] Variable ID
     $     tau_grv_id,          ! [enm] Variable ID
     $     tau_mbl_id,          ! [enm] Variable ID
     $     tau_pcp_id,          ! [enm] Variable ID
     $     tau_trb_id,          ! [enm] Variable ID
     $     tm_ttl_id            ! [enm] Variable ID

