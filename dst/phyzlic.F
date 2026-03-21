c     $Header$

c     Purpose: phyzlic.F acts like a timestepping routine which calls
c     all the dust source and sink tendencies

c     Requires dst.h for DST_MSS_BDG and DST_CHM
#include <dst.h> /* Dust preprocessor tokens */ 

      subroutine phyzlic(
     $     lat_idx)             ! I [idx] Latitude index
c     Purpose: Impersonate physics driving routine such as MATCH:src/physlic(),
c     and call all dust physics packages
c     phyzlic() is called by aer()
      use dstmbl
      use dstdpswet
      use dstdpsdry
      use precision
      use histout
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <aernvr.h> /* Aerosol environmental properties */
#include <dstnm.h> /* Nomenclature for outfld() */
#include <dstchm.h> /* Chemical properties of dust */
#include <dstctl.h> /* Model control variables */
c     Input
      integer lat_idx           ! L [idx] Latitude index
      integer i                 ! L [idx] Longitude index
      integer m                 ! L [idx] Constituent index
c     Local
      character fl_out*80       ! [sng] Name of netCDF output file
      real(r8)::cld_vlm(plond,plev) ! [frc] Cloud volume
      real(r8)::dst_trn_ttl(plond,plevp,plevp,bnd_nbr_LW) ! [frc] Total transmission through all size classes between interfaces k1 and k2
      real(r8)::mpc_dst(plond,dst_nbr) ! [kg m-2] Column mass path of dust
      real(r8)::mpc_dst_ttl(plond)   ! [kg m-2] Total column mass path of dust
      real(r8)::mpl_dst(plond,plev,dst_nbr) ! [kg m-2] Layer dust amount
      real(r8)::mpp_dst(plond,plevp,dst_nbr) ! [kg m-2] Dust path above kth interface level
      real(r8)::odxc_dst(plond,dst_nbr) ! [frc] Column dust optical depth
      real(r8)::odxc_dst_ttl(plond)  ! [frc] Total column dust optical depth
      real(r8)::q_dst_ttl(plond,plev) ! [kg kg-1] Total dust mixing ratio
      real(r8)::wnd_rfr(plond)       ! [m s-1] Wind speed at reference height

c     Main code

c     Convert advected tracers (except H2O vapor) from dry to moist mass mixing ratios
      call mmrd2mmrm(
     $     plev,                ! I [nbr] Number of levels
     $     plon,                ! I [nbr] Number of longitudes
     $     plond,               ! I [nbr] First dimension of arrays
     $     q_H2O_vpr,           ! I [kg H2O vapor kg-1 moist air] Specific humidity
     $     q_cst(1,1,dst_idx_srt), ! I [kg tracer kg-1 dry air] Tracer dry mass mixing ratio
     $     dst_nbr,             ! I [nbr] Number of tracers
     $     q_cst(1,1,dst_idx_srt)) ! O [kg tracer kg-1 moist air] Tracer moist mass mixing ratio

c     Driver for dust mobilization
c     dst_mbl() is called by CCM:physics/tphysbc(), MATCH:src/physlic()
      if (.true.) then         ! fxm fxm fxm
      call dst_mbl(obuf,
     $     doy,                 ! I [day] Day of year [1.0..366.0)
     $     hgt_mdp(1,plev),     ! I [m] Midpoint height above surface
     $     lat_idx,             ! I [idx] Model latitude index
     $     lat_rdn(lat_idx),    ! I [rdn] Latitude
     $     oro,                 ! I [frc] Orography
     $     prs_dlt(1,plev),     ! I [Pa] Pressure thickness
     $     prs_mdp(1,plev),     ! I [Pa] Pressure
     $     prs_ntf(1,plevp),    ! I [Pa] Surface pressure NB: plevp
     $     q_H2O_vpr(1,plev),   ! I [kg kg-1] Water vapor mixing ratio
     $     q_cst(:,plev,dst_idx_srt:dst_idx_end), ! I/O [kg kg-1] Dust mixing ratio
     $     snw_hgt_lqd,         ! I [m] Equivalent liquid water snow depth
     $     tm_adj,              ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     tpt_mdp(1,plev),     ! I [K] Temperature
     $     tpt_ptn(1,plev),     ! I [K] Potential temperature
     $     wnd_mrd_mdp,         ! I [m s-1] Meridional wind component
     $     wnd_znl_mdp)         ! I [m s-1] Zonal wind component
      endif                     ! fxm fxm fxm

#ifdef DST_MSS_BDG
c     Archive mass path of dust after mobilization
      call dst_add_lon_lev(q_cst(1,1,dst_idx_srt),q_dst_ttl)
      call bdg_gam_wet('mpc_dst_a_mbl',lat_idx,prs_dlt,q_dst_ttl) ! Tracer is moist mass mixing ratio
#endif /* not DST_MSS_BDG */

c      if (.false.) then         ! fxm fxm fxm
c     Driver for aerosol dry deposition
c     dst_dps_dry() is called by CCM:physics/tphysac(), MATCH:src/physlic()
      call dst_dps_dry(obuf,
     $     hgt_mdp,             ! I [m] Midpoint height above surface
     $     lat_idx,             ! I [idx] Model latitude index
     $     mno_lng,             ! O [m] Monin-Obukhov length
     $     oro,                 ! I [frc] Orography
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     prs_mdp,             ! I [Pa] Midlayer pressure
     $     q_H2O_vpr,           ! I [kg kg-1] Water vapor mixing ratio
     $     q_cst(1,1,dst_idx_srt), ! I/O [kg kg-1] Dust mixing ratio
     $     snw_hgt_lqd,         ! I [m] Equivalent liquid water snow depth
     $     tm_adj,              ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     tpt_mdp,             ! I [K] Temperature
     $     tpt_ptn(1,plev),     ! I [K] Potential temperature
     $     tpt_sfc,             ! I [K] Surface temperature
     $     wnd_frc,             ! O [m s-1] Surface friction velocity
     $     wnd_mrd_mdp,         ! I [m s-1] Meridional wind component
     $     wnd_rfr,             ! O [m s-1] Wind speed at reference height
     $     wnd_znl_mdp)         ! I [m s-1] Zonal wind component
c      endif                     ! fxm fxm fxm

#ifdef DST_MSS_BDG
c     Archive mass path of dust after dry deposition
      call dst_add_lon_lev(q_cst(1,1,dst_idx_srt),q_dst_ttl)
      call bdg_gam_wet('mpc_dst_a_dry',lat_idx,prs_dlt,q_dst_ttl) ! Tracer is moist mass mixing ratio
#endif /* not DST_MSS_BDG */

#ifdef DST_MSS_BDG
c     After mass fixer, before aphys(), called in CCM:dynamics/eul/linemsbc()
      call dst_add_lon_lev(q_cst(1,1,dst_idx_srt),q_dst_ttl) ! q3
      call bdg_gam_dry('mpc_dst_a_fxr',lat_idx,prs_dlt,q_cst,q_dst_ttl) ! q3
#endif /* not DST_MSS_BDG */

c     Cloud diagnostics needed for wet deposition
c      if (.false.) then         ! fxm: remove old wet dep routine soon
      call cld_dgn( 
     $     cld_frc,             ! I [frc] Local total cloud fraction
     $     cld_vlm,             ! O [frc] Cloud volume
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_cnd2pcp_tnd,   ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
     $     q_H2O_cnd_tnd,       ! I [kg kg-1 s-1] Net H2O condensate formation tendency
     $     q_H2O_pcp_lqd,       ! O [kg kg-1] Rainwater mixing ratio
     $     q_H2O_pcp2vpr_tnd,   ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
     $     q_H2O_vpr2pcp_cnv_tnd, ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
     $     tpt_mdp)             ! I [K] Temperature

c     Wet deposition
c     dst_dps_wet() called by CCM:physics/aphys(), MATCH:physlic()
      call dst_dps_wet(obuf,
     $     cld_frc,             ! I [frc] Local total cloud fraction
     $     cld_frc_cnv,         ! I/O [frc] Convective cloud fraction
     $     cld_vlm,             ! I [frc] Cloud volume
     $     lat_idx,             ! I [idx] Model latitude index
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_cnd,           ! I [kg kg-1] Condensed H2O mixing ratio
     $     q_H2O_cnd2pcp_tnd,   ! I/O [kg kg-1 s-1] Condensed H2O to precipitation tendency
     $     q_H2O_pcp2vpr_tnd,   ! I/O [kg kg-1 s-1] H2O precipitation to vapor tendency
     $     q_H2O_vpr2pcp_cnv_tnd, ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
     $     q_cst(1,1,dst_idx_srt), ! I/O [kg kg-1] Dust mixing ratio
     $     tm_adj,              ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     tpt_mdp)             ! I [K] Temperature
c      endif                     ! endif

      if (.false.) then         ! fxm: remove old wet dep routine soon
c     Old wet deposition routine
      call dst_dps_wet_old(obuf,
     $     cld_frc,             ! I [frc] Local total cloud fraction
     $     lat_idx,             ! I [idx] Model latitude index
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_cnd2pcp_tnd,   ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
     $     q_H2O_pcp2vpr_tnd,   ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
     $     q_H2O_vpr2pcp_cnv_tnd, ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
     $     q_cst(1,1,dst_idx_srt), ! I/O [kg kg-1] Dust mixing ratio
     $     tm_adj,              ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     tpt_vrt)             ! I [K] Virtual temperature
      endif                     ! endif false

#ifdef DST_MSS_BDG
c     Archive mass path of dust after wet deposition
      call dst_add_lon_lev(q_cst(1,1,dst_idx_srt),q_dst_ttl)
      call bdg_gam_wet('mpc_dst_a_wet',lat_idx,prs_dlt,q_dst_ttl) ! Tracer is moist mass mixing ratio
#endif /* not DST_MSS_BDG */

#ifdef DST_CHM
c     Chemistry reaction rates
c     dst_chm_rxr() is called by CCM:physics/aphys(), MATCH:physlic()
      call dst_chm_rxr(obuf,
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     lat_idx,             ! I [idx] Model latitude index
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_cst(1,1,dst_idx_srt), ! I [kg kg-1] Dust mixing ratio
     $     rxrc_chm(1,1,idx_HNO3_gas), ! O [s-1] Pseudo first order rate coefficient for HNO3
     $     rxrc_chm(1,1,idx_HO2_gas), ! O [s-1] Pseudo first order rate coefficient for HO2
     $     rxrc_chm(1,1,idx_N2O5_gas), ! O [s-1] Pseudo first order rate coefficient for N2O5
     $     rxrc_chm(1,1,idx_O3_gas), ! O [s-1] Pseudo first order rate coefficient for O3
     $     rxrc_chm(1,1,idx_SO2_gas), ! O [s-1] Pseudo first order rate coefficient for SO2
     $     tpt_mdp,             ! I [K] Temperature
     $     tpt_vrt,             ! I [K] Virtual temperature
     $     vmr_chm(1,1,idx_SO4_aer)) ! I [mlc mlc-1] Particulate SO4 volume mixing ratio

c     Chemistry solver and mixing ratio adjustment
c     dst_chm_slv() is a standin for MOZART's chemistry solver
      call dst_chm_slv(obuf,
     $     lat_idx,             ! I [idx] Model latitude index
     $     rxrc_chm(1,1,idx_HNO3_gas), ! I [s-1] Pseudo first order rate coefficient for HNO3
     $     rxrc_chm(1,1,idx_HO2_gas), ! I [s-1] Pseudo first order rate coefficient for HO2
     $     rxrc_chm(1,1,idx_N2O5_gas), ! I [s-1] Pseudo first order rate coefficient for N2O5
     $     rxrc_chm(1,1,idx_O3_gas), ! I [s-1] Pseudo first order rate coefficient for O3
     $     rxrc_chm(1,1,idx_SO2_gas), ! I [s-1] Pseudo first order rate coefficient for SO2
     $     tm_adj,              ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     vmr_chm(1,1,idx_HNO3_gas), ! I/O [mlc mlc-1] Gaseous HNO3 volume mixing ratio
     $     vmr_chm(1,1,idx_HO2_gas), ! I/O [mlc mlc-1] Gaseous HO2 volume mixing ratio
     $     vmr_chm(1,1,idx_N2O5_gas), ! I/O [mlc mlc-1] Gaseous N2O5 volume mixing ratio
     $     vmr_chm(1,1,idx_NO3_aer), ! I/O [mlc mlc-1] Particulate NO3 volume mixing ratio
     $     vmr_chm(1,1,idx_O3_gas), ! I/O [mlc mlc-1] Gaseous O3 volume mixing ratio
     $     vmr_chm(1,1,idx_SO2_gas), ! I/O [mlc mlc-1] Gaseous SO2 volume mixing ratio
     $     vmr_chm(1,1,idx_SO4_aer)) ! I/O [mlc mlc-1] Particulate SO4 volume mixing ratio
#endif /* not DST_CHM */
      
c     Precompute dust paths needed by shortwave, longwave, and outfld() routines
c     dst_pth() is called by CCM:physics/tphysbc() before radctl()
      call dst_pth(
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     q_cst(1,1,dst_idx_srt), ! I [kg kg-1] Dust mixing ratio
     $     mpc_dst,             ! O [kg m-2] Column mass path of dust
     $     mpc_dst_ttl,         ! O [kg m-2] Total column mass path of dust
     $     mpl_dst,             ! O [kg m-2] Layer dust amount
     $     mpp_dst,             ! O [kg m-2] Dust path above kth interface level
     $     odxc_dst,            ! O [frc] Column dust optical depth
     $     odxc_dst_ttl)        ! O [frc] Total column dust optical depth

c     Outfield column mass paths
c     call mmrd2mpc(sh1,q_cst(1,1,dst_idx_srt),prs_dlt,mpc_dst,mpc_dst_ttl) ! MATCH:physlic()
      call outfld('DSTMPC  ',mpc_dst_ttl,plond,lat_idx,obuf)
      do m=1,dst_nbr
         call outfld(mpc_dst_nm(m),mpc_dst(1,m),plond,lat_idx,obuf)
      end do                    ! end loop over cst

c     Diagnostic optical depths
c     call mpc2odxc(mpc_dst,odxc_dst,odxc_dst_ttl) ! MATCH:physlic()
      call outfld('DSTODXC',odxc_dst_ttl,plond,lat_idx,obuf)
      do m=1,dst_nbr
         call outfld(odxc_dst_nm(m),odxc_dst(1,m),plond,lat_idx,obuf)
      end do                    ! end loop over cst

c     Compute transmission through dust absorption continuum
c     dst_trn() is called by CCM:physics/radclw()
      call dst_trn(
     $     mpp_dst,             ! I [kg m-2] Dust path above kth interface level
     $     dst_trn_ttl)         ! O [frc] Total transmission through all size classes between interfaces k1 and k2

c     Update netCDF file
      fl_out='aer.nc' ! [sng] Name of netCDF output file
      call strnul(fl_out)
      call clm2nc(
     $     dst_trn_ttl,         ! I [frc] Total transmission through all size classes between interfaces k1 and k2
     $     fl_out,              ! I [sng] Name of netCDF output file
     $     lat_idx,             ! I [idx] Model latitude index
     $     mpc_dst,             ! I [kg m-2] Column mass path of dust
     $     mpc_dst_ttl,         ! I [kg m-2] Total column mass path of dust
     $     odxc_dst,            ! I [frc] Column dust optical depth
     $     odxc_dst_ttl,        ! I [frc] Total column dust optical depth
     $     q_cst(1,1,dst_idx_srt), ! I [kg kg-1] Dust mixing ratio
     $     q_dst_ttl)           ! I [kg kg-1] Total dust mixing ratio

c     Convert advected tracers (except H2O vapor) from moist to dry mass mixing ratios
      call mmrm2mmrd(
     $     plev,                ! I [nbr] Number of levels
     $     plon,                ! I [nbr] Number of longitudes
     $     plond,               ! I [nbr] First dimension of arrays
     $     q_H2O_vpr,           ! I [kg H2O vapor kg-1 moist air] Specific humidity
     $     q_cst(1,1,dst_idx_srt), ! I [kg tracer kg-1 moist air] Tracer moist mass mixing ratio
     $     dst_nbr,             ! I [nbr] Number of tracers
     $     q_cst(1,1,dst_idx_srt)) ! O [kg tracer kg-1 dry air] Tracer dry mass mixing ratio

#ifdef DST_MSS_BDG
c     Budgets after aphys(), time filter, called in CCM:dynamics/eul/linemsbc()
      call dst_add_lon_lev(q_cst(1,1,dst_idx_srt),q_dst_ttl) ! q3m1
      call bdg_gam_dry('mpc_dst_a_flt',lat_idx,prs_dlt,q_cst,q_dst_ttl) ! q3m1
#endif /* not DST_MSS_BDG */

c     Outfield total dust mixing ratio 
c     Called by CCM:dynamics/eul/linemsbc(), MATCH: physlic()
      call dst_add_lon_lev(q_cst(1,1,dst_idx_srt),q_dst_ttl) ! q3m1
      call outfld('DSTQ    ',q_dst_ttl,plond,lat_idx,obuf)
#ifdef DST_MSS_BDG
      call bdg_gam_dry('mpc_dst_mdl',lat_idx,prs_dlt,q_cst,q_dst_ttl) ! q3m1
      call bdg_gam_wet('mpc_H2O',lat_idx,prs_dlt,q_cst) ! q3m1 NB: Only H2O uses bdg_gam_wet()
#endif /* not DST_MSS_BDG */

      return
      end                       ! end phyzlic()
      
      subroutine clm2nc(            
     $     dst_trn_ttl,         ! I [frc] Total transmission through all size classes between interfaces k1 and k2
     $     fl_out,              ! I [sng] Name of netCDF output file
     $     lat_idx,             ! I [idx] Model latitude index
     $     mpc_dst,             ! I [kg m-2] Column mass path of dust
     $     mpc_dst_ttl,         ! I [kg m-2] Total column mass path of dust
     $     odxc_dst,            ! I [frc] Column dust optical depth
     $     odxc_dst_ttl,        ! I [frc] Total column dust optical depth
     $     q_dst,               ! I [kg kg-1] Dust mixing ratio
     $     q_dst_ttl)           ! I [kg kg-1] Total dust mixing ratio
c     Purpose: Output aerosol mobilization properties to netCDF file
      use precision
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstctl.h> /* Model control variables */
#include <dbg.com> /* Debugging constants */
#include <netcdf.inc>
c     Input
      character*(*) fl_out      ! I [sng] Name of netCDF output file
      integer    lat_idx        ! I [idx] Model latitude index
      real(r8)::dst_trn_ttl(plond,plevp,plevp,bnd_nbr_LW) ! I [frc] Total transmission through all size classes between interfaces k1 and k2
      real(r8)::mpc_dst(plond,dst_nbr) ! I [kg m-2] Column mass path of dust
      real(r8)::mpc_dst_ttl(plond)   ! [kI g m-2] Total column mass path of dust
      real(r8)::odxc_dst(plond,dst_nbr) ! I [frc] Column dust optical depth
      real(r8)::odxc_dst_ttl(plond)  ! I [frc] Total column dust optical depth
      real(r8)::q_dst(plond,plev,dst_nbr) ! I [kg kg-1] Dust mixing ratio
      real(r8)::q_dst_ttl(plond,plev) ! I [kg kg-1] Total dust mixing ratio
c     Output
c     Local
c     File metadata and dimension IDs
      integer cnt_lon_lev_sz_time(4) ! Count array
      integer cnt_lon_lev_time(3) ! Count array
      integer cnt_lon_sz_time(3) ! Count array
      integer cnt_lon_time(2)   ! Count array
      integer dim_lon_lev_sz_time(4) ! [enm] Dimension IDs
      integer dim_lon_lev_time(3) ! [enm] Dimension IDs
      integer dim_lon_sz_time(3) ! [enm] Dimension IDs
      integer dim_lon_time(2)   ! [enm] Dimension IDs
      integer fll_mode_old      ! Old fill mode
      integer lev_dim_id        ! [enm] Dimension ID for lev
      integer lon_dim_id        ! [enm] Dimension ID for lon
      integer nc_id             ! File handle
      integer rcd               ! [rcd] Return success code
      integer srt_lon_lev_sz_time(4) ! Starting index array
      integer srt_lon_lev_time(3) ! Starting index array
      integer srt_lon_sz_time(3) ! Starting index array
      integer srt_lon_time(2)   ! Starting index array
      integer sz_dim_id         ! [enm] Dimension ID for sz
      integer sz_grd_dim_id     ! [enm] Dimension ID for sz grid
      integer time_dim_id       ! [enm] Dimension ID for time
c     Variable IDs
      integer mpc_dst_id        ! [enm] Variable ID
      integer mpc_dst_ttl_id    ! [enm] Variable ID
      integer odxc_dst_id       ! [enm] Variable ID
      integer odxc_dst_ttl_id   ! [enm] Variable ID
      integer q_dst_id          ! [enm] Variable ID
      integer q_dst_ttl_id      ! [enm] Variable ID
c     Variable data
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
      
c     Begin netCDF output routines
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf_open(fl_out,NF_WRITE,nc_id)
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR opening ',fl_out(1:f_strlen(fl_out)),
     $     ' in clm2nc(), rcd = ',rcd
      rcd=rcd+nf_redef(nc_id)
      rcd=rcd+nf_inq_dimid(nc_id,'lon',lon_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'sz',sz_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'lev',lev_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'time',time_dim_id)
c     Add global attributes
c     Define dimension IDs
c     Assemble ID and count vectors for each multidimensional combination of dimensions
      dim_lon_time=(/lon_dim_id,time_dim_id/)
      cnt_lon_time=(/plon,1/)
      srt_lon_time=(/1,nstep/)

      dim_lon_sz_time=(/lon_dim_id,sz_dim_id,time_dim_id/)
      cnt_lon_sz_time=(/plon,dst_nbr,1/)
      srt_lon_sz_time=(/1,1,nstep/)

      dim_lon_lev_time=(/lon_dim_id,lev_dim_id,time_dim_id/)
      cnt_lon_lev_time=(/plon,plev,1/)
      srt_lon_lev_time=(/1,1,nstep/)

      dim_lon_lev_sz_time=(/lon_dim_id,lev_dim_id,sz_dim_id,time_dim_id/)
      cnt_lon_lev_sz_time=(/plon,plev,dst_nbr,1/)
      srt_lon_lev_sz_time=(/1,1,1,nstep/)

      if (nstep == 1) then
c     Variable definitions
         rcd=rcd+nf_def_var(nc_id,'mpc_dst',NF_DOUBLE,3,dim_lon_sz_time,mpc_dst_id)
         rcd=rcd+nf_def_var(nc_id,'mpc_dst_ttl',NF_DOUBLE,2,dim_lon_time,mpc_dst_ttl_id)
         rcd=rcd+nf_def_var(nc_id,'odxc_dst',NF_DOUBLE,3,dim_lon_sz_time,odxc_dst_id)
         rcd=rcd+nf_def_var(nc_id,'odxc_dst_ttl',NF_DOUBLE,2,dim_lon_time,odxc_dst_ttl_id)
         rcd=rcd+nf_def_var(nc_id,'q_dst',NF_DOUBLE,4,dim_lon_lev_sz_time,q_dst_id)
         rcd=rcd+nf_def_var(nc_id,'q_dst_ttl',NF_DOUBLE,3,dim_lon_lev_time,q_dst_ttl_id)
c     Add english text descriptions
         rcd=rcd+nf_put_att_text(nc_id,mpc_dst_id,'long_name',24,'Column mass path of dust')
         rcd=rcd+nf_put_att_text(nc_id,mpc_dst_ttl_id,'long_name',30,'Total column mass path of dust')
         rcd=rcd+nf_put_att_text(nc_id,odxc_dst_id,'long_name',25,'Column dust optical depth')
         rcd=rcd+nf_put_att_text(nc_id,odxc_dst_ttl_id,'long_name',31,'Total column dust optical depth')
         rcd=rcd+nf_put_att_text(nc_id,q_dst_id,'long_name',17,'Dust mixing ratio')
         rcd=rcd+nf_put_att_text(nc_id,q_dst_ttl_id,'long_name',23,'Total dust mixing ratio')
c     Add units
         rcd=rcd+nf_put_att_text(nc_id,mpc_dst_id,'units',16,'kilogram meter-2')
         rcd=rcd+nf_put_att_text(nc_id,mpc_dst_ttl_id,'units',16,'kilogram meter-2')
         rcd=rcd+nf_put_att_text(nc_id,odxc_dst_id,'units',8,'fraction')
         rcd=rcd+nf_put_att_text(nc_id,odxc_dst_ttl_id,'units',8,'fraction')
         rcd=rcd+nf_put_att_text(nc_id,q_dst_id,'units',19,'kilogram kilogram-1')
         rcd=rcd+nf_put_att_text(nc_id,q_dst_ttl_id,'units',19,'kilogram kilogram-1')
      else                      ! endif nstep == 1
c     Get variable IDs
         rcd=rcd+nf_inq_varid(nc_id,'mpc_dst',mpc_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'mpc_dst_ttl',mpc_dst_ttl_id)
         rcd=rcd+nf_inq_varid(nc_id,'odxc_dst',odxc_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'odxc_dst_ttl',odxc_dst_ttl_id)
         rcd=rcd+nf_inq_varid(nc_id,'q_dst',q_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'q_dst_ttl',q_dst_ttl_id)
      endif                     ! endif nstep /= 1
c     End define mode
      rcd=rcd+nf_enddef(nc_id)
c     Write data
      rcd=rcd+nf_put_vara_double(nc_id,mpc_dst_id,srt_lon_sz_time,cnt_lon_sz_time,mpc_dst)
      rcd=rcd+nf_put_vara_double(nc_id,mpc_dst_ttl_id,srt_lon_time,cnt_lon_time,mpc_dst_ttl)
      rcd=rcd+nf_put_vara_double(nc_id,odxc_dst_id,srt_lon_sz_time,cnt_lon_sz_time,odxc_dst)
      rcd=rcd+nf_put_vara_double(nc_id,odxc_dst_ttl_id,srt_lon_time,cnt_lon_time,odxc_dst_ttl)
      rcd=rcd+nf_put_vara_double(nc_id,q_dst_id,srt_lon_lev_sz_time,cnt_lon_lev_sz_time,q_dst)
      rcd=rcd+nf_put_vara_double(nc_id,q_dst_ttl_id,srt_lon_lev_time,cnt_lon_lev_time,q_dst_ttl)
c     Close output file
      rcd=rcd+nf_close(nc_id)
      if (nstep == 1) then
         write (6,'(a,a51,a)') prg_nm(1:f_strlen(prg_nm)),': Initialized column mass and optical path data in ',
     $     fl_out(1:f_strlen(fl_out))
      endif                     ! endif
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR writing ',fl_out(1:f_strlen(fl_out)),
     $     ' in clm2nc(), rcd = ',rcd
      return 
      end                       ! end clm2nc()
