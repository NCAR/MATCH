c     $Header$ -*-fortran-*- 
      
c     Purpose: aernvr.h contains common blocks which store environmental properties
c     needed for aerosol simulations
      
c     These variables are initialized in aer_nvr_cmn_ini()
c     aer_nvr_cmn_ini() is called by aer()
c     aernvr.h MUST have access to pmgrid.h to work
      
c     Usage:
c     #include <aernvr.h> /* Aerosol environmental properties */
      
c     aer_nvr_cmn contains environmental properties needed for simulations, initialized in aer_nvr_cmn_ini()
      
c     Coordinate grid
      real(r8) lat(plat)            ! [dgr] Latitude
      real(r8) lev(plev)            ! [dgr] Midlayer pressure
      real(r8) levp(plevp)          ! [dgr] Interface pressure
      real(r8) lon(plon)            ! [dgr] Longitude
      common /aer_nvr_grd_cmn/  ! aer_nvr_grd_cmn is initialized in aer_nvr_cmn_ini()
     $     lat,                 ! [rdn] Latitude
     $     lev,                 ! [dgr] Midlayer pressure
     $     levp,                ! [dgr] Interface pressure
     $     lon                  ! [dgr] Longitude
      
      integer sfc_typ(plond)    ! [idx] LSM surface type (0..28)
      common /aer_nvr_cmn_int/  ! aer_nvr_cmn is initialized in aer_nvr_cmn_ini()
     $     sfc_typ              ! [idx] LSM surface type (0..28)
      logical wrhstts           ! Flag to output buffer to file this timestep
      common /aer_nvr_cmn_lgc/  ! aer_nvr_cmn is initialized in aer_nvr_cmn_ini()
     $     wrhstts              ! Flag to output buffer to file this timestep
      
      real(r8) cld_frc(plond,plevp)  ! [frc] Cloud fraction
      real(r8) cld_frc_cnv(plond,plev) ! [frc] Convective cloud fraction
      real(r8) dns_mdp(plond,plev)  ! [kg m-3] Midlayer density
      real(r8) doy                  ! [day] Day of year [1.0..366.0)
      real(r8) flx_LW_dwn_sfc(plond) ! [W m-2] Longwave downwelling flux at surface
      real(r8) flx_SW_abs_sfc(plond) ! [W m-2] Solar flux absorbed by ground
      real(r8) hbuf(1)              ! History buffer
      real(r8) hgt_mdp(plond,plev)  ! [m] Midlayer height above surface
      real(r8) hgt_ntf(plond,plevp) ! [m] Interface height above surface
      real(r8) lat_dgr(plat)        ! [dgr] Latitude
      real(r8) lat_rdn(plat)        ! [rdn] Latitude
      real(r8) lnd_frc_dry(plond)   ! [frc] Dry land fraction
      real(r8) mbl_bsn_fct(plond)   ! [frc] Topographic source limitation
      real(r8) mno_lng(plond)       ! [m] Monin-Obukhov length
      real(r8) mpl_air(plond,plev)  ! [kg m-2] Air mass path in layer
      real(r8) mss_cnc_dst(plond,plev)  ! [kg m-3] Mass concentration of dust
      real(r8) mss_frc_CaCO3(plond) ! [frc] Mass fraction CaCO3 
      real(r8) mss_frc_cly(plond)   ! [frc] Mass fraction clay 
      real(r8) mss_frc_slt(plond)   ! [frc] Mass fraction silt 
      real(r8) mss_frc_snd(plond)   ! [frc] Mass fraction sand 
      real(r8) obuf(1)              ! Output buffer
      real(r8) oro(plond)           ! [frc] Orography
      real(r8) prs_dlt(plond,plev)  ! [Pa] Pressure thickness
      real(r8) prs_mdp(plond,plev)  ! [Pa] Midlayer pressure 
      real(r8) prs_ntf(plond,plevp) ! [Pa] Interface pressure
      real(r8) q_H2O_cnd(plond,plev) ! [kg kg-1] Condensed H2O mixing ratio
      real(r8) q_H2O_cnd2pcp_tnd(plond,plev) ! [kg kg-1 s-1] Condensed H2O to precipitation tendency
      real(r8) q_H2O_cnd_pcp(plond,plev) ! [kg kg-1] H2O precipitation mixing ratio
      real(r8) q_H2O_cnd_tnd(plond,plev) ! [kg kg-1 s-1] Net H2O condensate formation tendency
      real(r8) q_H2O_pcp_lqd(plond,plev) ! [kg kg-1] Rain water mixing ratio
      real(r8) q_H2O_pcp2vpr_tnd(plond,plev) ! [kg kg-1 s-1] H2O precipitation to vapor tendency
      real(r8) q_H2O_vpr(plond,plev) ! [kg kg-1] Water vapor mixing ratio
      real(r8) q_H2O_vpr2pcp_cnv_tnd(plond,plev) ! [kg kg-1 s-1] H2O vapor to convective precipitation tendency
      real(r8) q_cst(plond,plev,pcnst) ! [kg kg-1] Full constituent array
      real(r8) q_dst(plond,plev,dst_nbr) ! [kg kg-1] Dust mixing ratio
      real(r8) rxrc_chm(plond,plev,chm_nbr) ! [s-1] Pseudo first order rate coefficients
      real(r8) snw_hgt_lqd(plond)   ! [m] Equivalent liquid water snow depth
      real(r8) tm_adj               ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
      real(r8) tpt_gnd(plond)       ! [K] Ground temperature
      real(r8) tpt_ice(plond)       ! [K] Ice temperature
      real(r8) tpt_mdp(plond,plev)  ! [K] Temperature
      real(r8) tpt_ptn(plond,plev)  ! [K] Potential temperature
      real(r8) tpt_sfc(plond)       ! [K] Surface temperature
      real(r8) tpt_soi(plond)       ! [K] Soil temperature
      real(r8) tpt_sst(plond)       ! [K] Sea surface temperature
      real(r8) tpt_vrt(plond,plev)  ! [K] Virtual temperature
      real(r8) vmr_HNO3_gas(plond,plev) ! [mlc mlc-1] Gaseous HNO3 volume mixing ratio
      real(r8) vmr_NO3_aer(plond,plev) ! [mlc mlc-1] Particulate NO3 volume mixing ratio
      real(r8) vmr_SO4_aer(plond,plev) ! [mlc mlc-1] Particulate SO4 volume mixing ratio
      real(r8) vmr_chm(plond,plev,chm_nbr) ! [mlc mlc-1] Chemical species volume mixing ratios
      real(r8) vwc_sfc(plond)       ! [m3 m-3] Volumetric water content
      real(r8) wnd_frc(plond)       ! [m s-1] Friction velocity
      real(r8) wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
      real(r8) wnd_mrd_mdp(plond)   ! [m s-1] Meridional wind component
      real(r8) wnd_znl_mdp(plond)   ! [m s-1] Zonal wind component
      common /aer_nvr_cmn_real/ ! aer_nvr_cmn is initialized in aer_nvr_cmn_ini()
     $     cld_frc,             ! [frc] Cloud fraction
     $     cld_frc_cnv,         ! [frc] Convective cloud fraction
     $     dns_mdp,             ! [kg m-3] Midlayer density
     $     doy,                 ! [day] Day of year [1.0..366.0)
     $     flx_LW_dwn_sfc,      ! [W m-2] Longwave downwelling flux at surface
     $     flx_SW_abs_sfc,      ! [W m-2] Solar flux absorbed by ground
     $     hbuf,                ! History buffer
     $     hgt_mdp,             ! [m] Midlayer height above surface
     $     hgt_ntf,             ! [m] Interface height above surface
     $     lat_dgr,             ! [dgr] Latitude
     $     lat_rdn,             ! [rdn] Latitude
     $     lnd_frc_dry,         ! [frc] Dry land fraction
     $     mbl_bsn_fct,         ! [frc] Topographic source limitation
     $     mno_lng,             ! [m] Monin-Obukhov length
     $     mpl_air,             ! [kg m-2] Air mass path in layer
     $     mss_cnc_dst,         ! [kg m-3] Mass concentration of dust
     $     mss_frc_CaCO3,       ! [frc] Mass fraction CaCO3 
     $     mss_frc_cly,         ! [frc] Mass fraction clay 
     $     mss_frc_slt,         ! [frc] Mass fraction silt 
     $     mss_frc_snd,         ! [frc] Mass fraction sand 
     $     obuf,                ! Output buffer
     $     oro,                 ! [frc] Orography
     $     prs_dlt,             ! [Pa] Pressure thickness
     $     prs_mdp,             ! [Pa] Midlayer pressure 
     $     prs_ntf,             ! [Pa] Interface pressure
     $     q_H2O_cnd,           ! [kg kg-1] Condensed H2O mixing ratio
     $     q_H2O_cnd2pcp_tnd,   ! [kg kg-1 s-1] Condensed H2O to precipitation tendency
     $     q_H2O_cnd_pcp,       ! [kg kg-1] H2O precipitation mixing ratio
     $     q_H2O_cnd_tnd,       ! [kg kg-1 s-1] Net H2O condensate formation tendency
     $     q_H2O_pcp2vpr_tnd,   ! [kg kg-1 s-1] H2O precipitation to vapor tendency
     $     q_H2O_pcp_lqd,       ! [kg kg-1] Rain water mixing ratio
     $     q_H2O_vpr,           ! [kg kg-1] Water vapor mixing ratio
     $     q_H2O_vpr2pcp_cnv_tnd, ! [kg kg-1 s-1] H2O vapor to convective precipitation tendency
     $     q_cst,               ! [kg kg-1] Full constituent array
     $     q_dst,               ! [kg kg-1] Dust mixing ratio
     $     rxrc_chm,            ! [s-1] Pseudo first order rate coefficients
     $     snw_hgt_lqd,         ! [m] Equivalent liquid water snow depth
     $     tm_adj,              ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     tpt_gnd,             ! [K] Ground temperature
     $     tpt_ice,             ! [K] Ice temperature
     $     tpt_mdp,             ! [K] Temperature
     $     tpt_ptn,             ! [K] Potential temperature
     $     tpt_sfc,             ! [K] Surface temperature
     $     tpt_soi,             ! [K] Soil temperature
     $     tpt_sst,             ! [K] Sea surface temperature
     $     tpt_vrt,             ! [K] Virtual temperature
     $     vmr_HNO3_gas,        ! [mlc mlc-1] Gaseous HNO3 volume mixing ratio
     $     vmr_NO3_aer,         ! [mlc mlc-1] Particulate NO3 volume mixing ratio
     $     vmr_SO4_aer,         ! [mlc mlc-1] Particulate SO4 volume mixing ratio
     $     vmr_chm,             ! [mlc mlc-1] Chemical species volume mixing ratios
     $     vwc_sfc,             ! [m3 m-3] Volumetric water content
     $     wnd_frc,             ! [m s-1] Friction velocity
     $     wnd_mdp,             ! [m s-1] Surface layer mean wind speed
     $     wnd_mrd_mdp,         ! [m s-1] Meridional wind component
     $     wnd_znl_mdp          ! [m s-1] Zonal wind component
