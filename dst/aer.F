c     $Header$
      
#ifdef BXM
#include <params.h>
#include <dst.h> /* Dust preprocessor tokens */
      
      program aer
c     Purpose: Offline driver for aerosol microphysics
c     Usage:
c     aer -D 1   ! Turns on debugging output
c     aer -n 10  ! Runs for 10 timesteps
#ifndef CRAY
c     etags *.F *.h ~/c++/*.cc ~/c++/*.hh ~/mie/*.cc ~/mie/*.hh ~/map/*.F ~/map/*.h ~/map/*.com
c     */
#endif
c     Compilation:
c     cd ~/aer;make OPTS=D aer;cd -
c     cd ~/aer;make OPTS=D precision=double aer;cd -
      use precision
      use dstszdst
      implicit none
c     Parameters
      character*(*),parameter::CVS_Header='$Header$'
c     Physics parameters needed for call to dst_cst_cmn_ini()
      real(r8),parameter::gas_cst_unv=8.31441 ! [J mol-1 K-1] Universal gas constant
      real(r8),parameter::mmw_dry_air=28.9644e-3 ! [kg mol-1] (Source: radcsw.F in CCM2/3)
      real(r8),parameter::gas_cst_dry_air=gas_cst_unv/mmw_dry_air ! (287.05) [J kg-1 K-1] IrG81 p. 25, p. 245
      real(r8),parameter::grv_mean_sfc=9.80665 ! [m s-2] Mean gravitational acceleration at Earth's surace
      real(r8),parameter::rds_earth=6.370e+06 ! [m] Radius of sphere of same volume as Earth
c     Commons
#include <dbg.com> /* Debugging constants */
#include <xtr.com> /* Extrapolation constants */
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstaer.h> /* Aerosol microphysical properties */
#include <aernvr.h> /* Aerosol environmental properties */
#include <dstctl.h> /* Model control variables */
c     Input
c     Input/Output
c     Output
c     Local
      character arg_val*80      ! Command line argument value
      character cmd_ln*200      ! Command line
      character::dsh_key*2      ! Command line dash and switch
      character*26::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
      character prg_ID*200      ! Program ID
      character::CVS_Date*28='$Date$' ! Date of CVS checkin
      character::CVS_Id*80='$Id$' ! CVS ID
      character::CVS_Revision*80='$Revision$' ! CVS Revision
      character fl_out*80       ! [sng] Name of netCDF output file
      character pth_out*80      ! [sng] Name of netCDF output path
      integer arg_nbr           ! [nbr] Number of command line arguments
      integer arg_idx           ! [idx] Counting index
      integer::exit_status=0    ! Program exit status
      integer lat_idx           ! [idx] Counting index
      integer::time_nbr=1       ! [m] Number of timesteps to simulate
      real(r8) lat_wgt(plat)        ! [frc] Latitude weights (currently must sum to 2.0)
c     External
      integer iargc             ! [nbr] Command line argument count
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
c     External
      external dstlsmbd         ! Initialize LSM block data
      external dstszdstbd       ! Initialize aerosol block data
      external dstradbd         ! Initialize radiation block data
c     Main code
c     Set defaults
c     Set defaults for command line arguments
      dbg_lvl=dbg_off           ! dbg_lvl allocated in dbg.com, Option D
c     tm_adj must be set before bdg_cmn_ini() is called
      tm_adj=2400.0             ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
c     Derived fields
      
c     Retrieve command line arguments
      call date_time_get(lcl_date_time)
      call cmd_ln_sng(cmd_ln)
      call prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
      write (6,'(a)') prg_ID(1:f_strlen(prg_ID))
      arg_nbr=iargc()
      arg_idx=1
      do while (arg_idx <= arg_nbr)
         call getarg(arg_idx,arg_val)
c     Increment counter for next read
         arg_idx=arg_idx+1
         if (dbg_lvl >= dbg_io) write (6,'(2a,i2,2a)')prg_nm(1:f_strlen(prg_nm)),': DEBUG arg_idx = ',arg_idx,',arg_val = ',arg_val
         dsh_key=arg_val(1:2)
         if (dsh_key == '-D') then
            call getarg(arg_idx,arg_val)
            arg_idx=arg_idx+1
            read (arg_val,'(i4)') dbg_lvl
         else if (dsh_key == '-n') then
            call getarg(arg_idx,arg_val)
            arg_idx=arg_idx+1
            read (arg_val,'(i4)') time_nbr
         else if (dsh_key == '-d') then
            call getarg(arg_idx,arg_val)
            arg_idx=arg_idx+1
            read (arg_val,'(f9.3)') tm_adj
         else if (dsh_key == '-v') then
            goto 1000
         else
c     Option not recognized
            write (6,'(2a)') prg_nm(1:f_strlen(prg_nm)),': ERROR Unrecognized switch'
            exit_status=-1
            goto 1000
         endif                  ! endif arg_val
      end do                    ! end while (arg_idx <= arg_nbr)
      
c     Construct dummy Gaussian weights
      do lat_idx=1,plat
         lat_wgt(lat_idx)=2.0/plat ! [frc] Latitude weights (currently must sum to 2.0)
      end do                    ! end loop over lat

c     Initialize time-invariant physical constants
c     dst_cst_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
      call dst_cst_cmn_ini(gas_cst_dry_air,grv_mean_sfc,rds_Earth)
c     Initialize dust names
c     dst_nm_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
      call dst_nm_cmn_ini()
c     Initialize size grid
c     dst_szdst_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
      call dst_szdst_ini()
c     Initialize miscellaneous common blocks
c     dst_msc_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
      call dst_msc_cmn_ini()
#ifdef DST_MSS_BDG
c     Initialize mass budget common block
c     bdg_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
      call bdg_cmn_ini(grv_mean_sfc,lat_wgt,tm_adj)
#endif /* not DST_MSS_BDG */
      
c     Initialize aerosol environmental properties common block aer_nvr_cmn
      call aer_nvr_cmn_ini()    ! Not needed in global models
c     Set surface common block
      call dst_sfc_set(plond,   ! Not needed in global models
     $     flx_LW_dwn_sfc,      ! I [W m-2] Longwave downwelling flux at surface
     $     flx_SW_abs_sfc,      ! I [W m-2] Solar flux absorbed by ground
     $     lnd_frc_dry,         ! I [frc] Dry land fraction
     $     mbl_bsn_fct,         ! I [frc] Topographic source limitation
     $     mss_frc_CaCO3,       ! I [frc] Mass fraction of CaCO3
     $     mss_frc_cly,         ! I [frc] Mass fraction of clay
     $     mss_frc_snd,         ! I [frc] Mass fraction of sand
     $     sfc_typ,             ! I [idx] LSM surface type (0..28)
     $     tpt_gnd,             ! I [K] Ground temperature
     $     tpt_soi,             ! I [K] Soil temperature
     $     vwc_sfc)             ! I [m3 m-3] Volumetric water content
      
c     Run physics on each timestep on each latitude slice
      do nstep=1,time_nbr
         if (dbg_lvl > dbg_off) write (6,'(a1,$)') '.'
c     Simplest possible OpenMP parallelization
!$omp parallel do
         do lat_idx=1,plat
            call phyzlic(lat_idx)
         end do                 ! end loop over lat
!$omp end parallel do
#ifdef DST_MSS_BDG
c     Output global average mass budget for this timestep
c     bdg_update() is called by CCM:dynamics/eul/stepon(), MATCH:main()
         call bdg_update( nstep, mcdate, mcsec )
#endif /* not DST_MSS_BDG */

c     Output time-varying variables
         fl_out='aer.nc'        ! [sng] Name of netCDF output file
         call tm2nc(            ! Not needed in global models
     $        fl_out,           ! I [sng] Name of netCDF output file
     $        tm_adj)           ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
      end do                    ! end loop over time
      
#ifdef DST_MSS_BDG
c     bdg_close() is called by CCM:dynamics/eul/stepon(), MATCH:main()
      call bdg_close()
#endif /* not DST_MSS_BDG */
      
c     Add environmental variables to file
      fl_out='aer.nc'           ! [sng] Name of netCDF output file
      call nvr2nc(fl_out)       ! Not needed in global models
      
      if (dbg_lvl == dbg_fl) then
         write (6,'(a)') 'ncks -C -F -u -H aer.nc | m'
         write (6,'(a)') 'ncks -C -F -u -H -v nstep,mcdate,mcsec,nbdate,nbsec,ndcur,nscur,time,tm_adj aer.nc'
         write (6,'(a)') 'ncks -C -H -F -u -v dst_slt_flx_rat_ttl,flx_mss_hrz_slt_ttl,'//
     $        'flx_mss_vrt_dst_ttl,lnd_frc_mbl,lnd_frc_dry,mbl_bsn_fct,sfc_typ,lat_dgr,doy aer.nc'
         write (6,'(a)') 'ncks -C -H -F -u -v mno_lng_mbl,mno_lng_dps,'//
     $        'rgh_mmn_mbl,rgh_mmn_dps,hgt_zpd_mbl,hgt_zpd_dps,'//
     $        'wnd_frc_mbl,wnd_frc_slt,wnd_frc_dps,'//
     $        'wnd_frc_thr_slt,wnd_rfr_mbl,wnd_rfr_dps aer.nc'
         write (6,'(a)') 'ncks -C -H -F -u -v mpc_dst,mpc_dst_ttl,odxc_dst,odxc_dst_ttl,q_dst,q_dst_ttl,mss_cnc_dst aer.nc'
         write (6,'(a)') 'ncks -C -H -F -u -d sz,2 -v flx_mss_dry_sfc,flx_mss_dry_sfc_ttl,'//
     $        'flx_mss_grv_sfc,flx_mss_grv_sfc_ttl,flx_mss_trb_sfc,flx_mss_trb_sfc_ttl,'//
     $        'q_dst_tnd_dry,q_dst_tnd_dry_ttl aer.nc'
         write (6,'(a)') 'ncks -C -H -F -u -d sz,2 -v flx_mss_pcp_sfc,flx_mss_pcp_sfc_ttl,'//
     $        'spc_xsx_ncl_scv,spc_xsx_cll_scv,frc_trc_ist,pcp_flx_sfc,'//
     $        'q_dst_tnd_ncl,q_dst_tnd_wet,q_dst_tnd_evp,q_dst_tnd_evp_ttl,q_dst_tnd_wet,q_dst_tnd_pcp_ttl aer.nc'
         write (6,'(a)') 'ncks -C -H -F -u -d sz,2 -v cld_frc,cld_frc_cnv,'//
     $        'cld_vlm,q_H2O_cnd,q_H2O_cnd2pcp_tnd,q_H2O_cnd_pcp,'//
     $        'q_H2O_cnd_tnd,q_H2O_pcp2vpr_tnd,q_H2O_vpr2pcp_cnv_tnd aer.nc'
         write (6,'(a)') 
     $        'ncks -C -H -F -u -d sz,2 -v mss_frc_cly,flx_mss_hrz_slt_ttl,'//
     $        'flx_mss_vrt_dst_ttl,dst_slt_flx_rat_ttl,ovr_src_snk_mss,'//
     $        'ovr_src_snk_frc,ovr_src_snk_mss_ttl,flx_mss_vrt_dst aer.nc'
         write (6,'(a)') 
     $        'ncks -C -H -F -u -v dns_mdp,mss_frc_cly,mss_frc_snd,'//
     $        'wnd_frc_mbl,wnd_frc_slt,wnd_frc_thr_slt,wnd_rfr_mbl,wnd_mdp,'//
     $        'wnd_rfr_thr_slt,frc_thr_rgh_fct,frc_thr_wet_fct,'//
     $        'rgh_mmn_mbl,snw_frc,snw_hgt_lqd,vwc_sfc,hgt_zpd_mbl,hgt_mdp aer.nc'
         write (6,'(a)') 
     $        'ncks -C -H -F -u -v dns_mdp,mno_lng_mbl,prs_mdp,prs_ntf,q_H2O_vpr,'//
     $        'sfc_typ,tpt_mdp,tpt_vrt,wnd_frc_mbl,wnd_frc_slt,wnd_mdp,hgt_mdp,'//
     $        'rgh_mmn_mbl,hgt_zpd_mbl,wnd_rfr_mbl aer.nc'
         write (6,'(a)') 
     $        'ncks -C -H -F -u -v rss_lmn,shm_nbr,stk_nbr,vlc_grv,vlc_trb aer.nc'
         write (6,'(a)') 
     $        'ncks -C -H -F -u -d sz,2 -v nbr_spc_rsl,sfc_spc_rsl,vlm_spc_rsl aer.nc'
         write (6,'(a)') 
     $        'ncks -C -H -F -u -d sz,2 -v ovr_src_snk_mss aer.nc'
         write (6,'(a)') 'ncks -C -H -F -u -d sz,2 -v dmt_nwr,dmt_swr,dmt_vwr,'//
     $        'dmt_nmr,dmt_smr,dmt_vmr aer.nc'
         write (6,'(a)') 'ncks -C -H -F -u -v '//
     $        'rxrc_HNO3_dst,upt_cff_HNO3_dst,vmr_HNO3_gas,'//
     $        'rxrc_N2O5_dst,upt_cff_N2O5_dst,vmr_N2O5_gas,'//
     $        'rxrc_SO2_dst,upt_cff_SO2_dst,vmr_SO2_gas aer.nc'
         write (6,'(a)') 'ncks -C -H -F -u -v tau_dps,tau_mbl,tau_pcp,tau_dry,'//
     $        'tm_ttl,mpc_dst_mdl_avg,'//
     $        'tau_grv,tau_trb dst_mss_bdg.nc'
         write (6,'(a)') 'ncks -C -H -F -u -v mpc_err,mpc_err_frc,'//
     $        'mpc_err_tnd_frc dst_mss_bdg.nc'
      endif                     ! endif dbg
 1000 continue
      call exit(exit_status)
      end                       ! end aer()
      
      subroutine aer_nvr_cmn_ini()
c     Initialize aerosol environmental properties common block aer_nvr_cmn
c     aer_nvr_cmn_ini() is called by aer
      use precision
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <aernvr.h> /* Aerosol environmental properties */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstblm.h> /* Physical constants for boundary layer meteorology */
#include <dstchm.h> /* Chemical properties of dust */
c     Parameters
c     Local
      integer i                 ! [idx] Counting index for longitude
      integer j                 ! [idx] Counting index for latitude
      integer k                 ! [idx] Counting index for level
      integer m                 ! [idx] Counting index for constituent
      double precision pi       ! [frc] 3
c     Main code
c     Initialize defaults
      pi=4.0*atan(1.0d0)        ! [frc] 3
      
c     Set environmental pressure and temperature for settling speed
c     fxm: 19990901 coordinate doy with nbdate,ndcur, etc.
      doy=135.0                 ! [day] Day of year [1.0..366.0)
      do j=1,plat
         lat_dgr(j)=40.0        ! [dgr]
         lat(j)=lat_dgr(j)      ! [dgr]
         lat_rdn(j)=pi*lat_dgr(j)/180.0 ! [rdn]
      end do                    ! end loop over lat
      do i=1,plon
         lon(i)=15.0*(i-1)      ! [dgr]
      end do                    ! end loop over lon
      do k=1,plevp              ! Loop ends at plevp
         levp(k)=101325.0-1000.0*(plevp-k) ! [Pa]
      end do                    ! end loop over lev
      do k=1,plev
         lev(k)=0.5*(levp(k)+levp(k+1)) ! [Pa]
      end do                    ! end loop over lev
      
      do k=1,plev
         do i=1,plon
            q_H2O_vpr(i,k)=1.76721e-2 ! [kg kg-1] Water vapor mixing ratio
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      
c     Variables needed for wet deposition parameterization
      do k=1,plev
         do i=1,plon
            cld_frc(i,k)=1.0    ! [frc] Cloud fraction
            cld_frc_cnv(i,k)=0.0 ! [frc] Convective cloud fraction
            q_H2O_cnd(i,k)=0.0  ! [kg kg-1] Condensed H2O mixing ratio
            q_H2O_cnd2pcp_tnd(i,k)=1.0e-9 ! [kg kg-1 s-1] Condensed H2O to precipitation tendency
            q_H2O_cnd_pcp(i,k)=0.0 ! [kg kg-1] H2O precipitation mixing ratio
            q_H2O_cnd_tnd(i,k)=0.0 ! [kg kg-1 s-1] Net H2O condensate formation tendency
            q_H2O_pcp2vpr_tnd(i,k)=1.0e-9 ! [kg kg-1 s-1] H2O precipitation to vapor tendency
            q_H2O_vpr2pcp_cnv_tnd(i,k)=1.0e-9 ! [kg kg-1 s-1] H2O vapor to convective precipitation tendency
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      
      do i=1,plon
         flx_LW_dwn_sfc(i)=350.0 ! [W m-2] Longwave downwelling flux at surface
         flx_SW_abs_sfc(i)=450.0 ! [W m-2] Solar flux absorbed by ground
         oro(i)=1.0             ! [frc] Orography
         sfc_typ(i)=2           ! [idx] LSM surface type (0..28)
         snw_hgt_lqd(i)=0.0     ! [m] Equivalent liquid water snow depth
         tpt_gnd(i)=300.0       ! [K] Ground temperature
         tpt_soi(i)=297.0       ! [K] Soil temperature
         tpt_ice(i)=tpt_frz_pnt ! [K] Ice temperature
         tpt_sst(i)=297.0       ! [K] Sea surface temperature
         wnd_znl_mdp(i)=10.0    ! [m s-1] Zonal wind component
         wnd_mrd_mdp(i)=0.0     ! [m s-1] Meridional wind component
      end do                    ! end loop over lon
      do i=1,plon
         if (nint(oro(i)) == 0) then
            lnd_frc_dry(i)=0.0  ! [frc] Dry land fraction
            mbl_bsn_fct(i)=0.0  ! [frc] Topographic source limitation
            mss_frc_CaCO3(i)=0.0 ! [frc] Mass fraction of CaCO3
            mss_frc_cly(i)=0.0  ! [frc] Mass fraction of clay 
            mss_frc_snd(i)=0.0  ! [frc] Mass fraction of sand
            tpt_sfc(i)=tpt_sst(i) ! [K] Surface temperature
            vwc_sfc(i)=1.0e30   ! [m3 m-3] Volumetric water content
         else if (nint(oro(i)) == 1) then
            lnd_frc_dry(i)=1.0  ! [frc] Dry land fraction
            mbl_bsn_fct(i)=1.0  ! [frc] Topographic source limitation
            mss_frc_CaCO3(i)=0.05 ! [frc] Mass fraction of CaCO3
            mss_frc_cly(i)=0.19 ! [frc] Mass fraction of clay 
            mss_frc_snd(i)=0.777 ! [frc] Mass fraction of sand
            tpt_sfc(i)=tpt_gnd(i) ! [K] Surface temperature
            vwc_sfc(i)=0.03     ! [m3 m-3] Volumetric water content
         else if (nint(oro(i)) == 2) then
            lnd_frc_dry(i)=0.0  ! [frc] Dry land fraction
            mbl_bsn_fct(i)=0.0  ! [frc] Topographic source limitation
            mss_frc_CaCO3(i)=0.0 ! [frc] Mass fraction of CaCO3
            mss_frc_cly(i)=0.0  ! [frc] Mass fraction of clay 
            mss_frc_snd(i)=0.0  ! [frc] Mass fraction of sand
            tpt_sfc(i)=tpt_ice(i) ! [K] Surface temperature
            vwc_sfc(i)=1.0      ! [m3 m-3] Volumetric water content
         else 
            stop
         endif
      end do                    ! end loop over lon
      
      do i=1,plon
         mss_frc_slt(i)=1.0-mss_frc_cly(i)-mss_frc_snd(i) ! [frc] Mass fraction silt
         wnd_mdp(i)=sqrt(wnd_znl_mdp(i)*wnd_znl_mdp(i)+wnd_mrd_mdp(i)*wnd_mrd_mdp(i)) ! [m s-1] Surface layer mean wind speed
      end do                    ! end loop over lon
      
      do k=1,plevp              ! Loop ends at plevp
         do i=1,plon
            hgt_ntf(i,k)=100.0-10.0*(k-1) ! [m] Interface height above surface
            prs_ntf(i,k)=levp(k) ! [Pa] Interface pressure
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      
      do k=1,plev
         do i=1,plon
            prs_mdp(i,k)=lev(k) ! [Pa] Midlayer pressure 
            prs_dlt(i,k)=prs_ntf(i,k+1)-prs_ntf(i,k) ! [Pa] Pressure thickness
            tpt_mdp(i,k)=300.0+3.0*(k-1) ! [K] Temperature
            tpt_ptn(i,k)=tpt_mdp(i,k)*(prs_ntf(i,plevp)/prs_mdp(i,k))**kappa_dry_air ! [K] Virtual temperature
            tpt_vrt(i,k)=tpt_mdp(i,k)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(i,k)) ! [K] Virtual temperature
            dns_mdp(i,k)=prs_mdp(i,k)/(tpt_vrt(i,k)*gas_cst_dry_air) ! [kg m-3]
            hgt_mdp(i,k)=0.5*(hgt_ntf(i,k)+hgt_ntf(i,k+1)) ! [m] Midlayer height above surface
            mpl_air(i,k)=prs_dlt(i,k)*grv_sfc_rcp ! [kg m-2] Air mass path in layer
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      
      do m=1,dst_nbr
         do k=1,plev
            do i=1,plon
c     Initializing aerosol mixing ratio to zero ensures diagnostics in dst_mss_bdg() file are synchronized
               q_dst(i,k,m)=0.0 ! [kg kg-1] Dust mixing ratio
c               q_dst(i,k,m)=1.0e-9 ! [kg kg-1] Dust mixing ratio
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst

      do k=1,plev
         do i=1,plon
c     Initialize gaseous and particulate VMRs with PEM-WEST-A data from ZhC99 p. 358 Tbl. 3
c     Initial value = 0.0 ppbv for short-lived species
c     Gas VMRs
            vmr_chm(i,k,idx_H2O2_gas)=0.0e-09 ! [mlc mlc-1] Gaseous H2O2 volume mixing ratio
            vmr_chm(i,k,idx_HNO3_gas)=0.05e-09 ! [mlc mlc-1] Gaseous HNO3 volume mixing ratio
            vmr_chm(i,k,idx_HO2_gas)=0.0e-09 ! [mlc mlc-1] Gaseous HO2 volume mixing ratio
            vmr_chm(i,k,idx_N2O5_gas)=0.0e-09 ! [mlc mlc-1] Gaseous N2O5 volume mixing ratio
            vmr_chm(i,k,idx_NO3_gas)=0.0e-09 ! [mlc mlc-1] Gaseous NO3 volume mixing ratio
            vmr_chm(i,k,idx_O3_gas)=30.0e-09 ! [mlc mlc-1] Gaseous O3 volume mixing ratio
            vmr_chm(i,k,idx_OH_gas)=0.0e-09 ! [mlc mlc-1] Gaseous OH volume mixing ratio
            vmr_chm(i,k,idx_SO2_gas)=1.0e-09 ! [mlc mlc-1] Gaseous SO2 volume mixing ratio
c     Particulate VMRs
            vmr_SO4_aer(i,k)=0.6e-9 ! [mlc mlc-1] Particulate SO4 volume mixing ratio
            vmr_NO3_aer(i,k)=0.0e-9 ! [mlc mlc-1] Particulate NO3 volume mixing ratio
            vmr_chm(i,k,idx_SO4_aer)=vmr_SO4_aer(i,k) ! [mlc mlc-1] Particulate SO4 volume mixing ratio
            vmr_chm(i,k,idx_NO3_aer)=vmr_NO3_aer(i,k) ! [mlc mlc-1] Particulate NO3 volume mixing ratio
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      
      mss_cnc_dst=sum(q_dst,dim=3)*dns_mdp ! [kg m-3] Mass concentration of dust

      do m=1,pcnst
         do k=1,plev
            do i=1,plon
               if (m == 1) then 
c     CCM, MATCH adhere to convention that water vapor is stored as specific humidity [kg H2O vapor kg-1 moist air]
                  q_cst(i,k,m)=q_H2O_vpr(i,k) ! [kg kg-1] Water vapor mixing ratio
               else
c     CCM, MATCH adhere to convention that non-water vapor tracers are stored as dry mass mixing ratios [kg tracer kg-1 dry air]
c     Generally speaking, non-water tracers are carried as dry mass mixing ratios everywhere but the physics routines
                  q_cst(i,k,m)=q_dst(i,k,m-dst_idx_srt+1) ! [kg kg-1] Dust mixing ratio
               endif            ! endif
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst
      
      return
      end                       ! end aer_nvr_cmn_ini()
      
      subroutine nvr2nc(            
     $     fl_out               ! [sng] Name of netCDF output file
     $     )                    
c     Purpose: Output environmental properties to netCDF file
      use precision
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <aernvr.h> /* Aerosol environmental properties */
#include <dstctl.h> /* Model control variables */
#include <dbg.com> /* Debugging constants */
#include <netcdf.inc>
c     Input
      character*(*),intent(in)::fl_out      ! [sng] Name of netCDF output file
c     Output
c     Local
c     File metadata and dimension IDs
      integer dim_lon_lev_sz(3) ! [enm] Dimension IDs
      integer dim_lon_lev(2)    ! [enm] Dimension IDs
      integer dim_lon_levp(2)   ! [enm] Dimension IDs
      integer dim_lon_sz(2)     ! [enm] Dimension IDs
      integer dim_sz_src_sz(2)  ! [enm] Dimension IDs
      integer lev_dim_id        ! [enm] Dimension ID for lev
      integer levp_dim_id       ! [enm] Dimension ID for lev
      integer lon_dim_id        ! [enm] Dimension ID for lon
      integer lat_dim_id        ! [enm] Dimension ID for lat
      integer nc_id             ! File handle
      integer rcd               ! [rcd] Return success code
      integer sz_dim_id         ! [enm] Dimension ID for sz
      integer sz_grd_dim_id     ! [enm] Dimension ID for sz grid
      integer sz_src_dim_id     ! [enm] Dimension ID for sz_src
c     Variable IDs
      integer dns_mdp_id        ! [enm] Variable ID
      integer doy_id            ! [enm] Variable ID
      integer flx_LW_dwn_sfc_id ! [enm] Variable ID
      integer flx_SW_abs_sfc_id ! [enm] Variable ID
      integer hgt_mdp_id        ! [enm] Variable ID
      integer lat_dgr_id        ! [enm] Variable ID
      integer lat_id            ! Coordinate ID
      integer lat_rdn_id        ! [enm] Variable ID
      integer lev_id            ! Coordinate ID
      integer levp_id           ! Coordinate ID
      integer lnd_frc_dry_id    ! [enm] Variable ID
      integer mbl_bsn_fct_id    ! [enm] Variable ID
      integer lon_id            ! Coordinate ID
      integer mpl_air_id        ! [enm] Variable ID
      integer mss_cnc_dst_id    ! [enm] Variable ID
      integer mss_frc_CaCO3_id    ! [enm] Variable ID
      integer mss_frc_cly_id    ! [enm] Variable ID
      integer mss_frc_slt_id    ! [enm] Variable ID
      integer mss_frc_snd_id    ! [enm] Variable ID
      integer oro_id            ! [enm] Variable ID
      integer prs_dlt_id        ! [enm] Variable ID
      integer prs_mdp_id        ! [enm] Variable ID
      integer prs_ntf_id        ! [enm] Variable ID
      integer q_H2O_cnd_pcp_id  ! [enm] Variable ID
      integer q_H2O_cnd_tnd_id  ! [enm] Variable ID
      integer q_H2O_pcp_lqd_id  ! [enm] Variable ID
      integer q_H2O_vpr_id      ! [enm] Variable ID
      integer sfc_typ_id        ! [enm] Variable ID
      integer snw_hgt_lqd_id    ! [enm] Variable ID
      integer tm_adj_id         ! [enm] Variable ID
      integer tpt_gnd_id        ! [enm] Variable ID
      integer tpt_mdp_id        ! [enm] Variable ID
      integer tpt_ptn_id        ! [enm] Variable ID
      integer tpt_sfc_id        ! [enm] Variable ID
      integer tpt_soi_id        ! [enm] Variable ID
      integer tpt_vrt_id        ! [enm] Variable ID
      integer vwc_sfc_id        ! [enm] Variable ID
      integer wnd_mrd_mdp_id    ! [enm] Variable ID
      integer wnd_znl_mdp_id    ! [enm] Variable ID
c     Variable data
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
      
c     Begin netCDF output routines
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf_open(fl_out,NF_WRITE,nc_id)
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR opening ',fl_out(1:f_strlen(fl_out)),
     $     ' in nvr2nc(), rcd = ',rcd
      rcd=rcd+nf_redef(nc_id)
      rcd=rcd+nf_inq_dimid(nc_id,'levp',levp_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'lev',lev_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'lon',lon_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'lat',lat_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'sz',sz_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'sz_src',sz_src_dim_id)
c     Add global attributes
c     Define dimension IDs
c     Assemble ID and count vectors for each multidimensional combination of dimensions
      dim_lon_sz=(/lon_dim_id,sz_dim_id/)
      dim_lon_lev=(/lon_dim_id,lev_dim_id/)
      dim_lon_levp=(/lon_dim_id,levp_dim_id/)
      dim_lon_lev_sz=(/lon_dim_id,lev_dim_id,sz_dim_id/)
      dim_sz_src_sz=(/sz_src_dim_id,sz_dim_id/)
c     Variable definitions
      rcd=rcd+nf_def_var(nc_id,'dns_mdp',NF_DOUBLE,2,dim_lon_lev,dns_mdp_id)
      rcd=rcd+nf_def_var(nc_id,'doy',NF_DOUBLE,0,0,doy_id)
      rcd=rcd+nf_def_var(nc_id,'flx_LW_dwn_sfc',NF_DOUBLE,1,lon_dim_id,flx_LW_dwn_sfc_id)
      rcd=rcd+nf_def_var(nc_id,'flx_SW_abs_sfc',NF_DOUBLE,1,lon_dim_id,flx_SW_abs_sfc_id)
      rcd=rcd+nf_def_var(nc_id,'hgt_mdp',NF_DOUBLE,2,dim_lon_lev,hgt_mdp_id)
      rcd=rcd+nf_def_var(nc_id,'lat',NF_DOUBLE,1,lat_dim_id,lat_id)
      rcd=rcd+nf_def_var(nc_id,'lat_dgr',NF_DOUBLE,1,lat_dim_id,lat_dgr_id)
      rcd=rcd+nf_def_var(nc_id,'lat_rdn',NF_DOUBLE,1,lat_dim_id,lat_rdn_id)
      rcd=rcd+nf_def_var(nc_id,'lev',NF_DOUBLE,1,lev_dim_id,lev_id)
      rcd=rcd+nf_def_var(nc_id,'levp',NF_DOUBLE,1,levp_dim_id,levp_id)
      rcd=rcd+nf_def_var(nc_id,'lnd_frc_dry',NF_DOUBLE,1,lon_dim_id,lnd_frc_dry_id)
      rcd=rcd+nf_def_var(nc_id,'mbl_bsn_fct',NF_DOUBLE,1,lon_dim_id,mbl_bsn_fct_id)
      rcd=rcd+nf_def_var(nc_id,'lon',NF_DOUBLE,1,lon_dim_id,lon_id)
      rcd=rcd+nf_def_var(nc_id,'mpl_air',NF_DOUBLE,2,dim_lon_lev,mpl_air_id)
      rcd=rcd+nf_def_var(nc_id,'mss_cnc_dst',NF_DOUBLE,2,dim_lon_lev,mss_cnc_dst_id)
      rcd=rcd+nf_def_var(nc_id,'mss_frc_CaCO3',NF_DOUBLE,1,lon_dim_id,mss_frc_CaCO3_id)
      rcd=rcd+nf_def_var(nc_id,'mss_frc_cly',NF_DOUBLE,1,lon_dim_id,mss_frc_cly_id)
      rcd=rcd+nf_def_var(nc_id,'mss_frc_slt',NF_DOUBLE,1,lon_dim_id,mss_frc_slt_id)
      rcd=rcd+nf_def_var(nc_id,'mss_frc_snd',NF_DOUBLE,1,lon_dim_id,mss_frc_snd_id)
      rcd=rcd+nf_def_var(nc_id,'oro',NF_DOUBLE,1,lon_dim_id,oro_id)
      rcd=rcd+nf_def_var(nc_id,'prs_dlt',NF_DOUBLE,2,dim_lon_lev,prs_dlt_id)
      rcd=rcd+nf_def_var(nc_id,'prs_mdp',NF_DOUBLE,2,dim_lon_lev,prs_mdp_id)
      rcd=rcd+nf_def_var(nc_id,'prs_ntf',NF_DOUBLE,2,dim_lon_levp,prs_ntf_id)
      rcd=rcd+nf_def_var(nc_id,'q_H2O_cnd_pcp',NF_DOUBLE,2,dim_lon_lev,q_H2O_cnd_pcp_id)
      rcd=rcd+nf_def_var(nc_id,'q_H2O_cnd_tnd',NF_DOUBLE,2,dim_lon_lev,q_H2O_cnd_tnd_id)
      rcd=rcd+nf_def_var(nc_id,'q_H2O_pcp_lqd',NF_DOUBLE,2,dim_lon_lev,q_H2O_pcp_lqd_id)
      rcd=rcd+nf_def_var(nc_id,'q_H2O_vpr',NF_DOUBLE,2,dim_lon_lev,q_H2O_vpr_id)
      rcd=rcd+nf_def_var(nc_id,'sfc_typ',NF_INT,1,lon_dim_id,sfc_typ_id)
      rcd=rcd+nf_def_var(nc_id,'snw_hgt_lqd',NF_DOUBLE,1,lon_dim_id,snw_hgt_lqd_id)
      rcd=rcd+nf_def_var(nc_id,'tm_adj',NF_DOUBLE,0,0,tm_adj_id)
      rcd=rcd+nf_def_var(nc_id,'tpt_gnd',NF_DOUBLE,1,lon_dim_id,tpt_gnd_id)
      rcd=rcd+nf_def_var(nc_id,'tpt_mdp',NF_DOUBLE,2,dim_lon_lev,tpt_mdp_id)
      rcd=rcd+nf_def_var(nc_id,'tpt_ptn',NF_DOUBLE,2,dim_lon_lev,tpt_ptn_id)
      rcd=rcd+nf_def_var(nc_id,'tpt_sfc',NF_DOUBLE,1,lon_dim_id,tpt_sfc_id)
      rcd=rcd+nf_def_var(nc_id,'tpt_soi',NF_DOUBLE,1,lon_dim_id,tpt_soi_id)
      rcd=rcd+nf_def_var(nc_id,'tpt_vrt',NF_DOUBLE,2,dim_lon_lev,tpt_vrt_id)
      rcd=rcd+nf_def_var(nc_id,'vwc_sfc',NF_DOUBLE,1,lon_dim_id,vwc_sfc_id)
      rcd=rcd+nf_def_var(nc_id,'wnd_mrd_mdp',NF_DOUBLE,1,lon_dim_id,wnd_mrd_mdp_id)
      rcd=rcd+nf_def_var(nc_id,'wnd_znl_mdp',NF_DOUBLE,1,lon_dim_id,wnd_znl_mdp_id)
c     Add english text descriptions
      rcd=rcd+nf_put_att_text(nc_id,dns_mdp_id,'long_name',16,'Midlayer density')
      rcd=rcd+nf_put_att_text(nc_id,doy_id,'long_name',24,'Day of year [1.0..366.0)')
      rcd=rcd+nf_put_att_text(nc_id,flx_LW_dwn_sfc_id,'long_name',36,'Longwave downwelling flux at surface')
      rcd=rcd+nf_put_att_text(nc_id,flx_SW_abs_sfc_id,'long_name',29,'Solar flux absorbed by ground')
      rcd=rcd+nf_put_att_text(nc_id,hgt_mdp_id,'long_name',29,'Midlayer height above surface')
      rcd=rcd+nf_put_att_text(nc_id,lat_dgr_id,'long_name',8,'Latitude')
      rcd=rcd+nf_put_att_text(nc_id,lat_id,'long_name',8,'Latitude')
      rcd=rcd+nf_put_att_text(nc_id,lat_rdn_id,'long_name',8,'Latitude')
      rcd=rcd+nf_put_att_text(nc_id,lev_id,'long_name',17,'Midlayer pressure')
      rcd=rcd+nf_put_att_text(nc_id,levp_id,'long_name',18,'Interface pressure')
      rcd=rcd+nf_put_att_text(nc_id,lnd_frc_dry_id,'long_name',17,'Dry land fraction')
      rcd=rcd+nf_put_att_text(nc_id,mbl_bsn_fct_id,'long_name',53,'Topographic source limitation')
      rcd=rcd+nf_put_att_text(nc_id,lon_id,'long_name',9,'Longitude')
      rcd=rcd+nf_put_att_text(nc_id,mpl_air_id,'long_name',22,'Air mass path in layer')
      rcd=rcd+nf_put_att_text(nc_id,mss_cnc_dst_id,'long_name',26,'Mass concentration of dust')
      rcd=rcd+nf_put_att_text(nc_id,mss_frc_CaCO3_id,'long_name',19,'Mass fraction CaCO3')
      rcd=rcd+nf_put_att_text(nc_id,mss_frc_cly_id,'long_name',18,'Mass fraction clay')
      rcd=rcd+nf_put_att_text(nc_id,mss_frc_slt_id,'long_name',18,'Mass fraction silt')
      rcd=rcd+nf_put_att_text(nc_id,mss_frc_snd_id,'long_name',18,'Mass fraction sant')
      rcd=rcd+nf_put_att_text(nc_id,oro_id,'long_name',43,'Orography: ocean=0.0, land=1.0, sea ice=2.0')
      rcd=rcd+nf_put_att_text(nc_id,prs_dlt_id,'long_name',18,'Pressure thickness')
      rcd=rcd+nf_put_att_text(nc_id,prs_mdp_id,'long_name',17,'Midlayer pressure')
      rcd=rcd+nf_put_att_text(nc_id,prs_ntf_id,'long_name',18,'Interface pressure')
      rcd=rcd+nf_put_att_text(nc_id,q_H2O_cnd_pcp_id,'long_name',30,'H2O precipitation mixing ratio')
      rcd=rcd+nf_put_att_text(nc_id,q_H2O_cnd_tnd_id,'long_name',37,'Net H2O condensate formation tendency')
      rcd=rcd+nf_put_att_text(nc_id,q_H2O_pcp_lqd_id,'long_name',23,'Rain water mixing ratio')
      rcd=rcd+nf_put_att_text(nc_id,q_H2O_vpr_id,'long_name',24,'Water vapor mixing ratio')
      rcd=rcd+nf_put_att_text(nc_id,sfc_typ_id,'long_name',24,'LSM surface type (0..28)')
      rcd=rcd+nf_put_att_text(nc_id,snw_hgt_lqd_id,'long_name',34,'Equivalent liquid water snow depth')
      rcd=rcd+nf_put_att_text(nc_id,tm_adj_id,'long_name',42,'Adjustment timestep (CCM: 2*dt, MATCH: dt)')
      rcd=rcd+nf_put_att_text(nc_id,tpt_gnd_id,'long_name',18,'Ground temperature')
      rcd=rcd+nf_put_att_text(nc_id,tpt_mdp_id,'long_name',20,'Midlayer temperature')
      rcd=rcd+nf_put_att_text(nc_id,tpt_ptn_id,'long_name',21,'Potential temperature')
      rcd=rcd+nf_put_att_text(nc_id,tpt_sfc_id,'long_name',19,'Surface temperature')
      rcd=rcd+nf_put_att_text(nc_id,tpt_soi_id,'long_name',16,'Soil temperature')
      rcd=rcd+nf_put_att_text(nc_id,tpt_vrt_id,'long_name',19,'Virtual temperature')
      rcd=rcd+nf_put_att_text(nc_id,vwc_sfc_id,'long_name',24,'Volumetric water content')
      rcd=rcd+nf_put_att_text(nc_id,wnd_mrd_mdp_id,'long_name',25,'Meridional wind component')
      rcd=rcd+nf_put_att_text(nc_id,wnd_znl_mdp_id,'long_name',20,'Zonal wind component')
c     Add units
      rcd=rcd+nf_put_att_text(nc_id,dns_mdp_id,'units',16,'kilogram meter-3')
      rcd=rcd+nf_put_att_text(nc_id,doy_id,'units',3,'day')
      rcd=rcd+nf_put_att_text(nc_id,flx_LW_dwn_sfc_id,'units',12,'watt meter-2')
      rcd=rcd+nf_put_att_text(nc_id,flx_SW_abs_sfc_id,'units',12,'watt meter-2')
      rcd=rcd+nf_put_att_text(nc_id,hgt_mdp_id,'units',5,'meter')
      rcd=rcd+nf_put_att_text(nc_id,lat_dgr_id,'units',13,'degrees north')
      rcd=rcd+nf_put_att_text(nc_id,lat_id,'units',12,'radians north')
      rcd=rcd+nf_put_att_text(nc_id,lat_rdn_id,'units',13,'radians north')
      rcd=rcd+nf_put_att_text(nc_id,lev_id,'units',6,'Pascal')
      rcd=rcd+nf_put_att_text(nc_id,levp_id,'units',6,'Pascal')
      rcd=rcd+nf_put_att_text(nc_id,lnd_frc_dry_id,'units',8,'fraction')
      rcd=rcd+nf_put_att_text(nc_id,mbl_bsn_fct_id,'units',8,'fraction')
      rcd=rcd+nf_put_att_text(nc_id,lon_id,'units',12,'degrees east')
      rcd=rcd+nf_put_att_text(nc_id,mpl_air_id,'units',16,'kilogram meter-2')
      rcd=rcd+nf_put_att_text(nc_id,mss_cnc_dst_id,'units',16,'kilogram meter-3')
      rcd=rcd+nf_put_att_text(nc_id,mss_frc_CaCO3_id,'units',8,'fraction')
      rcd=rcd+nf_put_att_text(nc_id,mss_frc_cly_id,'units',8,'fraction')
      rcd=rcd+nf_put_att_text(nc_id,mss_frc_slt_id,'units',8,'fraction')
      rcd=rcd+nf_put_att_text(nc_id,mss_frc_snd_id,'units',8,'fraction')
      rcd=rcd+nf_put_att_text(nc_id,oro_id,'units',8,'fraction')
      rcd=rcd+nf_put_att_text(nc_id,prs_dlt_id,'units',6,'pascal')
      rcd=rcd+nf_put_att_text(nc_id,prs_mdp_id,'units',6,'pascal')
      rcd=rcd+nf_put_att_text(nc_id,prs_ntf_id,'units',6,'pascal')
      rcd=rcd+nf_put_att_text(nc_id,q_H2O_cnd_pcp_id,'units',19,'kilogram kilogram-1')
      rcd=rcd+nf_put_att_text(nc_id,q_H2O_cnd_tnd_id,'units',28,'kilogram kilogram-1 second-1')
      rcd=rcd+nf_put_att_text(nc_id,q_H2O_pcp_lqd_id,'units',19,'kilogram kilogram-1')
      rcd=rcd+nf_put_att_text(nc_id,q_H2O_vpr_id,'units',19,'kilogram kilogram-1')
      rcd=rcd+nf_put_att_text(nc_id,sfc_typ_id,'units',5,'index')
      rcd=rcd+nf_put_att_text(nc_id,snw_hgt_lqd_id,'units',5,'meter')
      rcd=rcd+nf_put_att_text(nc_id,tm_adj_id,'units',6,'second')
      rcd=rcd+nf_put_att_text(nc_id,tpt_gnd_id,'units',6,'kelvin')
      rcd=rcd+nf_put_att_text(nc_id,tpt_mdp_id,'units',6,'kelvin')
      rcd=rcd+nf_put_att_text(nc_id,tpt_ptn_id,'units',6,'kelvin')
      rcd=rcd+nf_put_att_text(nc_id,tpt_sfc_id,'units',6,'kelvin')
      rcd=rcd+nf_put_att_text(nc_id,tpt_soi_id,'units',6,'kelvin')
      rcd=rcd+nf_put_att_text(nc_id,tpt_vrt_id,'units',6,'kelvin')
      rcd=rcd+nf_put_att_text(nc_id,vwc_sfc_id,'units',14,'meter3 meter-3')
      rcd=rcd+nf_put_att_text(nc_id,wnd_mrd_mdp_id,'units',14,'meter second-1')
      rcd=rcd+nf_put_att_text(nc_id,wnd_znl_mdp_id,'units',14,'meter second-1')
c     End define mode
      rcd=rcd+nf_enddef(nc_id)
c     Write data
      rcd=rcd+nf_put_var_double(nc_id,dns_mdp_id,dns_mdp)
      rcd=rcd+nf_put_var_double(nc_id,doy_id,doy)
      rcd=rcd+nf_put_var_double(nc_id,flx_LW_dwn_sfc_id,flx_LW_dwn_sfc)
      rcd=rcd+nf_put_var_double(nc_id,flx_SW_abs_sfc_id,flx_SW_abs_sfc)
      rcd=rcd+nf_put_var_double(nc_id,hgt_mdp_id,hgt_mdp)
      rcd=rcd+nf_put_var_double(nc_id,lat_dgr_id,lat_dgr)
      rcd=rcd+nf_put_var_double(nc_id,lat_id,lat)
      rcd=rcd+nf_put_var_double(nc_id,lat_rdn_id,lat_rdn)
      rcd=rcd+nf_put_var_double(nc_id,lev_id,lev)
      rcd=rcd+nf_put_var_double(nc_id,levp_id,levp)
      rcd=rcd+nf_put_var_double(nc_id,lnd_frc_dry_id,lnd_frc_dry)
      rcd=rcd+nf_put_var_double(nc_id,mbl_bsn_fct_id,mbl_bsn_fct)
      rcd=rcd+nf_put_var_double(nc_id,lon_id,lon)
      rcd=rcd+nf_put_var_double(nc_id,mpl_air_id,mpl_air)
      rcd=rcd+nf_put_var_double(nc_id,mss_cnc_dst_id,mss_cnc_dst)
      rcd=rcd+nf_put_var_double(nc_id,mss_frc_CaCO3_id,mss_frc_CaCO3)
      rcd=rcd+nf_put_var_double(nc_id,mss_frc_cly_id,mss_frc_cly)
      rcd=rcd+nf_put_var_double(nc_id,mss_frc_slt_id,mss_frc_slt)
      rcd=rcd+nf_put_var_double(nc_id,mss_frc_snd_id,mss_frc_snd)
      rcd=rcd+nf_put_var_double(nc_id,oro_id,oro)
      rcd=rcd+nf_put_var_double(nc_id,prs_dlt_id,prs_dlt)
      rcd=rcd+nf_put_var_double(nc_id,prs_mdp_id,prs_mdp)
      rcd=rcd+nf_put_var_double(nc_id,prs_ntf_id,prs_ntf)
      rcd=rcd+nf_put_var_double(nc_id,q_H2O_cnd_pcp_id,q_H2O_cnd_pcp)
      rcd=rcd+nf_put_var_double(nc_id,q_H2O_cnd_tnd_id,q_H2O_cnd_tnd)
      rcd=rcd+nf_put_var_double(nc_id,q_H2O_pcp_lqd_id,q_H2O_pcp_lqd)
      rcd=rcd+nf_put_var_double(nc_id,q_H2O_vpr_id,q_H2O_vpr)
      rcd=rcd+nf_put_var_double(nc_id,snw_hgt_lqd_id,snw_hgt_lqd)
      rcd=rcd+nf_put_var_double(nc_id,tm_adj_id,tm_adj)
      rcd=rcd+nf_put_var_double(nc_id,tpt_gnd_id,tpt_gnd)
      rcd=rcd+nf_put_var_double(nc_id,tpt_mdp_id,tpt_mdp)
      rcd=rcd+nf_put_var_double(nc_id,tpt_ptn_id,tpt_ptn)
      rcd=rcd+nf_put_var_double(nc_id,tpt_sfc_id,tpt_sfc)
      rcd=rcd+nf_put_var_double(nc_id,tpt_soi_id,tpt_soi)
      rcd=rcd+nf_put_var_double(nc_id,tpt_vrt_id,tpt_vrt)
      rcd=rcd+nf_put_var_double(nc_id,vwc_sfc_id,vwc_sfc)
      rcd=rcd+nf_put_var_double(nc_id,wnd_mrd_mdp_id,wnd_mrd_mdp)
      rcd=rcd+nf_put_var_double(nc_id,wnd_znl_mdp_id,wnd_znl_mdp)
      rcd=rcd+nf_put_var_int(nc_id,sfc_typ_id,sfc_typ)
c     Close output file
      rcd=rcd+nf_close(nc_id)
      call strnul(fl_out)
      if (nstep == 1) then
         write (6,'(a,a45,a)') prg_nm(1:f_strlen(prg_nm)),': Initialized environmental data archive in',fl_out(1:f_strlen(fl_out))
      endif                     ! endif
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR writing ',fl_out(1:f_strlen(fl_out)),
     $     ' in nvr2nc(), rcd = ',rcd
      
      return 
      end                       ! end nvr2nc()
      
      subroutine tm2nc(            
     $     fl_out,              ! [sng] Name of netCDF output file
     $     tm_adj)              ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
c     Purpose: Output a single time step to netCDF file
      use precision
      implicit none
#include <dstctl.h> /* Model control variables */
#include <dbg.com> /* Debugging constants */
#include <netcdf.inc>
c     Input
      character*(*),intent(in)::fl_out      ! [sng] Name of netCDF output file
      real(r8),intent(in)::tm_adj               ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
c     Output
c     Local
      character*31 time_unit_sng ! [sng] Units attribute for time coordinate
      double precision time     ! [day] Days since simulation start
c     File metadata and dimension IDs
      integer time_dim_id       ! [enm] Dimension ID for time
      integer nc_id             ! File handle
      integer rcd               ! [rcd] Return success code
c     Variable IDs
      integer nstep_id          ! [enm] Variable ID
      integer mcdate_id         ! [enm] Variable ID
      integer mcsec_id          ! [enm] Variable ID
      integer nbdate_id         ! [enm] Variable ID
      integer nbsec_id          ! [enm] Variable ID
      integer ndcur_id          ! [enm] Variable ID
      integer nscur_id          ! [enm] Variable ID
      integer time_id           ! [enm] Variable ID
c     Variable data
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
      character*10 date2sng     ! [sng] Convert YYYYMMDD integer to YYYY-MM-DD string (libcsz_f77)
      character*8 sec2sng       ! [sng] Convert integer to 'HH:MM:SS' string (libcsz_f77)
      
c     Must be careful that variable "time" is defined in routine without access to F77 intrinsic time() function
      if (nstep == 1) then
         nbdate=640312          ! [day] Simulation start date in YYMMDD format
         nbsec=0                ! [s] Simulation start second relative to nbdate
         mcdate=nbdate          ! Current date in YYMMDD format
         mcsec=nbsec            ! Seconds past current date at 0Z
         ndcur=0                ! [day] Current day number of simulation[day] 
         nscur=0                ! [s] Seconds relative to ndcur
         call strnul(time_unit_sng) ! [sng] Units attribute for time coordinate
         time_unit_sng='days since '//date2sng(mcdate)//' '//sec2sng(mcsec)
      endif                     ! endif nstep == 1

c     Increment time coordinates
      mcsec=mcsec+tm_adj        ! [s] Seconds past current date at 0Z
      nscur=nscur+tm_adj        ! [s] Seconds relative to ndcur

c     Take care of odometer issues
      if (nscur >= 86400.0) then
         ndcur=ndcur+nscur/86400 ! [day] Current day number of simulation
         nscur=mod(nscur,86400) ! [s] Seconds relative to ndcur
      end if                    ! endif nscur >= 86400.0
      if (mcsec >= 86400.0) then
         mcdate=mcdate+mcsec/86400 ! [day] Current date in YYMMDD format
         mcsec=mod(mcsec,86400) ! [s] Seconds past current date at 0Z
      end if                    ! endif nscur >= 86400.0
c     Construct time from ndcur and nscur
      time=ndcur+nscur/86400.0  ! [day] Days since simulation start

c     Begin netCDF output routines
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf_open(fl_out,NF_WRITE,nc_id)
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR opening ',fl_out(1:f_strlen(fl_out)),
     $     ' in tm2nc(), rcd = ',rcd
      rcd=rcd+nf_redef(nc_id)
c     Get dimension IDs
      rcd=rcd+nf_inq_dimid(nc_id,'time',time_dim_id)
      if (nstep == 1) then
c     Variable definitions
         rcd=rcd+nf_def_var(nc_id,'nstep',NF_INT,1,time_dim_id,nstep_id)
         rcd=rcd+nf_def_var(nc_id,'mcdate',NF_INT,1,time_dim_id,mcdate_id)
         rcd=rcd+nf_def_var(nc_id,'mcsec',NF_INT,1,time_dim_id,mcsec_id)
         rcd=rcd+nf_def_var(nc_id,'nbdate',NF_INT,0,0,nbdate_id)
         rcd=rcd+nf_def_var(nc_id,'nbsec',NF_INT,0,0,nbsec_id)
         rcd=rcd+nf_def_var(nc_id,'ndcur',NF_INT,1,time_dim_id,ndcur_id)
         rcd=rcd+nf_def_var(nc_id,'nscur',NF_INT,1,time_dim_id,nscur_id)
         rcd=rcd+nf_def_var(nc_id,'time',NF_DOUBLE,1,time_dim_id,time_id)
c     Add english text descriptions
         rcd=rcd+nf_put_att_text(nc_id,nstep_id,'long_name',8,'Timestep')
         rcd=rcd+nf_put_att_text(nc_id,mcdate_id,'long_name',29,'Current date in YYMMDD format')
         rcd=rcd+nf_put_att_text(nc_id,mcsec_id,'long_name',31,'Seconds past current date at 0Z')
         rcd=rcd+nf_put_att_text(nc_id,nbdate_id,'long_name',38,'Simulation start date in YYMMDD format')
         rcd=rcd+nf_put_att_text(nc_id,nbsec_id,'long_name',42,'Simulation start second relative to nbdate')
         rcd=rcd+nf_put_att_text(nc_id,ndcur_id,'long_name',32,'Current day number of simulation')
         rcd=rcd+nf_put_att_text(nc_id,nscur_id,'long_name',25,'Seconds relative to ndcur')
         rcd=rcd+nf_put_att_text(nc_id,time_id,'long_name',27,'Days since simulation start')
c     Add units
         rcd=rcd+nf_put_att_text(nc_id,nstep_id,'units',5,'index')
         rcd=rcd+nf_put_att_text(nc_id,mcdate_id,'units',3,'day')
         rcd=rcd+nf_put_att_text(nc_id,mcsec_id,'units',6,'second')
         rcd=rcd+nf_put_att_text(nc_id,nbdate_id,'units',3,'day')
         rcd=rcd+nf_put_att_text(nc_id,nbsec_id,'units',6,'second')
         rcd=rcd+nf_put_att_text(nc_id,ndcur_id,'units',3,'day')
         rcd=rcd+nf_put_att_text(nc_id,nscur_id,'units',6,'second')
         rcd=rcd+nf_put_att_text(nc_id,time_id,'units',31,time_unit_sng)
      else                      ! endif nstep == 1
c     Get variable IDs
         rcd=rcd+nf_inq_varid(nc_id,'nstep',nstep_id)
         rcd=rcd+nf_inq_varid(nc_id,'mcdate',mcdate_id)
         rcd=rcd+nf_inq_varid(nc_id,'mcsec',mcsec_id)
         rcd=rcd+nf_inq_varid(nc_id,'nbdate',nbdate_id)
         rcd=rcd+nf_inq_varid(nc_id,'nbsec',nbsec_id)
         rcd=rcd+nf_inq_varid(nc_id,'ndcur',ndcur_id)
         rcd=rcd+nf_inq_varid(nc_id,'nscur',nscur_id)
         rcd=rcd+nf_inq_varid(nc_id,'time',time_id)
      endif                     ! endif nstep /= 1
c     End define mode
      rcd=rcd+nf_enddef(nc_id)
c     Write data
      if (nstep == 1) then
         rcd=rcd+nf_put_var_int(nc_id,nbdate_id,nbdate)
         rcd=rcd+nf_put_var_int(nc_id,nbsec_id,nbsec)
      endif                     ! endif nstep == 1
      rcd=rcd+nf_put_vara_int(nc_id,nstep_id,nstep,1,nstep)
      rcd=rcd+nf_put_vara_int(nc_id,mcdate_id,nstep,1,mcdate)
      rcd=rcd+nf_put_vara_int(nc_id,mcsec_id,nstep,1,mcsec)
      rcd=rcd+nf_put_vara_int(nc_id,ndcur_id,nstep,1,ndcur)
      rcd=rcd+nf_put_vara_int(nc_id,nscur_id,nstep,1,nscur)
      rcd=rcd+nf_put_vara_double(nc_id,time_id,nstep,1,time)
c     Close output file
      rcd=rcd+nf_close(nc_id)
      call strnul(fl_out)
      if (nstep == 1) then
         write (6,'(a,a38,a)') prg_nm(1:f_strlen(prg_nm)),': Initialized control data archive in ',fl_out(1:f_strlen(fl_out))
      endif                     ! endif
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR writing ',fl_out(1:f_strlen(fl_out)),
     $     ' in tm2nc(), rcd = ',rcd
      
      return 
      end                       ! end tm2nc()
      
      subroutine dst_sfc_set(lon_nbr_in, ! I Dimension size (normally plond)
     $     flx_LW_dwn_sfc_in,   ! I [W m-2] Longwave downwelling flux at surface
     $     flx_SW_abs_sfc_in,   ! I [W m-2] Solar flux absorbed by ground
     $     lnd_frc_dry_in,      ! I [frc] Dry land fraction
     $     mbl_bsn_fct_in,      ! I [frc] Topographic source limitation
     $     mss_frc_CaCO3_in,    ! I [frc] Mass fraction of CaCO3
     $     mss_frc_cly_in,      ! I [frc] Mass fraction of clay
     $     mss_frc_snd_in,      ! I [frc] Mass fraction of sand
     $     sfc_typ_in,          ! I [idx] LSM surface type (0..28)
     $     tpt_gnd_in,          ! I [K] Ground temperature
     $     tpt_soi_in,          ! I [K] Soil temperature
     $     vwc_sfc_in)          ! I [m3 m-3] Volumetric water content
c     Set up common block dst_sfc in one call
      use precision
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstsfc.h> /* Surface fields for dust parameterization */
c     Parameters
      integer,parameter::lat_nbr=PLAT
      integer,parameter::lon_nbr=PLON
c     Input
      integer,intent(in)::lon_nbr_in        ! Dimension size (normally plond)
      integer,intent(in)::sfc_typ_in(lon_nbr_in) ! I [idx] LSM surface type (0..28)
      real(r8),intent(in)::lnd_frc_dry_in(lon_nbr_in) ! I [frc] Dry land fraction
      real(r8),intent(in)::mbl_bsn_fct_in(lon_nbr_in) ! I [frc] Topographic source limitation
      real(r8),intent(in)::flx_LW_dwn_sfc_in(lon_nbr_in) ! I [W m-2] Longwave downwelling flux at surface
      real(r8),intent(in)::flx_SW_abs_sfc_in(lon_nbr_in) ! I [W m-2] Solar flux absorbed by ground
      real(r8),intent(in)::mss_frc_CaCO3_in(lon_nbr_in) ! [frc] Mass fraction of CaCO3
      real(r8),intent(in)::mss_frc_cly_in(lon_nbr_in) ! [frc] Mass fraction of clay
      real(r8),intent(in)::mss_frc_snd_in(lon_nbr_in) ! [frc] Mass fraction of sand
      real(r8),intent(in)::tpt_gnd_in(lon_nbr_in) ! I [m3 m-3] Volumetric water content
      real(r8),intent(in)::tpt_soi_in(lon_nbr_in) ! I [m3 m-3] Volumetric water content
      real(r8),intent(in)::vwc_sfc_in(lon_nbr_in) ! I [m3 m-3] Volumetric water content
c     Output
c     Local
      integer lat_idx           ! Latitude index
      integer lon_idx           ! [idx] Counting index for lon
      
      if (lon_nbr_in < lon_nbr) then
         write (6,'(a)') 'dst: ERROR dst_sfc_set() reports lon_nbr_in < lon_nbr'
         stop
      endif                     ! endif err
      
      do lat_idx=1,lat_nbr
         do lon_idx=1,lon_nbr
            flx_LW_dwn_sfc(lon_idx,lat_idx)=flx_LW_dwn_sfc_in(lon_idx)
            flx_SW_abs_sfc(lon_idx,lat_idx)=flx_SW_abs_sfc_in(lon_idx)
            lnd_frc_dry(lon_idx,lat_idx)=lnd_frc_dry_in(lon_idx)
            mbl_bsn_fct(lon_idx,lat_idx)=mbl_bsn_fct_in(lon_idx)
            mss_frc_CaCO3(lon_idx,lat_idx)=mss_frc_CaCO3_in(lon_idx)
            mss_frc_cly(lon_idx,lat_idx)=mss_frc_cly_in(lon_idx)
            mss_frc_snd(lon_idx,lat_idx)=mss_frc_snd_in(lon_idx)
            sfc_typ(lon_idx,lat_idx)=sfc_typ_in(lon_idx)
            tpt_gnd(lon_idx,lat_idx)=tpt_gnd_in(lon_idx)
            tpt_soi(lon_idx,lat_idx)=tpt_soi_in(lon_idx)
            vwc_sfc(lon_idx,lat_idx)=vwc_sfc_in(lon_idx)
         end do                 ! end loop over lon
      end do                    ! end loop over lat
      
      return
      end                       ! end dst_sfc_set()
      
#else /* not BXM */
c     The following is just to allow the compiler to compile something
c     when it attempts to compile this routine in nested mode.
c     This routine is not needed or used when Dust is run stand-alone
c     But, if the user tells the compiler to compile this code without
c     anything here -- the compiler complains.  A more sophisticated
c     build mechanism could decide which files should and should not
c     be built, but the simple solution is to do as follows...
      subroutine aer_stb
      stop 'aer_stb() error: BXM not set, this routine should not be called'
      end                       ! end aer_stb()
#endif  /* not BXM */
      
