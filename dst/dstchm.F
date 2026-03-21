c     $Header$
      
c     Purpose: dstchm.F contains chemistry subroutines for the dust model
      
      subroutine dst_chm_cmn_ini()
c     Purpose: Initialize chemistry common blocks
c     dst_chm_cmn_ini() is called by dst_msc_cmn_ini()
      use precision
      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstchm.h> /* Physical constants for dust routines */
c     Input
c     Output
c     Input/Output
c     Local
c     Main code
      
c     Sanity check
      if (idx_chm_end /= chm_nbr) stop 'dst_chm_cmn_ini() reports idx_chm_end /= chm_nbr'
      
c     H2O2 --> dust 
c     ZhC99 p. 357 Tbl. 2 use gamma(H2O2->dust)=1.0e-4
      upt_cff_H2O2_dst=1.0e-4   ! [frc] Uptake coefficient for H2O2 to dust
      
c     HNO3 --> dust 
c     Vicki Grassian (personal communication, 1999) suggests gamma(HNO3->dust) ~ 0.005
c     DCZ96 p. 22872 (9) use gamma(HNO3->dust)=0.1 when soil balance is alkaline
c     ZhC99 p. 357 Tbl. 2 use gamma(HNO3->dust)=0.01
      upt_cff_HNO3_dst=0.005    ! [frc] Uptake coefficient for HNO3 to dust
      
c     HO2 --> dust 
c     DCZ96 p. 22871 use gamma(HO2->dust)=0.1
c     ZhC99 p. 357 Tbl. 2 use gamma(HO2->dust)=0.1
      upt_cff_HO2_dst=0.1       ! [frc] Uptake coefficient for HO2 to dust
      
c     N2O5 --> dust 
c     DCZ96 p. 22871 use gamma(N2O5->dust)=0.1 but suggest ~0.001 for low RH
c     ZhC99 p. 357 Tbl. 2 use gamma(N2O5->dust)=0.1
      upt_cff_N2O5_dst=0.001    ! [frc] Uptake coefficient for N2O5 to dust
      
c     NO3 --> dust 
c     ZhC99 p. 357 Tbl. 2 use gamma(NO3->dust)=0.1
      upt_cff_NO3_dst=0.1       ! [frc] Uptake coefficient for NO3 to dust
      
c     O3 --> dust 
c     DCZ96 p. 22874 use gamma(O3->dust)=5.0e-5 based on inversion of deposition measurements
c     ZhC99 p. 357 Tbl. 2 use gamma(O3->dust)=1.0e-4
      upt_cff_O3_dst=5.0e-5     ! [frc] Uptake coefficient for O3 to dust
      
c     OH --> dust 
c     ZhC99 p. 357 Tbl. 2 use gamma(OH->dust)=0.1
      upt_cff_OH_dst=0.1        ! [frc] Uptake coefficient for OH to dust
      
c     SO2 --> dust 
c     DCZ96 p. 22872 use gamma(SO2->dust)=3.0e-4 when soil balance is alkaline
c     ZhC99 p. 357 Tbl. 2 use gamma(SO2->dust)=1.0e-4
      upt_cff_SO2_dst=3.0e-4    ! [frc] Uptake coefficient for SO2 to dust
      
      vmr_HNO3_gas_mxm=1.0e-3   ! [mlc mlc-1] Maximum gaseous HNO3 volume mixing ratio
      return
      end subroutine dst_chm_cmn_ini
      
      subroutine dst_chm_rxr(hbuf,obuf,wrhstts,
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     lat_idx,             ! I [idx] Model latitude index
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_dst,               ! I [kg kg-1] Dust mixing ratio
     $     rxrc_HNO3_dst,       ! O [s-1] Pseudo first order rate coefficient for HNO3
     $     rxrc_HO2_dst,        ! O [s-1] Pseudo first order rate coefficient for HO2
     $     rxrc_N2O5_dst,       ! O [s-1] Pseudo first order rate coefficient for N2O5
     $     rxrc_O3_dst,         ! O [s-1] Pseudo first order rate coefficient for O3
     $     rxrc_SO2_dst,        ! O [s-1] Pseudo first order rate coefficient for SO2
     $     tpt_mdp,             ! I [K] Temperature
     $     tpt_vrt,             ! I [K] Virtual temperature
     $     vmr_SO4_aer)         ! I [mlc mlc-1] Particulate SO4 volume mixing ratio
c     Heterogeneous chemical processes on dust
c     dst_chm_rxr() will be called by MATCH:src/physlic()
      use precision
      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstaer.h> /* Aerosol microphysical properties */
#include <dstchm.h> /* Chemical properties of dust */
c     Parameters
      real(r8),parameter::cst_Avagadro=6.022045e+23 ! [mlc mol-1] Avagadro's number
      real(r8),parameter::mmw_Ca=40.078e-03 ! [kg mol-1] Mean molecular weight of Ca
      real(r8),parameter::mmw_CaCO3=100.087e-03 ! [kg mol-1] Mean molecular weight of CaCO3
      real(r8),parameter::mmw_HNO3=62.995644e-03 ! [kg mol-1] Mean molecular weight of HNO3
      real(r8),parameter::mmw_HO2=33.0067e-03 ! [kg mol-1] Mean molecular weight of HO2
      real(r8),parameter::mmw_N2O5=108.01e-03 ! [kg mol-1] Mean molecular weight of N2O5
      real(r8),parameter::mmw_NO3=62.0049e-03 ! [kg mol-1] Mean molecular weight of NO3
      real(r8),parameter::mmw_O3=47.997832e-03 ! [kg mol-1] Mean molecular weight of O3
      real(r8),parameter::mmw_SO2=64.046674e-03 ! [kg mol-1] Mean molecular weight of SO2
      real(r8),parameter::mmw_SO4=96.0636e-03 ! [kg mol-1] Mean molecular weight of SO4
c     Global continental mean CaCO3 content from IGBP-DIS is 2.5%
c     ncwa -O -m lnd_msk -M 1.0 -o eq -a lat,lon -w gw -C -v mss_frc_CaCO3 dst_T42.nc foo.nc ; ncks -H foo.nc
c     Approximate mean over dust source regions is 5%
c     NB: DCZ96 quote sources saying Ca (not CaCO3) content is 3.6%, which clearly disagrees with IGBP-DIS
      real(r8),parameter::mss_frc_CaCO3=0.05 ! [frc] Soil calcium carbonate fraction 
c     Input
      integer,intent(in)::lat_idx ! I [idx] Model latitude index
      logical,intent(in)::wrhstts ! I [flg] Flag to output buffer to file this timestep
      real(r8),intent(in)::dns_mdp(plond,plev) ! I [kg m-3] Midlayer density
      real(r8),intent(in)::hbuf(*)  ! I [ptr] History buffer
      real(r8),intent(in)::obuf(*)  ! I [ptr] Output buffer
      real(r8),intent(in)::prs_dlt(plond,plev) ! I [Pa] Pressure thickness
      real(r8),intent(in)::prs_mdp(plond,plev) ! I [Pa] Midlayer pressure
      real(r8),intent(in)::tpt_mdp(plond,plev) ! I [K] Temperature
      real(r8),intent(in)::tpt_vrt(plond,plev) ! I [K] Virtual temperature
      real(r8),intent(in)::q_dst(plond,plev,dst_nbr) ! I [kg kg-1] Dust mixing ratio
      real(r8),intent(in)::vmr_SO4_aer(plond,plev) ! I [mlc mlc-1] Particulate SO4 volume mixing ratio
c     Input/Output
c     Output
      real(r8),intent(out)::rxrc_HNO3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HNO3
      real(r8),intent(out)::rxrc_HO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HO2
      real(r8),intent(out)::rxrc_N2O5_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for N2O5
      real(r8),intent(out)::rxrc_O3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for O3
      real(r8),intent(out)::rxrc_SO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for SO2
c     Local Output
c     Local
      character fl_out*80       ! [sng] Name of netCDF output file
      double precision pi       ! [frc] 3
      integer i                 ! [idx] Counting index over lon
      integer k                 ! [idx] Counting index over lev
      integer m                 ! [idx] Counting index over cst
      real(r8) mpl_air(plond,plev)  ! [kg m-2] Air mass path in layer
      real(r8) sfc_ttl(plond,plev)  ! [m2 kg-1] Total surface area
      real(r8) tm_dlt               ! [s] Chemistry timestep
      real(r8) q_dst_ttl(plond,plev) ! [kg kg-1] Total dust mixing ratio
      real(r8) rxrc_fct(plond,plev) ! [frc] Rate coefficient factor = (1/4)*area*sqrt(T)
      real(r8) vlc_MWB_HNO3_fct     ! [m s-1] Thermal speed of HNO3 factor
      real(r8) vlc_MWB_HO2_fct      ! [m s-1] Thermal speed of HO2 factor
      real(r8) vlc_MWB_N2O5_fct     ! [m s-1] Thermal speed of N2O5 factor
      real(r8) vlc_MWB_O3_fct       ! [m s-1] Thermal speed of O3 factor
      real(r8) vlc_MWB_SO2_fct      ! [m s-1] Thermal speed of SO2 factor
      real(r8) vmr_CaCO3_aer(plond,plev) ! [mlc mlc-1] Particulate CaCO3 volume mixing ratio
      real(r8) vmr_CaCO3_fct        ! [kg kg-1] Factor needed for CaCO3 
      real(r8) vmr_NO3_aer(plond,plev) ! [mlc mlc-1] Particulate NO3 volume mixing ratio
c     Local output
c     GCM diagnostics
      
c     Main Code
      
c     Initialize fluxes and tendencies
      pi=4.0*atan(1.0d0)        ! [frc] 3
      sfc_ttl=0.0               ! [m2 kg-1] Total surface area
      rxrc_HNO3_dst=0.0         ! [s-1] Pseudo first order rate coefficient for HNO3
      rxrc_HO2_dst=0.0          ! [s-1] Pseudo first order rate coefficient for HO2
      rxrc_N2O5_dst=0.0         ! [s-1] Pseudo first order rate coefficient for N2O5
      rxrc_O3_dst=0.0           ! [s-1] Pseudo first order rate coefficient for O3
      rxrc_SO2_dst=0.0          ! [s-1] Pseudo first order rate coefficient for SO2
      vmr_NO3_aer=0.0           ! [mlc mlc-1] Particulate NO3 volume mixing ratio
      
c     Time invariant scalar factor for thermal molecular speeds faciliates rate coefficient computation
      vlc_MWB_HNO3_fct=sqrt(8.0*gas_cst_unv/(pi*mmw_HNO3)) ! [m s-1] Thermal speed of HNO3 factor
      vlc_MWB_HO2_fct=sqrt(8.0*gas_cst_unv/(pi*mmw_HO2)) ! [m s-1] Thermal speed of HO2 factor
      vlc_MWB_N2O5_fct=sqrt(8.0*gas_cst_unv/(pi*mmw_N2O5)) ! [m s-1] Thermal speed of N2O5 factor
      vlc_MWB_O3_fct=sqrt(8.0*gas_cst_unv/(pi*mmw_O3)) ! [m s-1] Thermal speed of O3 factor
      vlc_MWB_SO2_fct=sqrt(8.0*gas_cst_unv/(pi*mmw_SO2)) ! [m s-1] Thermal speed of SO2 factor
      
c     Compute necessary derived fields
c     Mass of air currently in gridbox
      mpl_air=prs_dlt*grv_sfc_rcp ! [kg m-2]
      
c     Surface area and CaCO3 mixing ratio
      call dst_add_lon_lev(q_dst,q_dst_ttl) ! [kg kg-1] Total dust mixing ratio
      do m=1,dst_nbr
         do k=1,plev
            do i=1,plon
               sfc_ttl(i,k)=sfc_ttl(i,k) ! [m2 kg-1] Total surface area
     $              +q_dst(i,k,m)*mpl_air(i,k)*sfc_spc_rsl(m)
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst
      
c     Volume mixing ratios for electroneutrality condition
      vmr_CaCO3_fct=mmw_dry_air/mmw_CaCO3 ! [kg kg-1] Factor needed for CaCO3 volume mixing ratio
      vmr_CaCO3_aer=            ! [mlc mlc-1] Particulate CaCO3 volume mixing ratio
     $     q_dst_ttl*mss_frc_CaCO3 ! Converts dust to CaCO3 mass mixing ratio
     $     *vmr_CaCO3_fct       ! Converts CaCO3 mmr to CaCO3 vmr
      
c     Pseudo first order reaction rate coefficients asumme Maxwell-Boltzmann statistics and kinetic regime (surface area limited) uptake
c     k = ( Aerosol surface area * gamma * molecular velocity ) / 4.0
      rxrc_fct=0.25*sfc_ttl*sqrt(tpt_mdp) ! [frc] Rate coefficient factor = (1/4)*area*sqrt(T)
      rxrc_HO2_dst=rxrc_fct*upt_cff_HO2_dst*vlc_MWB_HO2_fct ! [s-1] Pseudo first order rate coefficient for HO2
      rxrc_N2O5_dst=rxrc_fct*upt_cff_N2O5_dst*vlc_MWB_N2O5_fct ! [s-1] Pseudo first order rate coefficient for N2O5
      rxrc_O3_dst=rxrc_fct*upt_cff_O3_dst*vlc_MWB_O3_fct ! [s-1] Pseudo first order rate coefficient for O3
      
c     Electroneutrality condition based on DCZ96 p. 22872:
c     Uptake of SO2 and of HNO3 proceeds only when aerosol contains enough basic (alkaline) material to neutralize the acidity of the added species
c     This stoichiometric relationship is based on in mineral dust samples from Asiac     Since the reaction only proceeds when a threshold is exceeded, the process causes some spikiness in the output
      where (vmr_NO3_aer+2.0*vmr_SO4_aer < 2.0*vmr_CaCO3_aer) 
         rxrc_HNO3_dst=rxrc_fct*upt_cff_HNO3_dst*vlc_MWB_HNO3_fct ! [s-1] Pseudo first order rate coefficient for HNO3
         rxrc_SO2_dst=rxrc_fct*upt_cff_SO2_dst*vlc_MWB_SO2_fct ! [s-1] Pseudo first order rate coefficient for SO2
      end where                 ! end where electroneutrality condition satisfied
      
      return
      end subroutine dst_chm_rxr
      
      subroutine dst_chm_slv(hbuf,obuf,wrhstts,
     $     lat_idx,             ! I [idx] Model latitude index
     $     rxrc_HNO3_dst,       ! I [s-1] Pseudo first order rate coefficient for HNO3
     $     rxrc_HO2_dst,        ! I [s-1] Pseudo first order rate coefficient for HO2
     $     rxrc_N2O5_dst,       ! I [s-1] Pseudo first order rate coefficient for N2O5
     $     rxrc_O3_dst,         ! I [s-1] Pseudo first order rate coefficient for O3
     $     rxrc_SO2_dst,        ! I [s-1] Pseudo first order rate coefficient for SO2
     $     tm_adj,              ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     vmr_HNO3_gas,        ! I/O [mlc mlc-1] Gaseous HNO3 volume mixing ratio
     $     vmr_HO2_gas,         ! I/O [mlc mlc-1] Gaseous HO2 volume mixing ratio
     $     vmr_N2O5_gas,        ! I/O [mlc mlc-1] Gaseous N2O5 volume mixing ratio
     $     vmr_NO3_aer,         ! I/O [mlc mlc-1] Particulate NO3 volume mixing ratio
     $     vmr_O3_gas,          ! I/O [mlc mlc-1] Gaseous O3 volume mixing ratio
     $     vmr_SO2_gas,         ! I/O [mlc mlc-1] Gaseous SO2 volume mixing ratio
     $     vmr_SO4_aer)         ! I/O [mlc mlc-1] Particulate SO4 volume mixing ratio
c     Purpose: Chemistry solver for heterogeneous chemical processes on dust
c     dst_chm_slv() will be called by MATCH:src/physlic()
      use precision
      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
c     Parameters
c     Input
      integer,intent(in)::lat_idx ! I [idx] Model latitude index
      logical,intent(in)::wrhstts ! I [flg] Flag to output buffer to file this timestep
      real(r8),intent(in)::hbuf(*)  ! I [ptr] History buffer
      real(r8),intent(in)::obuf(*)  ! I [ptr] Output buffer
      real(r8),intent(in)::tm_adj   ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
      real(r8),intent(in)::rxrc_HNO3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HNO3
      real(r8),intent(in)::rxrc_HO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HO2
      real(r8),intent(in)::rxrc_N2O5_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for N2O5
      real(r8),intent(in)::rxrc_O3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for O3
      real(r8),intent(in)::rxrc_SO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for SO2
c     Input/Output
      real(r8),intent(inout)::vmr_HNO3_gas(plond,plev) ! I [mlc mlc-1] Gaseous HNO3 volume mixing ratio
      real(r8),intent(inout)::vmr_HO2_gas(plond,plev) ! I [mlc mlc-1] Gaseous HO2 volume mixing ratio
      real(r8),intent(inout)::vmr_N2O5_gas(plond,plev) ! I [mlc mlc-1] Gaseous N2O5 volume mixing ratio
      real(r8),intent(inout)::vmr_O3_gas(plond,plev) ! I [mlc mlc-1] Gaseous O3 volume mixing ratio
      real(r8),intent(inout)::vmr_SO2_gas(plond,plev) ! I [mlc mlc-1] Gaseous SO2 volume mixing ratio
      real(r8),intent(inout)::vmr_NO3_aer(plond,plev) ! I/O [mlc mlc-1] Particulate NO3 volume mixing ratio
      real(r8),intent(inout)::vmr_SO4_aer(plond,plev) ! I/O [mlc mlc-1] Particulate SO4 volume mixing ratio
c     Output
c     Local Output
      real(r8) rxr_HNO3_gas_dst_vmr(plond,plev) ! O [mlc mlc-1 s-1] Mean rate of HNO3 removal by dust
      real(r8) rxr_HO2_gas_dst_vmr(plond,plev) ! O [mlc mlc-1 s-1] Mean rate of HO2 removal by dust
      real(r8) rxr_N2O5_gas_dst_vmr(plond,plev) ! O [mlc mlc-1 s-1] Mean rate of N2O5 removal by dust
      real(r8) rxr_SO2_gas_dst_vmr(plond,plev) ! O [mlc mlc-1 s-1] Mean rate of SO2 removal by dust
      real(r8) rxr_O3_gas_dst_vmr(plond,plev) ! O [mlc mlc-1 s-1] Mean rate of O3 removal by dust
c     Local
      character fl_out*80       ! [sng] Name of netCDF output file
      integer i                 ! [idx] Counting index over lon
      integer k                 ! [idx] Counting index over lev
      integer m                 ! [idx] Counting index over cst
      real(r8) tm_dlt               ! [s] Chemistry timestep
      real(r8) vmr_HNO3_gas_dlt(plond,plev) ! [mlc mlc-1] Change in gaseous HNO3 volume mixing ratio
      real(r8) vmr_SO2_gas_dlt(plond,plev) ! [mlc mlc-1] Change in gaseous SO2 volume mixing ratio
c     Local output
c     GCM diagnostics
      
c     Main Code
c     Timesplit if desired
      tm_dlt=tm_adj             ! [s] Chemistry timestep
      
c     Initialize fluxes and tendencies
      rxr_HNO3_gas_dst_vmr=0.0  ! [mlc mlc-1 s-1] Rate of HNO3 removal by dust
      rxr_HO2_gas_dst_vmr=0.0   ! [mlc mlc-1 s-1] Rate of HO2 removal by dust
      rxr_N2O5_gas_dst_vmr=0.0  ! [mlc mlc-1 s-1] Rate of N2O5 removal by dust
      rxr_O3_gas_dst_vmr=0.0    ! [mlc mlc-1 s-1] Rate of O3 removal by dust
      rxr_SO2_gas_dst_vmr=0.0   ! [mlc mlc-1 s-1] Rate of SO2 removal by dust
      
c     Compute heterogeneous removal rates and adjust mixing ratios
      
c     HNO3
      vmr_HNO3_gas_dlt=vmr_HNO3_gas*(1.0-exp(-rxrc_HNO3_dst*tm_dlt)) ! [mlc mlc-1] Change in gaseous HNO3 volume mixing ratio
      rxr_HNO3_gas_dst_vmr=vmr_HNO3_gas_dlt/tm_dlt ! [mlc mlc-1 s-1] Mean rate of HNO3 removal by dust
      vmr_HNO3_gas=vmr_HNO3_gas-vmr_HNO3_gas_dlt ! [mlc mlc-1] Gaseous HNO3 volume mixing ratio
      vmr_NO3_aer=vmr_NO3_aer+vmr_HNO3_gas_dlt ! [mlc mlc-1] Particulate NO3 volume mixing ratio
      
c     SO2
      vmr_SO2_gas_dlt=vmr_SO2_gas*(1.0-exp(-rxrc_SO2_dst*tm_dlt)) ! [mlc mlc-1] Change in gaseous SO2 volume mixing ratio
      vmr_SO2_gas=vmr_SO2_gas-vmr_SO2_gas_dlt ! [mlc mlc-1] Gaseous SO2 volume mixing ratio
      vmr_SO4_aer=vmr_SO4_aer+vmr_SO2_gas_dlt ! [mlc mlc-1] Particulate SO4 volume mixing ratio
      
#ifdef DST_DBG
      do k=1,plev
         do i=1,plon
            if (vmr_HNO3_gas(i,k) > vmr_HNO3_gas_mxm) write(6,'(a,i2,a,i3,a,i2,a,1p,e8.1,0p,a)') 
     $           'dst_chm_slv: lat = ',lat_idx,' vmr_HNO3_gas(',i,',',m,') = ',vmr_HNO3_gas(i,k),' mlc mlc-1'
         end do                 ! end loop over lon
      end do                    ! end loop over lev
#endif /* not DST_DBG */
      
#ifdef BXM
c     Update netCDF file
      fl_out='aer.nc'           ! [sng] Name of netCDF output file
      call strnul(fl_out)
      call chm2nc(            
     $     fl_out,              ! I [sng] Name of netCDF output file
     $     lat_idx,             ! I [idx] Model latitude index
     $     rxr_HNO3_gas_dst_vmr, ! I [mlc mlc-1 s-1] Mean rate of HNO3 removal by dust
     $     rxrc_HNO3_dst,       ! I [s-1] Pseudo first order rate coefficient for HNO3
     $     rxrc_HO2_dst,        ! I [s-1] Pseudo first order rate coefficient for HO2
     $     rxrc_N2O5_dst,       ! I [s-1] Pseudo first order rate coefficient for N2O5
     $     rxrc_O3_dst,         ! I [s-1] Pseudo first order rate coefficient for O3
     $     rxrc_SO2_dst,        ! I [s-1] Pseudo first order rate coefficient for SO2
     $     vmr_HNO3_gas,        ! I [mlc mlc-1] Gaseous HNO3 volume mixing ratio
     $     vmr_HO2_gas,         ! I [mlc mlc-1] Gaseous HO2 volume mixing ratio
     $     vmr_N2O5_gas,        ! I [mlc mlc-1] Gaseous N2O5 volume mixing ratio
     $     vmr_NO3_aer,         ! I [mlc mlc-1] Particulate NO3 volume mixing ratio
     $     vmr_O3_gas,          ! I [mlc mlc-1] Gaseous O3 volume mixing ratio
     $     vmr_SO2_gas,         ! I [mlc mlc-1] Gaseous SO2 volume mixing ratio
     $     vmr_SO4_aer)         ! I [mlc mlc-1] Particulate SO4 volume mixing ratio
#endif /* not BXM */
      
      return
      end subroutine dst_chm_slv
      
      subroutine chm2nc(            
     $     fl_out,              ! I [sng] Name of netCDF output file
     $     lat_idx,             ! I [idx] Model latitude index
     $     rxr_HNO3_gas_dst_vmr, ! I [mlc mlc-1 s-1] Mean rate of HNO3 removal by dust
     $     rxrc_HNO3_dst,       ! I [s-1] Pseudo first order rate coefficient for HNO3
     $     rxrc_HO2_dst,        ! I [s-1] Pseudo first order rate coefficient for HO2
     $     rxrc_N2O5_dst,       ! I [s-1] Pseudo first order rate coefficient for N2O5
     $     rxrc_O3_dst,         ! I [s-1] Pseudo first order rate coefficient for O3
     $     rxrc_SO2_dst,        ! I [s-1] Pseudo first order rate coefficient for SO2
     $     vmr_HNO3_gas,        ! I [mlc mlc-1] Gaseous HNO3 volume mixing ratio
     $     vmr_HO2_gas,         ! I [mlc mlc-1] Gaseous HO2 volume mixing ratio
     $     vmr_N2O5_gas,        ! I [mlc mlc-1] Gaseous N2O5 volume mixing ratio
     $     vmr_NO3_aer,         ! I [mlc mlc-1] Particulate NO3 volume mixing ratio
     $     vmr_O3_gas,          ! I [mlc mlc-1] Gaseous O3 volume mixing ratio
     $     vmr_SO2_gas,         ! I [mlc mlc-1] Gaseous SO2 volume mixing ratio
     $     vmr_SO4_aer          ! I [mlc mlc-1] Particulate SO4 volume mixing ratio
     $     )                    
c     Purpose: Output aerosol chemistry properties to netCDF file
      use precision
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstchm.h> /* Chemical properties of dust */
#include <dstctl.h> /* Model control variables */
#include <dbg.com> /* Debugging constants */
#include <netcdf.inc>
c     Input
      character*(*),intent(in)::fl_out ! [sng] Name of netCDF output file
      integer,intent(in)::lat_idx ! I [idx] Model latitude index
      real(r8),intent(in)::vmr_HNO3_gas(plond,plev) ! [s-1] Gaseous HNO3 volume mixing ratio
      real(r8),intent(in)::vmr_HO2_gas(plond,plev) ! [s-1] Gaseous HO2 volume mixing ratio
      real(r8),intent(in)::vmr_N2O5_gas(plond,plev) ! [s-1] Gaseous N2O5 volume mixing ratio
      real(r8),intent(in)::vmr_O3_gas(plond,plev) ! [s-1] Gaseous O3 volume mixing ratio
      real(r8),intent(in)::vmr_SO2_gas(plond,plev) ! [s-1] Gaseous SO2 volume mixing ratio
      real(r8),intent(in)::vmr_NO3_aer(plond,plev) ! [s-1] Particulate NO3 volume mixing ratio
      real(r8),intent(in)::vmr_SO4_aer(plond,plev) ! [s-1] Particulate SO4 volume mixing ratio
      real(r8),intent(in)::rxr_HNO3_gas_dst_vmr(plond,plev) ! [s-1] Mean rate of HNO3 removal by dust
      real(r8),intent(in)::rxrc_HNO3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HNO3
      real(r8),intent(in)::rxrc_HO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HO2
      real(r8),intent(in)::rxrc_N2O5_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for N2O5
      real(r8),intent(in)::rxrc_O3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for O3
      real(r8),intent(in)::rxrc_SO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for SO2
c     Output
c     Local
c     File metadata and dimension IDs
      integer cnt_lon_lev_sz_time(4) ! Count array
      integer cnt_lon_lev_time(3) ! Count array
      integer dim_lon_lev_sz_time(4) ! [enm] Dimension IDs
      integer dim_lon_lev_time(3) ! [enm] Dimension IDs
      integer fll_mode_old      ! Old fill mode
      integer lev_dim_id        ! [enm] Dimension ID for lev
      integer lon_dim_id        ! [enm] Dimension ID for lon
      integer nc_id             ! File handle
      integer rcd               ! [rcd] Return success code
      integer srt_lon_lev_sz_time(4) ! Starting index array
      integer srt_lon_lev_time(3) ! Starting index array
      integer sz_dim_id         ! [enm] Dimension ID for sz
      integer sz_grd_dim_id     ! [enm] Dimension ID for sz grid
      integer time_dim_id       ! [enm] Dimension ID for time
c     Variable IDs
      integer vmr_HNO3_gas_id   ! [enm] Variable ID
      integer vmr_HO2_gas_id    ! [enm] Variable ID
      integer vmr_N2O5_gas_id   ! [enm] Variable ID
      integer vmr_O3_gas_id     ! [enm] Variable ID
      integer vmr_SO2_gas_id    ! [enm] Variable ID
      integer vmr_NO3_aer_id    ! [enm] Variable ID
      integer vmr_SO4_aer_id    ! [enm] Variable ID
      integer rxr_HNO3_gas_dst_vmr_id ! [enm] Variable ID
      integer rxrc_HNO3_dst_id  ! [enm] Variable ID
      integer rxrc_HO2_dst_id   ! [enm] Variable ID
      integer rxrc_N2O5_dst_id  ! [enm] Variable ID
      integer rxrc_O3_dst_id    ! [enm] Variable ID
      integer rxrc_SO2_dst_id   ! [enm] Variable ID
      integer upt_cff_H2O2_dst_id ! [enm] Variable ID
      integer upt_cff_HNO3_dst_id ! [enm] Variable ID
      integer upt_cff_HO2_dst_id ! [enm] Variable ID
      integer upt_cff_N2O5_dst_id ! [enm] Variable ID
      integer upt_cff_NO3_dst_id ! [enm] Variable ID
      integer upt_cff_O3_dst_id ! [enm] Variable ID
      integer upt_cff_OH_dst_id ! [enm] Variable ID
      integer upt_cff_SO2_dst_id ! [enm] Variable ID
c     Variable data
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
      
c     Begin netCDF output routines
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf_open(fl_out,NF_WRITE,nc_id)
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR opening ',fl_out(1:f_strlen(fl_out)),
     $     ' in chm2nc(), rcd = ',rcd
      rcd=rcd+nf_redef(nc_id)
      rcd=rcd+nf_inq_dimid(nc_id,'lon',lon_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'lev',lev_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'time',time_dim_id)
      rcd=rcd+nf_inq_dimid(nc_id,'sz',sz_dim_id)
c     Add global attributes
c     Define dimension IDs
c     Assemble ID and count vectors for each multidimensional combination of dimensions
      dim_lon_lev_time=(/lon_dim_id,lev_dim_id,time_dim_id/)
      cnt_lon_lev_time=(/plon,plev,1/)
      srt_lon_lev_time=(/1,1,nstep/)

      dim_lon_lev_sz_time=(/lon_dim_id,lev_dim_id,sz_dim_id,time_dim_id/)
      cnt_lon_lev_sz_time=(/plon,plev,dst_nbr,1/)
      srt_lon_lev_sz_time=(/1,1,1,nstep/)
      
      if (nstep == 1) then
c     Variable definitions
         rcd=rcd+nf_def_var(nc_id,'vmr_HNO3_gas',NF_DOUBLE,3,dim_lon_lev_time,vmr_HNO3_gas_id)
         rcd=rcd+nf_def_var(nc_id,'vmr_HO2_gas',NF_DOUBLE,3,dim_lon_lev_time,vmr_HO2_gas_id)
         rcd=rcd+nf_def_var(nc_id,'vmr_N2O5_gas',NF_DOUBLE,3,dim_lon_lev_time,vmr_N2O5_gas_id)
         rcd=rcd+nf_def_var(nc_id,'vmr_O3_gas',NF_DOUBLE,3,dim_lon_lev_time,vmr_O3_gas_id)
         rcd=rcd+nf_def_var(nc_id,'vmr_SO2_gas',NF_DOUBLE,3,dim_lon_lev_time,vmr_SO2_gas_id)
         rcd=rcd+nf_def_var(nc_id,'vmr_NO3_aer',NF_DOUBLE,3,dim_lon_lev_time,vmr_NO3_aer_id)
         rcd=rcd+nf_def_var(nc_id,'vmr_SO4_aer',NF_DOUBLE,3,dim_lon_lev_time,vmr_SO4_aer_id)
         rcd=rcd+nf_def_var(nc_id,'rxr_HNO3_gas_dst_vmr',NF_DOUBLE,3,dim_lon_lev_time,rxr_HNO3_gas_dst_vmr_id)
         rcd=rcd+nf_def_var(nc_id,'rxrc_HNO3_dst',NF_DOUBLE,3,dim_lon_lev_time,rxrc_HNO3_dst_id)
         rcd=rcd+nf_def_var(nc_id,'rxrc_HO2_dst',NF_DOUBLE,3,dim_lon_lev_time,rxrc_HO2_dst_id)
         rcd=rcd+nf_def_var(nc_id,'rxrc_N2O5_dst',NF_DOUBLE,3,dim_lon_lev_time,rxrc_N2O5_dst_id)
         rcd=rcd+nf_def_var(nc_id,'rxrc_O3_dst',NF_DOUBLE,3,dim_lon_lev_time,rxrc_O3_dst_id)
         rcd=rcd+nf_def_var(nc_id,'rxrc_SO2_dst',NF_DOUBLE,3,dim_lon_lev_time,rxrc_SO2_dst_id)
         rcd=rcd+nf_def_var(nc_id,'upt_cff_H2O2_dst',NF_DOUBLE,0,0,upt_cff_H2O2_dst_id)
         rcd=rcd+nf_def_var(nc_id,'upt_cff_HNO3_dst',NF_DOUBLE,0,0,upt_cff_HNO3_dst_id)
         rcd=rcd+nf_def_var(nc_id,'upt_cff_HO2_dst',NF_DOUBLE,0,0,upt_cff_HO2_dst_id)
         rcd=rcd+nf_def_var(nc_id,'upt_cff_N2O5_dst',NF_DOUBLE,0,0,upt_cff_N2O5_dst_id)
         rcd=rcd+nf_def_var(nc_id,'upt_cff_NO3_dst',NF_DOUBLE,0,0,upt_cff_NO3_dst_id)
         rcd=rcd+nf_def_var(nc_id,'upt_cff_O3_dst',NF_DOUBLE,0,0,upt_cff_O3_dst_id)
         rcd=rcd+nf_def_var(nc_id,'upt_cff_OH_dst',NF_DOUBLE,0,0,upt_cff_OH_dst_id)
         rcd=rcd+nf_def_var(nc_id,'upt_cff_SO2_dst',NF_DOUBLE,0,0,upt_cff_SO2_dst_id)
c     Add english text descriptions
         rcd=rcd+nf_put_att_text(nc_id,vmr_HNO3_gas_id,'long_name',32,'Gaseous HNO3 volume mixing ratio')
         rcd=rcd+nf_put_att_text(nc_id,vmr_HO2_gas_id,'long_name',32,'Gaseous HO2 volume mixing ratio')
         rcd=rcd+nf_put_att_text(nc_id,vmr_N2O5_gas_id,'long_name',32,'Gaseous N2O5 volume mixing ratio')
         rcd=rcd+nf_put_att_text(nc_id,vmr_O3_gas_id,'long_name',32,'Gaseous O3 volume mixing ratio')
         rcd=rcd+nf_put_att_text(nc_id,vmr_SO2_gas_id,'long_name',32,'Gaseous SO2 volume mixing ratio')
         rcd=rcd+nf_put_att_text(nc_id,vmr_NO3_aer_id,'long_name',35,'Particulate NO3 volume mixing ratio')
         rcd=rcd+nf_put_att_text(nc_id,vmr_SO4_aer_id,'long_name',35,'Particulate SO4 volume mixing ratio')
         rcd=rcd+nf_put_att_text(nc_id,rxr_HNO3_gas_dst_vmr_id,'long_name',33,'Mean rate of HNO3 removal by dust')
         rcd=rcd+nf_put_att_text(nc_id,rxrc_HNO3_dst_id,'long_name',44,'Pseudo first order rate coefficient for HNO3')
         rcd=rcd+nf_put_att_text(nc_id,rxrc_HO2_dst_id,'long_name',43,'Pseudo first order rate coefficient for HO2')
         rcd=rcd+nf_put_att_text(nc_id,rxrc_N2O5_dst_id,'long_name',44,'Pseudo first order rate coefficient for N2O5')
         rcd=rcd+nf_put_att_text(nc_id,rxrc_O3_dst_id,'long_name',42,'Pseudo first order rate coefficient for O3')
         rcd=rcd+nf_put_att_text(nc_id,rxrc_SO2_dst_id,'long_name',43,'Pseudo first order rate coefficient for SO2')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_H2O2_dst_id,'long_name',35,'Uptake coefficient for H2O2 to dust')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_HNO3_dst_id,'long_name',35,'Uptake coefficient for HNO3 to dust')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_HO2_dst_id,'long_name',34,'Uptake coefficient for HO2 to dust')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_N2O5_dst_id,'long_name',35,'Uptake coefficient for N2O5 to dust')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_NO3_dst_id,'long_name',34,'Uptake coefficient for NO3 to dust')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_O3_dst_id,'long_name',33,'Uptake coefficient for O3 to dust')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_OH_dst_id,'long_name',33,'Uptake coefficient for OH to dust')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_SO2_dst_id,'long_name',34,'Uptake coefficient for SO2 to dust')
c     Add units
         rcd=rcd+nf_put_att_text(nc_id,vmr_HNO3_gas_id,'units',19,'molecule molecule-1')
         rcd=rcd+nf_put_att_text(nc_id,vmr_HO2_gas_id,'units',19,'molecule molecule-1')
         rcd=rcd+nf_put_att_text(nc_id,vmr_N2O5_gas_id,'units',19,'molecule molecule-1')
         rcd=rcd+nf_put_att_text(nc_id,vmr_O3_gas_id,'units',19,'molecule molecule-1')
         rcd=rcd+nf_put_att_text(nc_id,vmr_SO2_gas_id,'units',19,'molecule molecule-1')
         rcd=rcd+nf_put_att_text(nc_id,vmr_NO3_aer_id,'units',19,'molecule molecule-1')
         rcd=rcd+nf_put_att_text(nc_id,vmr_SO4_aer_id,'units',19,'molecule molecule-1')
         rcd=rcd+nf_put_att_text(nc_id,rxr_HNO3_gas_dst_vmr_id,'units',28,'molecule molecule-1 second-1')
         rcd=rcd+nf_put_att_text(nc_id,rxrc_HNO3_dst_id,'units',8,'second-1')
         rcd=rcd+nf_put_att_text(nc_id,rxrc_HO2_dst_id,'units',8,'second-1')
         rcd=rcd+nf_put_att_text(nc_id,rxrc_N2O5_dst_id,'units',8,'second-1')
         rcd=rcd+nf_put_att_text(nc_id,rxrc_O3_dst_id,'units',8,'second-1')
         rcd=rcd+nf_put_att_text(nc_id,rxrc_SO2_dst_id,'units',8,'second-1')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_H2O2_dst_id,'units',8,'fraction')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_HNO3_dst_id,'units',8,'fraction')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_HO2_dst_id,'units',8,'fraction')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_N2O5_dst_id,'units',8,'fraction')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_NO3_dst_id,'units',8,'fraction')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_O3_dst_id,'units',8,'fraction')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_OH_dst_id,'units',8,'fraction')
         rcd=rcd+nf_put_att_text(nc_id,upt_cff_SO2_dst_id,'units',8,'fraction')
      else                      ! endif nstep == 1
c     Get variable IDs
         rcd=rcd+nf_inq_varid(nc_id,'vmr_HNO3_gas',vmr_HNO3_gas_id)
         rcd=rcd+nf_inq_varid(nc_id,'vmr_HO2_gas',vmr_HO2_gas_id)
         rcd=rcd+nf_inq_varid(nc_id,'vmr_N2O5_gas',vmr_N2O5_gas_id)
         rcd=rcd+nf_inq_varid(nc_id,'vmr_O3_gas',vmr_O3_gas_id)
         rcd=rcd+nf_inq_varid(nc_id,'vmr_SO2_gas',vmr_SO2_gas_id)
         rcd=rcd+nf_inq_varid(nc_id,'vmr_NO3_aer',vmr_NO3_aer_id)
         rcd=rcd+nf_inq_varid(nc_id,'vmr_SO4_aer',vmr_SO4_aer_id)
         rcd=rcd+nf_inq_varid(nc_id,'rxr_HNO3_gas_dst_vmr',rxr_HNO3_gas_dst_vmr_id)
         rcd=rcd+nf_inq_varid(nc_id,'rxrc_HNO3_dst',rxrc_HNO3_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'rxrc_HO2_dst',rxrc_HO2_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'rxrc_N2O5_dst',rxrc_N2O5_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'rxrc_O3_dst',rxrc_O3_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'rxrc_SO2_dst',rxrc_SO2_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'upt_cff_H2O2_dst',upt_cff_H2O2_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'upt_cff_HNO3_dst',upt_cff_HNO3_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'upt_cff_HO2_dst',upt_cff_HO2_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'upt_cff_N2O5_dst',upt_cff_N2O5_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'upt_cff_NO3_dst',upt_cff_NO3_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'upt_cff_O3_dst',upt_cff_O3_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'upt_cff_OH_dst',upt_cff_OH_dst_id)
         rcd=rcd+nf_inq_varid(nc_id,'upt_cff_SO2_dst',upt_cff_SO2_dst_id)
      endif                     ! endif nstep /= 1
c     End define mode
      rcd=rcd+nf_enddef(nc_id)
c     Write data
      if (nstep == 1) then
         rcd=rcd+nf_put_var_double(nc_id,upt_cff_H2O2_dst_id,upt_cff_H2O2_dst)
         rcd=rcd+nf_put_var_double(nc_id,upt_cff_HNO3_dst_id,upt_cff_HNO3_dst)
         rcd=rcd+nf_put_var_double(nc_id,upt_cff_HO2_dst_id,upt_cff_HO2_dst)
         rcd=rcd+nf_put_var_double(nc_id,upt_cff_N2O5_dst_id,upt_cff_N2O5_dst)
         rcd=rcd+nf_put_var_double(nc_id,upt_cff_NO3_dst_id,upt_cff_NO3_dst)
         rcd=rcd+nf_put_var_double(nc_id,upt_cff_O3_dst_id,upt_cff_O3_dst)
         rcd=rcd+nf_put_var_double(nc_id,upt_cff_OH_dst_id,upt_cff_OH_dst)
         rcd=rcd+nf_put_var_double(nc_id,upt_cff_SO2_dst_id,upt_cff_SO2_dst)
      endif                     ! endif nstep == 1
      rcd=rcd+nf_put_vara_double(nc_id,vmr_HNO3_gas_id,srt_lon_lev_time,cnt_lon_lev_time,vmr_HNO3_gas)
      rcd=rcd+nf_put_vara_double(nc_id,vmr_HO2_gas_id,srt_lon_lev_time,cnt_lon_lev_time,vmr_HO2_gas)
      rcd=rcd+nf_put_vara_double(nc_id,vmr_N2O5_gas_id,srt_lon_lev_time,cnt_lon_lev_time,vmr_N2O5_gas)
      rcd=rcd+nf_put_vara_double(nc_id,vmr_O3_gas_id,srt_lon_lev_time,cnt_lon_lev_time,vmr_O3_gas)
      rcd=rcd+nf_put_vara_double(nc_id,vmr_SO2_gas_id,srt_lon_lev_time,cnt_lon_lev_time,vmr_SO2_gas)
      rcd=rcd+nf_put_vara_double(nc_id,vmr_NO3_aer_id,srt_lon_lev_time,cnt_lon_lev_time,vmr_NO3_aer)
      rcd=rcd+nf_put_vara_double(nc_id,vmr_SO4_aer_id,srt_lon_lev_time,cnt_lon_lev_time,vmr_SO4_aer)
      rcd=rcd+nf_put_vara_double(nc_id,rxr_HNO3_gas_dst_vmr_id,srt_lon_lev_time,cnt_lon_lev_time,rxr_HNO3_gas_dst_vmr)
      rcd=rcd+nf_put_vara_double(nc_id,rxrc_HNO3_dst_id,srt_lon_lev_time,cnt_lon_lev_time,rxrc_HNO3_dst)
      rcd=rcd+nf_put_vara_double(nc_id,rxrc_HO2_dst_id,srt_lon_lev_time,cnt_lon_lev_time,rxrc_HO2_dst)
      rcd=rcd+nf_put_vara_double(nc_id,rxrc_N2O5_dst_id,srt_lon_lev_time,cnt_lon_lev_time,rxrc_N2O5_dst)
      rcd=rcd+nf_put_vara_double(nc_id,rxrc_O3_dst_id,srt_lon_lev_time,cnt_lon_lev_time,rxrc_O3_dst)
      rcd=rcd+nf_put_vara_double(nc_id,rxrc_SO2_dst_id,srt_lon_lev_time,cnt_lon_lev_time,rxrc_SO2_dst)
c     Close output file
      rcd=rcd+nf_close(nc_id)
      if (nstep == 1) then
         write (6,'(a,a40,a)') prg_nm(1:f_strlen(prg_nm)),': Initialized chemistry data archive in ',fl_out(1:f_strlen(fl_out))
      endif                     ! endif
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR writing ',fl_out(1:f_strlen(fl_out)),
     $     ' in chm2nc(), rcd = ',rcd
      return 
      end subroutine chm2nc
      
