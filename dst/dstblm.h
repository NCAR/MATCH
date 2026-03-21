c     $Header$ -*-fortran-*-
      
c     Purpose: Common block dst_blm_cmn stores fundamental physical constants
c     needed for boundary layer meteorology
      
c     These variables are initialized in dstcmnini:dst_cst_cmn_ini()
c     dst_cst_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun() 
      
c     Usage: 
c     #include <dstblm.h> /* Physical constants for boundary layer meteorology */
      
      real(r8) cp_vpr_rcp_cp_dry_m1 ! 0.83745 [frc] Constant for moist specific heat IrG81 p. 77
      real(r8) cst_Stefan_Boltzmann ! (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant
      real(r8) dns_H2O_lqd_std      ! [kg m-3] Density of liquid water
      real(r8) ltn_heat_fsn_H2O_std ! [J kg-1] Latent heat of fusion of H2O at 0 C, standard
      real(r8) ltn_heat_sbl_H2O_std ! [J kg-1] Latent heat of sublimation of H2O, standard
      real(r8) ltn_heat_vpr_H2O_std ! [J kg-1] Latent heat of vaporization of H2O, standard
      real(r8) one_mns_eps_H2O          ! (0.378) [frc] Constant for saturation specific humidity
      real(r8) spc_heat_H2O_ice_vlm ! [J m-3 K-1] Volumetric specific heat capacity of ice water
      real(r8) spc_heat_H2O_lqd_vlm ! [J m-3 K-1] Volumetric specific heat capacity of liquid water
      real(r8) spc_heat_dry_air     ! (1005.0) [J kg-1 K-1] IrG81 p. 25
      real(r8) tpt_frz_pnt          ! [K] Kelvin--Celsius scale offset
      common /dst_blm_cmn/      ! dst_blm_cmn is initialized in dst_cst_cmn_ini()
     $     cp_vpr_rcp_cp_dry_m1, ! 0.83745 [frc] Constant for moist specific heat IrG81 p. 77
     $     cst_Stefan_Boltzmann, ! (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant
     $     dns_H2O_lqd_std,     ! [kg m-3] Density of liquid water
     $     ltn_heat_fsn_H2O_std, ! [J kg-1] Latent heat of fusion of H2O at 0 C, standard
     $     ltn_heat_sbl_H2O_std, ! [J kg-1] Latent heat of sublimation of H2O, standard
     $     ltn_heat_vpr_H2O_std, ! [J kg-1] Latent heat of vaporization of H2O, standard
     $     one_mns_eps_H2O,         ! (0.378) [frc] Constant for saturation specific humidity
     $     spc_heat_H2O_ice_vlm, ! [J m-3 K-1] Volumetric specific heat capacity of ice water
     $     spc_heat_H2O_lqd_vlm, ! [J m-3 K-1] Volumetric specific heat capacity of liquid water
     $     spc_heat_dry_air,    ! (1005.0) [J kg-1 K-1] IrG81 p. 25
     $     tpt_frz_pnt          ! [K] Kelvin--Celsius scale offset
      
      
