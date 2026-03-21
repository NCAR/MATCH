c     $Header$ -*-fortran-*-

c     Purpose: Common block dst_cst_cmn stores fundamental physical constants
c     needed by dust routines.

c     These variables are initialized in dstmss:dst_cst_cmn_ini()
c     dst_cst_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun() 
c     CCM:physics/inti() passes values from CCM:physics/comcon.h which were set by CCM:eul/initcom()

c     Usage: 
c     #include <dstcst.h> /* Physical constants for dust routines */

      real(r8) cst_von_krm          ! [frc] Von Karman's constant
      real(r8) cst_von_krm_rcp      ! [frc] Reciprocal of Von Karman's constant
      real(r8) eps_H2O              ! (0.622) [frc] mmw(H2O)/mmw(dry air)
      real(r8) eps_H2O_rcp_m1       ! [frc] (0.60777) Constant for virtual temperature
      real(r8) gas_cst_dry_air      ! [J kg-1 K-1] (287.05) Gas constant of dry air
      real(r8) gas_cst_unv          ! [J mol-1 K-1] Universal gas constant
      real(r8) grv_sfc              ! [m s-2] (9.80616) Gravity (mean surface)
      real(r8) grv_sfc_rcp          ! [s2 m-1] (0.101977) Reciprocal of gravity
      real(r8) kappa_dry_air        ! (0.286 = 2/7) [frc] Constant in potential temperature
      real(r8) mmw_H2O              ! [kg mol-1] Mean molecular weight of water
      real(r8) mmw_dry_air          ! [kg mol-1] Mean molecular weight of dry air
      real(r8) rds_Earth            ! [m] (6.37122e+6) Radius of sphere of same volume as Earth

      common /dst_cst_cmn/      ! dst_cst_cmn is initialized in dst_cst_cmn_ini()
     $     cst_von_krm,         ! [frc] Von Karman's constant
     $     cst_von_krm_rcp,     ! [frc] Reciprocal of Von Karman's constant
     $     eps_H2O,             ! (0.622) [frc] mmw(H2O)/mmw(dry air)
     $     eps_H2O_rcp_m1,      ! [frc] (0.60777) Constant for virtual temperature
     $     gas_cst_dry_air,     ! [J kg-1 K-1] (287.05) Gas constant of dry air
     $     gas_cst_unv,         ! [J mol-1 K-1] Universal gas constant
     $     grv_sfc,             ! [m s-2] (9.80616) Gravity (mean surface)
     $     grv_sfc_rcp,         ! [s2 m-1] (0.101977) Reciprocal of gravity
     $     kappa_dry_air,       ! (0.286 = 2/7) [frc] Constant in potential temperature
     $     mmw_H2O,             ! [kg mol-1] Mean molecular weight of water
     $     mmw_dry_air,         ! [kg mol-1] Mean molecular weight of dry air
     $     rds_Earth            ! [m] (6.37122e+6) Radius of sphere of same volume as Earth
