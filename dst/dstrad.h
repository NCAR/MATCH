c     $Header$ -*-fortran-*-

c     Purpose: Common block dstrad stores optical properties needed by dust parameterization
c     These common blocks are initialized in dstradbd() which is called by CCM:main()

c     Usage:
c#include <dstrad.h> /* Dust radiative properties */

c     Microphysical optical properties
      real(r8) abs_cff_mss_dst(bnd_nbr_LW,dst_nbr) ! [m2 kg-1] Mass absorption coefficient of dust
      real(r8) asm_prm_dst(bnd_nbr_SW,dst_nbr) ! [frc] Asymmetry parameter of dust
      real(r8) ext_cff_mss_dst(bnd_nbr_SW,dst_nbr) ! [m2 kg-1] Mass extinction coefficient of dust
      real(r8) ss_alb_dst(bnd_nbr_SW,dst_nbr) ! [frc] Single scattering albedo of dust
      common / dstrad /         ! dstrad is currently initialized in dstradbd()
     $     abs_cff_mss_dst,     ! [m2 kg-1] Mass absorption coefficient of dust
     $     asm_prm_dst,         ! [frc] Asymmetry parameter of dust
     $     ext_cff_mss_dst,     ! [m2 kg-1] Mass extinction coefficient of dust
     $     ss_alb_dst           ! [frc] Single scattering albedo of dust

      real(r8) asm_prm_dst_GSFC_LW(bnd_nbr_GSFC_LW,dst_nbr) ! [frc] Asymmetry parameter of dust
      real(r8) ext_cff_mss_dst_GSFC_LW(bnd_nbr_GSFC_LW,dst_nbr) ! [m2 kg-1] Mass extinction coefficient of dust
      real(r8) ss_alb_dst_GSFC_LW(bnd_nbr_GSFC_LW,dst_nbr) ! [frc] Single scattering albedo of dust
      common / dstrad_GSFC_LW / ! dstrad_GSFC_LW is currently initialized in dstradbd()
     $     asm_prm_dst_GSFC_LW, ! [frc] Asymmetry parameter of dust
     $     ext_cff_mss_dst_GSFC_LW, ! [m2 kg-1] Mass extinction coefficient of dust
     $     ss_alb_dst_GSFC_LW   ! [frc] Single scattering albedo of dust







