      module blmutl

      use precision
      use dstmblutl, only: snw_frc_get,wnd_rfr_get

      real(r8) :: half = 0.5
      real(r8) :: fifty = 50.
      real(r8) :: ten = 10.
      real(r8) :: one = 1.


      contains

c     $Header$
      
c     Purpose: blmutl.F contains boundary layer meteorology utilities
      
      subroutine rgh_zpd_get(
     $     hgt_zpd,             ! O [m] Zero plane displacement
     $     oro,                 ! I [frc] Orography
     $     rgh_mmn,             ! O [m] Roughness length momentum
     $     sfc_typ,             ! I [idx] LSM surface type (0..28)
     $     snw_frc,             ! I [frc] Fraction of surface covered by snow
     $     wnd_10m)             ! I [m s-1] 10 m wind speed
c     Purpose: Set roughness length and zero plane displacement 
c     NB: SeP97 apparently uses rlm but rlh may be more correct for aerosol
c     NB: Currently all lakes are treated as unfrozen
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstlsm.h> /* LSM block data */
c     Parameters
      real(r8),parameter::rgh_mmn_ice_lak=0.04 ! [m] Roughness length over frozen lakes Bon96 p. 59
      real(r8),parameter::rgh_mmn_ice_lnd=0.05 ! [m] Roughness length over ice, bare ground, wetlands Bon96 p. 59
      real(r8),parameter::rgh_mmn_ice_ocn=0.0005 ! [m] Roughness length over sea ice BKL97 p. F-3
      real(r8),parameter::rgh_mmn_lak_wrm=0.001 ! [m] Roughness length over unfrozen lakes Bon96 p. 59
      real(r8),parameter::rgh_mmn_snw=0.04 ! [m] Roughness length over snow Bon96 p. 59 CCM:lsm/snoconi.F
      real(r8),parameter::wnd_min_dps=1.0 ! [m s-1] Minimum windspeed for momentum exchange
c     Input
      integer sfc_typ(plond)    ! [idx] LSM surface type (0..28)
      real(r8) oro(plond)           ! [frc] Orography
      real(r8) snw_frc(plond)       ! [frc] Fraction of surface covered by snow
      real(r8) wnd_10m(plond)       ! [m s-1] 10 m wind speed
c     Output
      real(r8) hgt_zpd(plond)       ! [m] Zero plane displacement
      real(r8) rgh_mmn(plond)       ! [m] Roughness length momentum
c     Local
      integer ice_idx(plond)    ! [idx] Longitude index array (sea ice)
      integer ice_nbr           ! [nbr] Number of sea ice points
      integer idx_idx           ! [idx] Counting index
      integer lnd_idx(plond)    ! [idx] Longitude index array (land)
      integer lnd_nbr           ! [nbr] Number of land points
      integer lon_idx           ! [idx] Counting index
      integer ocn_idx(plond)    ! [idx] Longitude index array (ocean)
      integer ocn_nbr           ! [nbr] Number of ocean points
      integer pln_typ_idx       ! [idx] Plant type index
      integer sfc_typ_idx       ! [idx] Surface type index
      integer sgs_idx           ! [idx] Surface sub-gridscale index
      real(r8) rlm_crr              ! [m] Roughness length of current sub-gridscale
      real(r8) wnd_10m_bnd          ! [m s-1] Bounded wind speed at 10 m
      real(r8) wnd_10m_ntr          ! [m s-1] Wind speed at 10 m
      real(r8) xch_cff_mmn_ocn_ntr  ! [frc] Neutral 10 m drag coefficient over ocean
c     Statement functions
      logical oro_is_ice        ! [flg] True if > 50% sea ice
      logical oro_is_lnd        ! [flg] True if > 50% land
      logical oro_is_ocn        ! [flg] True if > 50% ocean
      real(r8) xch_cff_mmn_ocn_ntr_get ! [frc] Neutral 10 m drag coefficient over ocean
      real(r8) oro_val              ! [frc] Orography
      oro_is_ice(oro_val)=nint(oro_val) == 2
      oro_is_lnd(oro_val)=nint(oro_val) == 1
      oro_is_ocn(oro_val)=nint(oro_val) == 0
      xch_cff_mmn_ocn_ntr_get(wnd_10m_ntr)=0.0027/wnd_10m_ntr+0.000142+0.0000764*wnd_10m_ntr ! [frc] LaP82 CCM:dom/flxoce(), NOS97 p. I2
      
c     Main Code
      
c     Initialize array
      rgh_mmn=0.0               ! [m]
      hgt_zpd=0.0               ! [m]
      
c     Initialize ocean, sea-ice, and land vectors
      ocn_nbr=0
      do lon_idx=1,plon
         if (oro_is_ocn(oro(lon_idx))) then
            ocn_nbr=ocn_nbr+1
            ocn_idx(ocn_nbr)=lon_idx
         end if                 ! endif
      enddo                     ! end loop over lon
      ice_nbr=0
      do lon_idx=1,plon
         if (oro_is_ice(oro(lon_idx))) then
            ice_nbr=ice_nbr+1
            ice_idx(ice_nbr)=lon_idx
         end if                 ! endif
      enddo                     ! end loop over lon
      lnd_nbr=0
      do lon_idx=1,plon
         if (oro_is_lnd(oro(lon_idx))) then
            lnd_nbr=lnd_nbr+1
            lnd_idx(lnd_nbr)=lon_idx
         end if                 ! endif
      enddo                     ! end loop over lon
      
c     Ocean points
      do idx_idx=1,ocn_nbr
         lon_idx=ocn_idx(idx_idx)
c     Convert wind speed to roughness length over ocean
         wnd_10m_bnd=max(wnd_min_dps,wnd_10m(lon_idx)) ! [m s-1]
c     Approximation: neutral 10 m wind speed unavailable, use 10 m wind speed
         xch_cff_mmn_ocn_ntr=xch_cff_mmn_ocn_ntr_get(wnd_10m_bnd) ! [frc]
         rgh_mmn(lon_idx)=10.0*exp(-cst_von_krm/sqrt(xch_cff_mmn_ocn_ntr)) ! [m] BKL97 p. F-4, LaP81 p. 327 (14) 
         hgt_zpd(lon_idx)=0.0   ! [m]
      enddo                     ! end loop over ocn lon
      
c     Sea ice points
      do idx_idx=1,ice_nbr
         lon_idx=ice_idx(idx_idx)
         rgh_mmn(lon_idx)=snw_frc(lon_idx)*rgh_mmn_snw+(1.0-snw_frc(lon_idx))*rgh_mmn_ice_ocn ! [m] Bon96 p. 59
         hgt_zpd(lon_idx)=0.0   ! [m]
      enddo                     ! end loop over ice lon
      
c     Land points
      do idx_idx=1,lnd_nbr
         lon_idx=lnd_idx(idx_idx)
c     Store surface blend for current gridpoint
         sfc_typ_idx=sfc_typ(lon_idx)
         if (sfc_typ_idx == 0) then ! Inland lake
c     fxm: Need to add temperature input and so ability to discriminate warm from frozen lakes here
            rgh_mmn(lon_idx)=rgh_mmn_lak_wrm ! [m] Bon96 p. 59
            hgt_zpd(lon_idx)=0.0 ! [m]
         else if(sfc_typ_idx == 1) then ! Land ice
            rgh_mmn(lon_idx)=snw_frc(lon_idx)*rgh_mmn_snw+(1.0-snw_frc(lon_idx))*rgh_mmn_ice_lnd ! [m] Bon96 p. 59
            hgt_zpd(lon_idx)=0.0 ! [m]
         else                   ! Normal land
            do sgs_idx=1,3
c     Bare ground is pln_typ=14, ocean is pln_typ=0
               pln_typ_idx=pln_typ(sfc_typ_idx,sgs_idx)
               if (pln_typ_idx == 14) then ! Bare ground
                  rlm_crr=snw_frc(lon_idx)*rgh_mmn_snw+(1.0-snw_frc(lon_idx))*rgh_mmn_ice_lnd ! [m] Bon96 p. 59 (glacial ice is same as bare ground)
               else if (pln_typ_idx > 0) then ! Regular plant type
                  rlm_crr=snw_frc(lon_idx)*rgh_mmn_snw+(1.0-snw_frc(lon_idx))*z0mvt(pln_typ_idx) ! [m] Bon96 p. 59
               else             ! Presumably ocean snuck through
                  stop 'dst: ERROR: pln_typ_idx == 0 in rgh_zpd_get()'
               endif            ! endif
               rgh_mmn(lon_idx)= ! [m]
     $              rgh_mmn(lon_idx)+
     $              pln_frc(sfc_typ_idx,sgs_idx)* ! [frc]
     $              rlm_crr     ! [m]
               hgt_zpd(lon_idx)= ! [m]
     $              hgt_zpd(lon_idx)+
     $              pln_frc(sfc_typ_idx,sgs_idx)* ! [frc]
     $              zpdvt(pln_typ_idx) ! [m]
            enddo               ! end loop over number of possible plant types
         endif                  ! endif normal land
      end do                    ! end loop over lon
      return
      end  subroutine rgh_zpd_get                      ! end rgh_zpd_get()
      
      subroutine rgh_mmn_get(
     $     oro,                 ! I [frc] Orography
     $     rgh_mmn,             ! O [m] Roughness length momentum
     $     sfc_typ,             ! I [idx] LSM surface type (0..28)
     $     snw_frc,             ! I [frc] Fraction of surface covered by snow
     $     wnd_10m)             ! I [m s-1] 10 m wind speed
c     Purpose: Set roughness length
c     NB: SeP97 apparently uses rlm but rlh may be more correct for aerosol
c     NB: Currently all lakes are treated as unfrozen
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstlsm.h> /* LSM block data */
c     Parameters
      real(r8),parameter::rgh_mmn_ice_lak=0.04 ! [m] Roughness length over frozen lakes Bon96 p. 59
      real(r8),parameter::rgh_mmn_ice_lnd=0.05 ! [m] Roughness length over ice, bare ground, wetlands Bon96 p. 59
      real(r8),parameter::rgh_mmn_ice_ocn=0.0005 ! [m] Roughness length over sea ice BKL97 p. F-3
      real(r8),parameter::rgh_mmn_lak_wrm=0.001 ! [m] Roughness length over unfrozen lakes Bon96 p. 59
      real(r8),parameter::rgh_mmn_snw=0.04 ! [m] Roughness length over snow Bon96 p. 59 CCM:lsm/snoconi.F
      real(r8),parameter::wnd_min_dps=1.0 ! [m s-1] Minimum windspeed for momentum exchange
c     Input
      integer sfc_typ(plond)    ! [idx] LSM surface type (0..28)
      real(r8) oro(plond)           ! [frc] Orography
      real(r8) snw_frc(plond)       ! [frc] Fraction of surface covered by snow
      real(r8) wnd_10m(plond)       ! [m s-1] 10 m wind speed
c     Output
      real(r8) rgh_mmn(plond)       ! [m] Roughness length momentum
c     Local
      integer ice_idx(plond)    ! [idx] Longitude index array (sea ice)
      integer ice_nbr           ! [nbr] Number of sea ice points
      integer idx_idx           ! [idx] Counting index
      integer lnd_idx(plond)    ! [idx] Longitude index array (land)
      integer lnd_nbr           ! [nbr] Number of land points
      integer lon_idx           ! [idx] Counting index
      integer ocn_idx(plond)    ! [idx] Longitude index array (ocean)
      integer ocn_nbr           ! [nbr] Number of ocean points
      integer pln_typ_idx       ! [idx] Plant type index
      integer sfc_typ_idx       ! [idx] Surface type index
      integer sgs_idx           ! [idx] Surface sub-gridscale index
      real(r8) rlm_crr              ! [m] Roughness length of current sub-gridscale
      real(r8) wnd_10m_bnd          ! [m s-1] Bounded wind speed at 10 m
      real(r8) wnd_10m_ntr          ! [m s-1] Wind speed at 10 m
      real(r8) xch_cff_mmn_ocn_ntr  ! [frc] Neutral 10 m drag coefficient over ocean
c     Statement functions
      logical oro_is_ice        ! [flg] True if > 50% sea ice
      logical oro_is_lnd        ! [flg] True if > 50% land
      logical oro_is_ocn        ! [flg] True if > 50% ocean
      real(r8) xch_cff_mmn_ocn_ntr_get ! [frc] Neutral 10 m drag coefficient over ocean
      real(r8) oro_val              ! [frc] Orography
      oro_is_ice(oro_val)=nint(oro_val) == 2
      oro_is_lnd(oro_val)=nint(oro_val) == 1
      oro_is_ocn(oro_val)=nint(oro_val) == 0
      xch_cff_mmn_ocn_ntr_get(wnd_10m_ntr)=0.0027/wnd_10m_ntr+0.000142+0.0000764*wnd_10m_ntr ! [frc] LaP82 CCM:dom/flxoce(), NOS97 p. I2
      
c     Main Code
      
c     Initialize array
      rgh_mmn=0.0               ! [m]
      
c     Initialize ocean, sea-ice, and land vectors
      ocn_nbr=0
      do lon_idx=1,plon
         if (oro_is_ocn(oro(lon_idx))) then
            ocn_nbr=ocn_nbr+1
            ocn_idx(ocn_nbr)=lon_idx
         end if                 ! endif
      enddo                     ! end loop over lon
      ice_nbr=0
      do lon_idx=1,plon
         if (oro_is_ice(oro(lon_idx))) then
            ice_nbr=ice_nbr+1
            ice_idx(ice_nbr)=lon_idx
         end if                 ! endif
      enddo                     ! end loop over lon
      lnd_nbr=0
      do lon_idx=1,plon
         if (oro_is_lnd(oro(lon_idx))) then
            lnd_nbr=lnd_nbr+1
            lnd_idx(lnd_nbr)=lon_idx
         end if                 ! endif
      enddo                     ! end loop over lon
      
c     Ocean points
      do idx_idx=1,ocn_nbr
         lon_idx=ocn_idx(idx_idx)
c     Convert wind speed to roughness length over ocean
         wnd_10m_bnd=max(wnd_min_dps,wnd_10m(lon_idx)) ! [m s-1]
c     Approximation: neutral 10 m wind speed unavailable, use 10 m wind speed
         xch_cff_mmn_ocn_ntr=xch_cff_mmn_ocn_ntr_get(wnd_10m_bnd) ! [frc]
         rgh_mmn(lon_idx)=10.0*exp(-cst_von_krm/sqrt(xch_cff_mmn_ocn_ntr)) ! [m] BKL97 p. F-4, LaP81 p. 327 (14) 
      enddo                     ! end loop over ocn lon
      
c     Sea ice points
      do idx_idx=1,ice_nbr
         lon_idx=ice_idx(idx_idx)
         rgh_mmn(lon_idx)=snw_frc(lon_idx)*rgh_mmn_snw+(1.0-snw_frc(lon_idx))*rgh_mmn_ice_ocn ! [m] Bon96 p. 59
      enddo                     ! end loop over ice lon
      
c     Land points
      do idx_idx=1,lnd_nbr
         lon_idx=lnd_idx(idx_idx)
c     Store surface blend for current gridpoint
         sfc_typ_idx=sfc_typ(lon_idx)
         if (sfc_typ_idx == 0) then ! Inland lake
c     fxm: Add temperature input and so ability to discriminate warm from frozen lakes here
            rgh_mmn(lon_idx)=rgh_mmn_lak_wrm ! [m] Bon96 p. 59
         else if(sfc_typ_idx == 1) then ! Land ice
            rgh_mmn(lon_idx)=snw_frc(lon_idx)*rgh_mmn_snw+(1.0-snw_frc(lon_idx))*rgh_mmn_ice_lnd ! [m] Bon96 p. 59
         else                   ! Normal land
            do sgs_idx=1,3
c     Bare ground is pln_typ=14, ocean is pln_typ=0
               pln_typ_idx=pln_typ(sfc_typ_idx,sgs_idx)
               if (pln_typ_idx == 14) then ! Bare ground
                  rlm_crr=snw_frc(lon_idx)*rgh_mmn_snw+(1.0-snw_frc(lon_idx))*rgh_mmn_ice_lnd ! [m] Bon96 p. 59 (glacial ice is same as bare ground)
               else if (pln_typ_idx > 0) then ! Regular plant type
                  rlm_crr=snw_frc(lon_idx)*rgh_mmn_snw+(1.0-snw_frc(lon_idx))*z0mvt(pln_typ_idx) ! [m] Bon96 p. 59
               else             ! Presumably ocean snuck through
                  stop 'dst: ERROR: pln_typ_idx == 0 in rgh_mmn_get()'
               endif            ! endif
               rgh_mmn(lon_idx)= ! [m]
     $              rgh_mmn(lon_idx)+
     $              pln_frc(sfc_typ_idx,sgs_idx)* ! [frc]
     $              rlm_crr     ! [m]
            enddo               ! end loop over number of possible plant types
         endif                  ! endif normal land
      end do                    ! end loop over lon
      
c     Sanity check
c     do lon_idx=1,plon
c     if (rgh_mmn(lon_idx) <= 0.0) then
c     write(6,'(a,3(a,1p,e8.1,0p),a)') 
c     $           'rgh_mmn_get: ','oro = ',oro(lon_idx),', rgh_mmn = ',rgh_mmn(lon_idx),' m, wnd_10m = ',wnd_10m(lon_idx),' m s-1'
c     endif                  ! endif err
c     enddo                     ! end loop over ocn lon
      return
      end subroutine rgh_mmn_get                      ! end rgh_mmn_get()
      
      subroutine zpd_get(
     $     hgt_zpd,             ! O [m] Zero plane displacement
     $     oro,                 ! I [frc] Orography
     $     sfc_typ,             ! I [idx] LSM surface type (0..28)
     $     snw_frc)             ! I [frc] Fraction of surface covered by snow
c     Purpose: Set zero plane displacement 
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstlsm.h> /* LSM block data */
c     Parameters
c     Input
      integer sfc_typ(plond)    ! [idx] LSM surface type (0..28)
      real(r8) oro(plond)           ! [frc] Orography
      real(r8) snw_frc(plond)       ! [frc] Fraction of surface covered by snow
c     Output
      real(r8) hgt_zpd(plond)       ! [m] Zero plane displacement
c     Local
      integer idx_idx           ! [nbr] Counting index
      integer lnd_idx(plond)    ! [nbr] Longitude index array (land)
      integer lnd_nbr           ! [nbr] Number of land points
      integer lon_idx           ! [nbr] Counting index
      integer pln_typ_idx       ! [idx] Plant type index
      integer sfc_typ_idx       ! [idx] Surface type index
      integer sgs_idx           ! [nbr] Surface sub-gridscale index
c     Statement functions
      logical oro_is_lnd        ! [flg] True if > 50% land
      real(r8) oro_val              ! [frc] Orography
      oro_is_lnd(oro_val)=nint(oro_val) == 1
      
c     Main Code
      
c     Initialize array
c     Zero plane displacement is identically 0.0 everywhere except land
      hgt_zpd=0.0               ! [m]
      
c     Land ahoy!
      lnd_nbr=0
      do lon_idx=1,plon
         if (oro_is_lnd(oro(lon_idx))) then
            lnd_nbr=lnd_nbr+1
            lnd_idx(lnd_nbr)=lon_idx
         end if                 ! endif
      enddo                     ! end loop over lon
      
c     Land points
      do idx_idx=1,lnd_nbr
         lon_idx=lnd_idx(idx_idx)
c     Store surface blend for current gridpoint
         sfc_typ_idx=sfc_typ(lon_idx)
         if (sfc_typ_idx == 0) then ! Inland lake
            hgt_zpd(lon_idx)=0.0 ! [m]
         else if(sfc_typ_idx == 1) then ! Land ice
            hgt_zpd(lon_idx)=0.0 ! [m]
         else                   ! Normal land
            do sgs_idx=1,3
c     Bare ground is pln_typ=14, ocean is pln_typ=0
               pln_typ_idx=pln_typ(sfc_typ_idx,sgs_idx)
               hgt_zpd(lon_idx)= ! [m]
     $              hgt_zpd(lon_idx)+
     $              pln_frc(sfc_typ_idx,sgs_idx)* ! [frc]
     $              zpdvt(pln_typ_idx) ! [m]
            enddo               ! end loop over number of possible plant types
         endif                  ! endif normal land
      end do                    ! end loop over lon
      return
      end subroutine zpd_get                      ! end zpd_get()
      
      subroutine blm_ocn(
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     flg_ocn,             ! I [flg] Ocean flag
     $     hgt_mdp,             ! I [m] Midlayer height above surface
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_vpr,           ! I [kg kg-1] Specific humidity
     $     tpt_mdp,             ! I [K] Midlayer temperature
     $     tpt_ptn_mdp,         ! I [K] Midlayer local potential temperature
     $     tpt_sfc,             ! I [K] Surface temperature
     $     wnd_mrd_mdp,         ! I [m s-1] Surface layer meridional wind speed
     $     wnd_znl_mdp,         ! I [m s-1] Surface layer zonal wind speed
     $     mno_lng,             ! O [m] Monin-Obukhov length
     $     rgh_mmn,             ! O [m] Roughness length momentum
     $     wnd_frc)             ! O [m s-1] Surface friction velocit
c     Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
c     compute the boundary layer exchange properties over open ocean
c     blm_ocn() is called by blm_glb()
c     Routine uses specified surface temperature rather than solving energy balance equation for new Ts
c     Theory and algorithms: Large and Pond (1981,1982) CCM:dom/flxoce(), Bonan (1996) CCM:lsm/surtem()
c     Coding: Zender
      
c     Notes on algorithm:
c     In this routine surface temperature tpt_sfc = Sea surface temperature SST
c     Routine uses sfc suffix rather than SST simply for consistency with land and sea ice routines
c     Suffix mdp quantity evaluated at height hgt_mdp
c     Suffix sfc quantity evaluated at (sea) surface temperature
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstblm.h> /* Physical constants for boundary layer meteorology */
c     Parameters
      integer,parameter::itr_max=5 ! [nbr] Maximum number of iterations for tpt_gnd loop
      real(r8),parameter::dlt_nbr_ntr_10m_cff=0.0346 ! [frc] Coefficient for neutral 10 m Dalton number CCM:dom/flxoce() LaP82 p. 477
      real(r8),parameter::hgt_rfr_LaP81=10.0 ! [m] Reference height for turbulent flux parameterization of LaP81
      real(r8),parameter::ssh_H2O_sln_crc=0.98 ! [frc] Salinity correction to saturation specific humidity of H2O LaP81 p. 328 CCM:dom/flxoce()
      real(r8),parameter::stn_nbr_ntr_10m_cff_stb=0.0180 ! [frc] Coefficient for neutral, stable, 10 m Stanton number CCM:dom/flxoce() LaP82 p. 476
      real(r8),parameter::stn_nbr_ntr_10m_cff_uns=0.0327 ! [frc] Coefficient for neutral, unstable, 10 m Stanton number CCM:dom/flxoce() LaP82 p. 476
      real(r8),parameter::wnd_min_dps=1.0 ! [m s-1] Minimum windspeed used for mobilization
c     Input
      logical flg_ocn(plond)    ! I [flg] Ocean flag
      real(r8) dns_mdp(plond)       ! I [kg m-3] Midlayer density
      real(r8) hgt_mdp(plond)       ! I [m] Midlayer height above surface
      real(r8) prs_mdp(plond)       ! I [Pa] Pressure
      real(r8) q_H2O_vpr(plond)     ! I [kg kg-1] Specific humidity
      real(r8) tpt_mdp(plond)       ! I [K] Midlayer temperature
      real(r8) tpt_ptn_mdp(plond)   ! I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
      real(r8) tpt_sfc(plond)       ! I [K] Surface temperature
      real(r8) wnd_mrd_mdp(plond)   ! I [m s-1] Surface layer meridional wind speed
      real(r8) wnd_znl_mdp(plond)   ! I [m s-1] Surface layer zonal wind speed
c     Input/Output
c     Output
      real(r8) mno_lng(plond)       ! O [m] Monin-Obukhov length
      real(r8) rgh_mmn(plond)       ! O [m] Roughness length momentum
      real(r8) wnd_frc(plond)       ! O [m s-1] Surface friction velocity
c     Local
      integer itr_idx           ! [idx] Counting index
      integer lon_idx           ! [idx] Counting index for lon
      real(r8) hgt_rat_log          ! [frc] Log of ratio of local to reference height
      real(r8) ltn_heat_trn         ! [J kg-1] Latent heat of sublimation or evaporation
      real(r8) mno_stb_crc_heat     ! [frc] Monin-Obukhov stability correction heat
      real(r8) mno_stb_crc_mmn      ! [frc] Monin-Obukhov stability correction momentum
      real(r8) sml_fnc_mmn_uns_rcp  ! [frc] Reciprocal of similarity function for momentum, unstable atmosphere
      real(r8) mno_stb_prm(plond)   ! [frc] Monin-Obukhov stability parameter 
      real(r8) ntp_fct              ! [frc] Interpolation factor in reference height temperature calculation
      real(r8) q_H2O_vpr_dlt        ! [kg kg-1] Humidity change
      real(r8) spc_heat_mst_air     ! [J kg-1 K-1] Specific heat of moist air
      real(r8) ssh_H2O_sfc          ! [kg kg-1] Saturation specific humidity of H2O at surface
      real(r8) stb_val              ! [flg] 1.0 if stable, 0.0 if unstable
      real(r8) svp_H2O_sfc          ! [Pa] Saturation vapor pressure over planar condensed water at surface
      real(r8) tpt_bnd_cls          ! [C] Temperature bounded celsius
      real(r8) tpt_dlt              ! [K] Temperature change
      real(r8) tpt_ptn_vrt_mdp      ! [K] Midlayer virtual potential temperature
      real(r8) tpt_scl(plond)       ! [K] Temperature scale
      real(r8) vpr_scl(plond)       ! [kg kg-1] Moisture scale
      real(r8) wnd_frc_old          ! [m s-1] Surface friction velocity old
      real(r8) wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
      real(r8) wnd_mdp_bnd          ! [m s-1] Surface layer mean wind speed bounded
      real(r8) wnd_rfr_ntr(plond)   ! [m s-1] Neutral 10 m wind speed
      real(r8) xch_cff_heat_ntr_sqrt ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
      real(r8) xch_cff_heat_sqrt    ! [frc] Squareroot of mid-layer Stanton number for heat exchange
      real(r8) xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m drag coefficient
      real(r8) xch_cff_mmn_sqrt     ! [frc] Squareroot of mid-layer drag coefficient
      real(r8) xch_cff_vpr_ntr_sqrt ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
      real(r8) xch_cff_vpr_sqrt     ! [frc] Squareroot of mid-layer Dalton number for vapor exchange
      real(r8) xpn_heat_fct         ! [frc] Factor of heat exchange
      real(r8) xpn_mmn_fct          ! [frc] Factor of momentum exchange

#ifdef UNDEFINED
c     Variables used in diagnostics
      real(r8) flx_LW_up_sfc(plond) ! [W m-2] Longwave upwelling flux at surface
      real(r8) flx_ltn(plond)       ! [W m-2] Latent heat flux to atmosphere
      real(r8) flx_sns(plond)       ! [W m-2] Sensible heat flux to atmosphere
      real(r8) wnd_str(plond)       ! [kg m-1 s-2] Wind stress
      real(r8) wnd_str_mrd(plond)   ! [kg m-1 s-2] Meridional wind stress
      real(r8) wnd_str_znl(plond)   ! [kg m-1 s-2] Zonal wind stress
#endif /* not UNDEFINED */

cc     Next three variables needed only for 10 m wind speed and mid-layer drag coefficient
c      real(r8) hgt_zpd(plond)       ! [m] Zero plane displacement
c      real(r8) wnd_rfr(plond)       ! [m s-1] Wind speed at reference height LaP81 p. 327 (14)
c      real(r8) xch_cff_mmn(plond)   ! [frc] Drag coefficient at mid-layer
      
c     Inline statement functions for thermodynamics from Lowe and Ficke (1974) as reported in PrK78 p. 625. Range of validity is -50 C < T < 50 C.
      real(r8) a0,a1,a2,a3,a4,a5,a6 ! Coefficients for SVP over planar liquid water
      real(r8) svp_H2O_lqd_PrK78_fst_scl ! [Pa] Saturation vapor pressure over planar liquid water
      real(r8) tpt                  ! [K] Temperature 
      real(r8) tpt_cls              ! [C] Temperature in celsius
      real(r8) tpt_bnd_cls_get      ! [C] Temperature bounded celsius
      
c     Inline statement functions for boundary layer meteorology
      real(r8) mno_stb_crc_heat_uns_get ! [frc] Stability correction for heat, unstable case
      real(r8) mno_stb_crc_mmn_uns_get ! [frc] Stability correction for momentum, unstable case
      real(r8) wnd_10m_ntr          ! [m s-1] Wind speed at 10 m
      real(r8) xch_cff_mmn_ocn_ntr_get ! [frc] Neutral 10 m drag coefficient over ocean
      
      parameter(                ! Coefficients for SVP over planar liquid water
     $     a0=6.107799961,    a1=4.436518521e-01,
     $     a2=1.428945805e-02,a3=2.650648471e-04,
     $     a4=3.031240396e-06,a5=2.034080948e-08,
     $     a6=6.136820929e-11)
      tpt_bnd_cls_get(tpt)=min(fifty,max(-fifty,(tpt-tpt_frz_pnt))) ! [C] Temperature bounded celsius
c     Saturation vapor pressure over planar liquid water
      svp_H2O_lqd_PrK78_fst_scl(tpt_cls)= ! [Pa] 
     $     100.0*(a0+tpt_cls*(a1+tpt_cls*(a2+tpt_cls*(a3+tpt_cls*(a4+tpt_cls*(a5+tpt_cls*a6))))))
c     Stability correction for momentum, unstable case
      mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp)= ! [frc] Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325
     $     log((1.0+sml_fnc_mmn_uns_rcp*(2.0+sml_fnc_mmn_uns_rcp))*(1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/8.0)
     $     -2.0*atan(sml_fnc_mmn_uns_rcp)+1.571
c     Stability correction for heat, unstable case
      mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp)= ! [frc] Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325
     $     2.0*log((1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0)
c     Neutral 10 m drag coefficient over ocean
      xch_cff_mmn_ocn_ntr_get(wnd_10m_ntr)= ! [frc] LaP81 CCM:dom/flxoce(), NOS97 p. I2
     $     0.0027/wnd_10m_ntr+0.000142+0.0000764*wnd_10m_ntr
      
c     Main Code
      
c     Initialize variables which are independent of stability iteration
      do lon_idx=1,plon
         if (flg_ocn(lon_idx)) then
            
c     Midlayer wind speeds
            wnd_mdp(lon_idx)=   ! [m s-1] Surface layer mean wind speed
     $           sqrt(wnd_znl_mdp(lon_idx)*wnd_znl_mdp(lon_idx)+ 
     $           wnd_mrd_mdp(lon_idx)*wnd_mrd_mdp(lon_idx)) 
            wnd_mdp_bnd=max(wnd_mdp(lon_idx),wnd_min_dps) ! [m s-1] Surface layer mean wind speed bounded
            
            tpt_ptn_vrt_mdp=tpt_ptn_mdp(lon_idx)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(lon_idx)) ! [K] Midlayer Virtual potential temperature
            tpt_bnd_cls=tpt_bnd_cls_get(tpt_sfc(lon_idx)) ! [C]
c     Liquid sea water supercools befored freezing to sea ice (CCM:dom/parsst.h/tsice = -1.7999 C)
c     Assume gridpoint liquid sea water thus use heat of vaporization not sublimation
            svp_H2O_sfc=svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
            ltn_heat_trn=ltn_heat_vpr_H2O_std ! [J kg-1]
            ssh_H2O_sfc=eps_H2O*svp_H2O_sfc/(prs_mdp(lon_idx)-one_mns_eps_H2O*svp_H2O_sfc) ! [kg kg-1] Saturation specific humidity of H2O at surface
c     Correct freshwater saturation specific humidity for salinity effects
            ssh_H2O_sfc=ssh_H2O_sfc*ssh_H2O_sln_crc ! [kg kg-1] Saturation specific humidity of H2O at surface
            
            tpt_dlt=tpt_ptn_mdp(lon_idx)-tpt_sfc(lon_idx) ! [K] Temperature change
            q_H2O_vpr_dlt=q_H2O_vpr(lon_idx)-ssh_H2O_sfc ! [kg kg-1] Humidity change
            hgt_rat_log=log(hgt_mdp(lon_idx)/hgt_rfr_LaP81) ! [frc] Log of ratio of local to reference height
            spc_heat_mst_air=spc_heat_dry_air*(1.0+cp_vpr_rcp_cp_dry_m1*ssh_H2O_sfc) ! [J kg-1 K-1] Specific heat of moist air
            
c     Stable if tpt_ptn_mdp > tpt_sfc
            stb_val=half+sign(half,tpt_dlt) ! [flg] 1.0 if stable, 0.0 if unstable
c     Initial guess for roots of neutral exchange coefficients: z/L=0 and u10n=vmag
            xch_cff_mmn_ntr_sqrt=sqrt(xch_cff_mmn_ocn_ntr_get(wnd_mdp_bnd)) ! [frc] Squareroot of neutral 10 m drag coefficient
            xch_cff_heat_ntr_sqrt=(1.0-stb_val)*stn_nbr_ntr_10m_cff_uns+stb_val*stn_nbr_ntr_10m_cff_stb ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
            xch_cff_vpr_ntr_sqrt=dlt_nbr_ntr_10m_cff ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
            
      
c     Iteration loop
            do itr_idx=1,itr_max
               
c     Save old friction speed for convergence diagnostic
               if (itr_idx == 1) then
                  wnd_frc_old=0.0 ! [m s-1]
               else 
                  wnd_frc_old=wnd_frc(lon_idx) ! [m s-1]
               endif            ! endif
               
               if (itr_idx == 1) then
c     First iteration: Use estimated roots of neutral exchange coefficients
                  wnd_frc(lon_idx)=xch_cff_mmn_ntr_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
                  tpt_scl(lon_idx)=xch_cff_heat_ntr_sqrt*tpt_dlt ! [K] Temperature scale
                  vpr_scl(lon_idx)=xch_cff_vpr_ntr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
               else
c     Subsequently: Use current roots of (non-neutral) exchange coefficients at measurement height
                  wnd_frc(lon_idx)=xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
                  tpt_scl(lon_idx)=xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale
                  vpr_scl(lon_idx)=xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
               endif            ! endif first iteration
               
c     Compute stability parameter at midlayer and evaluate stability corrections  
c     Stable if tpt_ptn_mdp > tpt_sfc or if mno_stb_prm > 0
               mno_stb_prm(lon_idx)= ! [frc] Monin-Obukhov stability parameter 
     $              cst_von_krm*grv_sfc*hgt_mdp(lon_idx)
     $              *(tpt_scl(lon_idx)/tpt_ptn_vrt_mdp+vpr_scl(lon_idx)/(1.0/eps_H2O_rcp_m1+q_H2O_vpr(lon_idx)))
     $              /(wnd_frc(lon_idx)*wnd_frc(lon_idx))
c     Ensure |z/L| < 10.0 because similarity function asymptotes to (-z/L)^(-1/3) for z/L << -1, Ary88 p. 166
               mno_stb_prm(lon_idx)=sign(min(abs(mno_stb_prm(lon_idx)),ten),mno_stb_prm(lon_idx)) ! [frc] Monin-Obukhov stability parameter
               stb_val=half+sign(half,mno_stb_prm(lon_idx)) ! [flg] 1.0 if stable, 0.0 if unstable
               sml_fnc_mmn_uns_rcp=max(sqrt(abs(1.0-16.0*mno_stb_prm(lon_idx))),one) ! [frc] BKL97 p. F1, LaP81 p. 325
               sml_fnc_mmn_uns_rcp=sqrt(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
               mno_stb_crc_mmn=-5.0*mno_stb_prm(lon_idx)*stb_val+(1.0-stb_val)*mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
               mno_stb_crc_heat=-5.0*mno_stb_prm(lon_idx)*stb_val+(1.0-stb_val)*mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
               
c     Shift old neutral 10 m exchange coefficient to measurement height and stability 
               xch_cff_mmn_sqrt=xch_cff_mmn_ntr_sqrt/(1.0+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)) ! [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
               
c     Define neutral 10 m wind speed 
               wnd_rfr_ntr(lon_idx)=wnd_mdp_bnd*xch_cff_mmn_sqrt/xch_cff_mmn_ntr_sqrt ! [m s-1] Neutral 10 m wind speed
               
c     Update neutral 10 m exchange coefficients
               xch_cff_mmn_ntr_sqrt=sqrt(xch_cff_mmn_ocn_ntr_get(wnd_rfr_ntr(lon_idx))) ! [frc] Squareroot of neutral 10 m drag coefficient
               xch_cff_vpr_ntr_sqrt=dlt_nbr_ntr_10m_cff ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
               xch_cff_heat_ntr_sqrt=(1.0-stb_val)*stn_nbr_ntr_10m_cff_uns+stb_val*stn_nbr_ntr_10m_cff_stb ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
               
c     Shift old neutral 10 m exchange coefficients to measurement height and stability
               xch_cff_mmn_sqrt=xch_cff_mmn_ntr_sqrt/(1.0+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)) ! [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
               xch_cff_heat_sqrt=xch_cff_heat_ntr_sqrt/(1.0+xch_cff_heat_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Stanton number for heat exchange LaP82 p. 466 (10)
               xch_cff_vpr_sqrt=xch_cff_vpr_ntr_sqrt/(1.0+xch_cff_vpr_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Dalton number for vapor exchange LaP82 p. 466 (10)
               
#ifdef UNDEFINED
               if (.true.) then
                  mno_lng(lon_idx)=hgt_mdp(lon_idx)/mno_stb_prm(lon_idx) ! [m] Monin-Obukhov length
                  wnd_str(lon_idx)=dns_mdp(lon_idx)*wnd_frc(lon_idx)*wnd_frc(lon_idx) ! [kg m-1 s-2] Wind stress
                  flx_sns(lon_idx)=-spc_heat_mst_air*wnd_str(lon_idx)*tpt_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Sensible heat flux to atmosphere
                  flx_ltn(lon_idx)=-ltn_heat_trn*wnd_str(lon_idx)*vpr_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Latent heat flux to atmosphere
                  flx_LW_up_sfc(lon_idx)=cst_Stefan_Boltzmann*(tpt_sfc(lon_idx)**4.0) ! [W m-2] Longwave upwelling flux to atmosphere
                  write (6,'(i3,1x,f9.3,1x,f8.3,1x,f7.4,1x,f7.3,1x,f7.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3)')
     $                 itr_idx,mno_lng(lon_idx),mno_stb_prm(lon_idx),
     $                 wnd_frc(lon_idx),wnd_rfr_ntr(lon_idx),
     $                 1000.0*xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt,
     $                 1000.0*xch_cff_mmn_sqrt*xch_cff_mmn_sqrt,
     $                 flx_sns(lon_idx),flx_ltn(lon_idx),flx_LW_up_sfc(lon_idx)
               endif            ! endif dbg
#endif /* not UNDEFINED */
               
            enddo               ! end loop over itr
            

c     Update wnd_frc, tpt_scl, and vpr_scl using updated exchange coefficients at measurement height
            wnd_frc(lon_idx)=xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity 
            tpt_scl(lon_idx)=xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale 
            vpr_scl(lon_idx)=xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale 
            
c     Compute surface stress components
c            wnd_str(lon_idx)=dns_mdp(lon_idx)*wnd_frc(lon_idx)*wnd_frc(lon_idx) ! [kg m-1 s-2] Wind stress
c            wnd_str_znl(lon_idx)=-wnd_str(lon_idx)*wnd_znl_mdp(lon_idx)/wnd_mdp_bnd ! [kg m-1 s-2] Zonal wind stress
c            wnd_str_mrd(lon_idx)=-wnd_str(lon_idx)*wnd_mrd_mdp(lon_idx)/wnd_mdp_bnd ! [kg m-1 s-2] Meridional wind stress
            
c     Compute heat flux components at current surface temperature
c     Define positive latent and sensible heat as upwards into atmosphere 
c            flx_sns(lon_idx)=-spc_heat_mst_air*wnd_str(lon_idx)*tpt_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Sensible heat flux to atmosphere
c            flx_ltn(lon_idx)=-ltn_heat_trn*wnd_str(lon_idx)*vpr_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Latent heat flux to atmosphere
c            flx_LW_up_sfc(lon_idx)=cst_Stefan_Boltzmann*(tpt_sfc(lon_idx)**4.0) ! [W m-2] Longwave upwelling flux to atmosphere
            
c     Additional diagnostics
            mno_lng(lon_idx)=hgt_mdp(lon_idx)/mno_stb_prm(lon_idx) ! [m] Monin-Obukhov length
            rgh_mmn(lon_idx)=hgt_rfr_LaP81*exp(-cst_von_krm/xch_cff_mmn_ntr_sqrt) ! [m] BKL97 p. F-4, LaP81 p. 327 (14) 

c           xch_cff_mmn(lon_idx)=xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt ! [frc] Drag coefficient at mid-layer
            
         endif                  ! endif flg_ocn
      enddo                     ! end loop over lon
      
cc     Initialize hgt_zpd to 0.0 m for ocean
c      hgt_zpd=0.0               ! [m]
cc     Interpolate midlayer wind speed to 10 m
c      call wnd_rfr_get(
c     $     flg_ocn,             ! I [flg] Ocean flag
c     $     hgt_mdp,             ! I [m] Midpoint height above surface
c     $     hgt_rfr_LaP81,       ! I [m] Reference height for turbulent flux parameterization of LaP81
c     $     hgt_zpd,             ! I [m] Zero plane displacement
c     $     mno_lng,             ! I [m] Monin-Obukhov length
c     $     wnd_frc,             ! I [m s-1] Surface friction velocity
c     $     wnd_mdp,             ! I [m s-1] Surface layer mean wind speed
c     $     wnd_min_dps,         ! I [m s-1] Minimum windspeed
c     $     wnd_rfr)             ! O [m s-1] Wind speed at reference height

      return
      end subroutine blm_ocn                      ! end blm_ocn()
      
      subroutine blm_ice(
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     flg_ice,             ! I [flg] Sea ice flag
     $     hgt_mdp,             ! I [m] Midlayer height above surface
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_vpr,           ! I [kg kg-1] Specific humidity
     $     rgh_mmn,             ! I [m] Roughness length momentum
     $     tpt_mdp,             ! I [K] Midlayer temperature
     $     tpt_ptn_mdp,         ! I [K] Midlayer local potential temperature
     $     tpt_sfc,             ! I [K] Surface temperature
     $     wnd_mrd_mdp,         ! I [m s-1] Surface layer meridional wind speed
     $     wnd_znl_mdp,         ! I [m s-1] Surface layer zonal wind speed
     $     mno_lng,             ! O [m] Monin-Obukhov length
     $     wnd_frc)             ! O [m s-1] Surface friction velocit
c     Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
c     compute the boundary layer exchange properties over sea ice
c     blm_ice() is called by blm_glb()
c     Routine uses specified surface temperature rather than solving energy balance equation for new Ts
c     Theory and algorithms: Large and Pond (1981,1982) CCM:dom/flxoce(), Bonan (1996) CCM:lsm/surtem()
c     Coding: Zender
      
c     Notes on algorithm:
c     Suffix mdp quantity evaluated at height hgt_mdp
c     Suffix sfc quantity evaluated at (sea) surface temperature
c     The roughness length of sea ice depends on many factors
c     A globally uniform roughness length for sea ice is foolish but necessary
c     CCM1/2/3 used z0m=0.04 m (CCM:dom/parpbl), based on a glacial ice value
c     LSM also adopted this value for land ice
c     This value is more appropriate for ridged, multi-year ice
c     The NCAR Oceanography section used z0m=0.05 m (BKL97 p. F-3)
c     These values lead to large drag coefficients and excessive ice extent off of Antarctica
c     CSM later adopted z0m=0.0005 m, appropriate for new-formed, seasonal ice
c     This value greatly improved sea-ice dynamics in CSM, so we adopt it
c     This blm_ice() routine reads roughness length as an input parameter
c     for algorithmic reasons, and checks that its value agrees with CCM sea ice
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstblm.h> /* Physical constants for boundary layer meteorology */
c     Parameters
      integer,parameter::itr_max=5 ! [nbr] Maximum number of iterations for tpt_gnd loop
      real(r8),parameter::hgt_rfr_LaP81=10.0 ! [m] Reference height for turbulent flux parameterization of LaP81
      real(r8),parameter::rgh_mmn_ice_ocn=0.0005 ! [m] Roughness length over sea ice BKL97 p. F-3
      real(r8),parameter::wnd_min_dps=1.0 ! [m s-1] Minimum windspeed used for mobilization
c     Input
      logical flg_ice(plond)    ! I [flg] Sea ice flag
      real(r8) dns_mdp(plond)       ! I [kg m-3] Midlayer density
      real(r8) hgt_mdp(plond)       ! I [m] Midlayer height above surface
      real(r8) prs_mdp(plond)       ! I [Pa] Pressure
      real(r8) q_H2O_vpr(plond)     ! I [kg kg-1] Specific humidity
      real(r8) tpt_mdp(plond)       ! I [K] Midlayer temperature
      real(r8) tpt_ptn_mdp(plond)   ! I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
      real(r8) tpt_sfc(plond)       ! I [K] Surface temperature
      real(r8) wnd_mrd_mdp(plond)   ! I [m s-1] Surface layer meridional wind speed
      real(r8) wnd_znl_mdp(plond)   ! I [m s-1] Surface layer zonal wind speed
c     Input/Output
c     Output
      real(r8) mno_lng(plond)       ! O [m] Monin-Obukhov length
      real(r8) rgh_mmn(plond)       ! O [m] Roughness length momentum
      real(r8) wnd_frc(plond)       ! O [m s-1] Surface friction velocity
c     Local
      integer itr_idx           ! [idx] Counting index
      integer lon_idx           ! [idx] Counting index for lon
      real(r8) hgt_rat_log          ! [frc] Log of ratio of local to reference height
      real(r8) ltn_heat_trn         ! [J kg-1] Latent heat of sublimation or evaporation
      real(r8) mno_stb_crc_heat     ! [frc] Monin-Obukhov stability correction heat
      real(r8) mno_stb_crc_mmn      ! [frc] Monin-Obukhov stability correction momentum
      real(r8) sml_fnc_mmn_uns_rcp  ! [frc] Reciprocal of similarity function for momentum, unstable atmosphere
      real(r8) mno_stb_prm(plond)   ! [frc] Monin-Obukhov stability parameter 
      real(r8) ntp_fct              ! [frc] Interpolation factor in reference height temperature calculation
      real(r8) q_H2O_vpr_dlt        ! [kg kg-1] Humidity change
      real(r8) spc_heat_mst_air     ! [J kg-1 K-1] Specific heat of moist air
      real(r8) ssh_H2O_sfc          ! [kg kg-1] Saturation specific humidity of H2O at surface
      real(r8) stb_val              ! [flg] 1.0 if stable, 0.0 if unstable
      real(r8) svp_H2O_sfc          ! [Pa] Saturation vapor pressure over planar condensed water at surface
      real(r8) tpt_bnd_cls          ! [C] Temperature bounded celsius
      real(r8) tpt_dlt              ! [K] Temperature change
      real(r8) tpt_ptn_vrt_mdp      ! [K] Midlayer virtual potential temperature
      real(r8) tpt_scl(plond)       ! [K] Temperature scale
      real(r8) vpr_scl(plond)       ! [kg kg-1] Moisture scale
      real(r8) wnd_frc_old          ! [m s-1] Surface friction velocity old
      real(r8) wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
      real(r8) wnd_mdp_bnd          ! [m s-1] Surface layer mean wind speed bounded
      real(r8) wnd_rfr_ntr(plond)   ! [m s-1] Neutral 10 m wind speed
      real(r8) xch_cff_heat_ntr_sqrt ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
      real(r8) xch_cff_heat_sqrt    ! [frc] Squareroot of mid-layer Stanton number for heat exchange
      real(r8) xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m drag coefficient
      real(r8) xch_cff_mmn_sqrt     ! [frc] Squareroot of mid-layer drag coefficient
      real(r8) xch_cff_vpr_ntr_sqrt ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
      real(r8) xch_cff_vpr_sqrt     ! [frc] Squareroot of mid-layer Dalton number for vapor exchange
      real(r8) xpn_heat_fct         ! [frc] Factor of heat exchange
      real(r8) xpn_mmn_fct          ! [frc] Factor of momentum exchange

#ifdef UNDEFINED
c     Variables used in diagnostics
      real(r8) flx_LW_up_sfc(plond) ! [W m-2] Longwave upwelling flux at surface
      real(r8) flx_ltn(plond)       ! [W m-2] Latent heat flux to atmosphere
      real(r8) flx_sns(plond)       ! [W m-2] Sensible heat flux to atmosphere
      real(r8) wnd_str(plond)       ! [kg m-1 s-2] Wind stress
      real(r8) wnd_str_mrd(plond)   ! [kg m-1 s-2] Meridional wind stress
      real(r8) wnd_str_znl(plond)   ! [kg m-1 s-2] Zonal wind stress
#endif /* not UNDEFINED */

c     Inline statement functions for thermodynamics from Lowe and Ficke (1974) as reported in PrK78 p. 625. Range of validity is -50 C < T < 0 C.
      real(r8) b0,b1,b2,b3,b4,b5,b6 ! Coefficients for SVP over planar ice water
      real(r8) svp_H2O_ice_PrK78_fst_scl ! [Pa] Saturation vapor pressure over planar ice water
      real(r8) tpt                  ! [K] Temperature 
      real(r8) tpt_cls              ! [C] Temperature in celsius
      real(r8) tpt_bnd_cls_get      ! [C] Temperature bounded celsius
      
c     Inline statement functions for boundary layer meteorology
      real(r8) mno_stb_crc_heat_uns_get ! [frc] Stability correction for heat, unstable case
      real(r8) mno_stb_crc_mmn_uns_get ! [frc] Stability correction for momentum, unstable case
      parameter(                ! Coefficients for SVP over planar ice water
     $     b0=6.109177956,    b1=5.034698970e-01,
     $     b2=1.886013408e-02,b3=4.176223716e-04,
     $     b4=5.824720280e-06,b5=4.838803174e-08,
     $     b6=1.838826904e-10)
      tpt_bnd_cls_get(tpt)=min(fifty,max(-fifty,(tpt-tpt_frz_pnt))) ! [C] Temperature bounded celsius
c     Saturation vapor pressure over planar ice water
      svp_H2O_ice_PrK78_fst_scl(tpt_cls)= ! [Pa] 
     $     100.0*(b0+tpt_cls*(b1+tpt_cls*(b2+tpt_cls*(b3+tpt_cls*(b4+tpt_cls*(b5+tpt_cls*b6))))))
c     Stability correction for momentum, unstable case
      mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp)= ! [frc] Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325
     $     log((1.0+sml_fnc_mmn_uns_rcp*(2.0+sml_fnc_mmn_uns_rcp))*(1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/8.0)
     $     -2.0*atan(sml_fnc_mmn_uns_rcp)+1.571
c     Stability correction for heat, unstable case
      mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp)= ! [frc] Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325
     $     2.0*log((1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0)
      
c     Main Code
      
c     Initialize variables which are independent of stability iteration
      do lon_idx=1,plon
         if (flg_ice(lon_idx)) then
            
c     Midlayer wind speeds
            wnd_mdp(lon_idx)=   ! [m s-1] Surface layer mean wind speed
     $           sqrt(wnd_znl_mdp(lon_idx)*wnd_znl_mdp(lon_idx)+ 
     $           wnd_mrd_mdp(lon_idx)*wnd_mrd_mdp(lon_idx)) 
            wnd_mdp_bnd=max(wnd_mdp(lon_idx),wnd_min_dps) ! [m s-1] Surface layer mean wind speed bounded
            
            tpt_ptn_vrt_mdp=tpt_ptn_mdp(lon_idx)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(lon_idx)) ! [K] Midlayer Virtual potential temperature
            tpt_bnd_cls=tpt_bnd_cls_get(tpt_sfc(lon_idx)) ! [C]
c     Assume gridpoint is sea ice thus use heat of sublimation not vaporization (contrary to CSM)
c     NB: Prognostic models must sometimes use heat of vaporization everywhere to conserve energy
            svp_H2O_sfc=svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
            ltn_heat_trn=ltn_heat_sbl_H2O_std ! [J kg-1]
            ssh_H2O_sfc=eps_H2O*svp_H2O_sfc/(prs_mdp(lon_idx)-one_mns_eps_H2O*svp_H2O_sfc) ! [kg kg-1] Saturation specific humidity of H2O at surface
c     NB: No salinity correction for sea ice (contrary to CSM)
            
            tpt_dlt=tpt_ptn_mdp(lon_idx)-tpt_sfc(lon_idx) ! [K] Temperature change
            q_H2O_vpr_dlt=q_H2O_vpr(lon_idx)-ssh_H2O_sfc ! [kg kg-1] Humidity change
            hgt_rat_log=log(hgt_mdp(lon_idx)/hgt_rfr_LaP81) ! [frc] Log of ratio of local to reference height
            spc_heat_mst_air=spc_heat_dry_air*(1.0+cp_vpr_rcp_cp_dry_m1*ssh_H2O_sfc) ! [J kg-1 K-1] Specific heat of moist air
            
c     Iteration loop
            do itr_idx=1,itr_max
               
c     Save old friction speed for convergence diagnostic
               if (itr_idx == 1) then
                  wnd_frc_old=0.0 ! [m s-1]
               else 
                  wnd_frc_old=wnd_frc(lon_idx) ! [m s-1]
               endif            ! endif
               
c     Roots of neutral, 10 m exchange coefficients
               xch_cff_mmn_ntr_sqrt=cst_von_krm/log(hgt_rfr_LaP81/rgh_mmn(lon_idx)) ! [frc] Squareroot of neutral 10 m drag coefficient
               xch_cff_heat_ntr_sqrt=xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
               xch_cff_vpr_ntr_sqrt=xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange

               if (itr_idx == 1) then
c     First iteration: Use estimated roots of neutral exchange coefficients
                  wnd_frc(lon_idx)=xch_cff_mmn_ntr_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
                  tpt_scl(lon_idx)=xch_cff_heat_ntr_sqrt*tpt_dlt ! [K] Temperature scale
                  vpr_scl(lon_idx)=xch_cff_vpr_ntr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
               else
c     Subsequently: Use current roots of (non-neutral) exchange coefficients at measurement height
                  wnd_frc(lon_idx)=xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
                  tpt_scl(lon_idx)=xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale
                  vpr_scl(lon_idx)=xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
               endif            ! endif first iteration
               
c     Compute stability parameter at midlayer and evaluate stability corrections  
c     Stable if tpt_ptn_mdp > tpt_sfc or if mno_stb_prm > 0
               mno_stb_prm(lon_idx)= ! [frc] Monin-Obukhov stability parameter 
     $              cst_von_krm*grv_sfc*hgt_mdp(lon_idx)
     $              *(tpt_scl(lon_idx)/tpt_ptn_vrt_mdp+vpr_scl(lon_idx)/(1.0/eps_H2O_rcp_m1+q_H2O_vpr(lon_idx)))
     $              /(wnd_frc(lon_idx)*wnd_frc(lon_idx))
c     Ensure |z/L| < 10.0 because similarity function asymptotes to (-z/L)^(-1/3) for z/L << -1, Ary88 p. 166
               mno_stb_prm(lon_idx)=sign(min(abs(mno_stb_prm(lon_idx)),ten),mno_stb_prm(lon_idx)) ! [frc] Monin-Obukhov stability parameter
               stb_val=half+sign(half,mno_stb_prm(lon_idx)) ! [flg] 1.0 if stable, 0.0 if unstable
               sml_fnc_mmn_uns_rcp=max(sqrt(abs(1.0-16.0*mno_stb_prm(lon_idx))),one) ! [frc] BKL97 p. F1, LaP81 p. 325
               sml_fnc_mmn_uns_rcp=sqrt(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
               mno_stb_crc_mmn=-5.0*mno_stb_prm(lon_idx)*stb_val+(1.0-stb_val)*mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
               mno_stb_crc_heat=-5.0*mno_stb_prm(lon_idx)*stb_val+(1.0-stb_val)*mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
               
c     Shift old neutral 10 m exchange coefficients to measurement height and stability
               xch_cff_mmn_sqrt=xch_cff_mmn_ntr_sqrt/(1.0+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)) ! [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
               xch_cff_heat_sqrt=xch_cff_heat_ntr_sqrt/(1.0+xch_cff_heat_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Stanton number for heat exchange LaP82 p. 466 (10)
               xch_cff_vpr_sqrt=xch_cff_vpr_ntr_sqrt/(1.0+xch_cff_vpr_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Dalton number for vapor exchange LaP82 p. 466 (10)
               
#ifdef UNDEFINED
               if (.true.) then
                  mno_lng(lon_idx)=hgt_mdp(lon_idx)/mno_stb_prm(lon_idx) ! [m] Monin-Obukhov length
                  wnd_str(lon_idx)=dns_mdp(lon_idx)*wnd_frc(lon_idx)*wnd_frc(lon_idx) ! [kg m-1 s-2] Wind stress
                  flx_sns(lon_idx)=-spc_heat_mst_air*wnd_str(lon_idx)*tpt_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Sensible heat flux to atmosphere
                  flx_ltn(lon_idx)=-ltn_heat_trn*wnd_str(lon_idx)*vpr_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Latent heat flux to atmosphere
                  flx_LW_up_sfc(lon_idx)=cst_Stefan_Boltzmann*(tpt_sfc(lon_idx)**4.0) ! [W m-2] Longwave upwelling flux to atmosphere
                  write (6,'(i3,1x,f9.3,1x,f8.3,1x,f7.4,1x,f7.3,1x,f7.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3)')
     $                 itr_idx,mno_lng(lon_idx),mno_stb_prm(lon_idx),
     $                 wnd_frc(lon_idx),wnd_rfr_ntr(lon_idx),
     $                 1000.0*xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt,
     $                 1000.0*xch_cff_mmn_sqrt*xch_cff_mmn_sqrt,
     $                 flx_sns(lon_idx),flx_ltn(lon_idx),flx_LW_up_sfc(lon_idx)
               endif            ! endif dbg
#endif /* not UNDEFINED */
               
            enddo               ! end loop over itr
            

c     Update wnd_frc, tpt_scl, and vpr_scl using updated exchange coefficients at measurement height
            wnd_frc(lon_idx)=xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity 
            tpt_scl(lon_idx)=xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale 
            vpr_scl(lon_idx)=xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale 
            
c     Compute surface stress components
c            wnd_str(lon_idx)=dns_mdp(lon_idx)*wnd_frc(lon_idx)*wnd_frc(lon_idx) ! [kg m-1 s-2] Wind stress
c            wnd_str_znl(lon_idx)=-wnd_str(lon_idx)*wnd_znl_mdp(lon_idx)/wnd_mdp_bnd ! [kg m-1 s-2] Zonal wind stress
c            wnd_str_mrd(lon_idx)=-wnd_str(lon_idx)*wnd_mrd_mdp(lon_idx)/wnd_mdp_bnd ! [kg m-1 s-2] Meridional wind stress
            
c     Compute heat flux components at current surface temperature
c     Define positive latent and sensible heat as upwards into atmosphere 
c            flx_sns(lon_idx)=-spc_heat_mst_air*wnd_str(lon_idx)*tpt_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Sensible heat flux to atmosphere
c            flx_ltn(lon_idx)=-ltn_heat_trn*wnd_str(lon_idx)*vpr_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Latent heat flux to atmosphere
c            flx_LW_up_sfc(lon_idx)=cst_Stefan_Boltzmann*(tpt_sfc(lon_idx)**4.0) ! [W m-2] Longwave upwelling flux to atmosphere
            
c     Additional diagnostics
            mno_lng(lon_idx)=hgt_mdp(lon_idx)/mno_stb_prm(lon_idx) ! [m] Monin-Obukhov length
         endif                  ! endif flg_ice
      enddo                     ! end loop over lon
      
      return
      end subroutine blm_ice                       ! end blm_ice()

      subroutine blm_lnd(
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     flg_lnd,             ! I [flg] Land flag
     $     hgt_mdp,             ! I [m] Midlayer height above surface
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_vpr,           ! I [kg kg-1] Specific humidity
     $     rgh_mmn,             ! I [m] Roughness length momentum
     $     tpt_mdp,             ! I [K] Midlayer temperature
     $     tpt_ptn_mdp,         ! I [K] Midlayer local potential temperature
     $     tpt_sfc,             ! I [K] Surface temperature
     $     wnd_mrd_mdp,         ! I [m s-1] Surface layer meridional wind speed
     $     wnd_znl_mdp,         ! I [m s-1] Surface layer zonal wind speed
     $     mno_lng,             ! O [m] Monin-Obukhov length
     $     wnd_frc)             ! O [m s-1] Surface friction velocity

c     Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
c     compute the boundary layer exchange properties over land
c     blm_lnd() is called by blm_glb()
c     Routine uses specified surface temperature rather than solving energy balance equation for new Ts
c     Theory and algorithms: Large and Pond (1981,1982) CCM:dom/flxoce(), Bonan (1996) CCM:lsm/surtem()
c     Coding: Zender
      
c     Notes on algorithm:
c     Suffix mdp quantity evaluated at height hgt_mdp
c     Suffix sfc quantity evaluated at (sea) surface temperature
c     Currently, routine is virtually identical to blm_ice()
c     fxm: Using blm_ice()/blm_ocn()-type physics is unjustifiable over
c     heterogeneous land surface types because of varying transpiration
c     resistances over heterogeneous land surface types
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstblm.h> /* Physical constants for boundary layer meteorology */
c     Parameters
      integer,parameter::itr_max=5      ! [nbr] Maximum number of iterations for tpt_gnd loop
      real(r8),parameter::hgt_rfr_LaP81=10.0 ! [m] Reference height for turbulent flux parameterization of LaP81
      real(r8),parameter::wnd_min_dps=1.0 ! [m s-1] Minimum windspeed used for mobilization
c     Input
      logical flg_lnd(plond)    ! I [flg] Land flag
      real(r8) dns_mdp(plond)       ! I [kg m-3] Midlayer density
      real(r8) hgt_mdp(plond)       ! I [m] Midlayer height above surface
      real(r8) prs_mdp(plond)       ! I [Pa] Pressure
      real(r8) q_H2O_vpr(plond)     ! I [kg kg-1] Specific humidity
      real(r8) tpt_mdp(plond)       ! I [K] Midlayer temperature
      real(r8) tpt_ptn_mdp(plond)   ! I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
      real(r8) tpt_sfc(plond)       ! I [K] Surface temperature
      real(r8) wnd_mrd_mdp(plond)   ! I [m s-1] Surface layer meridional wind speed
      real(r8) wnd_znl_mdp(plond)   ! I [m s-1] Surface layer zonal wind speed
c     Input/Output
c     Output
      real(r8) mno_lng(plond)       ! O [m] Monin-Obukhov length
      real(r8) rgh_mmn(plond)       ! O [m] Roughness length momentum
      real(r8) wnd_frc(plond)       ! O [m s-1] Surface friction velocity
c     Local
      integer itr_idx           ! [idx] Counting index
      integer lon_idx           ! [idx] Counting index for lon

      real(r8) hgt_rat_log          ! [frc] Log of ratio of local to reference height
      real(r8) ltn_heat_trn         ! [J kg-1] Latent heat of sublimation or evaporation
      real(r8) mno_stb_crc_heat     ! [frc] Monin-Obukhov stability correction heat
      real(r8) mno_stb_crc_mmn      ! [frc] Monin-Obukhov stability correction momentum
      real(r8) sml_fnc_mmn_uns_rcp  ! [frc] Reciprocal of similarity function for momentum, unstable atmosphere
      real(r8) mno_stb_prm(plond)   ! [frc] Monin-Obukhov stability parameter 
      real(r8) ntp_fct              ! [frc] Interpolation factor in reference height temperature calculation
      real(r8) q_H2O_vpr_dlt        ! [kg kg-1] Humidity change
      real(r8) spc_heat_mst_air     ! [J kg-1 K-1] Specific heat of moist air
      real(r8) ssh_H2O_sfc          ! [kg kg-1] Saturation specific humidity of H2O at surface
      real(r8) stb_val              ! [flg] 1.0 if stable, 0.0 if unstable
      real(r8) svp_H2O_sfc          ! [Pa] Saturation vapor pressure over planar condensed water at surface
      real(r8) tpt_bnd_cls          ! [C] Temperature bounded celsius
      real(r8) tpt_dlt              ! [K] Temperature change
      real(r8) tpt_ptn_vrt_mdp      ! [K] Midlayer virtual potential temperature
      real(r8) tpt_scl(plond)       ! [K] Temperature scale
      real(r8) vpr_scl(plond)       ! [kg kg-1] Moisture scale
      real(r8) wnd_frc_old          ! [m s-1] Surface friction velocity old
      real(r8) wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
      real(r8) wnd_mdp_bnd          ! [m s-1] Surface layer mean wind speed bounded
      real(r8) wnd_rfr_ntr(plond)   ! [m s-1] Neutral 10 m wind speed
      real(r8) xch_cff_heat_ntr_sqrt ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
      real(r8) xch_cff_heat_sqrt    ! [frc] Squareroot of mid-layer Stanton number for heat exchange
      real(r8) xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m drag coefficient
      real(r8) xch_cff_mmn_sqrt     ! [frc] Squareroot of mid-layer drag coefficient
      real(r8) xch_cff_vpr_ntr_sqrt ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
      real(r8) xch_cff_vpr_sqrt     ! [frc] Squareroot of mid-layer Dalton number for vapor exchange
      real(r8) xpn_heat_fct         ! [frc] Factor of heat exchange
      real(r8) xpn_mmn_fct          ! [frc] Factor of momentum exchange


c- this is where to put #define UNDEFINED
#ifdef UNDEFINED
c     Variables used in diagnostics
      real(r8) flx_LW_up_sfc(plond) ! [W m-2] Longwave upwelling flux at surface
      real(r8) flx_ltn(plond)       ! [W m-2] Latent heat flux to atmosphere
      real(r8) flx_sns(plond)       ! [W m-2] Sensible heat flux to atmosphere
      real(r8) wnd_str(plond)       ! [kg m-1 s-2] Wind stress
      real(r8) wnd_str_mrd(plond)   ! [kg m-1 s-2] Meridional wind stress
      real(r8) wnd_str_znl(plond)   ! [kg m-1 s-2] Zonal wind stress
#endif /* not UNDEFINED */

c     Inline statement functions for thermodynamics from Lowe and Ficke (1974) as reported in PrK78 p. 625. Range of validity is -50 C < T < 50 C.
      real(r8) a0,a1,a2,a3,a4,a5,a6 ! Coefficients for SVP over planar liquid water
      real(r8) b0,b1,b2,b3,b4,b5,b6 ! Coefficients for SVP over planar ice water
      real(r8) svp_H2O_ice_PrK78_fst_scl ! [Pa] Saturation vapor pressure over planar ice water
      real(r8) svp_H2O_lqd_PrK78_fst_scl ! [Pa] Saturation vapor pressure over planar liquid water
      real(r8) tpt                  ! [K] Temperature 
      real(r8) tpt_cls              ! [C] Temperature in celsius
      real(r8) tpt_bnd_cls_get      ! [C] Temperature bounded celsius

c     Inline statement functions for boundary layer meteorology
      real(r8) mno_stb_crc_heat_uns_get ! [frc] Stability correction for heat, unstable case
      real(r8) mno_stb_crc_mmn_uns_get ! [frc] Stability correction for momentum, unstable case
      parameter(                ! Coefficients for SVP over planar liquid water
     $     a0=6.107799961,    a1=4.436518521e-01,
     $     a2=1.428945805e-02,a3=2.650648471e-04,
     $     a4=3.031240396e-06,a5=2.034080948e-08,
     $     a6=6.136820929e-11)
      parameter(                ! Coefficients for SVP over planar ice water
     $     b0=6.109177956,    b1=5.034698970e-01,
     $     b2=1.886013408e-02,b3=4.176223716e-04,
     $     b4=5.824720280e-06,b5=4.838803174e-08,
     $     b6=1.838826904e-10)
      tpt_bnd_cls_get(tpt)=min(fifty,max(-fifty,(tpt-tpt_frz_pnt))) ! [C] Temperature bounded celsius
c     Saturation vapor pressure over planar liquid water
      svp_H2O_lqd_PrK78_fst_scl(tpt_cls)= ! [Pa] 
     $     100.0*(a0+tpt_cls*(a1+tpt_cls*(a2+tpt_cls*(a3+tpt_cls*(a4+tpt_cls*(a5+tpt_cls*a6))))))
c     Saturation vapor pressure over planar ice water
      svp_H2O_ice_PrK78_fst_scl(tpt_cls)= ! [Pa] 
     $     100.0*(b0+tpt_cls*(b1+tpt_cls*(b2+tpt_cls*(b3+tpt_cls*(b4+tpt_cls*(b5+tpt_cls*b6))))))
c     Stability correction for momentum, unstable case
      mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp)= ! [frc] Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325
     $     log((1.0+sml_fnc_mmn_uns_rcp*(2.0+sml_fnc_mmn_uns_rcp))*(1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/8.0)
     $     -2.0*atan(sml_fnc_mmn_uns_rcp)+1.571
c     Stability correction for heat, unstable case
      mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp)= ! [frc] Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325
     $     2.0*log((1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0)
      
c     Main Code
      

c     Initialize variables which are independent of stability iteration


      do lon_idx=1,plon
         if (flg_lnd(lon_idx)) then
            
c     Midlayer wind speeds
            wnd_mdp(lon_idx)=   ! [m s-1] Surface layer mean wind speed
     $           sqrt(wnd_znl_mdp(lon_idx)*wnd_znl_mdp(lon_idx)+ 
     $           wnd_mrd_mdp(lon_idx)*wnd_mrd_mdp(lon_idx)) 
            wnd_mdp_bnd=max(wnd_mdp(lon_idx),wnd_min_dps) ! [m s-1] Surface layer mean wind speed bounded
            
            tpt_ptn_vrt_mdp=tpt_ptn_mdp(lon_idx)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(lon_idx)) ! [K] Midlayer Virtual potential temperature
c     Saturation vapor pressure of water at ground temperature
c     NB: Prognostic models must sometimes use heat of vaporization everywhere to conserve energy
            tpt_bnd_cls=tpt_bnd_cls_get(tpt_sfc(lon_idx)) ! [C]
            if (tpt_bnd_cls > 0.0) then 
               svp_H2O_sfc=svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
               ltn_heat_trn=ltn_heat_vpr_H2O_std ! [J kg-1]
            else
               svp_H2O_sfc=svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
               ltn_heat_trn=ltn_heat_sbl_H2O_std ! [J kg-1]
            endif               ! endif frozen
            ssh_H2O_sfc=eps_H2O*svp_H2O_sfc/(prs_mdp(lon_idx)-one_mns_eps_H2O*svp_H2O_sfc) ! [kg kg-1] Saturation specific humidity of H2O at surface
c     NB: No salinity correction for land
            
          tpt_dlt=tpt_ptn_mdp(lon_idx)-tpt_sfc(lon_idx) ! [K] Temperature change
            q_H2O_vpr_dlt=q_H2O_vpr(lon_idx)-ssh_H2O_sfc ! [kg kg-1] Humidity change
            hgt_rat_log=log(hgt_mdp(lon_idx)/hgt_rfr_LaP81) ! [frc] Log of ratio of local to reference height
            spc_heat_mst_air=spc_heat_dry_air*(1.0+cp_vpr_rcp_cp_dry_m1*ssh_H2O_sfc) ! [J kg-1 K-1] Specific heat of moist air
            
c     Iteration loop
            do itr_idx=1,itr_max

c     Save old friction speed for convergence diagnostic
               if (itr_idx == 1) then
                  wnd_frc_old=0.0 ! [m s-1]
               else 
                  wnd_frc_old=wnd_frc(lon_idx) ! [m s-1]
               endif            ! endif
               

c     Roots of neutral, 10 m exchange coefficients
               xch_cff_mmn_ntr_sqrt=cst_von_krm/log(hgt_rfr_LaP81/rgh_mmn(lon_idx)) ! [frc] Squareroot of neutral 10 m drag coefficient
               xch_cff_heat_ntr_sqrt=xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
               xch_cff_vpr_ntr_sqrt=xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange


               if (itr_idx == 1) then
c     First iteration: Use estimated roots of neutral exchange coefficients
                  wnd_frc(lon_idx)=xch_cff_mmn_ntr_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
                  tpt_scl(lon_idx)=xch_cff_heat_ntr_sqrt*tpt_dlt ! [K] Temperature scale
                  vpr_scl(lon_idx)=xch_cff_vpr_ntr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
               else
c     Subsequently: Use current roots of (non-neutral) exchange coefficients at measurement height
                  wnd_frc(lon_idx)=xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
                  tpt_scl(lon_idx)=xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale
                  vpr_scl(lon_idx)=xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
                  
               endif            ! endif first iteration

               
c     Compute stability parameter at midlayer and evaluate stability corrections  
c     Stable if tpt_ptn_mdp > tpt_sfc or if mno_stb_prm > 0
               mno_stb_prm(lon_idx)= ! [frc] Monin-Obukhov stability parameter 
     $              cst_von_krm*grv_sfc*hgt_mdp(lon_idx)
     $              *(tpt_scl(lon_idx)/tpt_ptn_vrt_mdp+vpr_scl(lon_idx)/(1.0/eps_H2O_rcp_m1+q_H2O_vpr(lon_idx)))
     $              /(wnd_frc(lon_idx)*wnd_frc(lon_idx))
c     Ensure |z/L| < 10.0 because similarity function asymptotes to (-z/L)^(-1/3) for z/L << -1, Ary88 p. 166

               mno_stb_prm(lon_idx)=sign(min(abs(mno_stb_prm(lon_idx)),ten),mno_stb_prm(lon_idx)) ! [frc] Monin-Obukhov stability parameter
               stb_val=half+sign(half,mno_stb_prm(lon_idx)) ! [flg] 1.0 if stable, 0.0 if unstable
               sml_fnc_mmn_uns_rcp=max(sqrt(abs(1.0-16.0*mno_stb_prm(lon_idx))),one) ! [frc] BKL97 p. F1, LaP81 p. 325
               sml_fnc_mmn_uns_rcp=sqrt(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
               mno_stb_crc_mmn=-5.0*mno_stb_prm(lon_idx)*stb_val+(1.0-stb_val)*mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
               mno_stb_crc_heat=-5.0*mno_stb_prm(lon_idx)*stb_val+(1.0-stb_val)*mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
               
c     Shift old neutral 10 m exchange coefficients to measurement height and stability
               xch_cff_mmn_sqrt=xch_cff_mmn_ntr_sqrt/(1.0+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)) ! [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
               xch_cff_heat_sqrt=xch_cff_heat_ntr_sqrt/(1.0+xch_cff_heat_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Stanton number for heat exchange LaP82 p. 466 (10)
               xch_cff_vpr_sqrt=xch_cff_vpr_ntr_sqrt/(1.0+xch_cff_vpr_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Dalton number for vapor exchange LaP82 p. 466 (10)
               

#ifdef UNDEFINED
               if (.true.) then
                  mno_lng(lon_idx)=hgt_mdp(lon_idx)/mno_stb_prm(lon_idx) ! [m] Monin-Obukhov length
                  wnd_str(lon_idx)=dns_mdp(lon_idx)*wnd_frc(lon_idx)*wnd_frc(lon_idx) ! [kg m-1 s-2] Wind stress
                  flx_sns(lon_idx)=-spc_heat_mst_air*wnd_str(lon_idx)*tpt_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Sensible heat flux to atmosphere
                  flx_ltn(lon_idx)=-ltn_heat_trn*wnd_str(lon_idx)*vpr_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Latent heat flux to atmosphere
                  flx_LW_up_sfc(lon_idx)=cst_Stefan_Boltzmann*(tpt_sfc(lon_idx)**4.0) ! [W m-2] Longwave upwelling flux to atmosphere
                  write (6,'(i3,1x,f12.6,1x,f12.6,1x,f12.6,1x,f12.6,1x,f12.6,1x,f12.6,1x,f12.6,1x,f12.6,1x,f12.6)')
     $                 itr_idx,mno_lng(lon_idx),mno_stb_prm(lon_idx),
     $                 wnd_frc(lon_idx),wnd_rfr_ntr(lon_idx),
     $                 1000.0*xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt,
     $                 1000.0*xch_cff_mmn_sqrt*xch_cff_mmn_sqrt,
     $                 flx_sns(lon_idx),flx_ltn(lon_idx),flx_LW_up_sfc(lon_idx)
               endif            ! endif dbg
#endif /* not UNDEFINED */

            enddo               ! end loop over itr


                  
            
c     Update wnd_frc, tpt_scl, and vpr_scl using updated exchange coefficients at measurement height
            wnd_frc(lon_idx)=xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity 
            tpt_scl(lon_idx)=xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale 
            vpr_scl(lon_idx)=xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale 

c      If wnd_frc comes out negative, something is wrong.
            if ( wnd_frc(lon_idx).lt.0.0 ) then
               write(6,*)'ABORT: blm_lnd: wnd_frc ',wnd_frc(lon_idx),' < 0 in blmutl.F at lon ',lon_idx
               call abort()
            end if


c     Compute surface stress components
c            wnd_str(lon_idx)=dns_mdp(lon_idx)*wnd_frc(lon_idx)*wnd_frc(lon_idx) ! [kg m-1 s-2] Wind stress
c            wnd_str_znl(lon_idx)=-wnd_str(lon_idx)*wnd_znl_mdp(lon_idx)/wnd_mdp_bnd ! [kg m-1 s-2] Zonal wind stress
c            wnd_str_mrd(lon_idx)=-wnd_str(lon_idx)*wnd_mrd_mdp(lon_idx)/wnd_mdp_bnd ! [kg m-1 s-2] Meridional wind stress
            
c     Compute heat flux components at current surface temperature
c     Define positive latent and sensible heat as upwards into atmosphere 
c            flx_sns(lon_idx)=-spc_heat_mst_air*wnd_str(lon_idx)*tpt_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Sensible heat flux to atmosphere
c            flx_ltn(lon_idx)=-ltn_heat_trn*wnd_str(lon_idx)*vpr_scl(lon_idx)/wnd_frc(lon_idx) ! [W m-2] Latent heat flux to atmosphere
c            flx_LW_up_sfc(lon_idx)=cst_Stefan_Boltzmann*(tpt_sfc(lon_idx)**4.0) ! [W m-2] Longwave upwelling flux to atmosphere
            
c     Additional diagnostics
            mno_lng(lon_idx)=hgt_mdp(lon_idx)/mno_stb_prm(lon_idx) ! [m] Monin-Obukhov length
         endif                  ! endif flg_lnd
      enddo                     ! end loop over lon
      

      return
      end subroutine blm_lnd                      ! end blm_lnd()

      subroutine blm_glb(       ! Solve boundary layer meteorology on global scale
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     hgt_mdp,             ! I [m] Midlayer height above surface
     $     oro,                 ! I [frc] Orography
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_vpr,           ! I [kg kg-1] Specific humidity
     $     sfc_typ,             ! I [idx] LSM surface type (0..28)
     $     snw_hgt_lqd,         ! I [m] Equivalent liquid water snow depth
     $     tpt_mdp,             ! I [K] Midlayer temperature
     $     tpt_ptn_mdp,         ! I [K] Potential temperature
     $     tpt_sfc,             ! I [K] Surface temperature
     $     wnd_mrd_mdp,         ! I [m s-1] Surface layer meridional wind speed
     $     wnd_znl_mdp,         ! I [m s-1] Surface layer zonal wind speed
     $     hgt_zpd,             ! O [m] Zero plane displacement
     $     mno_lng,             ! O [m] Monin-Obukhov length
     $     rgh_mmn,             ! O [m] Roughness length momentum
     $     snw_frc,             ! O [frc] Fraction of surface covered by snow
     $     wnd_frc,             ! O [m s-1] Surface friction velocity
     $     wnd_rfr)             ! O [m s-1] Wind speed at reference height

c     Purpose: Given vectors of state variables and surface type, 
c     divide vectors in land, ocean, and sea-ice components and call the routines 
c     to compute the boundary layer exchange properties
c     blm_glb() is called by dstdpsdry()
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
c     Parameters
      real(r8),parameter::hgt_rfr=10.0   ! [m] Reference height for deposition processes
      real(r8),parameter::wnd_min_dps=1.0 ! [m s-1] Minimum windspeed for momentum exchange
c     Input
      integer sfc_typ(plond)    ! I [idx] LSM surface type (0..28)
      real(r8) oro(plond)           ! I [frc] Orography
      real(r8) dns_mdp(plond)       ! I [kg m-3] Midlayer density
      real(r8) hgt_mdp(plond)       ! I [m] Midlayer height above surface
      real(r8) prs_mdp(plond)       ! I [Pa] Pressure
      real(r8) q_H2O_vpr(plond)     ! I [kg kg-1] Specific humidity
      real(r8) snw_hgt_lqd(plond)   ! I [m] Equivalent liquid water snow depth
      real(r8) tpt_mdp(plond)       ! I [K] Midlayer temperature
      real(r8) tpt_ptn_mdp(plond)   ! I [K] Midlayer local potential temperature
      real(r8) tpt_sfc(plond)       ! I [K] Surface temperature
      real(r8) wnd_mrd_mdp(plond)   ! I [m s-1] Surface layer meridional wind speed
      real(r8) wnd_znl_mdp(plond)   ! I [m s-1] Surface layer zonal wind speed
c     Input/Output
c     Output
      real(r8) hgt_zpd(plond)       ! O [m] Zero plane displacement
      real(r8) mno_lng(plond)       ! O [m] Monin-Obukhov length
      real(r8) rgh_mmn(plond)       ! O [m] Roughness length momentum
      real(r8) snw_frc(plond)       ! O [frc] Fraction of surface covered by snow
      real(r8) wnd_frc(plond)       ! O [m s-1] Surface friction velocity
      real(r8) wnd_rfr(plond)       ! O [m s-1] Wind speed at reference height
c     Local
      integer lon_idx           ! [idx] Counting index for lon

      logical flg_ice(plond)    ! [flg] Sea ice flag
      logical flg_lnd(plond)    ! [flg] Land flag
      logical flg_ocn(plond)    ! [flg] Ocean flag
      logical flg_true(plond)   ! [flg] True flag
      real(r8) wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
c     Statement functions
      logical oro_is_ice        ! [flg] True if > 50% sea ice
      logical oro_is_lnd        ! [flg] True if > 50% land
      logical oro_is_ocn        ! [flg] True if > 50% ocean
      real(r8) oro_val              ! [frc] Orography
      oro_is_ice(oro_val)=nint(oro_val) == 2
      oro_is_lnd(oro_val)=nint(oro_val) == 1
      oro_is_ocn(oro_val)=nint(oro_val) == 0

   
c     Main Code

c     Construct land, ocean, and sea-ice vectors
      do lon_idx=1,plon
         if (oro_is_lnd(oro(lon_idx))) then
            flg_lnd(lon_idx)=.true.
         else 
            flg_lnd(lon_idx)=.false.
         endif                  ! endif land
         if (oro_is_ocn(oro(lon_idx))) then
            flg_ocn(lon_idx)=.true.
         else 
            flg_ocn(lon_idx)=.false.
         endif                  ! endif ocean
         if (oro_is_ice(oro(lon_idx))) then
            flg_ice(lon_idx)=.true.
         else 
            flg_ice(lon_idx)=.false.
         endif                  ! endif ice
#ifdef DST_DBG
c     Sanity check
         if (.not.flg_lnd(lon_idx).and..not.flg_lnd(lon_idx).and..not.flg_lnd(lon_idx)) stop 'Invalid surface type in blm_glb()'
#endif /* not DST_DBG */
      enddo                     ! end loop over lon

c     Midlayer wind speeds
      do lon_idx=1,plon
         flg_true(lon_idx)=.true. ! [flg] True flag
         wnd_mdp(lon_idx)=      ! [m s-1] Surface layer mean wind speed
     $        sqrt(wnd_znl_mdp(lon_idx)*wnd_znl_mdp(lon_idx)+ 
     $        wnd_mrd_mdp(lon_idx)*wnd_mrd_mdp(lon_idx)) 
      end do                    ! end loop over lon

c     Fraction of surface covered by snow
      call snw_frc_get(
     $     snw_hgt_lqd,         ! I [m] Equivalent liquid water snow depth
     $     snw_frc)             ! O [frc] Fraction of surface covered by snow

c     Zero plane displacement
      call zpd_get(
     $     hgt_zpd,             ! O [m] Zero plane displacement
     $     oro,                 ! I [frc] Orography
     $     sfc_typ,             ! I [idx] LSM surface type (0..28)
     $     snw_frc)             ! I [frc] Fraction of surface covered by snow

c     Roughness length 
      call rgh_mmn_get(
     $     oro,                 ! I [frc] Orography
     $     rgh_mmn,             ! O [m] Roughness length momentum
     $     sfc_typ,             ! I [idx] LSM surface type (0..28)
     $     snw_frc,             ! I [frc] Fraction of surface covered by snow
     $     wnd_mdp)             ! I [m s-1] Surface layer mean wind speed
     

      call blm_lnd(
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     flg_lnd,             ! I [flg] Land flag
     $     hgt_mdp,             ! I [m] Midlayer height above surface
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_vpr,           ! I [kg kg-1] Specific humidity
     $     rgh_mmn,             ! I [m] Roughness length momentum
     $     tpt_mdp,             ! I [K] Midlayer temperature
     $     tpt_ptn_mdp,         ! I [K] Midlayer local potential temperature
     $     tpt_sfc,             ! I [K] Surface temperature
     $     wnd_mrd_mdp,         ! I [m s-1] Surface layer meridional wind speed
     $     wnd_znl_mdp,         ! I [m s-1] Surface layer zonal wind speed
     $     mno_lng,             ! O [m] Monin-Obukhov length
     $     wnd_frc)             ! O [m s-1] Surface friction velocit

      call blm_ocn(
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     flg_ocn,             ! I [flg] Ocean flag
     $     hgt_mdp,             ! I [m] Midlayer height above surface
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_vpr,           ! I [kg kg-1] Specific humidity
     $     tpt_mdp,             ! I [K] Midlayer temperature
     $     tpt_ptn_mdp,         ! I [K] Midlayer local potential temperature
     $     tpt_sfc,             ! I [K] Surface temperature
     $     wnd_mrd_mdp,         ! I [m s-1] Surface layer meridional wind speed
     $     wnd_znl_mdp,         ! I [m s-1] Surface layer zonal wind speed
     $     mno_lng,             ! O [m] Monin-Obukhov length
     $     rgh_mmn,             ! O [m] Roughness length momentum
     $     wnd_frc)             ! O [m s-1] Surface friction velocit

      call blm_ice(
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     flg_ice,             ! I [flg] Sea ice flag
     $     hgt_mdp,             ! I [m] Midlayer height above surface
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_vpr,           ! I [kg kg-1] Specific humidity
     $     rgh_mmn,             ! I [m] Roughness length momentum
     $     tpt_mdp,             ! I [K] Midlayer temperature
     $     tpt_ptn_mdp,         ! I [K] Midlayer local potential temperature
     $     tpt_sfc,             ! I [K] Surface temperature
     $     wnd_mrd_mdp,         ! I [m s-1] Surface layer meridional wind speed
     $     wnd_znl_mdp,         ! I [m s-1] Surface layer zonal wind speed
     $     mno_lng,             ! O [m] Monin-Obukhov length
     $     wnd_frc)             ! O [m s-1] Surface friction velocit

c     Interpolate midlayer wind speed to 10 m
      call wnd_rfr_get(
     $     flg_true,            ! I [flg] True flag
     $     hgt_mdp,             ! I [m] Midpoint height above surface
     $     hgt_rfr,             ! I [m] Reference height for deposition processes
     $     hgt_zpd,             ! I [m] Zero plane displacement
     $     mno_lng,             ! I [m] Monin-Obukhov length
     $     wnd_frc,             ! I [m s-1] Surface friction velocity
     $     wnd_mdp,             ! I [m s-1] Surface layer mean wind speed
     $     wnd_min_dps,         ! I [m s-1] Minimum windspeed used for deposition 
     $     wnd_rfr)             ! O [m s-1] Wind speed at reference height
      
      return
      end subroutine blm_glb                      ! end blm_glb()

      end module blmutl
