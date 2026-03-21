      module dstmblutl

      use precision

      real(r8) :: small = 5*epsilon(1._r8)
      real(r8) :: one = 1.


      contains

c     $Header$
      
c     Purpose: dstmblutl.F contains routines to compute mineral dust source tendencies
      
      subroutine wnd_frc_thr_slt_get(
     $     dmt_aer,             ! I [m] Particle diameter (Currently not used)
     $     dns_aer,             ! I [kg m-3] Particle density (Currently not used)
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     wnd_frc_thr_slt)     ! O [m s-1] Threshold friction velocity for saltation
c     Purpose: Compute dry threshold friction velocity for saltation
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
c     Parameters
      real(r8),parameter::dmt_slt_opt=75.0e-6 ! [m] Optimal diameter for saltation, IvW82 p. 117 Fgr. 8, Pye87 p. 31, MBA97 p. 4388, SRL96 (2)
      real(r8),parameter::dns_slt=2650.0 ! [kg m-3] Density of optimal saltation particles, MBA97 p. 4388 
c     Input
      real(r8),intent(in)::dns_aer(dst_nbr)     ! [kg m-3] Particle density
      real(r8),intent(in)::dmt_aer(dst_nbr)     ! [m] Particle diameter
      real(r8),intent(in)::dns_mdp(plond)       ! [kg m-3] Midlayer density
c     Output
      real(r8),intent(out)::wnd_frc_thr_slt(plond) ! [m s-1] Threshold friction velocity for saltation
c     Local
      integer lon_idx           ! [idx] Counting index
      real(r8) ryn_nbr_frc_thr_prx_opt ! [frc] Threshold friction Reynolds number approximation for optimal size
      real(r8) ryn_nbr_frc_thr_opt_fnc ! [frc] Threshold friction Reynolds factor for saltation calculation
      real(r8) dns_fct              ! Density ratio factor for saltation calculation
      real(r8) icf_fct              ! Interparticle cohesive forces factor for saltation calculation
      real(r8) tmp1                 ! Factor in saltation computation
c     Main Code
c     Initialize some variables
c     MaB95 pzn. for Re*t(D_opt) circumvents iterative solution
      ryn_nbr_frc_thr_prx_opt=0.38+1331.0*(100.0*dmt_slt_opt)**1.56 ! [frc] "B" MaB95 p. 16417 (5)
c     Given Re*t(D_opt), compute time independent factors contributing to u*t
      icf_fct=1.0+6.0e-07/(dns_slt*grv_sfc*(dmt_slt_opt**2.5)) ! [frc] IvW82 p. 115 (6) MaB95 p. 16417 (4) Interparticle cohesive forces
      dns_fct=dns_slt*grv_sfc*dmt_slt_opt ! IvW82 p. 115 (6) MaB95 p. 16417 (4)
      if (ryn_nbr_frc_thr_prx_opt < 0.03) then
         stop 'dst: wnd_frc_thr_slt_get() reports ryn_nbr_frc_thr_prx_opt < 0.03'
      else if (ryn_nbr_frc_thr_prx_opt < 10.0) then
         ryn_nbr_frc_thr_opt_fnc=-1.0+1.928*(ryn_nbr_frc_thr_prx_opt**0.0922) ! [frc] IvW82 p. 114 (3), MaB95 p. 16417 (6)
         ryn_nbr_frc_thr_opt_fnc=0.1291*0.1291/ryn_nbr_frc_thr_opt_fnc ! [frc] 
      else 
         ryn_nbr_frc_thr_opt_fnc=1.0-0.0858*exp(-0.0617*(ryn_nbr_frc_thr_prx_opt-10.0)) ! [frc] IvW82 p. 114 (3), MaB95 p. 16417 (7)
         ryn_nbr_frc_thr_opt_fnc=0.120*0.120*ryn_nbr_frc_thr_opt_fnc*ryn_nbr_frc_thr_opt_fnc! [frc]
      endif                     ! endif
c     This method minimizes the number of square root computations performed
      tmp1=sqrt(icf_fct*dns_fct*ryn_nbr_frc_thr_opt_fnc)
      do lon_idx=1,plon
         wnd_frc_thr_slt(lon_idx)=tmp1/sqrt(dns_mdp(lon_idx)) ! [m s-1] Threshold friction velocity for saltation dry ground 
      end do                    ! end loop over lon
      return
      end  subroutine wnd_frc_thr_slt_get                     ! end wnd_frc_thr_slt_get()
      
      subroutine wnd_rfr_thr_slt_get(
     $     wnd_frc,             ! I [m s-1] Surface friction velocity
     $     wnd_frc_thr_slt,     ! I [m s-1] Threshold friction velocity for saltation
     $     wnd_mdp,             ! I [m s-1] Surface layer mean wind speed
     $     wnd_rfr,             ! I [m s-1] Wind speed at reference height
     $     wnd_rfr_thr_slt)     ! O [m s-1] Threshold 10 m wind speed for saltation
c     Purpose: Compute threshold 10 m wind speed
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
c     Parameters
c     Input
      real(r8),intent(in)::wnd_frc(plond)       ! [m s-1] Friction velocity
      real(r8),intent(in)::wnd_frc_thr_slt(plond) ! [m s-1] Threshold friction velocity for saltation
      real(r8),intent(in)::wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
      real(r8),intent(in)::wnd_rfr(plond)       ! [m s-1] Wind speed at reference height
c     Output
      real(r8),intent(out)::wnd_rfr_thr_slt(plond) ! [m s-1] Threshold 10 m wind speed for saltation
c     Local
      integer lon_idx           ! [idx] Counting index
c     Main Code
c     Compute threshold horizontal wind speed at reference height
      do lon_idx=1,plon
c     A more complicated procedure would recompute mno_lng for wnd_frc_thr,
c     and then integrate vertically from rgh_mmn+hgt_zpd to hgt_rfr
c     wnd_crc_fct is (1/k)*[ln(z-D)/z0 - psi(zeta2) + psi(zeta1)]
         wnd_rfr_thr_slt(lon_idx)=wnd_frc_thr_slt(lon_idx)*wnd_rfr(lon_idx)/wnd_frc(lon_idx) ! [m s-1]
      end do                    ! end loop over lon
      return
      end subroutine wnd_rfr_thr_slt_get                      ! end wnd_rfr_thr_slt_get()
      
      subroutine frc_thr_wet_fct_get(
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     frc_thr_wet_fct,     ! O [frc] Factor by which moisture increases threshold friction velocity
     $     mss_frc_cly,         ! I [frc] Mass fraction of clay
     $     vwc_sfc)             ! I [frc] Volumetric water content
c     Purpose: Compute factor by which soil moisture increases threshold friction velocity 
c     This parameterization is based on FMB99
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
c     Parameters
c     Input
      logical,intent(in)::flg_mbl(plond) ! I [flg] Mobilization candidate flag
      real(r8),intent(in)::mss_frc_cly(plond) ! [frc] Mass fraction of clay
      real(r8),intent(in)::vwc_sfc(plond) ! [frc] Volumetric water content
c     Output
      real(r8),intent(out)::frc_thr_wet_fct(plond) ! [frc] Factor by which moisture increases threshold friction velocity
c     Local
      integer lon_idx           ! [idx] Counting index
      real(r8) vwc_thr(plond)       ! [frc] Threshold volumetric water content
c     Main Code
c     Initialize output
      frc_thr_wet_fct=1.0       ! [frc]
c     Compute threshold soil moisture based on clay content
c     vwc_thr=mss_frc_cly*(0.17+0.14*mss_frc_cly) ! [m3 m-3] FMB99 p. 155 (14)
c     fxm: 19991105 remove factor of mss_frc_cly from vwc_thr to improve large scale behavior
      vwc_thr(1:plon)=0.17+0.14*mss_frc_cly(1:plon) ! [m3 m-3] 
c     Adjust threshold velocity for inhibition by moisture
c     frc_thr_wet_fct(lon_idx)=exp(22.7*vwc_sfc(lon_idx)) ! [frc] SRL96
      where (vwc_sfc(1:plon) > vwc_thr(1:plon).and.flg_mbl(1:plon))
     $     frc_thr_wet_fct(1:plon)=sqrt(1.0+1.21*(100.0*(vwc_sfc(1:plon)-vwc_thr(1:plon)))**0.68) ! [frc] FMB99 p. 155 (15)
c     Uncomment the following line to remove all dependence on vwc_sfc
c     frc_thr_wet_fct(lon_idx)=1.0 ! [frc]
      return
      end subroutine frc_thr_wet_fct_get                      ! end frc_thr_wet_fct_get()
      
      subroutine frc_thr_rgh_fct_get(
     $     frc_thr_rgh_fct,     ! O [frc] Factor by which surface roughness increases threshold friction velocity
     $     rgh_mmn_mbl,         ! I [m] Roughness length momentum for erodible surfaces
     $     rgh_mmn_smt)         ! I [m] Smooth roughness length
c     Purpose: Compute factor by which surface roughness increases threshold friction velocity
c     This parameterization is based on MaB95 and GMB98
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
c     Parameters
c     Input
      real(r8),intent(in)::rgh_mmn_mbl ! [m] Roughness length momentum for erodible surfaces
      real(r8),intent(in)::rgh_mmn_smt ! [m] Smooth roughness length
c     Output
      real(r8),intent(out)::frc_thr_rgh_fct(plond) ! [frc] Factor by which roughness increases threshold friction velocity
c     Local
      integer lon_idx           ! [idx] Counting index
      real(r8) wnd_frc_fsh_frc      ! [frc] Efficient fraction of wind friction
      real(r8) wnd_frc_fsh_frc_rcp  ! [frc] Reciprocal of wnd_frc_fsh_frc
c     Main Code
c     Adjust threshold velocity for inhibition by roughness elements
      wnd_frc_fsh_frc=          ! [frc] MaB95 p. 16420, GMB98 p. 6207
     $     +1.0-log(rgh_mmn_mbl/rgh_mmn_smt)/log(0.35*((0.1/rgh_mmn_smt)**0.8))
      if (wnd_frc_fsh_frc <= 0.0.or.wnd_frc_fsh_frc > 1.0) stop 
     $     'dst: ERROR frc_thr_rgh_fct_get() reports wnd_frc_fsh_frc out of range'
      wnd_frc_fsh_frc_rcp=1.0/wnd_frc_fsh_frc ! [frc] 
      frc_thr_rgh_fct=wnd_frc_fsh_frc_rcp ! [frc]
c     fxm: 19991012 Set frc_thr_rgh_fct=1.0, equivalent to assuming mobilization takes place at smooth roughness length
      frc_thr_rgh_fct=1.0
      return
      end  subroutine frc_thr_rgh_fct_get                     ! end frc_thr_rgh_fct_get()

      subroutine wnd_frc_slt_get(
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     wnd_frc,             ! I [m s-1] Surface friction velocity
     $     wnd_frc_slt,         ! O [m s-1] Saltating friction velocity
     $     wnd_rfr,             ! I [m s-1] Wind speed at reference height
     $     wnd_rfr_thr_slt)     ! I [m s-1] Threshold 10 m wind speed for saltation
c     Purpose: Compute the saltating friction velocity
c     Saltation increases friction speed by roughening surface, AKA "Owen's effect"
c     This acts as a positive feedback to the friction speed
c     GMB98 parameterized this feedback in terms of 10 m windspeeds
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
c     Parameters
c     Input
      logical,intent(in)::flg_mbl(plond)    ! I [flg] Mobilization candidate flag
      real(r8),intent(in)::wnd_frc(plond)       ! I [m s-1] Surface friction velocity
      real(r8),intent(in)::wnd_rfr(plond)       ! I [m s-1] Wind speed at reference height
      real(r8),intent(in)::wnd_rfr_thr_slt(plond) ! I [m s-1] Threshold 10 m wind speed for saltation
c     Output
      real(r8),intent(out)::wnd_frc_slt(plond)   ! O [m s-1] Saltating friction velocity
c     Local
      integer lon_idx           ! [idx] Counting index
      real(r8) wnd_rfr_dlt          ! [m s-1] Reference windspeed excess over threshold
      real(r8) wnd_frc_slt_dlt      ! [m s-1] Friction velocity increase from saltation
c     Main Code
c     Compute saltating friction velocity, accounting for "Owen's effect"
      wnd_frc_slt=wnd_frc       ! [m s-1] Saltating friction velocity
      do lon_idx=1,plon
         if (flg_mbl(lon_idx).and.wnd_rfr(lon_idx) >= wnd_rfr_thr_slt(lon_idx)) then
c     Saltation roughens the boundary layer, AKA "Owen's effect"
c     GMB98 p. 6206 Fig. 1 shows observed/computed u* dependence on observed U(1 m)
c     GMB98 eqn(12) p. 6209 has u* in cm s-1 and U, Ut in m s-1, personal communication, D. Gillette, 19990529
c     With everything in MKS, the 0.3 coefficient in GMB98 (12) becomes 0.003 
c     Increase in friction velocity due to saltation varies as square of 
c     difference between reference wind speed and reference threshold speed 
            wnd_rfr_dlt=wnd_rfr(lon_idx)-wnd_rfr_thr_slt(lon_idx)
            wnd_frc_slt_dlt=0.003*wnd_rfr_dlt*wnd_rfr_dlt ! [m s-1] Friction velocity increase from saltation GMB98 p. 6209
            wnd_frc_slt(lon_idx)=wnd_frc(lon_idx)+wnd_frc_slt_dlt ! [m s-1] Saltating friction velocity
         endif                  ! endif wnd_frc_mbl > wnd_frc_thr_slt
      end do                    ! end loop over lon
      return
      end subroutine wnd_frc_slt_get                      ! end wnd_frc_slt_get()

      subroutine wnd_rfr_get(
     $     flg_oro,             ! I [flg] Orography flag
     $     hgt_mdp,             ! I [m] Midpoint height above surface
     $     hgt_rfr,             ! I [m] Reference height
     $     hgt_zpd,             ! I [m] Zero plane displacement
     $     mno_lng,             ! I [m] Monin-Obukhov length
     $     wnd_frc,             ! I [m s-1] Surface friction velocity
     $     wnd_mdp,             ! I [m s-1] Surface layer mean wind speed
     $     wnd_min,             ! I [m s-1] Minimum windspeed
     $     wnd_rfr)             ! O [m s-1] Wind speed at reference height
c     Purpose: Interpolate wind speed at given height to wind speed at reference height
      implicit none
#include <pmgrid.h>
#include <dstcst.h> /* Physical constants for dust routines */
c     Parameters
      integer,parameter::rfr_hgt_idx=1 ! Named index for lower (target) height
      integer,parameter::gcm_hgt_idx=2 ! Named index for upper (known) height
c     Input
      logical,intent(in)::flg_oro(plond)    ! I [flg] Orography flag
      real(r8),intent(in)::hgt_mdp(plond)       ! I [m] Midlayer height above surface
      real(r8),intent(in)::hgt_rfr              ! I [m] Reference height
      real(r8),intent(in)::hgt_zpd(plond)       ! I [m] Zero plane displacement
      real(r8),intent(in)::mno_lng(plond)       ! I [m] Monin-Obukhov length
      real(r8),intent(in)::wnd_frc(plond)       ! I [m s-1] Surface friction velocity
      real(r8),intent(in)::wnd_mdp(plond)       ! I [m s-1] Surface layer mean wind speed
      real(r8),intent(in)::wnd_min              ! I [m s-1] Minimum windspeed
c     Output
      real(r8),intent(out)::wnd_rfr(plond)       ! O [m s-1] Wind speed at reference height
c     Local
      integer idx_idx           ! [idx] Counting index
      integer lon_idx           ! [idx] Counting index for lon
      integer lvl_idx           ! Stability computation loop index
      integer vld_idx(plond)    ! Valid indices
      integer vld_nbr           ! [nbr] Number of valid indices
      real(r8) mno_stb_crc_mmn(plond,2) ! [frc] Monin-Obukhov stability correction momentum
      real(r8) mno_stb_prm(plond,2) ! [frc] Monin-Obukhov stability parameter 
      real(r8) sml_fnc_mmn_uns_rcp  ! [frc] Reciprocal of similarity function for momentum, unstable atmosphere
      real(r8) tmp2                 ! Term in stability correction computation
      real(r8) tmp3                 ! Term in stability correction computation
      real(r8) tmp4                 ! Term in stability correction computation
      real(r8) wnd_crc_fct(plond)   ! [frc] Wind correction factor
      real(r8) hgt_rfr_rcp      ! [m-1] Reciprocal of reference height
c     Initialize
      hgt_rfr_rcp=1.0/hgt_rfr ! [m-1]
c     Initialize output
      wnd_rfr=wnd_min           ! [m s-1]
      
c     Compute horizontal wind speed at reference height
      do lon_idx=1,plon
         if (flg_oro(lon_idx).and.hgt_zpd(lon_idx) < hgt_rfr) then
c     Code uses notation of Bon96 p. 50, where lvl_idx=1 is 10 m ref. hgt, lvl_idx=2 is atm. hgt
            mno_stb_prm(lon_idx,rfr_hgt_idx)=min((hgt_rfr-hgt_zpd(lon_idx))/mno_lng(lon_idx),one) ! [frc]
            mno_stb_prm(lon_idx,gcm_hgt_idx)=min((hgt_mdp(lon_idx)-hgt_zpd(lon_idx))/mno_lng(lon_idx),one) ! [frc]
            do lvl_idx=1,2
               if (mno_stb_prm(lon_idx,lvl_idx) < 0.0) then
                  sml_fnc_mmn_uns_rcp=(1.0-16.0*mno_stb_prm(lon_idx,lvl_idx))**0.25
                  tmp2=log((1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0)
                  tmp3=log((1.0+sml_fnc_mmn_uns_rcp)/2.0)
                  mno_stb_crc_mmn(lon_idx,lvl_idx)=2.0*tmp3+tmp2-2.0*atan(sml_fnc_mmn_uns_rcp)+1.5707963
               else             ! not stable
                  mno_stb_crc_mmn(lon_idx,lvl_idx)=-5.0*mno_stb_prm(lon_idx,lvl_idx)
               endif            ! stable
            enddo               ! end loop over lvl_idx
            tmp4=log((hgt_mdp(lon_idx)-hgt_zpd(lon_idx))/(hgt_rfr-hgt_zpd(lon_idx)))
c     Correct neutral stability assumption
            wnd_crc_fct(lon_idx)=tmp4-mno_stb_crc_mmn(lon_idx,gcm_hgt_idx)+mno_stb_crc_mmn(lon_idx,rfr_hgt_idx) ! [frc]
            wnd_rfr(lon_idx)=wnd_mdp(lon_idx)-wnd_frc(lon_idx)*cst_von_krm_rcp*wnd_crc_fct(lon_idx) ! [m s-1]
            wnd_rfr(lon_idx)=max(wnd_rfr(lon_idx),wnd_min) ! [m s-1]
c     Use neutral stability assumption
c     wnd_rfr(lon_idx)=wnd_mdp(lon_idx)-wnd_frc(lon_idx)*cst_von_krm_rcp*tmp4
         endif                  ! endif (hgt_zpd(lon_idx) < hgt_rfr.and.flg_oro(lon_idx))
      end do                    ! end loop over lon
      
      return
      end subroutine wnd_rfr_get                      ! end wnd_rfr_get()
      
      subroutine snw_frc_get(
     $     snw_hgt_lqd,         ! I [m] Equivalent liquid water snow depth
     $     snw_frc              ! O [frc] Fraction of surface covered by snow
     $     )
c     Convert equivalent liquid water snow depth to fractional snow cover
c     Use algorithm of Bon96
      implicit none
#include <pmgrid.h>
c     Parameters
      real(r8),parameter::snw_dns=250.0 ! [kg m-3] Bulk density of snow, CCM:lsm/snoconi.F (NB: CCM:physics/tsinti() uses 100.0)
      real(r8),parameter::snw_hgt_thr=0.05 ! [m] Geometric snow thickness for 100% coverage, CCM:lsm/snoconi.F 
c     Input
      real(r8),intent(in)::snw_hgt_lqd(plond)   ! [m] Equivalent liquid water snow depth
c     Output
      real(r8),intent(out)::snw_frc(plond)       ! [frc] Fraction of surface covered by snow
c     Local
      integer lon_idx           ! [idx] Counting index for lon
      real(r8) snw_hgt(plond)       ! [m] Geometric bulk thickness of snow
c     Main code
c     Initialize some variables
c     Fractional snow cover
      do lon_idx=1,plon
         snw_hgt(lon_idx)=snw_hgt_lqd(lon_idx)*4.0 ! [m] NB: CCM and LSM seem to disagree on this
         snw_frc(lon_idx)=min(snw_hgt(lon_idx)/snw_hgt_thr,one) ! [frc]
      end do                    ! end loop over lon
      return
      end subroutine snw_frc_get                      ! end snw_frc_get()
      
      subroutine flx_mss_CaCO3_msk(
     $     dmt_vwr,             ! I [m] Mass weighted diameter resolved 
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     flx_mss_vrt_dst_CaCO3, ! I/O [kg m-2 s-1] Vertical mass flux of dust on input, CaCO3 on output
     $     mss_frc_CaCO3,       ! I [frc] Mass fraction of CaCO3
     $     mss_frc_cly,         ! I [frc] Mass fraction of clay
     $     mss_frc_snd)         ! I [frc] Mass fraction of sand
c     Purpose: Mask dust mass flux by CaCO3 mass fraction at source
c     Theory: Uses soil CaCO3 mass fraction from Global Soil Data Task, 1999 (Sch99)
c     Uses size dependent apportionment of CaCO3 from Claquin et al, 1999 (CSB99)
c     Dependencies: <pmgrid.h> <dstgrd.h>
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
c     Parameters
      real(r8),parameter::dmt_cly_max=2.0e-6 ! [m] Maximum diameter of Clay soil texture CSB99 p. 22250
      real(r8),parameter::dmt_slt_max=50.0e-6 ! [m] Maximum diameter of Silt soil texture CSB99 p. 22250
      real(r8),parameter::dns_CaCO3=2950.0 ! [kg m-3] Density of CaCO3 http://www.ssc.on.ca/mandm/calcit.htm
c     Input
      logical,intent(in)::flg_mbl(plond) ! I [flg] Mobilization candidate flag
      real(r8),intent(in)::dmt_vwr(dst_nbr) ! I [m] Mass weighted diameter resolved
      real(r8),intent(in)::mss_frc_CaCO3(plond) ! I [frc] Mass fraction of CaCO3
      real(r8),intent(in)::mss_frc_cly(plond) ! I [frc] Mass fraction of clay
      real(r8),intent(in)::mss_frc_snd(plond) ! I [frc] Mass fraction of sand
c     Output
c     Input/Output
c     Fluxes refer to dust on input and to CaCO3 on output
      real(r8),intent(inout)::flx_mss_vrt_dst_CaCO3(plond,dst_nbr) ! I/O [kg m-2 s-1] Vertical mass flux of dust on input, CaCO3 on output
      integer rcd               ! [rcd] Return success code
c     Local
      integer m                 ! [idx] Counting index
      integer lon_idx           ! [idx] Counting index for lon
      real(r8) mss_frc_slt(plond)   ! [frc] Mass fraction of silt
      real(r8) mss_frc_CaCO3_sz_crr ! [frc] Fraction of soil CaCO3 in size bin
      real(r8) mss_frc_CaCO3_cly    ! [frc] Fraction of CaCO3 in clay
      real(r8) mss_frc_CaCO3_slt    ! [frc] Fraction of CaCO3 in silt
      real(r8) mss_frc_CaCO3_snd    ! [frc] Fraction of CaCO3 in sand
c     Main Code
c     Initialize output
      do m=1,dst_nbr
         do lon_idx=1,plon
c     Simple technique is to mask dust mass by tracer mass fraction
c     Model transports (hence conserves) CaCO3 rather than total dust itself
c     Method assumes source, transport, and removal processes are linear with tracer mass
            if (flg_mbl(lon_idx)) then
c     20000320: Currently this is only process  in dust model requiring mss_frc_slt
               mss_frc_slt(lon_idx)= ! [frc] Mass fraction of silt
     $              max(0.0_r8,1.0-mss_frc_cly(lon_idx)-mss_frc_snd(lon_idx))
c     CSB99 showed that CaCO3 is not uniformly distributed across sizes
c     There is more CaCO3 per unit mass of silt than per unit mass of clay
               mss_frc_CaCO3_cly=max(0.0_r8,-0.045+0.5*min(0.5_r8,mss_frc_cly(lon_idx))) ! [frc] Fraction of CaCO3 in clay CSB99 p. 22249 Figure 1b
               mss_frc_CaCO3_slt=max(0.0_r8,-0.175+1.4*min(0.5_r8,mss_frc_slt(lon_idx))) ! [frc] Fraction of CaCO3 in silt CSB99 p. 22249 Figure 1a
               mss_frc_CaCO3_snd=1.0-mss_frc_CaCO3_cly-mss_frc_CaCO3_snd ! [frc] Fraction of CaCO3 in sand CSB99 p. 22249 Figure 1a
c     Set CaCO3 fraction of total CaCO3 for each transport bin 
               if (dmt_vwr(m) < dmt_cly_max) then
c     Transport bin carries Clay
                  mss_frc_CaCO3_sz_crr=mss_frc_CaCO3_cly ! [frc] Fraction of soil CaCO3 in size bin
               else if (dmt_vwr(m) < dmt_slt_max) then ! endif size bin carries clay
c     Transport bin carries Silt
                  mss_frc_CaCO3_sz_crr=mss_frc_CaCO3_slt ! [frc] Fraction of soil CaCO3 in size bin
               else             ! endif size bin carries silt
c     Transport bin carries Sand
                  mss_frc_CaCO3_sz_crr=mss_frc_CaCO3_snd ! [frc] Fraction of soil CaCO3 in size bin
               endif            ! endif size bin carries sand
c     Sanity check
               if (mss_frc_CaCO3_sz_crr < 0.0.or.mss_frc_CaCO3_sz_crr > 1.0) stop
     $              'dst: flx_mss_CaCO3_msk() mss_frc_CaCO3_sz_crr < 0.0.or.mss_frc_CaCO3_sz_crr > 1.0'
               if (mss_frc_CaCO3(lon_idx) < 0.0.or.mss_frc_CaCO3(lon_idx) > 1.0) stop
     $              'dst: flx_mss_CaCO3_msk() mss_frc_CaCO3 < 0.0.or.mss_frc_CaCO3 > 1.0'

c     Convert dust flux to CaCO3 flux
               flx_mss_vrt_dst_CaCO3(lon_idx,m)=flx_mss_vrt_dst_CaCO3(lon_idx,m) ! [kg m-2 s-1]
     $              *mss_frc_CaCO3(lon_idx) ! [frc] Mass fraction of CaCO3 (at this location)
     $              *mss_frc_CaCO3_sz_crr ! [frc] Fraction of soil CaCO3 in size bin
            endif               ! endif 
         enddo                  ! end loop over lon
      enddo                     ! end loop over sz

      return
      end subroutine  flx_mss_CaCO3_msk                       ! end flx_mss_CaCO3_msk()
      
      subroutine flx_mss_hrz_slt_ttl_Whi79_get(
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     flx_mss_hrz_slt_ttl, ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
     $     wnd_frc,             ! I [m s-1] Surface friction velocity
     $     wnd_frc_thr_slt)     ! I [m s-1] Threshold friction speed for saltation
c     Purpose: Compute vertically integrated streamwise mass flux of particles
c     Theory: Uses method proposed by White (1979)
c     Dependencies: <pmgrid.h> <dstcst.h>
c     fxm: use surface air density not midlayer density
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstcst.h> /* Physical constants for dust routines */
c     Parameters
      real(r8),parameter::cst_slt=2.61 ! [frc] Saltation constant Whi79 p. 4648, MaB97 p. 16422 
c     Input
      logical,intent(in)::flg_mbl(plond)    ! I [flg] Mobilization candidate flag
      real(r8),intent(in)::dns_mdp(plond)       ! I [kg m-3] Midlayer density
      real(r8),intent(in)::wnd_frc(plond)       ! I [m s-1] Surface friction velocity
      real(r8),intent(in)::wnd_frc_thr_slt(plond) ! I [m s-1] Threshold friction speed for saltation
c     Output
      real(r8),intent(out)::flx_mss_hrz_slt_ttl(plond) ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
      integer rcd               ! [rcd] Return success code
c     Local
      real(r8) wnd_frc_rat          ! [frc] Ratio of wind friction threshold to wind friction
      integer lon_idx           ! [idx] Counting index for lon
c     Main Code
c     Initialize output
      flx_mss_hrz_slt_ttl=0.0   ! [kg m-1 s-1]
      
      do lon_idx=1,plon
         if (flg_mbl(lon_idx).and.wnd_frc(lon_idx) > wnd_frc_thr_slt(lon_idx)) then
            wnd_frc_rat=wnd_frc_thr_slt(lon_idx)/wnd_frc(lon_idx) ! [frc]
            flx_mss_hrz_slt_ttl(lon_idx)= ! [kg m-1 s-1] 
     $           cst_slt*dns_mdp(lon_idx)*(wnd_frc(lon_idx)**3.0)*
     $           (1.0-wnd_frc_rat)*(1.0+wnd_frc_rat)*(1.0+wnd_frc_rat)/grv_sfc ! Whi79 p. 4648 (19), MaB97 p. 16422 (28)
         endif                  ! endif 
      enddo                     ! end loop over lon
      
      return
      end subroutine  flx_mss_hrz_slt_ttl_Whi79_get                      ! end flx_mss_hrz_slt_ttl_Whi79_get()
      
      subroutine flx_mss_vrt_dst_ttl_MaB95_get(
     $     dst_slt_flx_rat_ttl, ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     flx_mss_hrz_slt_ttl, ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
     $     flx_mss_vrt_dst_ttl, ! O [kg m-2 s-1] Total vertical mass flux of dust
     $     mss_frc_cly)         ! I [frc] Mass fraction clay 
c     Purpose: Diagnose total vertical mass flux of dust from vertically integrated streamwise mass flux
c     Theory: Uses clay-based method proposed by Marticorena & Bergametti (1995)
c     Their parameterization is based only on data for mss_frc_cly < 0.20
c     For clayier soils, dst_slt_flx_rat_ttl may behave dramatically differently
c     Whether this behavior changes when mss_frc_cly > 0.20 is unknown
c     Anecdotal evidence suggests vertical flux decreases for mss_frc_cly > 0.20
c     Thus we use min[mss_frc_cly,0.20] in MaB95 parameterization
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
c     Parameters
c     Input
      logical,intent(in)::flg_mbl(plond)    ! I [flg] Mobilization candidate flag
      real(r8),intent(in)::flx_mss_hrz_slt_ttl(plond) ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
      real(r8),intent(in)::mss_frc_cly(plond)   ! I [frc] Mass fraction clay 
c     Output
      real(r8),intent(out)::dst_slt_flx_rat_ttl(plond) ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
      real(r8),intent(out)::flx_mss_vrt_dst_ttl(plond) ! O [kg m-2 s-1] Total vertical mass flux of dust
      integer rcd               ! [rcd] Return success code
c     Local
      integer lon_idx           ! [idx] Counting index for lon
      real(r8) mss_frc_cly_vld      ! [frc] Mass fraction clay limited to 0.20
      real(r8) ln10                 ! [frc] Natural log of 10
      
c     Initialize some variables
      ln10=log(10.0)            ! [frc] Natural log of 10
      do lon_idx=1,plon
         if (flg_mbl(lon_idx)) then
c     19990603: fxm: Dust production is EXTREMELY sensitive to this parameter, which changes flux by 3 orders of magnitude in 0.0 < mss_frc_cly < 0.20
            mss_frc_cly_vld=min(mss_frc_cly(lon_idx),0.2_r8) ! [frc]
            dst_slt_flx_rat_ttl(lon_idx)= ! [m-1]
     $           100.0*exp(ln10*(13.4*mss_frc_cly_vld-6.0)) ! MaB95 p. 16423 (47)
            flx_mss_vrt_dst_ttl(lon_idx)=flx_mss_hrz_slt_ttl(lon_idx)*dst_slt_flx_rat_ttl(lon_idx) ! [kg m-1 s-1] 
         endif                  ! endif flg_mbl
      enddo                     ! end loop over lon
      
      return
      end subroutine  flx_mss_vrt_dst_ttl_MaB95_get                      ! end flx_mss_vrt_dst_ttl_MaB95_get()
      
      subroutine flx_mss_vrt_dst_prt(
     $     flg_mbl,             ! I [flg] Mobilization candidate flag
     $     flx_mss_vrt_dst,     ! O [kg m-2 s-1] Vertical mass flux of dust
     $     flx_mss_vrt_dst_ttl) ! I [kg m-2 s-1] Total vertical mass flux of dust
c     Purpose: Partition total vertical mass flux of dust into transport bins
c     Theory: 
c     Dependencies: <pmgrid.h>,<dstaer.h>
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstaer.h> /* Aerosol microphysical properties */
c     Input
      logical,intent(in)::flg_mbl(plond)    ! I [flg] Mobilization candidate flag
      real(r8),intent(in)::flx_mss_vrt_dst_ttl(plond) ! I [kg m-2 s-1] Total vertical mass flux of dust
c     Output
      real(r8),intent(out)::flx_mss_vrt_dst(plond,dst_nbr) ! O [kg m-2 s-1] Vertical mass flux of dust
      integer rcd               ! [rcd] Return success code
c     Local
      integer lon_idx           ! [idx] Counting index for lon
      integer src_idx           ! [idx] Counting index for src
      integer snk_idx           ! [idx] Counting index for snk
      integer snk_nbr           ! [nbr] Dimension size
c++mahow
      real(r8) src_fraction(4)
      src_fraction(1)=.1
      src_fraction(2)=.3
      src_fraction(3)=.3
      src_fraction(4)=.3
c--mahow

c     Initialize some variables
      flx_mss_vrt_dst=0.0       ! [frc]
c     Recall that snk_nbr=dst_nbr
      do lon_idx=1,plon         ! NB: Inefficient loop order
         if (flg_mbl(lon_idx)) then
            do snk_idx=1,dst_nbr
!pjr               do src_idx=1,dst_src_nbr
                  flx_mss_vrt_dst(lon_idx,snk_idx)= ! [kg m-2 s-1] 
     $                 flx_mss_vrt_dst(lon_idx,snk_idx)+
c-mahow     $                 ovr_src_snk_mss(src_idx,snk_idx)*
c++mahow
     $                 src_fraction(snk_idx)*
c--mahow
     $                 flx_mss_vrt_dst_ttl(lon_idx) 
!pjr               enddo            ! end loop over src
            enddo               ! end loop over snk
         endif                  ! endif flg_mbl
      enddo                     ! end loop over lon
      
      return
      end subroutine flx_mss_vrt_dst_prt                      ! end flx_mss_vrt_dst_prt()
      
      subroutine lai_lsm_get(
     $     doy,                 ! I [day] Day of year [1.0..366.0)
     $     lat_dgr,             ! I [dgr] Latitude
     $     sfc_typ,             ! I [idx] LSM surface type (0..28)
     $     lai_lsm)             ! O [m2 m-2] Leaf area index, one-sided
c     Purpose: Update time-varying surface fields to current timestep
c     Routine currently computes latitude slice of LAI based on surface type, date, and latitude
c     Algorithm is from CCM:lsm/phenol()
c     lai_lsm_get() is not currently called, but would be called by CCM:dynamics/advnce(), MATCH:src/main()
      implicit none
#include <pmgrid.h>
#include <dstlsm.h> /* LSM block data */
c     Parameters
      integer,parameter::mpy=12 ! Resolution of LAI data
      integer,parameter::dpy=365 ! Model days per year
c     Input
      integer,intent(in)::sfc_typ(plond)    ! [idx] LSM surface type (0..28)
      real(r8),intent(in)::lat_dgr              ! [dgr] Latitude
      real(r8),intent(in)::doy                  ! [day] Day of year [1.0..366.0)
c     Output
      real(r8),intent(out)::lai_lsm(plond)        ! [m2 m-2] Leaf area index, one-sided
c     Input/Output
c     Local
      integer iday_NH           ! [day] Day number Northern Hemisphere [1..365]
      integer iday_SH           ! [day] iday_NH shifted 6 mth for SH [1..365]
      integer idoy              ! [day] Current day of year [1..365]
      integer idx_mth_glb       ! Interpolation month, future
      integer idx_mth_lub       ! Interpolation month, past
      integer imoy              ! [mth] Current month of year [1..12]
      integer lon_idx           ! [idx] Counting index
      integer pln_typ_idx       ! [idx] Plant type index
      integer sfc_typ_idx       ! [idx] Surface type index
      integer sgs_idx           ! Surface sub-gridscale index
      real(r8) day_grow             ! [day] Days since Jan 1 in NH or Jul 1 in SH (1.0..365.0)
      real(r8) dom                  ! [day] Current day of month [1.0..32.0)
      real(r8) mth_grow             ! [mth] Months since Jan (NH) or Jul (SH) (1.0..12.0)
      real(r8) wgt_glb              ! [frc] Interpolation weight
      real(r8) wgt_lub              ! [frc] Interpolation weight
      
      integer dpm(mpy)          ! [day] Days per month
      integer fdom_doy(mpy)     ! [day] Day of year of first day of month 
      integer ldom_doy(mpy)     ! [day] Day of year of last day of month
c     NB: Initializing named data statements seems to be non-standard in F90
      data dpm     / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data fdom_doy/  1, 32, 60, 91,121,152,182,213,244,274,305,335/
      data ldom_doy/ 31, 59, 90,120,151,181,212,243,273,304,334,365/
      save dpm
      save fdom_doy
      save ldom_doy
      
c     Main Code
c     Sanity check
      if (doy < 1.0.or.doy >= 366.0) stop 'dst: ERROR doy < 1.0.or.doy >= 366.0 in lai_lsm_get()'                  
      
c     Initialize array
      lai_lsm=0.0               ! [m2 m-2]
      
c     Convert doy to moy and dom
      imoy=1
      idoy=int(doy)
      do while(ldom_doy(imoy) < idoy.and.imoy < mpy)
         imoy=imoy+1
      end do                    ! end while
      dom=doy-fdom_doy(imoy)+1.0
c     Sanity check
      if (imoy < 1.or.imoy > mpy) stop 'dst: ERROR imoy < 1.or.imoy > mpy in lai_lsm_get()'
      if (dom < 1.0.or.dom > 32.0) stop 'dst: ERROR dom < 1.0.or.dom > 32.0 in lai_lsm_get()'
      
c     Compute iday_NH and iday_SH:
c     iday_NH = day of current year since Jan 0 [1..365]
c     iday_SH = iday_NH shifted 6 mth for SH: 1->183, 183->365, 184->1, 365->182
      iday_NH=int(doy)          ! Day of year [1..365]
      iday_SH=mod(iday_NH-1+dpy/2,dpy)+1
      
c     Compute cyclic indices into monthly data idx_mth_glb,idx_mth_lub
      if (lat_dgr >= 0.0) then
         day_grow=iday_NH 
      else
         day_grow=iday_SH 
      end if                    ! endif
      mth_grow=12.0*(day_grow-0.5)/dpy
      idx_mth_glb=mth_grow+0.5
      idx_mth_lub=idx_mth_glb+1
      wgt_glb=(idx_mth_glb+0.5)-mth_grow
      wgt_lub=1.0-wgt_glb
      if (idx_mth_glb < 1) idx_mth_glb=mpy
      if (idx_mth_lub > mpy) idx_mth_lub=1
      
      do lon_idx=1,plon
c     Store surface blend of current gridpoint
         sfc_typ_idx=sfc_typ(lon_idx)
         do sgs_idx=1,3
c     Non-vegetated surfaces have pln_typ=14
            pln_typ_idx=pln_typ(sfc_typ_idx,sgs_idx)
c     Current LAI is weighted past LAI...
            lai_lsm(lon_idx)=   ! [m2 m-2]
     $           lai_lsm(lon_idx)+ ! [m2 m-2]
     $           pln_frc(sfc_typ_idx,sgs_idx)* ! [frc]
     $           gai(pln_typ_idx,idx_mth_glb)* ! [m2 m-2]
     $           wgt_glb        ! [frc] 
c     ...plus weighted future LAI
            lai_lsm(lon_idx)=   ! [m2 m-2]
     $           lai_lsm(lon_idx)+ ! [m2 m-2]
     $           pln_frc(sfc_typ_idx,sgs_idx)* ! [frc] 
     $           gai(pln_typ_idx,idx_mth_lub)* ! [m2 m-2] 
     $           wgt_lub        ! [frc]
         enddo                  ! end loop over number of possible plant types
      end do                    ! end loop over lon
      return 
      end subroutine lai_lsm_get          ! end lai_lsm_get()
      
      subroutine dst_lsm_ini(
     $     pln_frc,             ! I
     $     pln_typ,             ! I
     $     lnd_frc_dry,         ! I
     $     sfc_typ,             ! I
     $     sfc_dst_mbl_frc)     ! O
c     Purpose: LSM vegetation blends and surface type fields determine time-invariant fraction of each gridcell available for deflation
c     dst_lsm_ini() is not currently called, but used to be called by CCM:control/initext(), MATCH:src/inirun()
      implicit none
#include <pmgrid.h>
c     Input
      integer,intent(in)::pln_typ(0:28,3)   ! [idx] LSM plant type (1..14 = nbr_LSM_pln_typ)
      integer,intent(in)::sfc_typ(plon,plat) ! [idx] LSM surface type 
      
      real(r8),intent(in)::pln_frc(0:28,3)      ! [frc] LSM fractional coverage of corresponding plant type (sums to 1)
      real(r8),intent(in)::lnd_frc_dry(plon,plat)   ! [frc] Land surface that is not lake or wetland (by area)
c     Output
      real(r8),intent(out)::sfc_dst_mbl_frc(plon,plat) ! [frc] Non-lake, non-wetland area
c     Local
      integer sgs_idx           ! Surface sub-gridscale index
      integer lat_idx           ! Latitude index
      integer lon_idx           ! Longitude index
      integer sfc_typ_idx       ! Surface type index
      
c     Set the fraction of the surface available for mobilization by dust
      do lat_idx=1,plat
         do lon_idx=1,plon
c     LSM surface type described in Bon96 Table 5, p. 12
c     Exclude mobilization from beneath tall canopies, lakes, wetlands, etc
            
c     Get rid of ocean and sea ice points
c     if (nint(oro(i)) /= 1) flx_mss_mbl_sfc(i,m)=0.0
            
c     Store surface blend of current gridpoint
            sfc_typ_idx=sfc_typ(lon_idx,lat_idx)
c     NB: plant type (pln_typ) and plant fraction (pln_frc) arrays MUST be dimensioned
c     (0:28,3) in order to be correctly indexed by sfc_typ.
c     Dimensioning them (29,3) will fail unless you index them by sfc_typ+1
c     This is because sfc_typ=0 (Ocean) is allowed by LSM
            
c     Initialize the array
            sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0 
            
c     Include shrubs, grasses, crops, and bare ground (TeF94 p. 22900)
            do sgs_idx=1,3
               if (pln_typ(sfc_typ_idx,sgs_idx) >= 6.and.
     $              pln_typ(sfc_typ_idx,sgs_idx) <= 14) then
                  sfc_dst_mbl_frc(lon_idx,lat_idx)= ! [frc]
     $                 sfc_dst_mbl_frc(lon_idx,lat_idx)+
     $                 pln_frc(sfc_typ_idx,sgs_idx)
               endif            ! endif bare ground
            enddo               ! end loop over number of possible plant types
            
c     Exclude ocean
            if (sfc_typ_idx == 0) sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0 ! [frc]
            
c     Exclude ice and sea ice
            if (sfc_typ_idx == 1) sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0 ! [frc]
            
c     Exclude warm crops (mostly Europe)
            if (sfc_typ_idx == 26) sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0 ! [frc]

c     Exclude wetlands
            if (sfc_typ_idx == 27) sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0 ! [frc]
            if (sfc_typ_idx == 28) sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0 ! [frc]
            
c     Land fraction suitable for dust mobilization is known for each surface gridbox
c     However, many gridboxes are part land and part water
c     Renormalize sfc_dst_mbl_frc to account for fractional lakes, wetlands, and ocean
c     This removes as candidate source regions any areas which are wet year-round
c     GBB uses same procedure to normalize land fluxes to actual grid-cell area in CCM:lsm/lsmini.F
            sfc_dst_mbl_frc(lon_idx,lat_idx)= ! [frc]
     $           sfc_dst_mbl_frc(lon_idx,lat_idx)* ! [frc]
     $           lnd_frc_dry(lon_idx,lat_idx) ! [frc]
c     sfc_dst_mbl_frc now contains the (time-invariant) surface fraction of each gridcell considered suitable for dust mobilization. 
            
c     Sanity check
            if ((sfc_dst_mbl_frc(lon_idx,lat_idx) > 1.0).or.
     $           (sfc_dst_mbl_frc(lon_idx,lat_idx) < 0.0)) then
               write(6,'(a,i3,a,i2,a,f6.4)') 
     $              'dst: dst_lsm_ini(): sfc_dst_mbl_frc(',
     $              lon_idx,',',lat_idx,') = ',
     $              sfc_dst_mbl_frc(lon_idx,lat_idx)
            endif               ! endif out of bounds error
            
         end do                 ! end loop over lon
      end do                    ! end loop over lat
      
      return
      end subroutine dst_lsm_ini                  ! end dst_lsm_ini()

      subroutine lnd_frc_mbl_get(
     $     doy,                 ! I [day] Day of year [1.0..366.0)
     $     flg_mbl,             ! O [flg] Mobilization candidate flag
     $     lat_rdn,             ! I [rdn] Latitude
     $     lnd_frc_dry,         ! I [frc] Dry land fraction
     $     lnd_frc_mbl,         ! O [frc] Bare ground fraction
     $     mbl_nbr,             ! O [flg] Number of mobilization candidates
     $     oro,                 ! I [frc] Orography
     $     sfc_typ,             ! I [idx] LSM surface type (0..28)
c-mahow     $     snw_frc)             ! I [frc] Fraction of surface covered by snow
c++mahow
     $     snw_frc             ! I [frc] Fraction of surface covered by snow
     $ ,   lat)                ! I latitidue number
c--mahow

c     Purpose: Return fraction of each gridcell suitable for dust mobilization
c     lnd_frc_mbl_get() is called by dst_mbl()
c     The date is used to obtain the time-varying vegetation cover
c     Routine currently computes latitude slice of LAI based on surface type, date, and latitude
c     LAI algorithm is from CCM:lsm/phenol() Bon96
c     The LSM data are mid-month values, i.e., valid on the 15th of the month
      implicit none
#include <pmgrid.h>
#include <dstlsm.h> /* LSM block data */
#include <dbg.com> /* Debugging constants */
c     Parameters
      integer,parameter::mpy=12 ! [nbr] Resolution of LAI data
      integer,parameter::dpy=365 ! [day] Model days per year
      real(r8),parameter::vai_mbl_thr=0.2 ! [m2 m-2] VAI threshold quenching dust mobilization
c     Input
      integer,intent(in)::sfc_typ(plond)    ! I [idx] LSM surface type (0..28)
      real(r8),intent(in)::doy                  ! I [day] Day of year [1.0..366.0)
      real(r8),intent(in)::lat_rdn              ! I [rdn] Latitude
      real(r8),intent(in)::lnd_frc_dry(plond)   ! I [frc] Dry land fraction
      real(r8),intent(in)::oro(plond)           ! I [frc] Orography
      real(r8),intent(in)::snw_frc(plond)       ! I [frc] Fraction of surface covered by snow
c     Output
      integer,intent(out)::mbl_nbr ! O [nbr] Number of mobilization candidates
      logical,intent(out)::flg_mbl(plond) ! O [flg] Mobilization candidates
      real(r8),intent(out)::lnd_frc_mbl(plond) ! O [frc] Bare ground fraction
c     Input/Output
c     Local
      integer idx_idx           ! [idx] Counting index
      integer lnd_idx(plond)    ! [idx] Longitude index array (land)
      integer lnd_nbr           ! [nbr] Number of land points
      integer iday_NH           ! [day] Day number Northern Hemisphere [1..365]
      integer iday_SH           ! [day] iday_NH shifted 6 mth for SH [1..365]
      integer idoy              ! [day] Current day of year [1..365]
      integer idx_mth_glb       ! [idx] Interpolation month, future
      integer idx_mth_lub       ! [idx] Interpolation month, past
      integer imoy              ! [mth] Current month of year [1..12]
      integer lon_idx           ! [idx] Counting index for longitude
      integer pln_typ_idx       ! [idx] Plant type index
      integer sfc_typ_idx       ! [idx] Surface type index
      integer sgs_idx           ! [idx] Surface sub-gridscale index
      real(r8) day_grow             ! [day] Days since Jan 1 in NH or Jul 1 in SH (1.0..365.0)
      real(r8) dom                  ! [day] Current day of month [1.0..32.0)
      real(r8) lat_dgr              ! [dgr] Latitude
      real(r8) mth_grow             ! [mth] Months since Jan (NH) or Jul (SH) (1.0..12.0)
      double precision pi       ! [frc] 3
      real(r8) vai_sgs              ! [m2 m-2] Leaf + stem area index, one-sided
      real(r8) wgt_glb              ! [frc] Interpolation weight
      real(r8) wgt_lub              ! [frc] Interpolation weight
      real(r8) pln_frc_mbl          ! [frc] "Bare ground" fraction of sub-gridscale cell
      real(r8) pln_frc_sgs          ! [frc] Plant fraction of current sub-gridscale cell
      
      integer dpm(mpy)          ! [day] Days per month
      integer fdom_doy(mpy)     ! [day] Day of year of first day of month 
      integer ldom_doy(mpy)     ! [day] Day of year of last day of month
c     NB: Initializing named data statements seems to be non-standard in F90
      data dpm     / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data fdom_doy/  1, 32, 60, 91,121,152,182,213,244,274,305,335/
      data ldom_doy/ 31, 59, 90,120,151,181,212,243,273,304,334,365/
c++mahow
      real(r8) source(plon,plat)
      common /source/ source
      integer lat,i,j
c--mahow

      save dpm
      save fdom_doy
      save ldom_doy
c     Statement functions
      logical oro_is_lnd        ! True if > 50% land
      real(r8) oro_val              ! [frc] Orography
      oro_is_lnd(oro_val)=nint(oro_val) == 1
c     Main Code
c     Sanity check
      if (doy < 1.0.or.doy >= 366.0) stop 'dst: ERROR doy < 1.0.or.doy >= 366.0 in lnd_frc_mbl_get()'                  
      if (vai_mbl_thr == 0.0) stop 'dst: ERROR vai_mbl_thr == 0.0 in lnd_frc_mbl_get()'                  
c     Initialize outputs
      pi=4.0*atan(1.0d0)        ! [frc] 3
      lat_dgr=180.0*lat_rdn/pi  ! [dgr] Latitude
      mbl_nbr=0                 ! [nbr] Number of mobilization candidates
      do lon_idx=1,plon         ! NB: plond
         flg_mbl(lon_idx)=.false. ! [flg] Mobilization candidates
      enddo                     ! end loop over lon
      lnd_frc_mbl=0.0           ! [frc]

c     Land ahoy!
      lnd_nbr=0
      do lon_idx=1,plon
         if (oro_is_lnd(oro(lon_idx))) then
            lnd_nbr=lnd_nbr+1
            lnd_idx(lnd_nbr)=lon_idx
         end if                 ! endif
      enddo                     ! end loop over lon
      
c     Much ado about nothing
      if (lnd_nbr == 0) return

c     Convert doy to moy and dom
      imoy=1
      idoy=int(doy)
      do while(ldom_doy(imoy) < idoy.and.imoy < mpy)
         imoy=imoy+1
      end do                    ! end while
      dom=doy-fdom_doy(imoy)+1.0
c     Sanity check
      if (imoy < 1.or.imoy > mpy) stop 'dst: ERROR imoy < 1.or.imoy > mpy in lnd_frc_mbl_get()'
      if (dom < 1.0.or.dom > 32.0) stop 'dst: ERROR dom < 1.0.or.dom > 32.0 in lnd_frc_mbl_get()'
      
c     Compute iday_NH and iday_SH:
c     iday_NH = day of current year since Jan 0 [1..365]
c     iday_SH = iday_NH shifted 6 mth for SH: 1->183, 183->365, 184->1, 365->182
      iday_NH=int(doy)          ! Day of year [1..365]
      iday_SH=mod(iday_NH-1+dpy/2,dpy)+1
      
c     Compute cyclic indices into monthly data idx_mth_glb,idx_mth_lub
      if (lat_dgr >= 0.0) then
         day_grow=iday_NH 
      else
         day_grow=iday_SH 
      end if                    ! endif
      mth_grow=12.0*(day_grow-0.5)/dpy
      idx_mth_glb=mth_grow+0.5
      idx_mth_lub=idx_mth_glb+1
      wgt_glb=(idx_mth_glb+0.5)-mth_grow
      wgt_lub=1.0-wgt_glb
      if (idx_mth_glb < 1) idx_mth_glb=mpy
      if (idx_mth_lub > mpy) idx_mth_lub=1
      
c     Land points
      do idx_idx=1,lnd_nbr
         lon_idx=lnd_idx(idx_idx)
c++mahow
         lnd_frc_mbl(lon_idx)=source(lon_idx,lat)
c--mahow
#ifdef OFFCHARLIE
c     Store surface blend of current gridpoint
         sfc_typ_idx=sfc_typ(lon_idx)
         if (sfc_typ_idx <= 1.or. ! Inland lakes and land ice
     $        sfc_typ_idx >= 27) then ! Wetlands
c??mahow    lnd_frc_mbl(lon_idx)=0.0 ! [frc] Bare ground fraction
            lnd_frc_mbl(lon_idx)=0.0 ! [frc] Bare ground fraction
         else                   ! Normal land
            do sgs_idx=1,3
c     Non-vegetated surfaces have pln_typ=14
               pln_typ_idx=pln_typ(sfc_typ_idx,sgs_idx) ! [idx]
               pln_frc_sgs=pln_frc(sfc_typ_idx,sgs_idx) ! [frc]
c     Current LAI is weighted past LAI plus weighted future LAI
               vai_sgs= ! [m2 m-2] Leaf + stem area index, one-sided
     $              +tai(pln_typ_idx,idx_mth_glb)*wgt_glb
     $              +tai(pln_typ_idx,idx_mth_lub)*wgt_lub
c     The "bare ground" fraction of the current sub-gridscale cell decreases
c     linearly from 1.0 to 0.0 as VAI increases from 0.0 to vai_mbl_thr
               pln_frc_mbl=     ! [frc] "Bare ground" fraction of sub-gridscale cell
     $              1.0-min(vai_sgs,vai_mbl_thr)/vai_mbl_thr
c               if (dbg_lvl == dbg_crr) then
c                  write(6,'(a,3(a,f6.4,a))') 'lnd_frc_mbl_get(): ',
c     $                 'vai_sgs = ',vai_sgs,' m-2 m-2, ',
c     $                 'pln_frc = ',pln_frc_sgs,', ',
c     $                 'pln_frc_mbl = ',pln_frc_mbl,''
c               endif            ! endif dbg
c??mahow       lnd_frc_mbl(lon_idx)=lnd_frc_mbl(lon_idx)+pln_frc_sgs*pln_frc_mbl ! [frc]
               lnd_frc_mbl(lon_idx)=lnd_frc_mbl(lon_idx)+pln_frc_sgs*pln_frc_mbl ! [frc]
            enddo               ! end loop over number of possible plant types
         endif                  ! endif normal land
#endif

c     Adjust for factors which constrain entire gridcell
         lnd_frc_mbl(lon_idx)=  
     $        lnd_frc_mbl(lon_idx) ! [frc] Bare ground fraction
     $        *lnd_frc_dry(lon_idx) ! Account for wetlands, lake, ocean, ice
     $        *(1.0-snw_frc(lon_idx)) ! Account for snow coverage

         if ((lnd_frc_mbl(lon_idx) > 1.0+small).or.(lnd_frc_mbl(lon_idx) < -small)) then
            write(6,'(a,i3,a,1p,e21.12)') 
     $           'dst: lnd_frc_mbl_get(): lnd_frc_mbl(',
     $           lon_idx,') = ',lnd_frc_mbl(lon_idx)
!pjr            stop
         endif                  ! endif out of bounds error

         lnd_frc_mbl(lon_idx) = max(0._r8, min(one, lnd_frc_mbl(lon_idx)))

         if (lnd_frc_mbl(lon_idx) > 0.0) then
            flg_mbl(lon_idx)=.true. ! [flg] Mobilization candidates
            mbl_nbr=mbl_nbr+1   ! [nbr] Number of mobilization candidates
         endif                  ! endif mobilization candidate
            
      end do                    ! end loop over land
      return 
      end subroutine lnd_frc_mbl_get                      ! end lnd_frc_mbl_get()

      end module dstmblutl
