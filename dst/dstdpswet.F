      module dstdpswet
      use precision

      real(r8) :: one = 1.0
      real(r8) :: zero = 0.0


      contains
c     $Header$

c     Purpose: dstdpswet.F contains routines which represent mineral dust 
c     wet deposition tendencies
      
c     Implementation is designed with two goals in mind:
c     1. Make scheme work with CCM3 clouds, since Paleo-CSM will use these
c     2. Make scheme work with PJR's clouds, since Neo-CSM will use these
c     PJR's prognostic scheme contains a more complete set of wet deposition predictors than the CCM3 cloud scheme 
c     In addition to q_H2O_cnd2pcp_tnd, which is available in CCM3, PJR predicts q_H2O_wtr2pcp_ttl_tnd and q_H2O_pcp2vpr_tnd
c     These are necessary for "exact" representation of wet deposition
c     Assuming a simple scavenging ratio is sufficient to parameterize the 
c     wet deposition process, the q_H2O_wtr2pcp_ttl_tnd, q_H2O_pcp2vpr_tnd, and q_H2O_cnd2pcp_tnd tendencies suffice to predict dust rain-out rates.
c     However, vanilla CCM3 only provides dqcond, the net vapor tendency due to condensation and evaporation
c     CCM3 does not provide condensation and evaporation tendencies separately
c     Thus we map available CCM3 tendencies onto tendencies of same name that PJR provides
c     Thus q_H2O_wtr2pcp_ttl_tnd and q_H2O_pcp2vpr_tnd in CCM3 do not mean quite the same thing as they do in CCM4
c     However, this strategy has the advantage that once inputs are defined in either
c     version of the model (CCM3 or CCM4), dust prediction can use the same code.
      
c     Cloud nomenclature:
c     cld_frc = [frc] Local total cloud fraction
c     = Horizontal fraction of gridcell covered by cloud
c     = MATCH:cldt
c     cld_frc_cnv = [frc] Convective cloud fraction
c     = Horizontal fraction of gridcell covered by convective cloud
c     = CCM/MATCH:cldc
c     cld_vlm = [frc] Cloud volume
c     = Fraction of gridcell volume occupied by rain or cloud water
c     = Is larger than cld_frc when rain from above passes through gridcell
c     = MATCH:cldv

c     Tracer nomenclature:
c     q_H2O_vpr = [kg kg-1] Water vapor mixing ratio
c     = CCM/MATCH:q

c     q_H2O_cnd = [kg kg-1] Condensed H2O mixing ratio
c     = Cloud water amount
c     = Sum of liquid and ice water
c     = MATCH:cwat

c     q_H2O_pcp_lqd = ! [kg kg-1] Rainwater mixing ratio
c     = Mixing ratio of liquid precipitation (for wet deposition purposes)
c     = MATCH:rain

c     q_H2O_vpr2cnd_tnd = [kg kg-1 s-1] H2O vapor to condensate tendency
c     q_H2O_cnd_tnd = [kg kg-1 s-1] Net H2O condensate formation tendency
c     = Condensation - evaporation 
c     = This is net condensation minus evaporation, so may be positive or negative 
c     = Positive when condensation exceeds evaporation
c     = Negative when evaporation exceeds condensation
c     = Evaporation here refers to in-cloud evaporation, and does not include
c     = evaporation from falling rain which is counted separately in q_H2O_pcp2vpr
c     = MATCH:cme but MATCH sometimes uses this as an in-cloud rate
c     = MATCH:conds

c     q_H2O_pcp2vpr_tnd = [kg kg-1 s-1] H2O precipitation to vapor tendency
c     = H2O stratiform precipitation to vapor tendency
c     = Evaporation of non-local, stratiform rain into local gridcell
c     = Rate of evaporation of falling stratiform precipitation
c     = Does not include any evaporation from in-cloud water
c     = Evaporation from convective rain (convective downdrafts) is not included (it may be accounted for directly in convection routines, I do not know, but it is not included in this term)
c     = Positive definite
c     = MATCH:evapr,evaps

c     q_H2O_cnd2pcp_tnd = [kg kg-1 s-1] Condensed H2O to precipitation tendency
c     = Local water to stratiform precipitation tendency
c     = Stratiform precipitation tendency
c     = CCM 3.6: vpr2pcp as predicted by physics/cond()
c     = CCM 3.7+: stratiform precipitation is taken directly from condensate
c     = Production tendency only, does not include local evaporation of falling rain
c     = Positive definite
c     = MATCH:prain,precs

c     q_H2O_vpr2pcp_cnv_tnd = [kg kg-1 s-1] H2O vapor to convective precipitation tendency
c     = Precipitation formation tendency from convection
c     = Production tendency only, does not include local evaporation of falling rain
c     = Positive definite
c     = CCM/MATCH:cmfdqr

c     q_H2O_wtr2pcp_ttl_tnd = [kg kg-1 s-1] H2O vapor + condensate conversion to precipitation tendency
c     = Local water (vapor + condensate) to total (convective + stratiform) precipitation tendency
c     = Sum of convective and stratiform precipitation production tendencies
c     = Positive definite
c     = q_H2O_cnd2pcp_tnd + q_H2O_vpr2pcp_cnv_tnd

c     Precipitation flux nomenclature:
c     pcp_flx_sfc = [kg m-2 s-1] Total precipitation reaching surface
c     = CCM:prect (but in mass flux units not m s-1)

c     pcp_flx_str = [kg m-2 s-1] Stratiform precipitation flux
c     = CCM:precl (but in mass flux units not m s-1)

c     pcp_flx_cnv = [kg m-2 s-1] Convective precipitation flux
c     = CCM:precc (but in mass flux units not m s-1)

c     wtr2pcp_flx = [kg m-2 s-1] Local conversion of water to precip (q_H2O_wtr2pcp_ttl_tnd*mpl_air)
c     = Total gridcell precipitation formation mass flux
c     = Sum of vapor->precipitation and condensate->precipitation
c     = Until 20000413 wtr2pcp_flx was mislabeled as vpr2pcp_flx but was defined correctly (=vpr2pcp_flx+cnd2pcp_flx)
c     = CCM/MATCH: q_H2O_wtr2pcp_ttl_tnd*mpl_air

c     Total tendency equations in terms of the above:
c     tpt_mdp = tpt_mdp + 
c     $     latvap/cpair*(q_H2O_vpr2cnd_tnd-q_H2O_pcp2vpr)*tm_dlt +
c     $     rmelt(i,k)*tm_dlt
c     q_H2O_vpr = q_H2O_vpr - 
c     $     (q_H2O_vpr2cnd_tnd - q_H2O_pcp2vpr)*tm_dlt
c     q_H2O_cnd = q_H2O_cnd + 
c     $     (q_H2O_vpr2cnd_tnd - q_H2O_cnd2pcp_tnd(i,k))*tm_dlt
c     pcp_flx_str(i) = pcp_flx_str(i) +
c     $     (q_H2O_cnd2pcp_tnd(i,k)-q_H2O_pcp2vpr)*prs_dlt/grv_mean_sfc

c     Requires dst.h for DST_MSS_BDG
#include <dst.h> /* Dust preprocessor tokens */
      
c      subroutine dst_dps_wet(hbuf,obuf,
      subroutine dst_dps_wet(obuf,
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
     $     q_dst,               ! I/O [kg kg-1] Dust mixing ratio
     $     tm_adj,              ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     tpt_mdp)             ! I [K] Temperature
c     $     q_H2O_vpr2cnd_tnd,   ! I [kg kg-1 s-1] H2O vapor to condensate tendency
c     Purpose: Simulate wet deposition processes 
c     dst_dps_wet() is called by CCM:physics/aphys(), MATCH:src/physlic()
c     Algorithm based on Barth, Rasch, and Kiehl (2000) (BRK00)
c++match
      use histout, only: outfld
c--match
      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstblm.h> /* Physical constants for boundary layer meteorology */
#include <dstmsc.h> /* Miscellaneous variables for dust parameterization */
#include <dstnm.h> /* Nomenclature for outfld() */
c     Parameters
      real(r8),parameter::pcp_flx_str_min=1.0e-12 ! [kg m-2 s-1] Minimum stratiform precipitation flux for determining evaporated fraction
      real(r8),parameter::q_H2O_cnd_min=1.0e-12 ! [kg kg-1] Minimum condensed water for determining condensation to precipitation conversion fraction
      real(r8),parameter::cld_frc_cll_min=1.0e-5 ! [frc] Minimum cloud fraction for collision scavenging for determining potential scavenging
      real(r8),parameter::q_dst_cld_min=1.0e-30 ! [kg kg-1] Epsilon to prevent divide by zero for interstitial fraction
      real(r8),parameter::q_trc_tnd_scv_eps=1.0e-30 ! [kg kg-1 s-1] Epsilon to prevent to prevent divide by zero in denominators of scavenged fractions

c     Input
      integer,intent(in)::lat_idx ! I [idx] Model latitude index
c     real(r8),intent(in)::hbuf(*)  ! I [ptr] History buffer
      real(r8),intent(inout)::obuf(*)  ! I [ptr] Output buffer
      real(r8),intent(in)::cld_frc(plond,plevp) ! I [frc] Local total cloud fraction
      real(r8),intent(in)::cld_vlm(plond,plev) ! I [frc] Cloud volume
      real(r8),intent(in)::q_H2O_cnd(plond,plev) ! I [kg kg-1] Condensed H2O mixing ratio
      real(r8),intent(in)::q_H2O_vpr2pcp_cnv_tnd(plond,plev) ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
      real(r8),intent(in)::prs_dlt(plond,plev) ! I [Pa] Pressure thickness
      real(r8),intent(in)::prs_mdp(plond,plev) ! I [Pa] Midlayer pressure
      real(r8),intent(in)::tpt_mdp(plond,plev) ! I [K] Temperature
      real(r8),intent(in)::tm_adj   ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
c     Input/Output
c     Fields besides q_dst are only altered when they are physically inconsistent
      real(r8),intent(inout)::cld_frc_cnv(plond,plev) ! I/O [frc] Convective cloud fraction
      real(r8),intent(inout)::q_H2O_pcp2vpr_tnd(plond,plev) ! I/O [kg kg-1 s-1] H2O precipitation to vapor tendency
      real(r8),intent(inout)::q_H2O_cnd2pcp_tnd(plond,plev) ! I/O [kg kg-1 s-1] Condensed H2O to precipitation tendency
      real(r8),intent(inout)::q_dst(plond,plev,dst_nbr) ! I/O [kg dust/kg moist air]
c     Local Output
      real(r8) flx_mss_pcp_sfc(plond,dst_nbr) ! [kg m-2 s-1] Dust reaching surface in precipitation
      real(r8) flx_mss_pcp_sfc_ttl(plond) ! [kg m-2 s-1] Total dust reaching surface in precipitation
      real(r8) frc_trc_ist(plond,plev,dst_nbr) ! O [frc] Interstitial tracer fraction
      real(r8) pcp_flx_sfc(plond)   ! [kg m-2 s-1] Total precipitation reaching surface
c      real(r8) prect2(plond)        ! [m s-1] Total precipitation reaching surface (diagnostic, should be same as PRECT)
      real(r8) q_dst_tnd_evp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Evaporation tendency
      real(r8) q_dst_tnd_evp_ttl(plond,plev) ! [kg kg-1 s-1] Total evaporation tendency
      real(r8) q_dst_tnd_ncl(plond,plev,dst_nbr) ! O [kg kg-1 s-1] Nucleation scavenging tendency (positive definite)
      real(r8) q_dst_tnd_pcp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Scavenging tendency
      real(r8) q_dst_tnd_pcp_ttl(plond,plev) ! [kg kg-1 s-1] Total scavenging tendency
      real(r8) q_dst_tnd_wet(plond,plev,dst_nbr) ! [kg kg-1 s-1] Wet deposition (evaporation minus scavenging) tendency
      real(r8) spc_xsx_ncl_scv(dst_nbr) ! [m2 kg-1] Specific cross section for nucleation scavenging, GBS98, BJG93
      real(r8) spc_xsx_cll_scv(dst_nbr) ! [m2 kg-1] Specific cross section for collision scavenging, GBS98, BJG93
c     Local
      character fl_out*80       ! [sng] Name of netCDF output file
      integer i                 ! [idx] lon index
      integer k                 ! [idx] lev index
      integer m                 ! [idx] cst index 
      real(r8) cld_frc_cll_cnv      ! [frc] Cloud fraction for collision scavenging by convective rain
      real(r8) cld_frc_cll_str      ! [frc] Cloud fraction for collision scavenging by stratiform rain
      real(r8) cld_frc_eff          ! [frc] Fraction of dust removable by precip
      real(r8) cld_frc_ncl_cnv      ! [frc] Cloud fraction for nucleation scavenging in convective clouds
      real(r8) cld_frc_ncl_str      ! [frc] Cloud fraction for nucleation scavenging in stratiform clouds
c      real(r8) cnd_frc_ice          ! [frc] Ice fraction of condensate
      real(r8) dst2pcp_flx          ! [kg m-2 s-1] Local removal of dust by precip (q_H2O_wtr2pcp_ttl_tnd*q_dst*z_scv)
      real(r8) flx_mss_evp          ! [kg m-2 s-1] Local source of dust due to evaporation of precipitation from above (pcpfrc_evp*flx_mss_pcp_sfc)
      real(r8) flx_mss_scv_abv_cnv(plond,dst_nbr) ! [kg m-2 s-1] Tracer mass flux scavenged from above, convective
      real(r8) flx_mss_scv_abv_str(plond,dst_nbr) ! [kg m-2 s-1] Tracer mass flux scavenged from above, stratiform
      real(r8) mpl_air(plond,plev)  ! [kg m-2] Air mass path in layer
      real(r8) pcp2vpr_flx          ! [kg m-2 s-1] Local conversion of precip to vapor (q_H2O_pcp2vpr_tnd*mpl_air)
      real(r8) pcp_flx_cnv(plond)   ! [kg m-2 s-1] Convective precipitation flux
      real(r8) pcp_flx_str(plond)   ! [kg m-2 s-1] Stratiform precipitation flux
      real(r8) pcp_str_frc_evp      ! [frc] Stratiform precipitation from above that evaporates locally (pcp2vpr_flx/pcp_flx_sfc)
      real(r8) ptn_frc_scv_cll_cnv  ! [frc] Potential tracer fraction removed by collision scavenging in convective precipitation
      real(r8) ptn_frc_scv_cll_str  ! [frc] Potential tracer fraction removed by collisions scavenging in stratiform precipitation
      real(r8) ptn_frc_scv_ncl_cnv  ! [frc] Potential tracer fraction removed by nucleation scavenging in convective precipitation
      real(r8) ptn_frc_scv_ncl_str  ! [frc] Potential tracer fraction removed by nucleation scavenging in stratiform precipitation
      real(r8) q_H2O_cnd2pcp_frc    ! [frc] Fraction of condensate converted to (stratiform) precipitation
      real(r8) q_H2O_wtr2pcp_ttl_tnd(plond,plev) ! [kg kg-1 s-1] Local water to total precipitation tendency
      real(r8) q_trc_tnd_scv_cll_cnv ! [kg kg-1 s-1] Convective rain-collision scavenging tracer flux
      real(r8) q_trc_tnd_scv_cll_str ! [kg kg-1 s-1] Stratiform rain-collision scavenging tracer tendency
      real(r8) q_trc_tnd_scv_cnv    ! [kg kg-1 s-1] Total convective scavenging tracer tendency
      real(r8) q_trc_tnd_scv_ncl_cnv ! [kg kg-1 s-1] Convective nucleation scavenging tracer tendency
      real(r8) q_trc_tnd_scv_ncl_str ! [kg kg-1 s-1] Stratiform nucleation scavenging tracer tendency
      real(r8) q_trc_tnd_scv_str    ! [kg kg-1 s-1] Total stratitform scavenging tracer tendency
      real(r8) scv_ncl_frc_cnv      ! [frc] Nucleation scavenging fraction in convective clouds
      real(r8) scv_ncl_frc_str      ! [frc] Nucleation scavenging fraction in stratiform clouds
      real(r8) scv_tnd_fct          ! Factor limiting scavenging tendencies to amount present [frc]
      real(r8) tm_dlt               ! [s] Wet deposition timestep
      real(r8) tpt_cls              ! [C] Temperature
      real(r8) wtr2pcp_flx          ! [kg m-2 s-1] Local conversion of water to precip (q_H2O_wtr2pcp_ttl_tnd*mpl_air)
c      real(r8) cld_frc_max_abv_cnv  ! [frc] Maximum convective cloud fraction at or above this level
c      real(r8) cld_frc_max_abv_str  ! [frc] Maximum stratiform cloud fraction at or above this level
c      real(r8) mpc_trc_abv          ! [kg m-2] Column integrated tracer path

c     Main Code
c     Timesplit if desired
      tm_dlt=tm_adj             ! [s] (default CCM: 2*dt, MATCH: dt)
      
c     Initialize fluxes and tendencies
      q_dst_tnd_pcp=0.0         ! [kg kg-1 s-1] Scavenging tendency
      q_dst_tnd_evp=0.0         ! [kg kg-1 s-1] Evaporation tendency
      flx_mss_pcp_sfc=0.0       ! [kg m-2 s-1]

c     Specific cross-sections adopted were defined empirically by Balkanski
c     fxm: specific cross-sections should use rain and aerosol size information
c     Set size and species dependent scavenging cross-sections
      do m=1,dst_nbr
c         spc_xsx_ncl_scv(m)=0.1 ! [m2 kg-1] Specific cross section for nucleation scavenging, GBS98, BJG93
         spc_xsx_ncl_scv(m)=0.05 ! [m2 kg-1] Specific cross section for nucleation scavenging, csz 50% of value used for hygrophilic aerosols
c        spc_xsx_cll_scv(m)=0.1 ! [m2 kg-1] Specific cross section for collision scavenging, GBS98, BJG93
         spc_xsx_cll_scv(m)=0.05 ! [m2 kg-1] Specific cross section for collision scavenging, csz 50% of value used for hygrophilic aerosols
      end do                    ! end loop over cst

c     Sanity checks
c     Alter fields that are physically inconsistent
c     fxm: Put this loop in an #ifdef DST_DBG block once the routines responsible are fixed
      do k=1,plev
         do i=1,plon
c     q_H2O_pcp2vpr_tnd, and q_H2O_cnd2pcp_tnd may be < 0.0, though I do not understand why, so just make them positive definite
            q_H2O_pcp2vpr_tnd(i,k)=max(q_H2O_pcp2vpr_tnd(i,k),zero) ! [kg kg-1 s-1] 
            q_H2O_cnd2pcp_tnd(i,k)=max(q_H2O_cnd2pcp_tnd(i,k),zero) ! [kg kg-1 s-1] 
c     if (q_H2O_pcp2vpr_tnd(i,k) < 0.0) write(6,100) 'dst: dst_dps_wet: lat = ',lat_idx,' q_H2O_pcp2vpr_tnd(',i,',',k,') = ',q_H2O_pcp2vpr_tnd(i,k)
c     if (q_H2O_cnd2pcp_tnd(i,k) < 0.0) write(6,100) 'dst: dst_dps_wet: lat = ',lat_idx,' q_H2O_cnd2pcp_tnd(',i,',',k,') = ',q_H2O_cnd2pcp_tnd(i,k)
c     20000612: MATCH sets max(cld_frc=0.999) but allows cld_frc_cnv=1.0
            cld_frc_cnv(i,k)=min(cld_frc_cnv(i,k),cld_frc(i,k))
c            if (cld_frc(i,k) < cld_frc_cnv(i,k)) then
c               write(6,'(a)') 'dst: ERROR dst_dps_wet() reports cloud fraction error:'
c               write(6,'(a,i2,a,i3,a,i2,a,1p,e15.9,0p,a,1p,e15.9,0p)')
c     $              'dst_pcp: lat = ',lat_idx,' cld_frc(',i,',',k,') = ',cld_frc(i,k)
c     $              ,' < cld_frc_cnv = ',cld_frc_cnv(i,k)
c            endif               ! endif
         end do                 ! end loop over lon
      end do                    ! end loop over lev
 100  format(a,i2,a,i3,a,i2,a,1p,e8.1,0p) 
      
c     Compute necessary derived fields
      do k=1,plev
         do i=1,plon
c     Mass of air currently in gridbox
            mpl_air(i,k)=prs_dlt(i,k)*grv_sfc_rcp ! [kg m-2]
c     Total precipitation formation
            q_H2O_wtr2pcp_ttl_tnd(i,k)= ! [kg kg-1 s-1] H2O vapor + condensate conversion to precipitation tendency
     $           q_H2O_cnd2pcp_tnd(i,k)+q_H2O_vpr2pcp_cnv_tnd(i,k)
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      
c     Keep track of amount of stratiform precip and tracer falling from above
      pcp_flx_str=0.0           ! [kg m-2 s-1] Stratiform precipitation flux
      pcp_flx_cnv=0.0           ! [kg m-2 s-1] Convective precipitation flux
      flx_mss_scv_abv_str=0.0   ! [kg m-2 s-1] Tracer mass flux scavenged from above, stratiform
      flx_mss_scv_abv_cnv=0.0   ! [kg m-2 s-1] Tracer mass flux scavenged from above, convective
      
c     Notes on algorithm:
c     This algorithm originated by PJR and MKB for highly soluble aerosols (sulfate)
c     For some species, e.g., HNO3, any tracer within cloud volume is in cloud water
c     For other species, e.g., mineral dust, amount in cloud water varies
c     Convective clouds: We do not know cloud water amount 
c     Assume all convective cloud water (and thus all tracer in convective cloud water) falls out each time step
c     Stratiform clouds: Fraction of cloud water converted to precipitation defines amount of tracer which is scavenged

c     Invert loops (lon is inner dimension) process from top down on each longitude
c     Surface fluxes accumulate as local scavenging and evaporation modify falling precipitation
c     Top-down structure of algorithm based on CCM2 cond() routine of J. J. Hack
      do k=1,plev               ! NB: Inefficient loop order
         do i=1,plon

            tpt_cls=tpt_mdp(i,k)-tpt_frz_pnt ! [C] Temperature

c     Ice fraction is currently not used in dust wet deposition
c     PJR and MCB use cnd_frc_ice to disallow snow from contributing to nucleation scavenging (of, e.g., SO4)
c     CSZ allows dust to nucleate ice, thus allows snow to remove dust
c            cnd_frc_ice=max(zero,min(-tpt_cls*0.05,one)) ! [frc] Ice fraction of condensate
            
c     Convert local condensation and evaporation tendencies to fluxes
            wtr2pcp_flx=q_H2O_wtr2pcp_ttl_tnd(i,k)*mpl_air(i,k) ! [kg m-2 s-1] Local conversion of water to precip
            pcp2vpr_flx=q_H2O_pcp2vpr_tnd(i,k)*mpl_air(i,k) ! [kg m-2 s-1] Local conversion of precip to vapor
            
c     Compute fraction of stratiform rain from above that evaporates in this gridcell
c     This fraction will be applied to tracer in stratitform rain from above to compute
c     tracer released into this gridcell during evaporation of stratiform rain.
c     Use safety limits because pcp_str_frc_evp often equals, e.g., 1.0000000000000007
            pcp_str_frc_evp=min(pcp2vpr_flx/max(pcp_flx_str_min,pcp_flx_str(i)),one) ! [frc]
            
c     Compute fraction of condensate converted to (stratiform) precipitation
c     This fraction determines what fraction of tracer in stratiform clouds is potentially nucleation scavenged
c     [frc] Fraction of condensate converted to (stratiform) precipitation
            q_H2O_cnd2pcp_frc=max(zero,min(one,
     $           q_H2O_cnd2pcp_tnd(i,k)*tm_dlt/max(q_H2O_cnd(i,k),q_H2O_cnd_min)))
c     fxm: PJR and MCB decided not to add condensate formed in current timestep to denominator, why?
c     $              q_H2O_cnd2pcp_tnd(i,k)*tm_dlt/max(q_H2O_cnd(i,k)+q_H2O_vpr2cnd_tnd(i,k)*tm_dlt,q_H2O_cnd_min)))
            
c     Cloud fraction properties for each process are independent of aerosol:
c     MCB and PJR considered multiple candidates for each process, but chose
c     to use cld_vlm for all SO4 processes

c     Cloud used in convective nucleation scavenging should be local gridcell amount
c     MCB and PJR considered cld_frc, cld_frc_cnv, and cld_vlm for SO4
            cld_frc_ncl_cnv=cld_frc_cnv(i,k) ! [frc] Cloud fraction for nucleation scavenging in convective clouds
c     Cloud in convective rain-collision scavenging depends on local and higher cloud
c     MCB and PJR considered max(cld_frc,cld_frc_abv_cnv), max(cld_frc_cnv,cld_frc_abv_cnv), max(cld_vlm,cld_frc_abv_cnv), and simply cld_vlm for SO4
            cld_frc_cll_cnv=cld_vlm(i,k) ! [frc] Cloud fraction for collision scavenging by convective rain
c     Cloud used in stratiform nucleation scavenging should be local gridcell amount
c     MCB and PJR considered cld_frc-cld_frc_cnv, cld_frc, and cld_vlm for SO4
            cld_frc_ncl_str=cld_frc(i,k)-cld_frc_cnv(i,k) ! [frc] Cloud fraction for nucleation scavenging in stratiform clouds
c     Cloud in stratiform rain-collision scavenging depends on local and higher cloud
c     MCB and PJR considered max(cld_frc,cld_frc_abv_str), max(cld_vlm,cld_frc_abv_str), and simply cld_vlm for SO4
            cld_frc_cll_str=cld_vlm(i,k) ! [frc] Cloud fraction for collision scavenging by stratiform rain
            
c     Scavenging processes may be species and size dependent
            do m=1,dst_nbr
               
c     Scavenging in convective clouds:
c     Nucleation scavenging (SeP97 p. 1027):
c     In-cloud scavenging occurs through nucleation scavenging and interstitial aerosol collection by cloud droplets (as opposed to raindrops)
c     Nucleation scavenging is the growth (activation) of CCN to cloud drops
c     The following representation of nucleation scavenging involves a number of assumptions:
c     1. Nucleation scavenging permanently removes aerosol from cloud, in-cloud re-evaporation (de-activation) is not allowed
c     2. Nucleation scavenging in convective clouds is treated here, but the interstitial fraction is not. fxm: this would be somewhat bizarre if it were true
c

c     Interstitial aerosol collection by cloud droplets is inefficient unless the aerosol or the droplet is larger than about 20 microns in diameter
c     However, removal of interstitial aerosol is likely to be much more efficient in convective updrafts than in the stratiform region
c     PJR's wet deposition and convection scheme apparently handle interstitial coalescence separately, but based on the interstitial fraction computed here
c     fxm: need to understand whether some of interstitial fraction is removed in convection routines

c     For dust, interstitial collection by droplets is rare because large dust particles settle out quickly, before precipitation has a chance to act

c     Convective precipitation formation during timestep times specific cross-section
c     for nucleation scavenging yields fraction of homogeneously distributed tracer
c     in liquid cloud that is potentially nucleation scavenged this timestep
c     Species-dependence is incorporated through spc_xsx_ncl_scv factor
c     [frc] Potential tracer fraction removed by nucleation scavenging in convective precipitation
               ptn_frc_scv_ncl_cnv=max(zero,min(one,
     $              q_H2O_vpr2pcp_cnv_tnd(i,k)*mpl_air(i,k)*spc_xsx_ncl_scv(m)*tm_dlt))
               
c     [kg kg-1 s-1] Convective nucleation scavenging tracer tendency
c     Allow dust to nucleate ice, thus allow snow to remove dust
               q_trc_tnd_scv_ncl_cnv=cld_frc_ncl_cnv*ptn_frc_scv_ncl_cnv*q_dst(i,k,m)/tm_dlt ! csz
c     Following line disallows snow from contributing to nucleation scavenging
c              q_trc_tnd_scv_ncl_cnv=cld_frc_ncl_cnv*ptn_frc_scv_ncl_cnv*q_dst(i,k,m)*(1.0-cnd_frc_ice)/tm_dlt
               
c     Convective rain-collision scavenging (SeP97 p. 1021):
c     Proportional to cross-sectional area of gridcell swept by precipitation 
c     Collision efficiency E is probability that collision occurs within the
c     cross-sectional area swept out by the path of the falling raindrop.
c     For rain and hygroscopic aerosol the collision efficiency is the same
c     as the collection efficiency because the sticking efficiency is ~ 1.0. 
c     E is usually << 1 for atmospheric aerosol due to streamlines diverging around  around the falling raindrop
c     Very small (D < 0.1 um) particles are only collected by Brownian diffusion
c     Somewhere near D ~ 1.0 um particles exceed the critical Stokes number and 
c     inertial impaction collection becomes very efficient 
c     Particles D > 10.0 um are also collected by interception
c     In between (0.1 < D < 1.0 um) is the "Greenfield gap" where E << 1
c     A "reasonable" raindrop radius is ~ 1 mm RoY94 p. 179
               
c     Area flux = pcp_flx [kg m-2 s-1] * pi rds^2       [m2 H2O] * dt [s]
c                 ------- ------------   -------------- --------  
c                 dns_H2O [kg m-3 H2O]   (4/3) pi rds^3 [m3 H2O]
               
c     = 3*rds*pcp_flx*dt [m2 H2O m-2]
c       ----------------
c       4
               
c     [frc] Potential tracer fraction removed by collision scavenging in convective precipitation
               ptn_frc_scv_cll_cnv=max(zero,min(one,
     $              pcp_flx_cnv(i)*spc_xsx_cll_scv(m)*tm_dlt/max(cld_frc_cll_cnv,cld_frc_cll_min)))
               
c     [kg kg-1 s-1] Convective rain-collision scavenging tracer tendency
               q_trc_tnd_scv_cll_cnv=cld_frc_cll_cnv*ptn_frc_scv_cll_cnv*q_dst(i,k,m)/tm_dlt
               
c     [kg kg-1 s-1] Total convective scavenging tracer tendency
               q_trc_tnd_scv_cnv=q_trc_tnd_scv_ncl_cnv+q_trc_tnd_scv_cll_cnv
               
c     [frc] Nucleation scavenging fraction in convective clouds
               scv_ncl_frc_cnv=q_trc_tnd_scv_ncl_cnv/(q_trc_tnd_scv_cnv+q_trc_tnd_scv_eps) 
               
c     Scavenging in stratiform clouds:
               
c     Nucleation scavenging (SeP97 p. 1027):
c     This is handled slightly differently than in the convective case
c     First we determine fraction of local condensate that precipitates
c     We assume same fraction of tracer is removed from liquid portion of cloud
c     Species-dependence is currently not incorporated, i.e., there is no
c     factor to account for the scavenging efficiency for in-cloud nucleation
               
c     [kg kg-1 s-1] Stratiform nucleation scavenging tracer tendency
c     Allow dust to nucleate ice, thus allow snow to remove dust
               q_trc_tnd_scv_ncl_str=cld_frc_ncl_str*q_H2O_cnd2pcp_frc*q_dst(i,k,m)/tm_dlt ! csz
c     Following line disallows snow from contributing to nucleation scavenging
c               q_trc_tnd_scv_ncl_str=cld_frc_ncl_str*q_H2O_cnd2pcp_frc*q_dst(i,k,m)*(1.0-cnd_frc_ice)/tm_dlt
               
c     Stratiform rain-collision scavenging (SeP97 p. 1021):
               
c     [frc] Potential tracer fraction removed by collisions scavenging in stratiform precipitation
               ptn_frc_scv_cll_str=max(zero,min(one,
     $              pcp_flx_str(i)*spc_xsx_cll_scv(m)*tm_dlt/max(cld_frc_cll_str,cld_frc_cll_min)))
               
c     [kg kg-1 s-1] Stratiform rain-collision scavenging tracer tendency
               q_trc_tnd_scv_cll_str=cld_frc_cll_str*ptn_frc_scv_cll_str*q_dst(i,k,m)/tm_dlt
               
c     [kg kg-1 s-1] Total stratiform scavenging tracer tendency
               q_trc_tnd_scv_str=q_trc_tnd_scv_ncl_str+q_trc_tnd_scv_cll_str
               
c     [frc] Nucleation scavenging fraction in stratiform clouds
               scv_ncl_frc_str=q_trc_tnd_scv_ncl_str/(q_trc_tnd_scv_str+q_trc_tnd_scv_eps) 
               
c     Ensure sum of convective and stratiform scavenging tendencies does not exceed amount present
               scv_tnd_fct=q_dst(i,k,m)/max(tm_dlt*(q_trc_tnd_scv_cnv+q_trc_tnd_scv_str),q_trc_tnd_scv_eps) ! [frc]
               if (scv_tnd_fct < 1.0) then
c     Total stratiform scavenging tracer tendency
                  q_trc_tnd_scv_str=q_trc_tnd_scv_str*scv_tnd_fct ! [kg kg-1 s-1]
c     Total convective scavenging tracer tendency
                  q_trc_tnd_scv_cnv=q_trc_tnd_scv_cnv*scv_tnd_fct ! [kg kg-1 s-1]
c     It is now "safe" to use these tendencies to adjust the tracer mixing ratio
c     However, it is only "safe" to machine precision since order of operations can still cause problems
               endif            ! endif
               
c     Nomenclature: 
c     Scavenging tendency refers to removal by nucleation and collision scavenging [positive definite]
c     Evaporation tendency refers to source by evaporation of stratiform precipitation from above [positive definite]
c     Wet deposition tendency is evaporation tendency (source) minus scavenging tendency (sink) [positive when source > sink]
c     [kg kg-1 s-1] Scavenging tendency
               q_dst_tnd_pcp(i,k,m)=q_trc_tnd_scv_str+q_trc_tnd_scv_cnv 
               
c     [kg kg-1 s-1] Evaporation tendency
               q_dst_tnd_evp(i,k,m)=pcp_str_frc_evp*flx_mss_scv_abv_str(i,m)/mpl_air(i,k)
               
c     Tracer not removed within cloud is interstitial and subject to convective transport
               frc_trc_ist(i,k,m)= ! [frc] Interstitial tracer fraction
     $              1.0-max(zero,min(one,
     $              q_dst_tnd_pcp(i,k,m)*tm_dlt/max(cld_vlm(i,k)*q_dst(i,k,m),q_dst_cld_min)))
               
c     Diagnose nucleation scavenging tendency for possible archival (superfluous, could be deleted)
               q_dst_tnd_ncl(i,k,m)= ! [kg kg-1 s-1] Nucleation scavenging tendency
     $              q_trc_tnd_scv_cnv*scv_ncl_frc_cnv+q_trc_tnd_scv_str*scv_ncl_frc_str
               
c     Diagnose wet deposition (evaporation minus scavenging) tendency for possible archival (superflous, could be deleted)
               q_dst_tnd_wet(i,k,m)=q_dst_tnd_evp(i,k,m)-q_dst_tnd_pcp(i,k,m) ! [kg kg-1 s-1]
               
c     Tendencies are known, so adjust mixing ratios
c     max() clause is necessary because result cannot be guarranteed positive otherwise
               q_dst(i,k,m)=max(zero,q_dst(i,k,m)+tm_dlt*q_dst_tnd_wet(i,k,m)) ! [kg kg-1] Dust mixing ratio

c     Accumulate aerosol-dependent vertical integrals
c     All tracer in stratiform rain may potentially be evaporated at every level
               flx_mss_scv_abv_str(i,m)=flx_mss_scv_abv_str(i,m)*(1.0-pcp_str_frc_evp)+q_trc_tnd_scv_str*mpl_air(i,k) ! [kg m-2 s-1] Tracer mass flux scavenged from above, stratiform
c     Tracer in convective rain increases monotonically on way down--convective evaporation is not accounted for
               flx_mss_scv_abv_cnv(i,m)=flx_mss_scv_abv_cnv(i,m)+q_trc_tnd_scv_cnv*mpl_air(i,k) ! [kg m-2 s-1] Tracer mass flux scavenged from above, convective
               
            end do              ! end loop over cst
            
c     Accumulate aerosol-independent vertical integrals
c     Precipitation falling into cell below is precipitation that entered this cell from above, 
c     plus local source of precipitation, minus what was evaporated. 
            pcp_flx_str(i)=pcp_flx_str(i)+(q_H2O_cnd2pcp_tnd(i,k)-q_H2O_pcp2vpr_tnd(i,k))*mpl_air(i,k) ! [kg m-2 s-1] Stratiform precipitation flux
c     Evaporation of convective precipitation is not allowed here
            pcp_flx_cnv(i)=pcp_flx_cnv(i)+q_H2O_vpr2pcp_cnv_tnd(i,k)*mpl_air(i,k) ! [kg m-2 s-1] Convective precipitation flux

         end do                 ! end loop over lon
      end do                    ! end loop over lev
      
#ifdef DST_DBG
c     Sanity check
      do m=1,dst_nbr
         do k=1,plev
            do i=1,plon
               if (q_dst(i,k,m) < 0.0 .or. q_dst(i,k,m) > q_dst_mxm) then
                  write(6,'(a)') 'dst: ERROR dst_dps_wet() reports mixing ratio too high or too low:'
                  write(6,'(a,i2,a,i3,a,i2,a,i2,a,1p,e8.1,0p,a,
     $                 2(a,1p,e8.1,0p,a))') 
     $                 'dst_pcp: lat = ',lat_idx,' q_dst(',i,',',k,',',m,') = ',q_dst(i,k,m),' kg kg-1'
     $                 ,', q_dst_tnd_evp = ',q_dst_tnd_evp(i,k,m),' kg kg-1 s-1'
     $                 ,', q_dst_tnd_pcp = ',q_dst_tnd_pcp(i,k,m),' kg kg-1 s-1'
                  stop
               endif            ! endif err
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst
#endif /* not DST_DBG */

c     Diagnostic totals
      do i=1,plon
         pcp_flx_sfc(i)=pcp_flx_str(i)+pcp_flx_cnv(i) ! [kg m-2 s-1] Total precipitation reaching surface
cc     Convert [kg m-2 s-1] to [m s-1]
c         prect2(i)=pcp_flx_sfc(i)*0.001 ! [m s-1] Total precipitation reaching surface (same as PRECT)
      end do                    ! end loop over lon
      do m=1,dst_nbr
         do i=1,plon
c     Tracer reaching surface in precipitation
            flx_mss_pcp_sfc(i,m)=flx_mss_scv_abv_str(i,m)+flx_mss_scv_abv_cnv(i,m) ! [kg m-2 s-1]
         end do                 ! end loop over cst
      end do                    ! end loop over lon
      
c     Integrate over all size categories and levels for diagnostic output
      call dst_add_lon_lev(q_dst_tnd_evp,q_dst_tnd_evp_ttl)
      call dst_add_lon_lev(q_dst_tnd_pcp,q_dst_tnd_pcp_ttl)
      call dst_add_lon(flx_mss_pcp_sfc,flx_mss_pcp_sfc_ttl)
      
#ifdef BXM
c     Update netCDF file
      fl_out='aer.nc'           ! [sng] Name of netCDF output file
      call strnul(fl_out)
      call dpswet2nc(            
     $     cld_frc,             ! I [frc] Local total cloud fraction
     $     cld_frc_cnv,         ! I [frc] Convective cloud fraction
     $     cld_vlm,             ! I [frc] Cloud volume
     $     fl_out,              ! I [sng] Name of netCDF output file
     $     flx_mss_pcp_sfc,     ! I [kg m-2 s-1] Dust reaching surface in precipitation
     $     flx_mss_pcp_sfc_ttl, ! I [kg m-2 s-1] Total dust reaching surface in precipitation
     $     frc_trc_ist,         ! I [frc] Interstitial tracer fraction
     $     lat_idx,             ! I [idx] Model latitude index
     $     pcp_flx_sfc,         ! I [kg m-2 s-1] Total precipitation reaching surface
     $     q_H2O_cnd,           ! I [kg kg-1] Condensed H2O mixing ratio
     $     q_H2O_cnd2pcp_tnd,   ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
     $     q_H2O_pcp2vpr_tnd,   ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
     $     q_H2O_vpr2pcp_cnv_tnd, ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
     $     q_dst_tnd_evp,       ! I [kg kg-1 s-1] Evaporation tendency
     $     q_dst_tnd_evp_ttl,   ! I [kg kg-1 s-1] Total evaporation tendency
     $     q_dst_tnd_ncl,       ! I [kg kg-1 s-1] Nucleation scavenging tendency
     $     q_dst_tnd_pcp,       ! I [kg kg-1 s-1] Scavenging tendency
     $     q_dst_tnd_pcp_ttl,   ! I [kg kg-1 s-1] Total scavenging tendency
     $     q_dst_tnd_wet,       ! I [kg kg-1 s-1] Wet deposition (evaporation minus scavenging) tendency
     $     spc_xsx_cll_scv,     ! I [m2 kg-1] Specific cross section for collision scavenging
     $     spc_xsx_ncl_scv)     ! I [m2 kg-1] Specific cross section for nucleation scavenging
#endif /* not BXM */
#ifdef CCM
      call outfld('DSTSFPCP',flx_mss_pcp_sfc_ttl,plond,lat_idx,hbuf)
      call outfld('DSTSSPCP',q_dst_tnd_pcp_ttl,plond,lat_idx,hbuf)
c-drb      call outfld('DSTSSEVP',q_dst_tnd_evp_ttl,plond,lat_idx,hbuf)
      do m=1,dst_nbr
         call outfld(flx_mss_pcp_sfc_nm(m),flx_mss_pcp_sfc(1,m),plond,lat_idx,hbuf)
      end do                    ! end loop over cst
#else /* not CCM */ 
c      call outfld('PRECT2',prect2,plond,lat_idx,obuf)
      call outfld('DSTSFPCP',flx_mss_pcp_sfc_ttl,plond,lat_idx,obuf)
      call outfld('DSTSSPCP',q_dst_tnd_pcp_ttl,plond,lat_idx,obuf)
c-drb      call outfld('DSTSSEVP',q_dst_tnd_evp_ttl,plond,lat_idx,obuf)
      do m=1,dst_nbr
         call outfld(flx_mss_pcp_sfc_nm(m),flx_mss_pcp_sfc(1,m),plond,lat_idx,obuf)
      end do                    ! end loop over cst
#endif /* not CCM */
#ifdef DST_MSS_BDG
      call bdg_aal('dst_sf_pcp',lat_idx,flx_mss_pcp_sfc_ttl)
      call bdg_gam_wet('dst_ss_pcp',lat_idx,prs_dlt,q_dst_tnd_pcp_ttl)
      call bdg_gam_wet('dst_ss_evp',lat_idx,prs_dlt,q_dst_tnd_evp_ttl)
#endif /* not DST_MSS_BDG */
      return
      end subroutine dst_dps_wet                      ! end dst_dps_wet()

c++mahow
c-db      subroutine dst_dps_wet_old(hbuf,obuf,wrhstts,
      subroutine dst_dps_wet_old(obuf,
     $     cld_frc,             ! I [frc] Local total cloud fraction
     $     lat_idx,             ! I [idx] Model latitude index
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     prs_mdp,             ! I [Pa] Pressure
     $     q_H2O_cnd2pcp_tnd,   ! I/O [kg kg-1 s-1] Condensed H2O to precipitation tendency
     $     q_H2O_pcp2vpr_tnd,   ! I/O [kg kg-1 s-1] H2O precipitation to vapor tendency
     $     q_H2O_vpr2pcp_cnv_tnd, ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
     $     q_dst,               ! I/O [kg kg-1] Dust mixing ratio
     $     tm_adj,              ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     $     tpt_mdp)             ! I [K] Temperature
c     Simulate wet deposition processes 
c     dst_dps_wet_old() is called by CCM:physics/aphys(), MATCH:src/physlic()
c     Current processes considered are:
c     1. Nucleation scavenging
c     2. Evaporation of dust in stratiform rain evaporation
c     Algorithm based on Tegen and Fung 1995
c     Assumes dust is homogeneously distributed throughout a gridbox
c     Scavenging, however, should only be allowed to occur in that fraction of a gridbox which is covered by clouds
c     fxm: 20000525 Slight differences in diagnosed precipitation between this routine and host model
c     might be due to applying stratiform evaporation tendency to both stratiform and convective rain
c++db
      use histout, only: outfld
c--db

      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstmsc.h> /* Miscellaneous variables for dust parameterization */
#include <dstnm.h> /* Nomenclature for outfld() */
c     Input
      integer,intent(in)::lat_idx           ! Model latitude index
c-db      logical,intent(in)::wrhstts           ! Flag to output buffer to file this timestep
c-db      real(r8),intent(in)::hbuf(*)              ! History buffer
      real(r8),intent(in)::obuf(*)              ! Output buffer
      real(r8),intent(in)::cld_frc(plond,plevp) ! [frc] Local total cloud fraction
      real(r8),intent(in)::q_H2O_vpr2pcp_cnv_tnd(plond,plev) ! [kg kg-1 s-1] H2O vapor to convective precipitation tendency
      real(r8),intent(in)::prs_dlt(plond,plev)  ! [Pa] Pressure thickness
      real(r8),intent(in)::prs_mdp(plond,plev)  ! [Pa] Midlayer pressure
      real(r8),intent(in)::tpt_mdp(plond,plev)  ! [K] Temperature
      real(r8),intent(in)::tm_adj               ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
c     Input/Output
      real(r8),intent(inout)::q_H2O_pcp2vpr_tnd(plond,plev) ! [kg kg-1 s-1] H2O precipitation to vapor tendency
      real(r8),intent(inout)::q_H2O_cnd2pcp_tnd(plond,plev) ! [kg kg-1 s-1] Condensed H2O to precipitation tendency
      real(r8),intent(inout)::q_dst(plond,plev,dst_nbr) ! [kg dust/kg moist air]
c     Local Output
c     Local
      integer i                 ! [idx] lon index
      integer k                 ! [idx] lev index
      integer m                 ! [idx] cst index 
      real(r8) mpl_air(plond,plev)  ! [kg m-2] Air mass path in layer
      real(r8) tm_dlt               ! [s] Wet deposition timestep
      real(r8) q_H2O_wtr2pcp_ttl_tnd(plond,plev) ! [kg kg-1 s-1] Local water to total precipitation tendency
c     GCM diagnostics
      real(r8) cld_frc_eff          ! [frc] Fraction of dust removable by precip
      real(r8) dst2pcp_flx          ! [kg m-2 s-1] Local removal of dust by precip (q_H2O_wtr2pcp_ttl_tnd*q_dst*z_scv)
      real(r8) flx_mss_evp          ! [kg m-2 s-1] Local source of dust due to evaporation of precipitation from above (pcpfrc_evp*flx_mss_pcp_sfc)
      real(r8) flx_mss_pcp_sfc(plond,dst_nbr) ! [kg m-2 s-1] Dust reaching surface in precipitation
      real(r8) flx_mss_pcp_sfc_ttl(plond) ! [kg m-2 s-1] Total dust reaching surface in precipitation
      real(r8) pcp2vpr_flx          ! [kg m-2 s-1] Local conversion of precip to vapor (q_H2O_pcp2vpr_tnd*mpl_air)
      real(r8) pcp_flx_sfc(plond)   ! [kg m-2 s-1] Total precipitation reaching surface
c      real(r8) prect2(plond)        ! [m s-1] Total precipitation reaching surface (diagnostic, should be same as PRECT)
      real(r8) pcp_frc_evp          ! [frc] Precipitation from above that evaporates locally (pcp2vpr_flx/pcp_flx_sfc)
      real(r8) q_dst_tnd_evp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Evaporation tendency
      real(r8) q_dst_tnd_evp_ttl(plond,plev) ! [kg kg-1 s-1] Total evaporation tendency
      real(r8) q_dst_tnd_pcp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Scavenging tendency
      real(r8) q_dst_tnd_pcp_ttl(plond,plev) ! [kg kg-1 s-1] Total scavenging tendency
      real(r8) wtr2pcp_flx          ! [kg m-2 s-1] Local conversion of water to precip (q_H2O_wtr2pcp_ttl_tnd*mpl_air)
      
c     Main Code
c     Timesplit if desired
      tm_dlt=tm_adj             ! [s] (default CCM: 2*dt, MATCH: dt)
      
c     Initialize fluxes and tendencies
      q_dst_tnd_pcp=0.0         ! [kg kg-1 s-1]
      q_dst_tnd_evp=0.0         ! [kg kg-1 s-1]
      flx_mss_pcp_sfc=0.0       ! [kg m-2 s-1]

c     Sanity check
      do k=1,plev
         do i=1,plon
c     q_H2O_pcp2vpr_tnd, and q_H2O_cnd2pcp_tnd may be < 0.0, though I do not understand why
            q_H2O_pcp2vpr_tnd(i,k)=max(q_H2O_pcp2vpr_tnd(i,k),zero) ! [kg kg-1 s-1] 
            q_H2O_cnd2pcp_tnd(i,k)=max(q_H2O_cnd2pcp_tnd(i,k),zero) ! [kg kg-1 s-1] 
c     if (q_H2O_pcp2vpr_tnd(i,k) < 0.0) write(6,100) 'dst: dst_dps_wet: lat = ',lat_idx,' q_H2O_pcp2vpr_tnd(',i,',',k,') = ',q_H2O_pcp2vpr_tnd(i,k)
c     if (q_H2O_cnd2pcp_tnd(i,k) < 0.0) write(6,100) 'dst: dst_dps_wet: lat = ',lat_idx,' q_H2O_cnd2pcp_tnd(',i,',',k,') = ',q_H2O_cnd2pcp_tnd(i,k)
         end do                 ! end loop over lon
      end do                    ! end loop over lev
 100  format(a,i2,a,i3,a,i2,a,1p,e8.1,0p) 
      
c     Compute necessary derived fields
      do k=1,plev
         do i=1,plon
c     Mass of air currently in gridbox
            mpl_air(i,k)=prs_dlt(i,k)*grv_sfc_rcp ! [kg m-2]
c     Total precipitation formation
            q_H2O_wtr2pcp_ttl_tnd(i,k)= ! [kg kg-1 s-1] H2O vapor + condensate conversion to precipitation tendency
     $           q_H2O_cnd2pcp_tnd(i,k)+q_H2O_vpr2pcp_cnv_tnd(i,k)
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      
c     Scavenging of dust by precipitation (wet deposition)
c     Scavenging ratio relates precip tendency to dust removal tendency
      do m=1,dst_nbr
         do k=1,plev
            do i=1,plon

c     Compute gridpoint dust tendency for wet deposition
c     Currently, application of cloud fraction to dust in rain is rather haphazard
c     Needs to be consistent with physical picture of a cloud and rainshaft
c     fxm: Dust should be last thing to evaporate in a raindrop since it is a CCN
c     Thus, dust should only be evaporated from virga
c     Using cld_frc_eff = cld_frc(i,k) (dynamic cloud fraction) results in lifetime of clay of 30 days
c     In CCM, using cld_frc_eff = 0.85 results in precipitation lifetime of clay of 12-13 days for dynamic source scheme, 19-20 days for TOMS scheme 
c     In MATCH, using cld_frc_eff = 0.85 results in precipitation lifetime of clay of 4--5 days for dynamic source scheme, ??? days for TOMS scheme 
c     For a given model, The DSS/TOMS timescales differ because TOMS scheme injects more dust into subtropical easterlies which undergo fewer precipitation events (??)
c     Actual lifetime of clay should be 14--21 days
c     fxm: should use cld_frc_eff=f(cld_frc), should be fixed by PJR wet dps routine
               cld_frc_eff=0.85 ! [frc] Until dstmch18
c     cld_frc_eff=cld_frc(i,k) ! [frc] As of dstmch18
               q_dst_tnd_pcp(i,k,m)= ! [kg kg-1 s-1]
     $              cld_frc_eff*z_scv(m)*q_dst(i,k,m)*q_H2O_wtr2pcp_ttl_tnd(i,k)
c     Do not scavenge more than available dust
               if (tm_dlt*q_dst_tnd_pcp(i,k,m) < q_dst(i,k,m)) then
c     Tendencies are known, so adjust mixing ratios
                  q_dst(i,k,m)=q_dst(i,k,m)-
     $                 tm_dlt*q_dst_tnd_pcp(i,k,m) ! [kg kg-1]
               else             ! endif physical
c     Only bother printing a warning when it concerns significant amounts dust
                  if (q_dst(i,k,m) > q_dst_sgn) write (6,'(a,i2,a,i3,a,i2,a,i2,a,1p,e8.1,0p)') 
     $                 'dst: pcp_fxr: lat = ',lat_idx,' q_dst_tnd_pcp(',i,',',k,',',m,') = ',q_dst_tnd_pcp(i,k,m)
                  if (q_dst(i,k,m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,i2,a,1p,e8.1,0p)') 
     $                 'dst: pcp_fxr, set to 0.: lat = ',lat_idx,' q_dst(',i,',',k,',',m,') = ',q_dst(i,k,m)
c     Reset to physical value
                  q_dst_tnd_pcp(i,k,m)=q_dst(i,k,m)/tm_dlt ! [kg kg-1 s-1]
                  q_dst(i,k,m)=0.0 ! [kg kg-1]
               endif            ! endelse unphysical
#ifdef DST_DBG
               if (q_dst(i,k,m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,i2,a,1p,e8.1,0p)') 
     $              'dst_pcp: lat = ',lat_idx,' q_dst(',i,',',k,',',m,') = ',q_dst(i,k,m)
#endif /* not DST_DBG */ 
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst

c     Initialize surface precipitation with top level precipitation
c     Assume no evaporation in top level
      do i=1,plon
         pcp_flx_sfc(i)=q_H2O_wtr2pcp_ttl_tnd(i,1)*mpl_air(i,1) ! [kg m-2 s-1]
      end do                    ! end loop over lon
      do m=1,dst_nbr
         do i=1,plon
            flx_mss_pcp_sfc(i,m)=q_dst_tnd_pcp(i,1,m)*mpl_air(i,1) ! [kg m-2 s-1]
            q_dst_tnd_evp(i,1,m)=0.0 ! [kg kg-1 s-1]
         end do                 ! end loop over lon
      end do                    ! end loop over cst
c     Invert loops (lon is inner dimension) process from top down on each longitude
c     Surface fluxes accumulate as local washout and evaporation modify precipitation on its way down
c     Algorithm is based on JJH's cond() routine
      do k=2,plev               ! NB: level starts with 2
         do i=1,plon
c     Convert local condensation and evaporation tendencies to fluxes
            wtr2pcp_flx=q_H2O_wtr2pcp_ttl_tnd(i,k)*mpl_air(i,k) ! [kg m-2 s-1]
            pcp2vpr_flx=q_H2O_pcp2vpr_tnd(i,k)*mpl_air(i,k) ! [kg m-2 s-1]

c     Compute fraction of rain from above that evaporates in this gridcell
c     This fraction will be applied to dust in rain from above to compute
c     dust released into this cell due to evaporation.
            if (pcp_flx_sfc(i) > 0.0) then
c     Use safety limiter because pcp_frc_evp often equals, e.g., 1.0000000000000007
               pcp_frc_evp=min(pcp2vpr_flx/pcp_flx_sfc(i),one) ! [frc]
#ifdef DST_DBG
               if (pcp_frc_evp > 1.0.or.pcp_frc_evp < 0.0) write(6,'(a,i2,a,i3,a,i2,a,1p,e8.1,0p)') 
     $              'dst: pcp_evp: lat = ',lat_idx,' pcp_frc_evp(',i,',',k,') = ',pcp_frc_evp
#endif /* not DST_DBG */ 
            else
               pcp_frc_evp=0.0  ! [frc]
            endif               ! endif pcp > 0
            
c     Precipitation falling into cell below is precipitation entering cell from above, 
c     plus local source of precipitation, minus what was evaporated. 
            pcp_flx_sfc(i)=max( ! [kg m-2 s-1]
     $           pcp_flx_sfc(i)+wtr2pcp_flx-
     $           pcp2vpr_flx
     $           ,zero)
            do m=1,dst_nbr
               dst2pcp_flx=q_dst_tnd_pcp(i,k,m)*mpl_air(i,k) ! [kg m-2 s-1]
c     Evaporate dust into cell based on its concentration in rain from above
               flx_mss_evp=pcp_frc_evp*flx_mss_pcp_sfc(i,m) ! [kg m-2 s-1]
               q_dst_tnd_evp(i,k,m)=flx_mss_evp/mpl_air(i,k) ! [kg kg-1 s-1]
c     Do not let more dust be released by evaporation than is falling through gridcell
               if (q_dst_tnd_evp(i,k,m) > flx_mss_pcp_sfc(i,m)/mpl_air(i,k)) then


c     Reset to physical value
                  q_dst_tnd_evp(i,k,m)=flx_mss_pcp_sfc(i,m)/(mpl_air(i,k)*tm_dlt) ! [kg kg-1 s-1]
               endif            ! endif unphysical
c     Tendencies are known, so adjust mixing ratios
               q_dst(i,k,m)=q_dst(i,k,m)+
     $              tm_dlt*q_dst_tnd_evp(i,k,m) ! [kg kg-1]
#ifdef DST_DBG
               if (q_dst(i,k,m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,i2,a,1p,e8.1,0p)') 
     $              'dst_evp: lat = ',lat_idx,' q_dst(',i,',',k,',',m,') = ',q_dst(i,k,m)
#endif /* not DST_DBG */ 
               
c     Wet dust falling into next cell down is wet dust entering cell from above, 
c     plus local source of wet dust, minus what was evaporated. 
               flx_mss_pcp_sfc(i,m)=max(
     $              flx_mss_pcp_sfc(i,m)+dst2pcp_flx-
     $              flx_mss_evp
     $              ,zero)
            end do              ! end loop over cst
         end do                 ! end loop over lon
      end do                    ! end loop over lev

c     Get rid of unphysical values
c     20000415: Better to allow qneg() routines to perform these checks and fixes
      if (.false.) then
         do m=1,dst_nbr
            do k=1,plev
               do i=1,plon
                  q_dst(i,k,m)=max(q_dst(i,k,m),zero) ! [kg kg-1]
                  if (q_dst(i,k,m) < q_dst_sgn) q_dst(i,k,m)=0.0 ! [kg kg-1]
               end do           ! end loop over lon
            end do              ! end loop over lev
         end do                 ! end loop over cst
      endif                     ! endif false

c     Integrate over all size categories and levels for diagnostic output
      call dst_add_lon_lev(q_dst_tnd_evp,q_dst_tnd_evp_ttl)
      call dst_add_lon_lev(q_dst_tnd_pcp,q_dst_tnd_pcp_ttl)
      call dst_add_lon(flx_mss_pcp_sfc,flx_mss_pcp_sfc_ttl)
      
#ifdef CCM
      call outfld('DSTSFPCP',flx_mss_pcp_sfc_ttl,plond,lat_idx,hbuf)
      call outfld('DSTSSPCP',q_dst_tnd_pcp_ttl,plond,lat_idx,hbuf)
c-drb      call outfld('DSTSSEVP',q_dst_tnd_evp_ttl,plond,lat_idx,hbuf)
      do m=1,dst_nbr
         call outfld(flx_mss_pcp_sfc_nm(m),flx_mss_pcp_sfc(1,m),plond,lat_idx,hbuf)
      end do                    ! end loop over cst
#else /* not CCM */ 
c      do i=1,plon
cc     Convert [kg m-2 s-1] to [m s-1]
c         prect2(i)=pcp_flx_sfc(i)*0.001 ! [m s-1] Total precipitation reaching surface (same as PRECT)
c      end do                    ! end loop over lon
c      call outfld('PRECT2',prect2,plond,wrhstts,lat_idx,hbuf,obuf)
c-db      call outfld('DSTSFPCP',flx_mss_pcp_sfc_ttl,plond,wrhstts,lat_idx,hbuf,obuf)
c-db      call outfld('DSTSSPCP',q_dst_tnd_pcp_ttl,plond,wrhstts,lat_idx,hbuf,obuf)
c-db      call outfld('DSTSSEVP',q_dst_tnd_evp_ttl,plond,wrhstts,lat_idx,hbuf,obuf)
      call outfld('DSTSFPCP',flx_mss_pcp_sfc_ttl,plond,lat_idx,obuf)
      call outfld('DSTSSPCP',q_dst_tnd_pcp_ttl,plond,lat_idx,obuf)
c-drb      call outfld('DSTSSEVP',q_dst_tnd_evp_ttl,plond,lat_idx,obuf)
      do m=1,dst_nbr
c-db         call outfld(flx_mss_pcp_sfc_nm(m),flx_mss_pcp_sfc(1,m),plond,wrhstts,lat_idx,hbuf,obuf)
         call outfld(flx_mss_pcp_sfc_nm(m),flx_mss_pcp_sfc(1,m),plond,lat_idx,obuf)
      end do                    ! end loop over cst
#endif /* not CCM */
#ifdef DST_MSS_BDG
      call bdg_aal('dst_sf_pcp',lat_idx,flx_mss_pcp_sfc_ttl)
      call bdg_gam_wet('dst_ss_pcp',lat_idx,prs_dlt,q_dst_tnd_pcp_ttl)
      call bdg_gam_wet('dst_ss_evp',lat_idx,prs_dlt,q_dst_tnd_evp_ttl)
#endif /* not DST_MSS_BDG */
      return
      end subroutine  dst_dps_wet_old                 ! end dst_dps_wet_old()

c--mahow

      subroutine cld_dgn( 
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
c     Purpose: Diagnose statistics of cloud properties
c     Same as MATCH:src/cloud/clddiag()
c     cld_dgn() is called by MATCH:physlic() and CCM:aphys()
      implicit none
#include <pmgrid.h>
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstblm.h> /* Physical constants for boundary layer meteorology */
c     Parameters
      real(r8),parameter::rds_pcp=1000.0 ! [m] Radius of precipitation RoY94 p. 179
      real(r8),parameter::q_H2O_pcp_lqd_min=1.0e-14 ! [kg kg-1] Minimum rainwater mixing ratio
      real(r8),parameter::pcp_flx_sfc_min=1.0e-30 ! [kg m-2 s-1] Minimum precipitation rate to prevent divide by zero
c     Input
      real(r8),intent(in)::cld_frc(plond,plevp) ! I [frc] Local total cloud fraction
      real(r8),intent(in)::prs_dlt(plond,plev)  ! I [Pa] Pressure thickness
      real(r8),intent(in)::prs_mdp(plond,plev)  ! I [Pa] Midlayer pressure
      real(r8),intent(in)::q_H2O_cnd2pcp_tnd(plond,plev) ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
      real(r8),intent(in)::q_H2O_cnd_tnd(plond,plev) ! I [kg kg-1 s-1] Net H2O condensate formation tendency
      real(r8),intent(in)::q_H2O_pcp2vpr_tnd(plond,plev) ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
      real(r8),intent(in)::q_H2O_vpr2pcp_cnv_tnd(plond,plev) ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
      real(r8),intent(in)::tpt_mdp(plond,plev)  ! I [K] Temperature
c     Output
      real(r8),intent(out)::cld_vlm(plond,plev) ! O [frc] Cloud volume
      real(r8),intent(out)::q_H2O_pcp_lqd(plond,plev) ! O [kg kg-1] Rainwater mixing ratio
c     Local
      integer i                 ! [idx] Longitude index
      integer k                 ! [idx] Level index
      real(r8) cld_frc_pcp_wgt(plond) ! [kg m-2 s-1] Integral of precipitation-weighted cloud fraction
      real(r8) dns_mdp              ! [kg m-3] Midlayer density
      real(r8) pcp_flx_sfc(plond)   ! [kg m-2 s-1] Precipitation rate 
      real(r8) pcp_flx_pst_sfc(plond) ! [kg m-2 s-1] Sum of positive definite precipitation production fluxes
      real(r8) pcp_flx_lcl          ! [kg m-2 s-1] Local net precipitation production rate (+'ve when formation exceeds evaporation)
      real(r8) pcp_flx_lcl_pst      ! [kg m-2 s-1] Positive definite local net precipitation production rate
      real(r8) vlc_grv_fct          ! [kg1/2 m-1/2 s-1] Factor in terminal velocity formula SeP97 p. 467 (8.44)
      real(r8) vlc_grv_pcp          ! [m s-1] Fallspeed of precipitation

c     Collision scavenging depends on horizontal fraction of local gridcell through which rain from above falls
c     This ought to be true for collision scavenging by both stratiform and convective rain
c     Nucleation scavenging should depend solely on local cloud fraction
c     This ought to be true for nucleation scavenging by both stratiform and convective clouds
c     Note that pre-CCM4, no information about is known about vertical correlation functions between stratiform cloud at different levels

c     In developing the SO4 wet deposition scheme, PJR and MKB settled on using a single cloud fraction,
c     called cld_vlm, for all wet scavenging processes (nucleation and collision)
c     They define this factor as the maximum of local cloud amount and
c
c     sum above of (cloud*positive precipitation production)   sum real precipitation from above
c     ------------------------------------------------------ x --------------------------------------------
c     sum above of (      positive precipitation production)   sum positive precipitation from above
c
c     According to this definition, cld_vlm 
c     1. increases monotonically from the upper troposphere to the lower troposphere
c     2. is unity in and everywhere beneath the first socked-in layer

      vlc_grv_fct=1.94*2.13*sqrt(dns_H2O_lqd_std*rds_pcp*grv_sfc*2.7e-4) ! [kg1/2 m-1/2 s-1] Factor in terminal velocity formula SeP97 p. 467 (8.44)
      do i=1,plon
         pcp_flx_sfc(i)=0.0     ! [kg m-2 s-1] Precipitation rate
         cld_frc_pcp_wgt(i)=0.0 ! [kg m-2 s-1] Integral of precipitation-weighted cloud fraction
         pcp_flx_pst_sfc(i)=pcp_flx_sfc_min ! [kg m-2 s-1] Sum of positive definite precipitation production fluxes
      end do                    ! end loop over lon

      do k=1,plev
         do i=1,plon
            cld_vlm(i,k)=       ! [frc] Cloud volume
     $           max(min(one,
     $           cld_frc_pcp_wgt(i)/pcp_flx_pst_sfc(i)
     $           )*pcp_flx_sfc(i)/pcp_flx_pst_sfc(i),
     $           cld_frc(i,k))
            pcp_flx_lcl=        ! [kg m-2 s-1] Local net precipitation production rate (+'ve when formation exceeds evaporation)
     $           prs_dlt(i,k)*grv_sfc_rcp*
     $           (q_H2O_cnd2pcp_tnd(i,k)+q_H2O_vpr2pcp_cnv_tnd(i,k)-q_H2O_pcp2vpr_tnd(i,k))
            pcp_flx_sfc(i)=pcp_flx_sfc(i)+pcp_flx_lcl ! [kg m-2 s-1] Precipitation rate

c     Positive definite version of pcp_flx_lcl used to weight contribution of layer cloud to precipitation-weighted cloud fraction
            pcp_flx_lcl_pst=max(pcp_flx_lcl,zero) ! [kg m-2 s-1] Positive definite local net precipitation production rate
            cld_frc_pcp_wgt(i)=cld_frc_pcp_wgt(i)+cld_frc(i,k)*pcp_flx_lcl_pst ! [kg m-2 s-1] Integral of precipitation-weighted cloud fraction
            pcp_flx_pst_sfc(i)=pcp_flx_pst_sfc(i)+pcp_flx_lcl_pst ! [kg m-2 s-1] Sum of positive definite precipitation production fluxes

            q_H2O_pcp_lqd(i,k)=0.0 ! [kg kg-1] Rainwater mixing ratio
            if (tpt_mdp(i,k) > tpt_frz_pnt) then
c     fxm: Should use virtual temperature in density calculation
               dns_mdp=prs_mdp(i,k)/(gas_cst_dry_air*tpt_mdp(i,k)) ! [kg m-3] Midlayer density
               vlc_grv_pcp=vlc_grv_fct/sqrt(dns_mdp) ! [m s-1] Fallspeed of precipitation
               q_H2O_pcp_lqd(i,k)=pcp_flx_sfc(i)/(dns_mdp*vlc_grv_pcp) ! [kg kg-1] Rainwater mixing ratio
               if (q_H2O_pcp_lqd(i,k) < q_H2O_pcp_lqd_min) q_H2O_pcp_lqd(i,k)=0.0 ! [kg kg-1] Rainwater mixing ratio
            endif               ! endif not freezing
         end do                 ! end loop over lon
      end do                    ! end loop over lev

      return
      end subroutine cld_dgn

      subroutine dpswet2nc(            
     $     cld_frc,             ! I [frc] Local total cloud fraction
     $     cld_frc_cnv,         ! I [frc] Convective cloud fraction
     $     cld_vlm,             ! I [frc] Cloud volume
     $     fl_out,              ! I [sng] Name of netCDF output file
     $     flx_mss_pcp_sfc,     ! I [kg m-2 s-1] Dust reaching surface in precipitation
     $     flx_mss_pcp_sfc_ttl, ! I [kg m-2 s-1] Total dust reaching surface in precipitation
     $     frc_trc_ist,         ! I [frc] Interstitial tracer fraction
     $     lat_idx,             ! I [idx] Model latitude index
     $     pcp_flx_sfc,         ! I [kg m-2 s-1] Total precipitation reaching surface
     $     q_H2O_cnd,           ! I [kg kg-1] Condensed H2O mixing ratio
     $     q_H2O_cnd2pcp_tnd,   ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
     $     q_H2O_pcp2vpr_tnd,   ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
     $     q_H2O_vpr2pcp_cnv_tnd, ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
     $     q_dst_tnd_evp,       ! I [kg kg-1 s-1] Evaporation tendency
     $     q_dst_tnd_evp_ttl,   ! I [kg kg-1 s-1] Total evaporation tendency
     $     q_dst_tnd_ncl,       ! I [kg kg-1 s-1] Nucleation scavenging tendency
     $     q_dst_tnd_pcp,       ! I [kg kg-1 s-1] Scavenging tendency
     $     q_dst_tnd_pcp_ttl,   ! I [kg kg-1 s-1] Total scavenging tendency
     $     q_dst_tnd_wet,       ! I [kg kg-1 s-1] Wet deposition (evaporation minus scavenging) tendency
     $     spc_xsx_cll_scv,     ! I [m2 kg-1] Specific cross section for collision scavenging
     $     spc_xsx_ncl_scv)     ! I [m2 kg-1] Specific cross section for nucleation scavenging

      use netcdf

c     Purpose: Output time-varying wet deposition properties to netCDF file
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstctl.h> /* Model control variables */
#include <dbg.com> /* Debugging constants */
#include <netcdf.inc>
c     Input
      character*(*),intent(in)::fl_out ! [sng] Name of netCDF output file
      integer,intent(in)::lat_idx ! I [idx] Model latitude index
      real(r8),intent(in)::cld_frc(plond,plev) ! [frc] Cloud fraction
      real(r8),intent(in)::cld_frc_cnv(plond,plev) ! [frc] Convective cloud fraction
      real(r8),intent(in)::cld_vlm(plond,plev) ! [frc] Cloud volume
      real(r8),intent(in)::q_H2O_vpr2pcp_cnv_tnd(plond,plev) ! [kg kg-1 s-1] H2O vapor to convective precipitation tendency
      real(r8),intent(in)::q_H2O_pcp2vpr_tnd(plond,plev) ! [kg kg-1 s-1] H2O precipitation to vapor tendency
      real(r8),intent(in)::q_H2O_cnd(plond,plev) ! [kg kg-1] Condensed H2O mixing ratio
      real(r8),intent(in)::q_H2O_cnd2pcp_tnd(plond,plev) ! [kg kg-1 s-1] Condensed H2O to precipitation tendency
      real(r8),intent(in)::spc_xsx_ncl_scv(dst_nbr) ! [m2 kg-1] Specific cross section for nucleation scavenging
      real(r8),intent(in)::spc_xsx_cll_scv(dst_nbr) ! [m2 kg-1] Specific cross section for collision scavenging
      real(r8),intent(in)::flx_mss_pcp_sfc(plond,dst_nbr) ! [kg m-2 s-1] Dust reaching surface in precipitation
      real(r8),intent(in)::flx_mss_pcp_sfc_ttl(plond) ! [kg m-2 s-1] Total dust reaching surface in precipitation
      real(r8),intent(in)::frc_trc_ist(plond,plev,dst_nbr) ! O [frc] Interstitial tracer fraction
      real(r8),intent(in)::pcp_flx_sfc(plond) ! [kg m-2 s-1] Total precipitation reaching surface
      real(r8),intent(in)::q_dst_tnd_evp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Evaporation tendency
      real(r8),intent(in)::q_dst_tnd_evp_ttl(plond,plev) ! [kg kg-1 s-1] Total evaporation tendency
      real(r8),intent(in)::q_dst_tnd_ncl(plond,plev,dst_nbr) ! O [kg kg-1 s-1] Nucleation scavenging tendency
      real(r8),intent(in)::q_dst_tnd_pcp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Scavenging tendency
      real(r8),intent(in)::q_dst_tnd_pcp_ttl(plond,plev) ! [kg kg-1 s-1] Total scavenging tendency
      real(r8),intent(in)::q_dst_tnd_wet(plond,plev,dst_nbr) ! [kg kg-1 s-1] Wet deposition (evaporation minus scavenging) tendency
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
      integer time_dim_id       ! [enm] Dimension ID for time
c     Variable IDs
      integer cld_frc_cnv_id    ! [enm] Variable ID
      integer cld_frc_id        ! [enm] Variable ID
      integer cld_vlm_id        ! [enm] Variable ID
      integer q_H2O_cnd_id      ! [enm] Variable ID
      integer q_H2O_cnd2pcp_tnd_id ! [enm] Variable ID
      integer q_H2O_pcp2vpr_tnd_id ! [enm] Variable ID
      integer q_H2O_vpr2pcp_cnv_tnd_id ! [enm] Variable ID
      integer spc_xsx_ncl_scv_id ! [enm] Variable ID
      integer spc_xsx_cll_scv_id ! [enm] Variable ID
      integer flx_mss_pcp_sfc_id ! [enm] Variable ID
      integer flx_mss_pcp_sfc_ttl_id ! [enm] Variable ID
      integer frc_trc_ist_id    ! [enm] Variable ID
      integer pcp_flx_sfc_id    ! [enm] Variable ID
      integer q_dst_tnd_evp_id  ! [enm] Variable ID
      integer q_dst_tnd_evp_ttl_id ! [enm] Variable ID
      integer q_dst_tnd_ncl_id  ! [enm] Variable ID
      integer q_dst_tnd_pcp_id  ! [enm] Variable ID
      integer q_dst_tnd_pcp_ttl_id ! [enm] Variable ID
      integer q_dst_tnd_wet_id  ! [enm] Variable ID
c     Variable data
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
      
c     Begin netCDF output routines
      rcd=NF_NOERR              ! NF_NOERR == 0
      rcd=rcd+nf90_open(fl_out,NF_WRITE,nc_id)
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR opening ',fl_out(1:f_strlen(fl_out)),
     $     ' in dpswet2nc(), rcd = ',rcd
      rcd=rcd+nf90_redef(nc_id)
      rcd=rcd+nf90_inq_dimid(nc_id,'lev',lev_dim_id)
      rcd=rcd+nf90_inq_dimid(nc_id,'lon',lon_dim_id)
      rcd=rcd+nf90_inq_dimid(nc_id,'time',time_dim_id)
      rcd=rcd+nf90_inq_dimid(nc_id,'sz',sz_dim_id)
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
         rcd=rcd+nf90_def_var(nc_id,'cld_frc',NF_DOUBLE,dim_lon_lev_time,cld_frc_id)
         rcd=rcd+nf90_def_var(nc_id,'cld_frc_cnv',NF_DOUBLE,dim_lon_lev_time,cld_frc_cnv_id)
         rcd=rcd+nf90_def_var(nc_id,'cld_vlm',NF_DOUBLE,dim_lon_lev_time,cld_vlm_id)
         rcd=rcd+nf90_def_var(nc_id,'flx_mss_pcp_sfc',NF_DOUBLE,dim_lon_sz_time,flx_mss_pcp_sfc_id)
         rcd=rcd+nf90_def_var(nc_id,'flx_mss_pcp_sfc_ttl',NF_DOUBLE,dim_lon_time,flx_mss_pcp_sfc_ttl_id)
         rcd=rcd+nf90_def_var(nc_id,'frc_trc_ist',NF_DOUBLE,dim_lon_lev_sz_time,frc_trc_ist_id)
         rcd=rcd+nf90_def_var(nc_id,'pcp_flx_sfc',NF_DOUBLE,dim_lon_time,pcp_flx_sfc_id)
         rcd=rcd+nf90_def_var(nc_id,'q_H2O_cnd',NF_DOUBLE,dim_lon_lev_time,q_H2O_cnd_id)
         rcd=rcd+nf90_def_var(nc_id,'q_H2O_cnd2pcp_tnd',NF_DOUBLE,dim_lon_lev_time,q_H2O_cnd2pcp_tnd_id)
         rcd=rcd+nf90_def_var(nc_id,'q_H2O_pcp2vpr_tnd',NF_DOUBLE,dim_lon_lev_time,q_H2O_pcp2vpr_tnd_id)
         rcd=rcd+nf90_def_var(nc_id,'q_H2O_vpr2pcp_cnv_tnd',NF_DOUBLE,dim_lon_lev_time,q_H2O_vpr2pcp_cnv_tnd_id)
         rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_evp',NF_DOUBLE,dim_lon_lev_sz_time,q_dst_tnd_evp_id)
         rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_evp_ttl',NF_DOUBLE,dim_lon_lev_time,q_dst_tnd_evp_ttl_id)
         rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_ncl',NF_DOUBLE,dim_lon_lev_sz_time,q_dst_tnd_ncl_id)
         rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_pcp',NF_DOUBLE,dim_lon_lev_sz_time,q_dst_tnd_pcp_id)
         rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_pcp_ttl',NF_DOUBLE,dim_lon_lev_time,q_dst_tnd_pcp_ttl_id)
         rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_wet',NF_DOUBLE,dim_lon_lev_sz_time,q_dst_tnd_wet_id)
         rcd=rcd+nf90_def_var(nc_id,'spc_xsx_cll_scv',NF_DOUBLE,sz_dim_id,spc_xsx_cll_scv_id)
         rcd=rcd+nf90_def_var(nc_id,'spc_xsx_ncl_scv',NF_DOUBLE,sz_dim_id,spc_xsx_ncl_scv_id)
c     Add english text descriptions
         rcd=rcd+nf90_put_att(nc_id,cld_frc_cnv_id,'long_name','Convective cloud fraction')
         rcd=rcd+nf90_put_att(nc_id,cld_frc_id,'long_name','Cloud fraction')
         rcd=rcd+nf90_put_att(nc_id,cld_vlm_id,'long_name','Cloud volume')
         rcd=rcd+nf90_put_att(nc_id,flx_mss_pcp_sfc_id,'long_name','Dust reaching surface in precipitation')
         rcd=rcd+nf90_put_att(nc_id,flx_mss_pcp_sfc_ttl_id,'long_name','Total dust reaching surface in precipitation')
         rcd=rcd+nf90_put_att(nc_id,frc_trc_ist_id,'long_name','Interstitial tracer fraction')
         rcd=rcd+nf90_put_att(nc_id,pcp_flx_sfc_id,'long_name','Total precipitation reaching surface')
         rcd=rcd+nf90_put_att(nc_id,q_H2O_cnd2pcp_tnd_id,'long_name','Condensed H2O to precipitation tendency')
         rcd=rcd+nf90_put_att(nc_id,q_H2O_cnd_id,'long_name','Condensed H2O mixing ratio')
         rcd=rcd+nf90_put_att(nc_id,q_H2O_pcp2vpr_tnd_id,'long_name','H2O precipitation to vapor tendency')
         rcd=rcd+nf90_put_att(nc_id,q_H2O_vpr2pcp_cnv_tnd_id,'long_name','H2O vapor to convective precipitation tendency')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_evp_id,'long_name','Evaporation tendency')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_evp_ttl_id,'long_name','Total evaporation tendency')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_ncl_id,'long_name','Nucleation scavenging tendency')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_pcp_id,'long_name','Scavenging tendency')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_pcp_ttl_id,'long_name','Total scavenging tendency')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_wet_id,'long_name','Wet deposition (evaporation minus scavenging) tendency')
         rcd=rcd+nf90_put_att(nc_id,spc_xsx_cll_scv_id,'long_name','Specific cross section for collision scavenging')
         rcd=rcd+nf90_put_att(nc_id,spc_xsx_ncl_scv_id,'long_name','Specific cross section for nucleation scavenging')
c     Add units
         rcd=rcd+nf90_put_att(nc_id,cld_frc_cnv_id,'units','fraction')
         rcd=rcd+nf90_put_att(nc_id,cld_frc_id,'units','fraction')
         rcd=rcd+nf90_put_att(nc_id,cld_vlm_id,'units','fraction')
         rcd=rcd+nf90_put_att(nc_id,flx_mss_pcp_sfc_id,'units','kilogram meter-2 second-1')
         rcd=rcd+nf90_put_att(nc_id,flx_mss_pcp_sfc_ttl_id,'units','kilogram meter-2 second-1')
         rcd=rcd+nf90_put_att(nc_id,frc_trc_ist_id,'units','fraction')
         rcd=rcd+nf90_put_att(nc_id,pcp_flx_sfc_id,'units','kilogram meter-2 second-1')
         rcd=rcd+nf90_put_att(nc_id,q_H2O_cnd2pcp_tnd_id,'units','kilogram kilogram-1 second-1')
         rcd=rcd+nf90_put_att(nc_id,q_H2O_cnd_id,'units','kilogram kilogram-1')
         rcd=rcd+nf90_put_att(nc_id,q_H2O_pcp2vpr_tnd_id,'units','kilogram kilogram-1 second-1')
         rcd=rcd+nf90_put_att(nc_id,q_H2O_vpr2pcp_cnv_tnd_id,'units','kilogram kilogram-1 second-1')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_evp_id,'units','kilogram kilogram-1 second-1')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_evp_ttl_id,'units','kilogram kilogram-1 second-1')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_ncl_id,'units','kilogram kilogram-1 second-1')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_pcp_id,'units','kilogram kilogram-1 second-1')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_pcp_ttl_id,'units','kilogram kilogram-1 second-1')
         rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_wet_id,'units','kilogram kilogram-1 second-1')
         rcd=rcd+nf90_put_att(nc_id,spc_xsx_cll_scv_id,'units','meter2 kilogram-1')
         rcd=rcd+nf90_put_att(nc_id,spc_xsx_ncl_scv_id,'units','meter2 kilogram-1')
      else                      ! endif nstep == 1
c     Get variable IDs
         rcd=rcd+nf90_inq_varid(nc_id,'cld_frc',cld_frc_id)
         rcd=rcd+nf90_inq_varid(nc_id,'cld_frc_cnv',cld_frc_cnv_id)
         rcd=rcd+nf90_inq_varid(nc_id,'cld_vlm',cld_vlm_id)
         rcd=rcd+nf90_inq_varid(nc_id,'flx_mss_pcp_sfc',flx_mss_pcp_sfc_id)
         rcd=rcd+nf90_inq_varid(nc_id,'flx_mss_pcp_sfc_ttl',flx_mss_pcp_sfc_ttl_id)
         rcd=rcd+nf90_inq_varid(nc_id,'frc_trc_ist',frc_trc_ist_id)
         rcd=rcd+nf90_inq_varid(nc_id,'pcp_flx_sfc',pcp_flx_sfc_id)
         rcd=rcd+nf90_inq_varid(nc_id,'q_H2O_cnd',q_H2O_cnd_id)
         rcd=rcd+nf90_inq_varid(nc_id,'q_H2O_cnd2pcp_tnd',q_H2O_cnd2pcp_tnd_id)
         rcd=rcd+nf90_inq_varid(nc_id,'q_H2O_pcp2vpr_tnd',q_H2O_pcp2vpr_tnd_id)
         rcd=rcd+nf90_inq_varid(nc_id,'q_H2O_vpr2pcp_cnv_tnd',q_H2O_vpr2pcp_cnv_tnd_id)
         rcd=rcd+nf90_inq_varid(nc_id,'q_dst_tnd_evp',q_dst_tnd_evp_id)
         rcd=rcd+nf90_inq_varid(nc_id,'q_dst_tnd_evp_ttl',q_dst_tnd_evp_ttl_id)
         rcd=rcd+nf90_inq_varid(nc_id,'q_dst_tnd_ncl',q_dst_tnd_ncl_id)
         rcd=rcd+nf90_inq_varid(nc_id,'q_dst_tnd_pcp',q_dst_tnd_pcp_id)
         rcd=rcd+nf90_inq_varid(nc_id,'q_dst_tnd_pcp_ttl',q_dst_tnd_pcp_ttl_id)
         rcd=rcd+nf90_inq_varid(nc_id,'q_dst_tnd_wet',q_dst_tnd_wet_id)
      endif                     ! endif nstep /= 1
c     End define mode
      rcd=rcd+nf90_enddef(nc_id)
c     Write data
      if (nstep == 1) then
         rcd=rcd+nf90_put_var(nc_id,spc_xsx_ncl_scv_id,spc_xsx_ncl_scv)
         rcd=rcd+nf90_put_var(nc_id,spc_xsx_cll_scv_id,spc_xsx_cll_scv)
      endif                     ! endif nstep == 1
      rcd=rcd+nf90_put_var(nc_id,cld_frc_cnv_id,cld_frc_cnv,srt_lon_lev_time,cnt_lon_lev_time)
      rcd=rcd+nf90_put_var(nc_id,cld_frc_id,cld_frc,srt_lon_lev_time,cnt_lon_lev_time)
      rcd=rcd+nf90_put_var(nc_id,cld_vlm_id,cld_vlm,srt_lon_lev_time,cnt_lon_lev_time)
      rcd=rcd+nf90_put_var(nc_id,flx_mss_pcp_sfc_id,flx_mss_pcp_sfc,srt_lon_sz_time,cnt_lon_sz_time)
      rcd=rcd+nf90_put_var(nc_id,flx_mss_pcp_sfc_ttl_id,flx_mss_pcp_sfc_ttl,srt_lon_time,cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id,frc_trc_ist_id,frc_trc_ist,srt_lon_lev_sz_time,cnt_lon_lev_sz_time)
      rcd=rcd+nf90_put_var(nc_id,pcp_flx_sfc_id,pcp_flx_sfc,srt_lon_time,cnt_lon_time)
      rcd=rcd+nf90_put_var(nc_id,q_H2O_cnd2pcp_tnd_id,q_H2O_cnd2pcp_tnd,srt_lon_lev_time,cnt_lon_lev_time)
      rcd=rcd+nf90_put_var(nc_id,q_H2O_cnd_id,q_H2O_cnd,srt_lon_lev_time,cnt_lon_lev_time)
      rcd=rcd+nf90_put_var(nc_id,q_H2O_pcp2vpr_tnd_id,q_H2O_pcp2vpr_tnd,srt_lon_lev_time,cnt_lon_lev_time)
      rcd=rcd+nf90_put_var(nc_id,q_H2O_vpr2pcp_cnv_tnd_id,q_H2O_vpr2pcp_cnv_tnd,srt_lon_lev_time,cnt_lon_lev_time)
      rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_evp_id,q_dst_tnd_evp,srt_lon_lev_sz_time,cnt_lon_lev_sz_time)
      rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_evp_ttl_id,q_dst_tnd_evp_ttl,srt_lon_lev_time,cnt_lon_lev_time)
      rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_ncl_id,q_dst_tnd_ncl,srt_lon_lev_sz_time,cnt_lon_lev_sz_time)
      rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_pcp_id,q_dst_tnd_pcp,srt_lon_lev_sz_time,cnt_lon_lev_sz_time)
      rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_pcp_ttl_id,q_dst_tnd_pcp_ttl,srt_lon_lev_time,cnt_lon_lev_time)
      rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_wet_id,q_dst_tnd_wet,srt_lon_lev_sz_time,cnt_lon_lev_sz_time)
c     Close output file
      rcd=rcd+nf90_close(nc_id)
      if (nstep == 1) then
         write (6,'(a,a45,a)') prg_nm(1:f_strlen(prg_nm)),': Initialized wet deposition data archive in ',fl_out(1:f_strlen(fl_out))
      endif                     ! endif
      if (rcd /= NF_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR writing ',fl_out(1:f_strlen(fl_out)),
     $     ' in dpswet2nc(), rcd = ',rcd
      
      return 
      end subroutine dpswet2nc                      ! end dpswet2nc()
      end module dstdpswet
