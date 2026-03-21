c     $Header$
      
c     Purpose: dpsdryutl.F contains routines to compute 
c     dry deposition tendencies
      
      subroutine rss_aer_get(
     $     hgt_mdp,             ! I [m] Midlayer height above surface
     $     hgt_zpd,             ! I [m] Zero plane displacement height
     $     mno_lng,             ! I [m] Monin-Obukhov length
     $     rgh_mmn,             ! I [m] Roughness length momentum
     $     rss_aer,             ! O [s m-1] Aerodynamic resistance
     $     wnd_frc)             ! I [m s-1] Surface friction velocity

c     Purpose: Given the surface layer structure and dynamics properties,
c     compute and return the aerodynamic resistance

      use precision

      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstcst.h> /* Physical constants for dust routines */
c     Input
      real(r8) mno_lng(plond)       ! [m] Monin-Obukhov length
      real(r8) rgh_mmn(plond)       ! [m] Roughness length momentum
      real(r8) wnd_frc(plond)       ! [m s-1] Surface friction velocity from pbldif()
      real(r8) hgt_mdp(plond)       ! [m] Midlayer height above surface
      real(r8) hgt_zpd(plond)       ! [m] Zero plane displacement
c     Output
      real(r8) rss_aer(plond)       ! [s m-1] Aerodynamic resistance
c     Local
      integer lon_idx           ! [idx] Counting index
      real(r8) eta_sqr_rlm          ! Eta squared at roughness height
      real(r8) eta_sqr_gcm          ! Eta squared at GCM layer height
      real(r8) eta_rlm              ! Eta at roughness height
      real(r8) eta_gcm              ! Eta at GCM layer height
      real(r8) nmr_rlm              ! Numerator
      real(r8) dnm_gcm              ! Denominator
      real(r8) tmp4                 ! Correction to neutral atmosphere factor
      real(r8) tmp5                 ! Neutral atmosphere factor
      real(r8) mno_prm_rlm          ! [frc] Monin-Obukhov parameter at roughness height
      real(r8) mno_prm_gcm          ! [frc] Monin-Obukhov parameter at GCM layer height
c     Main Code
      
c     Sanity check
      do lon_idx=1,plon
         if (rgh_mmn(lon_idx) <= 0.0.or.hgt_mdp(lon_idx) <= hgt_zpd(lon_idx)) then
            write(6,'(a,3(a,1p,e8.1,0p),a)') 'rss_aer_get: ',
     $           'rgh_mmn = ',rgh_mmn(lon_idx),', hgt_mdp = ',hgt_mdp(lon_idx),' m, hgt_zpd = ',hgt_zpd(lon_idx),' m s-1'
         endif                  ! endif err
      enddo                     ! end loop over lon

c     Compute stability parameter
      do lon_idx=1,plon
c     Maximum value is 1 because stability correction function is valid only for zeta < 1, e.g., Bon96 p. 52, Bru82 p. 71, SeP97 p. 963
         mno_prm_rlm=min(rgh_mmn(lon_idx)/mno_lng(lon_idx),1.0_r8) ! [frc]
         mno_prm_gcm=min((hgt_mdp(lon_idx)-hgt_zpd(lon_idx))/mno_lng(lon_idx),1.0_r8) ! [frc]
         if (mno_lng(lon_idx) < 0.0) then
c     Difference between unstable corrections
            eta_sqr_rlm=sqrt(1.0-16.0*mno_prm_rlm)
            eta_sqr_gcm=sqrt(1.0-16.0*mno_prm_gcm)
            eta_rlm=sqrt(eta_sqr_rlm)
            eta_gcm=sqrt(eta_sqr_gcm)
            nmr_rlm=(eta_sqr_rlm+1.0)*(eta_rlm+1.0)*(eta_rlm+1.0)
            dnm_gcm=(eta_sqr_gcm+1.0)*(eta_gcm+1.0)*(eta_gcm+1.0)
            tmp4=log(nmr_rlm/dnm_gcm)+2.0*(atan(eta_gcm)-atan(eta_rlm)) ! [frc]
         else                   ! not stable
c     Difference between stable corrections
            tmp4=5.0*(mno_prm_gcm-mno_prm_rlm) ! [frc]
         endif                  ! not stable
         tmp5=log((hgt_mdp(lon_idx)-hgt_zpd(lon_idx))/rgh_mmn(lon_idx))
         rss_aer(lon_idx)=(tmp4+tmp5)/(cst_von_krm*wnd_frc(lon_idx)) ! [s m-1] Bon96 p. 54

      enddo                     ! end loop over lon
      return
      end                       ! end rss_aer_get()
      
      subroutine stk_crc_get(
     $     dmt_prt,             ! I [m] Particle diameter
     $     dns_prt,             ! I [kg m-3] Particle density
     $     stk_crc,             ! O [frc] Correction to Stokes settling velocity
     $     sz_nbr)              ! I [nbr] Number of particle sizes
c     Purpose: Compute factor which corrects Stokes settling velocity to actual
c     settling velocity when Reynolds' number is greater than 0.1
c     stk_crc_get() is called by dst_szdst_ini()
c     This implementation computes a time-invariant, size-dependent ratio which
c     corrects the solution of the (linear) Stokes settling velocity to the full
c     (non-linear, iterative solution) for the settling velocity.

      use precision

      implicit none
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
c     Parameters
c     Apparently Boltzmann's constant is not set anywhere by CCM or MATCH
      real(r8),parameter::cst_Boltzmann=1.38063e-23 ! [J K-1] Boltzmann's constant
      real(r8),parameter::prs_mdp=100000.0 ! [Pa] Pressure 
      real(r8),parameter::tpt_mdp=295.0 ! [K] Temperature
      real(r8),parameter::tpt_vrt=295.0 ! [K] Virtual temperature
c     Input
      integer,intent(in)::sz_nbr ! I [nbr] Number of particle sizes
      real(r8),intent(in)::dmt_prt(sz_nbr) ! I [m] Particle diameter
      real(r8),intent(in)::dns_prt(sz_nbr) ! I [kg m-3] Particle density
c     Output
      real(r8),intent(out)::stk_crc(sz_nbr) ! O [frc] Correction to Stokes settling velocity
c     Local
      double precision pi       ! [frc] 3
      integer itr_idx           ! [idx] Counting index
      integer sz_idx            ! [idx] Counting index
      real(r8) cff_drg_grv(sz_nbr)  ! [frc] Drag coefficient at terminal velocity
      real(r8) dns_mdp              ! [kg m-3] Midlayer density
      real(r8) eps_crr              ! [frc] Current relative accuracy
      real(r8) eps_max              ! [frc] Relative accuracy for convergence
      real(r8) mfp_atm              ! [m] Mean free path of air
      real(r8) ryn_nbr_grv(sz_nbr)  ! [frc] Reynolds number at terminal velocity
      real(r8) slp_crc(sz_nbr)      ! [frc] Slip correction factor
      real(r8) vlc_grv(sz_nbr)      ! [m s-1] Settling velocity
      real(r8) vlc_grv_old          ! [m s-1] Previous gravitational settling velocity
      real(r8) vlc_stk(sz_nbr)      ! [m s-1] Stokes settling velocity
      real(r8) vsc_dyn_atm          ! [kg m-1 s-1] Dynamic viscosity of air
      real(r8) vsc_knm_atm          ! [m2 s-1] Kinematic viscosity of atmosphere 
c     Initialize some variables
      pi=4.0*atan(1.0d0)        ! [frc] 3
      eps_max=1.0e-4            ! [frc] Relative accuracy for convergence
      dns_mdp=prs_mdp/(tpt_vrt*gas_cst_dry_air) ! [kg m-3]
c     Size-independent thermokinetic properties
      vsc_dyn_atm=1.72e-5*((tpt_mdp/273.0)**1.5)*393.0/(tpt_mdp+120.0) ! [kg m-1 s-1] RoY94  p. 102
      mfp_atm=2.0*vsc_dyn_atm/(prs_mdp*sqrt(8.0*mmw_dry_air/(pi*gas_cst_unv*tpt_mdp))) ! [m] Mean free path of air SeP97 p. 455
      vsc_knm_atm=vsc_dyn_atm/dns_mdp ! [m2 s-1] Kinematic viscosity of air
      
      do sz_idx=1,sz_nbr
         slp_crc(sz_idx)=1.0+2.0*mfp_atm*(1.257+0.4*exp(-1.1*dmt_prt(sz_idx)/(2.0*mfp_atm)))/dmt_prt(sz_idx) ! [frc] Slip correction factor SeP97 p. 464
         vlc_stk(sz_idx)=(1.0/18.0)*dmt_prt(sz_idx)*dmt_prt(sz_idx)*dns_prt(sz_idx)*grv_sfc*slp_crc(sz_idx)/vsc_dyn_atm ! [m s-1] Stokes settling velocity SeP97 p. 466
      enddo                     ! end loop over size
      
c     For Reynolds number flows Re < 0.1 Stokes' velocity is valid for vlc_grv SeP97 p. 466 (8.42)
c     For larger Re, inertial effects become important and empirical drag coefficients must be employed
c     Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
c     Using Stokes' velocity rather than iterative solution with empirical drag coefficient causes 60% errors for D = 200 um SeP97 p. 468
      
c     Iterative solution for drag coefficient, Reynolds number, and terminal velocity
      do sz_idx=1,sz_nbr
c     Initialize accuracy and counter
         eps_crr=eps_max+1.0    ! [frc] Current relative accuracy
         itr_idx=0              ! [idx] Counting index
c     Initial guess for vlc_grv is exact for Re < 0.1
         vlc_grv(sz_idx)=vlc_stk(sz_idx)  ! [m s-1]
         do while(eps_crr > eps_max)
c     Save terminal velocity for convergence test
            vlc_grv_old=vlc_grv(sz_idx) ! [m s-1] 
            ryn_nbr_grv(sz_idx)=vlc_grv(sz_idx)*dmt_prt(sz_idx)/vsc_knm_atm ! [frc] SeP97 p. 460
c     Update drag coefficient based on new Reynolds number
            if (ryn_nbr_grv(sz_idx) < 0.1) then
               cff_drg_grv(sz_idx)=24.0/ryn_nbr_grv(sz_idx) ! Stokes' law Sep97 p. 463 (8.32)
            else if (ryn_nbr_grv(sz_idx) < 2.0) then
               cff_drg_grv(sz_idx)=(24.0/ryn_nbr_grv(sz_idx))*
     $              (1.0+3.0*ryn_nbr_grv(sz_idx)/16.0+
     $              9.0*ryn_nbr_grv(sz_idx)*ryn_nbr_grv(sz_idx)*log(2.0*ryn_nbr_grv(sz_idx))/160.0) ! Sep97 p. 463 (8.32)
            else if (ryn_nbr_grv(sz_idx) < 500.0) then
               cff_drg_grv(sz_idx)=(24.0/ryn_nbr_grv(sz_idx))*(1.0+0.15*ryn_nbr_grv(sz_idx)**0.687) ! Sep97 p. 463 (8.32)
            else if (ryn_nbr_grv(sz_idx) < 2.0e5) then
               cff_drg_grv(sz_idx)=0.44 ! Sep97 p. 463 (8.32)
            else
               stop 'Reynolds number too large in stk_crc_get()'
            endif               ! end else
c     Update terminal velocity based on new Reynolds number and drag coefficient
            vlc_grv(sz_idx)=sqrt(4.0*grv_sfc*dmt_prt(sz_idx)*slp_crc(sz_idx)*dns_prt(sz_idx)/(3.0*cff_drg_grv(sz_idx)*dns_mdp)) ! [m s-1] Terminal velocity SeP97 p. 467 (8.44)
            eps_crr=abs((vlc_grv(sz_idx)-vlc_grv_old)/vlc_grv(sz_idx)) ! Relative convergence
            if (itr_idx == 12) then
c     Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0 due to discontinuities in derivative of drag coefficient
               vlc_grv(sz_idx)=0.5*(vlc_grv(sz_idx)+vlc_grv_old) ! [m s-1]
            endif               ! endif
            if (itr_idx > 20) then
               write (6,'(a)') 'dst: Terminal velocity not converging in stk_crc_get(), breaking loop...'
               goto 100         ! Jump to next iteration
            endif               ! endif
            itr_idx=itr_idx+1
         end do                 ! end while
 100     continue               ! Label to jump to when iteratio does not converge
      enddo                     ! end loop over size
      
c     Compute factors to convert Stokes' settling velocities to actual settling velocities
      do sz_idx=1,sz_nbr
         stk_crc(sz_idx)=vlc_grv(sz_idx)/vlc_stk(sz_idx) ! [frc] Correction to Stokes settling velocity
      enddo                     ! end loop over size
      return
      end                       ! end stk_crc_get()
      
      subroutine vlc_grv_get(
     $     dmt_prt,             ! I [m] Particle diameter
     $     dns_prt,             ! I [kg m-3] Particle density
     $     prs_mdp,             ! I [Pa] Pressure
     $     stk_crc,             ! I [frc] Correction to Stokes settling velocity
     $     tpt_mdp,             ! I [K] Temperature
     $     vlc_grv)             ! O [m s-1] Settling velocity
c     Purpose: Given size and density of particle, and column thermodynamic profile,
c     compute terminal fall speed

      use precision

      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
c     Parameters
c     Apparently Boltzmann's constant is not set anywhere by CCM or MATCH
      real(r8) cst_Boltzmann        ! [J K-1] Boltzmann's constant
      parameter(cst_Boltzmann=1.38063e-23) ! [J K-1] Boltzmann's constant
c     Input
      real(r8) dns_prt(dst_nbr)     ! I [kg m-3] Particle density
      real(r8) dmt_prt(dst_nbr)     ! I [m] Particle diameter
      real(r8) prs_mdp(plond,plev)  ! I [Pa] Pressure 
      real(r8) tpt_mdp(plond,plev)  ! I [K] Temperature
      real(r8) stk_crc(dst_nbr)     ! I [frc] Correction to Stokes settling velocity
c     Output
      real(r8) vlc_grv(plond,plev,dst_nbr) ! O [m s-1] Settling velocity
c     Local
      double precision pi       ! [frc] 3
      integer i                 ! [idx] Counting index
      integer k                 ! [idx] Counting index
      integer m                 ! [idx] Counting index
      real(r8) mfp_atm(plond,plev)  ! [m] Mean free path of air
      real(r8) slp_crc              ! [frc] Slip correction factor
      real(r8) vsc_dyn_atm(plond,plev) ! [kg m-1 s-1] Dynamic viscosity of air
c     Main Code
c     Initialize some variables
      pi=4.0*atan(1.0d0)        ! [frc] 3
      

c     Size-independent thermokinetic properties
      do k=1,plev
         do i=1,plon
            vsc_dyn_atm(i,k)=1.72e-5*((tpt_mdp(i,k)/273.0)**1.5)*393.0/(tpt_mdp(i,k)+120.0) ! [kg m-1 s-1] RoY94  p. 102
            mfp_atm(i,k)=2.0*vsc_dyn_atm(i,k)/(prs_mdp(i,k)*sqrt(8.0*mmw_dry_air/(pi*gas_cst_unv*tpt_mdp(i,k)))) ! [m] Mean free path of air SeP97 p. 455
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      
      do m=1,dst_nbr
         do k=1,plev
            do i=1,plon
               slp_crc=1.0+2.0*mfp_atm(i,k)*(1.257+0.4*exp(-1.1*dmt_prt(m)/(2.0*mfp_atm(i,k))))/dmt_prt(m) ! [frc] Slip correction factor SeP97 p. 464
               vlc_grv(i,k,m)=(1.0/18.0)*dmt_prt(m)*dmt_prt(m)*dns_prt(m)*grv_sfc*slp_crc/vsc_dyn_atm(i,k) ! [m s-1] Stokes settling velocity SeP97 p. 466
               vlc_grv(i,k,m)=vlc_grv(i,k,m)*stk_crc(m) ! [m s-1] Corrected Stokes settling velocity
               if (vlc_grv(i,k,m).lt.0.0 ) then
                  write(6,*)'ABORTING: vlc_grv ',vlc_grv(i,k,m),' < 0 in vlc_grv_get at i ',i,' k ',k,' m ',m
                  call abort()
               end if
            end do              ! end loop over lon
         end do                 ! end loop over lev
      enddo                     ! end loop over size

      return
      end                       ! end vlc_grv_get()

      subroutine rss_lmn_get(
     $     dmt_aer,             ! I [m] Particle diameter
     $     dns_aer,             ! I [kg m-3] Particle density
     $     dns_mdp,             ! I [kg m-3] Midlayer density
     $     prs_mdp,             ! I [Pa] Pressure
     $     rss_lmn,             ! O [s m-1] Quasi-laminar layer resistance
     $     shm_nbr,             ! O [frc] Schmidt number
     $     stk_crc,             ! I [frc] Correction to Stokes settling velocity
     $     stk_nbr,             ! O [frc] Stokes number
     $     tpt_mdp,             ! I [K] Temperature
     $     wnd_frc)             ! I [m s-1] Friction velocity

c     Purpose: Given size and density of aerosol, and surface environmental properties, compute quasi-laminar layer resistance
c     rss_lmn_get() describes physics occuring within a few meters (at most) of the zero plane displacement height
c     Thus rss_lmn_get() requires the aerodynamic pressure and temperature which is closer to lowest interface than to lowest midlayer
c     Currently we ignore this subtlety and use lowest midlayer quantities

      use precision

      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
c     Parameters
c     Apparently Boltzmann's constant is not set anywhere by CCM or MATCH
      real(r8) cst_Boltzmann        ! [J K-1] Boltzmann's constant
      parameter(cst_Boltzmann=1.38063e-23) ! [J K-1] Boltzmann's constant
c     Input
      real(r8) dns_aer(dst_nbr)     ! I [kg m-3] Particle density
      real(r8) dmt_aer(dst_nbr)     ! I [m] Particle diameter
      real(r8) stk_crc(dst_nbr)     ! I [frc] Correction to Stokes settling velocity
      real(r8) dns_mdp(plond)       ! I [kg m-3] Midlayer density
      real(r8) prs_mdp(plond)       ! I [Pa] Pressure 
      real(r8) tpt_mdp(plond)       ! I [K] Temperature
      real(r8) wnd_frc(plond)       ! I [m s-1] Friction velocity
c     Output
      real(r8) rss_lmn(plond,dst_nbr) ! O [s m-1] Quasi-laminar layer resistance
      real(r8) stk_nbr(plond,dst_nbr) ! O [frc] Stokes number
      real(r8) shm_nbr(plond,dst_nbr) ! O [frc] Schmidt number
c     Local
      double precision pi       ! [frc] 3
      integer i                 ! [idx] Counting index
      integer k                 ! [idx] Counting index
      integer m                 ! [idx] Counting index
      real(r8) dff_aer(plond,dst_nbr) ! [m2 s-1] Brownian diffusivity of particle
      real(r8) mfp_atm(plond)       ! [m] Mean free path of air
      real(r8) slp_crc(plond,dst_nbr) ! [frc] Slip correction factor
      real(r8) tmp                  ! [frc] Factor in rss_lmn computation
      real(r8) vlc_grv(plond,dst_nbr) ! [m s-1] Settling velocity
      real(r8) vsc_dyn_atm(plond)   ! [kg m-1 s-1] Dynamic viscosity of air
      real(r8) vsc_knm_atm(plond)   ! [m2 s-1] Kinematic viscosity of air

c     Main Code
c     Initialize some variables

      pi=4.0*atan(1.0d0)        ! [frc] 3
c     Size-independent thermokinetic properties
      do i=1,plon
         vsc_dyn_atm(i)=1.72e-5*((tpt_mdp(i)/273.0)**1.5)*393.0/(tpt_mdp(i)+120.0) ! [kg m-1 s-1] RoY94  p. 102
         mfp_atm(i)=2.0*vsc_dyn_atm(i)/(prs_mdp(i)*sqrt(8.0*mmw_dry_air/(pi*gas_cst_unv*tpt_mdp(i)))) ! [m] Mean free path of air SeP97 p. 455
         vsc_knm_atm(i)=vsc_dyn_atm(i)/dns_mdp(i) ! [m2 s-1] Kinematic viscosity of air
      end do                    ! end loop over lon

      do m=1,dst_nbr
         do i=1,plon
            slp_crc(i,m)=1.0+2.0*mfp_atm(i)*(1.257+0.4*exp(-1.1*dmt_aer(m)/(2.0*mfp_atm(i))))/dmt_aer(m) ! [frc] Slip correction factor SeP97 p. 464
            vlc_grv(i,m)=(1.0/18.0)*dmt_aer(m)*dmt_aer(m)*dns_aer(m)*grv_sfc*slp_crc(i,m)/vsc_dyn_atm(i) ! [m s-1] Stokes' settling velocity SeP97 p. 466
            vlc_grv(i,m)=vlc_grv(i,m)*stk_crc(m) ! [m s-1] Correction to Stokes settling velocity
         end do                 ! end loop over lon
      enddo                     ! end loop over size

      do m=1,dst_nbr
         do i=1,plon
            stk_nbr(i,m)=vlc_grv(i,m)*wnd_frc(i)*wnd_frc(i)/(grv_sfc*vsc_knm_atm(i)) ! [frc] Stokes number SeP97 p. 965
            dff_aer(i,m)=cst_Boltzmann*tpt_mdp(i)*slp_crc(i,m)/(3.0*pi*vsc_dyn_atm(i)*dmt_aer(m)) ! [m2 s-1] Brownian diffusivity of particle SeP97 p. 474
            shm_nbr(i,m)=vsc_knm_atm(i)/dff_aer(i,m) ! [frc] Schmidt number SeP97 p. 972
            tmp=shm_nbr(i,m)**(-0.6666667)+10.0**(-3.0/stk_nbr(i,m))
            rss_lmn(i,m)=1.0/(tmp*wnd_frc(i)) ! [s m-1] SeP97 p. 972, 965
         end do                 ! end loop over lon
      enddo                     ! end loop over size

      return
      end                       ! end rss_lmn_get()
      

