      module dstszdst

      use precision

      contains

c     $Header$
      
c     Purpose: dstszdst.F contains routines to initialize mineral dust size distributions
      
      subroutine dst_szdst_ini()
c     Purpose: Get dust size grid from data file and set other
c     important statistics of size distributions
c     dst_szdst_ini() is called by CCM:physics/inti(), MATCH:inirun()
      implicit none
c     Parameters
      integer,parameter::grd_typ_lgr=0 ! [enm] Logarithmic grid
      integer,parameter::grd_typ_lnr=1 ! [enm] Linear grid
c     Size range to consider--computed bins in dmt_grd will fill this size range
c     Range set here has no effect if preset dmt_grd configuration is used
!pjr cant do this cause grd_mk may reset it
!pjr      real(r8),parameter::dmt_min_min=0.1e-6; ! [m] Minimum grid diameter
!pjr      real(r8),parameter::dmt_max_max=100.0e-6; ! [m] Maximum grid diameter
      real(r8) ::dmt_min_min=0.1e-6; ! [m] Minimum grid diameter
      real(r8) ::dmt_max_max=100.0e-6; ! [m] Maximum grid diameter
c     Preset boundaries for common bin configurations
c     Values are used when a preset configurations exists for dst_nbr bins
c     Bin boundaries at 1.0, 2.5, and 10.0 microns are nice for comparison to observations
      real(r8),dimension(5),parameter::dmt_grd_04=(/ 0.1e-6 ,  1.0e-6,  2.5e-6,  5.0e-6, 10.0e-6 /) ! [m] Particle diameter grid
c     real(r8),dimension(5),parameter::dmt_grd_04=(/ 0.1e-6 ,  0.5e-6,  1.0e-6,  5.0e-6, 10.0e-6 /) ! [m] Particle diameter grid

c     Geometric standard deviation (GSD)
c     She84 p. 75 Table 1, BSM96 p. 73 Table 2, PaG77 p. 2080 Table 1, SBG98 p. 10581 Table 1
c     She84, BSM96, PaG77, and SBG98 show GSD is always near 2.0
c     SBG98 shows GSD >> 2.0 is unsupportable for long range transport
c     SBG98 shows GSD =~ 2.0 is optimal for long range transport mode
      real(r8),parameter::gsd_anl_dfl=2.0 ! [frc] Geometric standard deviation
c     real(r8),parameter::gsd_anl_dfl=2.0       ! [frc] Geometric standard deviation SBG98 p. 75 Table 1
c     real(r8),parameter::gsd_anl_dfl=1.90       ! [frc] Geometric standard deviation BSM96 p. 73 Table 2
c     real(r8),parameter::gsd_anl_dfl=2.2       ! [frc] Geometric standard deviation PaG77 p. 2080 Table 1

c     Aerosol density 
c     PaG77 p. 2076 say dns_aer = 2500 kg m-3, DKS91 p. 118 say dns_aer = 1600 kg m-3
      real(r8),parameter::dns_aer_dfl=2.5e+3 ! [kg m-3] Aerosol density

c     Commons
#include <dstgrd.h> /* Dust grid sizes */
#include <dstaer.h> /* Aerosol microphysical properties */
#include <dstszdst.h> /* Diagnostic aerosol microphysical properties */
c     Input
c     Output
c     Local
      character fl_out*80       ! [sng] Name of netCDF output file
      integer src_idx           ! [idx] Counting index for src
      integer snk_idx           ! [idx] Counting index for snk
      integer idx               ! [idx] Counting index
      integer sz_grd_typ        ! [enm] Size grid type
      real(r8) dmt_ctr(dst_nbr)     ! [m] Diameter at bin center
      real(r8) dmt_dlt(dst_nbr)     ! [m] Width of size bin
c     Main Code

c     Set grid type
      sz_grd_typ=grd_typ_lgr    ! [enm] Size grid type

c     Create size grid
      call grd_mk(sz_grd_typ,   ! I
     $     dmt_min_min,dmt_max_max,dst_nbr, ! I
     $     dmt_ctr,dmt_dlt,dmt_min,dmt_max,dmt_grd) ! O

c     Override automatic grid with preset grid if available
      if (dst_nbr == 4) then
         do idx=1,dst_nbr+1     ! NB: dst_nbr+1
            dmt_grd(idx)=dmt_grd_04(idx) ! [m] Particle diameter grid
         enddo                  ! end loop over grd
         do idx=1,dst_nbr
            dmt_min(idx)=dmt_grd_04(idx) ! [m] Maximum diameter in bin
            dmt_max(idx)=dmt_grd_04(idx+1) ! [m] Minimum diameter in bin
            dmt_ctr(idx)=0.5*(dmt_min(idx)+dmt_max(idx)) ! [m] Diameter at bin center
            dmt_dlt(idx)=dmt_max(idx)-dmt_min(idx) ! [m] Width of size bin
         enddo                  ! end loop over grd
      endif ! endif dst_nbr == 4

c     Bin physical properties
      do idx=1,dst_nbr
         gsd_anl(idx)=gsd_anl_dfl ! [frc] Geometric standard deviation
         dns_aer(idx)=dns_aer_dfl ! [kg m-3] Aerosol density
      enddo                     ! end loop over grd

c     Set a fundamental statistic for each bin
c     For a large number of bins, statistic is irrelevent since distribution is resolved anyway
c     Choosing "right" statistic is key to value of method for relatively small number of bins (less than about 10) 
c     Setting bin center equal to number median diameter may be too simplistic
c     Setting dmt_vma to observed size distribution improves quadrature points when modeled distribution approaches observations
c     Fit to long-lived transport mode in order to capture long-range effects better
      do idx=1,dst_nbr            
         dmt_vma(idx)=2.524e-6  ! [m] Mass median diameter analytic She84 p. 75 Table 1
c         dmt_vma(idx)=4.82e-6  ! [m] Mass median diameter analytic BSM96 p. 73 Table 2
c         dmt_vma(idx)=5.6e-6  ! [m] Mass median diameter analytic PaG77 p. 2080 Table 1
c     rds_nma = 0.5*dmt_nma is dst_a input to mie() program
c         dmt_nma(idx)=dmt_ctr(idx) ! [m] Number median particle size
      enddo                     ! end loop over grd
      
c     Compute analytic size statistics

c     Convert mass median diameter to number median diameter
      call vma2nma(dst_nbr,dmt_vma,gsd_anl, ! I
     $     dmt_nma)             ! O [m] Number median particle size

cc     Convert number median diameter to mass median diameter
c      call nma2vma(dst_nbr,dmt_nma,gsd_anl, ! I
c     $     dmt_vma)             ! O [m] Mass median diameter
      
c     Compute resolved size statistics for each size distribution
      do idx=1,dst_nbr
         call dst_sz_rsl(
     $        dmt_min(idx),dmt_max(idx), ! I
     $        dns_aer(idx),dmt_nma(idx),gsd_anl(idx), ! I
     $        dmt_nwr(idx),dmt_swr(idx),dmt_vwr(idx), ! O
     $        dmt_nmr(idx),dmt_smr(idx),dmt_vmr(idx), ! O
     $        nbr_spc_rsl(idx),sfc_spc_rsl(idx),vlm_spc_rsl(idx) ! O
     $        )
      enddo                     ! end loop over sz
      
c     Source/sink fractional overlap factors 
      call ovr_src_snk_frc_get(dst_src_nbr,dmt_vma_src,gsd_anl_src, ! I
     $     dst_nbr,dmt_min,dmt_max, ! I
     $     ovr_src_snk_frc)     ! O
      
c     Multiply ovr_src_snk_frc(src_idx,*) by mss_frc(src_idx) to obtain
c     absolute mass fraction mapping from source dists. to sink bins
c     \sum_{i=1}^{i=src_nbr} \sum_{j=1}^{j=snk_nbr} ovr_src_snk_mss(i,j) <= 1
c     with equality holding iff sources and sinks are completely overlapped
c     Sum of ovr_src_snk_mss(*,*) is fraction of vertical dust flux which is transported
      ovr_src_snk_mss_ttl=0.0   ! [frc]
      do snk_idx=1,dst_nbr
         do src_idx=1,dst_src_nbr
            ovr_src_snk_mss(src_idx,snk_idx)= ! [frc]
     $           ovr_src_snk_frc(src_idx,snk_idx)*mss_frc_src(src_idx) ! [frc]
            ovr_src_snk_mss_ttl= ! [frc]
     $           ovr_src_snk_mss_ttl+ovr_src_snk_mss(src_idx,snk_idx)
         enddo                  ! end loop over src
      enddo                     ! end loop over snk
      
c     Corrections from Stokes' settling velocities to terminal velocities
      call stk_crc_get(
     $     dmt_vwr,             ! I [m] Mass-weighted mean diameter resolved
     $     dns_aer,             ! I [kg m-3] Aerosol density
     $     stk_crc,             ! O [frc] Correction to Stokes settling velocity
     $     dst_nbr)             ! I [nbr] Number of dust constituents
      
c     Size distribution diagnostics
c     Convert number median diameter to number mean diameter
      call nma2naa(dst_nbr,dmt_nma,gsd_anl, ! I
     $     dmt_naa)             ! O [m] Number mean diameter analytic
c     Convert number median diameter to surface median diameter
      call nma2sma(dst_nbr,dmt_nma,gsd_anl, ! I
     $     dmt_sma)             ! O [m] Surface median diameter analytic
c     Convert number median diameter to surface area weighted mean diameter
      call nma2swa(dst_nbr,dmt_nma,gsd_anl, ! I
     $     dmt_swa)             ! O [m] Surface weighted mean diameter analytic
c     Convert number median diameter to mass weighted mean diameter
      call nma2vwa(dst_nbr,dmt_nma,gsd_anl, ! I
     $     dmt_vwa)             ! O [m] Mass-weighted mean diameter analytic
      
c     Sanity checks
      do idx=1,dst_nbr
         write (6,*) ' xxx: ', idx, dmt_min(idx), dmt_nmr(idx), dmt_max(idx), dmt_nmr(idx)
         if (dmt_min(idx) > dmt_nmr(idx).or.dmt_max(idx) < dmt_nmr(idx)) stop
     $        'dst: dst_szdst_ini() reports dmt_min(idx) > dmt_nmr(idx).or.dmt_max(idx) < dmt_nmr(idx)'
         if (dmt_min(idx) > dmt_vwr(idx).or.dmt_max(idx) < dmt_vwr(idx)) stop
     $        'dst: dst_szdst_ini() reports dmt_min(idx) > dmt_vwr(idx).or.dmt_max(idx) < dmt_vwr(idx)'
         if (gsd_anl(idx) < 1.0.or.gsd_anl(idx) > 3.0) stop
     $        'dst: dst_szdst_ini() reports gsd_anl(idx) < 1.0.or.gsd_anl(idx) > 3.0'
         if (dns_aer(idx) < 500.0.or.dns_aer(idx) > 3000.0) stop
     $        'dst: dst_szdst_ini() reports dns_aer(idx) < 500.0.or.dns_aer(idx) > 3000.0'
      enddo                     ! end loop over cst
      
c     Display banner
      call dst_bnr()

c     Update netCDF file
      fl_out='aer.nc' ! [sng] Name of netCDF output file
      call aersz2nc(
     $     fl_out               ! [sng] Name of netCDF output file
     $     )                    
      
      write (6,'(a,f4.3)') 'Total transported mass fraction of dust flux = ',ovr_src_snk_mss_ttl
      write (6,'(a)') 'Transported size classes:'
      write (6,'(a3,7(a10,1x))') 'idx','dmt_min','dmt_max','dmt_ctr','dmt_nmr','dmt_nwr','dmt_vmr','dmt_vwr'
      write (6,'(a3,7(a10,1x))') '','um','um','um','um','um','um','um'
      do idx=1,dst_nbr
         write (6,'(1p,i3,7(e10.3,1x)),0p') idx,
     $        dmt_min(idx)*1.0e6,dmt_max(idx)*1.0e6,dmt_ctr(idx)*1.0e6,
     $        dmt_nmr(idx)*1.0e6,dmt_nwr(idx)*1.0e6,
     $        dmt_vmr(idx)*1.0e6,dmt_vwr(idx)*1.0e6
      enddo                     ! end loop over cst

      write (6,'(a)') 'Size distribution information:'
      write (6,'(a3,4(a10,1x))') 'idx','dmt_nma','dmt_vma','gsd','dns'
      write (6,'(a3,4(a10,1x))') '','um','um','frc','kg m-3'
      do idx=1,dst_nbr
         write (6,'(1p,i3,4(e10.3,1x)),0p') idx,
     $        dmt_nma(idx)*1.0e6,dmt_vma(idx)*1.0e6,
     $        gsd_anl(idx),dns_aer(idx)
      enddo                     ! end loop over cst

      write (6,'(a)') 'Chemistry information:'
      write (6,'(a3,3(a10,1x))') 'idx','stk_crc','nbr_spc','sfc_spc'
      write (6,'(a3,3(a10,1x))') '','frc','# kg-1','m2 kg-1'
      do idx=1,dst_nbr
         write (6,'(1p,i3,3(e10.3,1x)),0p') idx,
     $        stk_crc(idx),nbr_spc_rsl(idx),sfc_spc_rsl(idx)
      enddo                     ! end loop over cst

      write (6,'(a)') '------------------------------------------------------------'
      return
      end subroutine dst_szdst_ini       ! end dst_szdst_ini()
      
      subroutine dst_sz_rsl(
     $     dmt_min,dmt_max,     ! I
     $     dns_aer,dmt_nma,gsd, ! I
     $     dmt_nwr,dmt_swr,dmt_vwr, ! O
     $     dmt_nmr,dmt_smr,dmt_vmr, ! O
     $     nbr_spc_rsl,sfc_spc_rsl,vlm_spc_rsl ! O
     $     )
c     Compute resolved size statistics for a single size distribution
      implicit none
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
c     Input
      integer grd_typ           ! Type of grid
      real(r8) dns_aer              ! [kg m-3] Particle density
      real(r8) dmt_nma              ! [m] Number median particle size
      real(r8) dmt_max              ! [m] Size grid maximum
      real(r8) dmt_min              ! [m] Size grid minimum
      real(r8) gsd                  ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_nmr              ! [m] Number median diameter resolved
      real(r8) dmt_smr              ! [m] Surface area median diameter resolved
      real(r8) dmt_vmr              ! [m] Mass median diameter resolved
      real(r8) dmt_nwr              ! [m] Number weighted diameter resolved
      real(r8) dmt_swr              ! [m] Surface area weighted diameter resolved
      real(r8) dmt_vwr              ! [m] Mass weighted diameter resolved
      real(r8) nbr_spc_rsl          ! [# kg-1] Specific concentration resolved 
      real(r8) sfc_spc_rsl          ! [m2 kg-1] Specific surface area resolved 
      real(r8) vlm_spc_rsl          ! [m3 kg-1] Specific volume resolved 
c     Local
      double precision pi       ! [frc] 3
      integer sz_grd_typ        ! Size grid type
      integer grd_typ_lgr       ! Logarithmic
      integer grd_typ_lnr       ! Linear
      integer idx               ! [idx] Counting index
      real(r8) cnc(sgs_nbr)         ! [# m-3] Lognormal concentration at sz_ctr
      real(r8) cnc_nbr_rsl          ! [# m-3] Number concentration resolved
      real(r8) dst(sgs_nbr)         ! [# m-3 m-1] Lognormal distribution at sz_ctr
      real(r8) sz_ctr(sgs_nbr)      ! [m] Size Bin centers
      real(r8) sz_dlt(sgs_nbr)      ! [m] Size Bin widths
      real(r8) sz_grd(sgs_nbr+1)    ! [m] Size Grid
      real(r8) sz_max(sgs_nbr)      ! [m] Size Bin maxima
      real(r8) sz_min(sgs_nbr)      ! [m] Size Bin minima
      real(r8) vlm(sgs_nbr)         ! [m3] Volume
      real(r8) vlm_rsl              ! [m3 m-3] Volume concentration resolved
      real(r8) sfc(sgs_nbr)         ! [m2] Surface area
      real(r8) sfc_rsl              ! [m2 m-3] Surface area concentration resolved
      real(r8) nbr_prt_rsl(sgs_nbr) ! [# m-3] Number concentration of smaller particles
      real(r8) nbr_prt_rsl_frc(sgs_nbr) ! [frc] Fraction of number concentration from smaller particles
      real(r8) sfc_prt_rsl(sgs_nbr) ! [m2 m-3] Surface area concentration of smaller particles
      real(r8) sfc_prt_rsl_frc(sgs_nbr) ! [frc] Fraction of surface area concentration from smaller particles
      real(r8) vlm_prt_rsl(sgs_nbr) ! [m3 m-3] Volume concentration of smaller particles
      real(r8) vlm_prt_rsl_frc(sgs_nbr) ! [frc] Fraction of volume concentration from smaller particles
      real(r8) half(1)
c     External
      real(r8) ntp_vec_one          ! Scalar interpolation (libcsz_f77)
      
c     Main Code
c     Initialize
      pi=4.0*atan(1.0d0)        ! [frc] 3
c     Set grid types
      grd_typ_lgr=0             ! Logarithmic
      grd_typ_lnr=1             ! Linear
      sz_grd_typ=grd_typ_lgr
      
c     Create size grid
      call grd_mk(sz_grd_typ,   ! I
     $     dmt_min,dmt_max,sgs_nbr, ! I
     $     sz_ctr,sz_dlt,sz_min,sz_max,sz_grd) ! O

c     Evaluate lognormal distribution for these sizes
      call lgn_evl(sgs_nbr,sz_ctr,dmt_nma,gsd, ! I
     $     dst)                 ! O
      
c     Integrate moments of size distribution
      cnc_nbr_rsl=0.0           ! [# m-3] Number concentration resolved
      sfc_rsl=0.0               ! [m2 m-3] Surface area concentration resolved
      vlm_rsl=0.0               ! [m3 m-3] Volume concentration resolved
      dmt_nwr=0.0               ! [m] Number weighted diameter resolved
      dmt_swr=0.0               ! [m] Surface area weighted diameter resolved
      dmt_vwr=0.0               ! [m] Mass weighted diameter resolved
      do idx=1,sgs_nbr
         cnc(idx)=dst(idx)*sz_dlt(idx) ! [# m-3] Number concentration
         cnc_nbr_rsl=cnc_nbr_rsl+cnc(idx) ! [# m-3] Number concentration resolved
         nbr_prt_rsl(idx)=cnc_nbr_rsl ! [# m-3] Number concentration of smaller particles
         sfc(idx)=pi*sz_ctr(idx)*sz_ctr(idx) ! [m2] Surface area 
         sfc_rsl=sfc_rsl+sfc(idx)*cnc(idx) ! [m2 m-3] Surface area distribution resolved so far
         sfc_prt_rsl(idx)=sfc_rsl ! [m2 m-3] Surface area of smaller particles
         vlm(idx)=(1.0/6.0)*pi*(sz_ctr(idx)**3.0) ! [m3] Volume 
         vlm_rsl=vlm_rsl+vlm(idx)*cnc(idx) ! [m3 m-3] Volume distribution resolved so far
         vlm_prt_rsl(idx)=vlm_rsl ! [m3 m-3] Volume of smaller particles
         dmt_nwr=dmt_nwr+sz_ctr(idx)*cnc(idx) ! [m] Number weighted diameter resolved
         dmt_swr=dmt_swr+sz_ctr(idx)*sfc(idx)*cnc(idx) ! [m] Surface area weighted diameter resolved
         dmt_vwr=dmt_vwr+sz_ctr(idx)*vlm(idx)*cnc(idx) ! [m] Mass weighted diameter resolved
      enddo                     ! end loop over size
      
c     Normalize moment weighted sizes by moment
      dmt_nwr=dmt_nwr/cnc_nbr_rsl ! [m] Number weighted diameter resolved
      dmt_swr=dmt_swr/sfc_rsl   ! [m] Surface area weighted diameter resolved
      dmt_vwr=dmt_vwr/vlm_rsl   ! [m] Mass weighted diameter resolved
      
c     Statistics per unit mass
      nbr_spc_rsl=cnc_nbr_rsl/(vlm_rsl*dns_aer) ! [# kg-1] Specific concentration resolved 
      sfc_spc_rsl=sfc_rsl/(vlm_rsl*dns_aer) ! [m2 kg-1] Specific surface area resolved 
      vlm_spc_rsl=vlm_rsl/(vlm_rsl*dns_aer) ! [m3 kg-1] Specific volume resolved 
      
      do idx=1,sgs_nbr
         nbr_prt_rsl_frc(idx)=nbr_prt_rsl(idx)/cnc_nbr_rsl ! [frc] Fraction of number concentration from smaller particles
         sfc_prt_rsl_frc(idx)=sfc_prt_rsl(idx)/sfc_rsl ! [frc] Fraction of surface area concentration from smaller particles
         vlm_prt_rsl_frc(idx)=vlm_prt_rsl(idx)/vlm_rsl ! [frc] Fraction of volume concentration from smaller particles
      enddo                     ! end loop over size
      
c     Find resolved median sizes by interpolating inverted size -> [number,surface,volume] relationships
      half(1) = 0.5
      if (sgs_nbr > 1) then
         dmt_nmr=ntp_vec_one(sgs_nbr,nbr_prt_rsl_frc,sz_ctr,half) ! [m] Number median diameter resolved
         dmt_smr=ntp_vec_one(sgs_nbr,sfc_prt_rsl_frc,sz_ctr,half) ! [m] Surface area median diameter resolved
         dmt_vmr=ntp_vec_one(sgs_nbr,vlm_prt_rsl_frc,sz_ctr,half) ! [m] Volume median diameter resolved
      else
         dmt_nmr=sz_ctr(1)      ! [m] Number median diameter resolved
         dmt_smr=sz_ctr(1)      ! [m] Surface area median diameter resolved
         dmt_vmr=sz_ctr(1)      ! [m] Volume median diameter resolved
      endif                     ! endif
      
      return
      end subroutine dst_sz_rsl                      ! end dst_sz_rsl()
      
      subroutine grd_mk(sz_grd_typ, ! I
     $     sz_grd_min,sz_grd_max,sz_nbr, ! I
     $     sz_ctr,sz_dlt,sz_min,sz_max,sz_grd) ! O
c     Purpose: Create a size grid
c     Routine uses sz nomenclature but works for any (linear) grid
      implicit none
c     Input
      integer sz_grd_typ        ! Type of grid
      integer sz_nbr            ! [nbr] Number of bins
      real(r8) sz_grd_min           ! Grid minimum
      real(r8) sz_grd_max           ! Grid maximum
c     Output
      real(r8) sz_grd(sz_nbr+1)     ! Grid
      real(r8) sz_ctr(sz_nbr)       ! Bin centers
      real(r8) sz_dlt(sz_nbr)       ! Bin widths
      real(r8) sz_min(sz_nbr)       ! Bin minima
      real(r8) sz_max(sz_nbr)       ! Bin maxima
c     Local
      integer idx               ! [idx] Counting index
      real(r8) sz_grd_rsn           ! Resolution for linear grid
      real(r8) max_min_ratio        ! Factor for logarithmic grid
      real(r8) series_ratio         ! Factor for logarithmic grid
      integer grd_typ_lgr       ! Logarithmic
      integer grd_typ_lnr       ! Linear
c     Main Code
c     Initialize
      grd_typ_lgr=0             ! Logarithmic
      grd_typ_lnr=1             ! Linear
      
c     write (6,*) ' sz_nbr ', sz_nbr

      if (sz_grd_typ == grd_typ_lnr) then
         sz_grd_rsn=(sz_grd_max-sz_grd_min)/sz_nbr
         do idx=1,sz_nbr
            sz_min(idx)=sz_grd_min+(idx-1)*sz_grd_rsn ! [m]
         enddo                  ! end loop over grd
      else if (sz_grd_typ == grd_typ_lgr) then
         max_min_ratio=sz_grd_max/sz_grd_min
         series_ratio=max_min_ratio**(1.0/sz_nbr)
         if (sz_grd_min == 0.0) stop 'dst: grd_mk() reports sz_grd_min = 0.0'
         sz_min(1)=sz_grd_min
         do idx=2,sz_nbr        ! Loop starts at 2
            sz_min(idx)=sz_min(idx-1)*series_ratio
         enddo                  ! end loop over grd
      else
         stop 'dst: grd_mk() reports unknown grid type'
      endif                     ! endif
c      write (6,*) ' grd_mk: sz_grd_min ', sz_grd_min
c      write (6,*) ' grd_mk: sz_min ', sz_min
      
c     Derived grid values 
!pjr moved before the do loop to get past a compiler problem with xlf90_r
!    the following line was not executed if after the do loop
      sz_max(sz_nbr)=sz_grd_max ! [m]
      do idx=1,sz_nbr-1         ! Loop ends at sz_nbr-1
         sz_max(idx)=sz_min(idx+1) ! [m]
      enddo                     ! end loop over grd
c     sz_max(sz_nbr)=sz_grd_max ! [m]
c     write (6,*) ' grd_mk: sz_max ', sz_max
      do idx=1,sz_nbr
         sz_grd(idx)=sz_min(idx) ! [m]
      enddo                     ! end loop over grd
      sz_grd(sz_nbr+1)=sz_max(sz_nbr)
c      write (6,*) ' grd_mk: sz_grd ', sz_grd
      sz_grd_max=sz_max(sz_nbr) ! [m]
      sz_grd_min=sz_min(1)      ! [m]
      
c     Final derived grid values 
      do idx=1,sz_nbr
         sz_ctr(idx)=0.5*(sz_min(idx)+sz_max(idx))
         sz_dlt(idx)=sz_max(idx)-sz_min(idx)
      enddo                     ! end loop over grd
      
      return
      end subroutine grd_mk                      ! end grd_mk()
      
      subroutine aersz2nc(            
     $     fl_out               ! [sng] Name of netCDF output file
     $     )                    
c     Purpose: Output time-constant aerosol properties to netCDF file
      use netcdf
      implicit none
#include <pmgrid.h> /* Spatial resolution parameters */
#include <dstgrd.h> /* Dust grid sizes */
#include <dstaer.h> /* Aerosol microphysical properties */
#include <dstszdst.h> /* Diagnostic aerosol microphysical properties */
#include <dstctl.h> /* Model control variables */
#include <dbg.com> /* Debugging constants */
#include <netcdf.inc>
c     Input
      character*(*) fl_out      ! [sng] Name of netCDF output file
c     Output
c     Local
      integer idx               ! [idx] Counting index
c     File metadata and dimension IDs
      integer dim_lon_lev_sz(3) ! [enm] Dimension IDs
      integer dim_lon_sz(2)     ! [enm] Dimension IDs
      integer dim_sz_src_sz(2)  ! [enm] Dimension IDs
      integer fll_mode_old      ! Old fill mode
      integer lev_dim_id        ! [enm] Dimension ID for lev
      integer levp_dim_id       ! [enm] Dimension ID for levp
      integer lon_dim_id        ! [enm] Dimension ID for lon
      integer lat_dim_id        ! [enm] Dimension ID for lat
      integer time_dim_id       ! [enm] Dimension ID for time
      integer nc_id             ! File handle
      integer rcd               ! [rcd] Return success code
      integer sz_dim_id         ! [enm] Dimension ID for sz
      integer sz_grd_dim_id     ! [enm] Dimension ID for sz grid
      integer sz_src_dim_id     ! [enm] Dimension ID for sz_src
c     Variable IDs
      integer dmt_naa_id        ! [enm] Variable ID
      integer dmt_vma_id        ! [enm] Variable ID
      integer dmt_vma_src_id    ! [enm] Variable ID
      integer dmt_nma_id        ! [enm] Variable ID
      integer dmt_sma_id        ! [enm] Variable ID
      integer dns_aer_id        ! [enm] Variable ID
      integer dmt_vwa_id        ! [enm] Variable ID
      integer dmt_swa_id        ! [enm] Variable ID
      integer gsd_anl_id        ! [enm] Variable ID
      integer gsd_anl_src_id    ! [enm] Variable ID
      integer dmt_max_id        ! [enm] Variable ID
      integer dmt_min_id        ! [enm] Variable ID
      integer dmt_nwr_id        ! [enm] Variable ID
      integer dmt_swr_id        ! [enm] Variable ID
      integer dmt_vwr_id        ! [enm] Variable ID
      integer dmt_nmr_id        ! [enm] Variable ID
      integer dmt_smr_id        ! [enm] Variable ID
      integer dmt_vmr_id        ! [enm] Variable ID
      integer mss_frc_src_id    ! [enm] Variable ID
      integer ovr_src_snk_frc_id ! [enm] Variable ID
      integer ovr_src_snk_mss_id ! [enm] Variable ID
      integer ovr_src_snk_mss_ttl_id ! [enm] Variable ID
      integer nbr_spc_rsl_id    ! [enm] Variable ID
      integer sfc_spc_rsl_id    ! [enm] Variable ID
      integer vlm_spc_rsl_id    ! [enm] Variable ID
      integer sz_grd_id         ! Coordinate ID
      integer sz_id             ! Coordinate ID
      integer sz_src_id         ! Coordinate ID
c     Variable data
      real(r8) sz(dst_nbr)          ! [m] Coordinate variable
! pjr modified on charlies advice
!      real(r8) sz_src(dst_nbr)      ! [m] Coordinate variable
      real(r8) sz_src(dst_src_nbr)      ! [m] Coordinate variable
      real(r8) sz_grd(dst_nbr+1)    ! [m] Coordinate variable
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
c     Initialize coordinates
      do idx=1,dst_src_nbr
         sz_src(idx)=dmt_vma_src(idx)
      enddo                     ! end loop over sz
      do idx=1,dst_nbr
         sz(idx)=0.5*(dmt_min(idx)+dmt_max(idx))
         sz_grd(idx)=dmt_min(idx)
      enddo                     ! end loop over sz
      sz_grd(dst_nbr+1)=dmt_max(dst_nbr)
c     Begin netCDF output routines
      rcd=NF90_NOERR              ! NF90_NOERR == 0
      rcd=rcd+nf90_create(fl_out,NF90_CLOBBER,nc_id)
      if (rcd /= NF90_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),
     $     ': ERROR creating ',fl_out(1:f_strlen(fl_out)),
     $     ' in aersz2nc(), rcd = ',rcd
c     Add global attributes
c     Define dimension IDs
      rcd=rcd+nf90_def_dim(nc_id,'time',NF90_UNLIMITED,time_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'sz_src',dst_src_nbr,sz_src_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'sz',dst_nbr,sz_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'lat',plat,lat_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'lon',plon,lon_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'lev',plev,lev_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'levp',plevp,levp_dim_id)
      rcd=rcd+nf90_def_dim(nc_id,'sz_grd',dst_nbr+1,sz_grd_dim_id)
c     Assemble ID and count vectors for each multidimensional combination of dimensions
      dim_lon_sz(1)=lon_dim_id
      dim_lon_sz(2)=sz_dim_id
      dim_lon_lev_sz(1)=lon_dim_id
      dim_lon_lev_sz(2)=lev_dim_id
      dim_lon_lev_sz(3)=sz_dim_id
      dim_sz_src_sz(1)=sz_src_dim_id
      dim_sz_src_sz(2)=sz_dim_id
c     Variable definitions
      rcd=rcd+nf90_def_var(nc_id,'dmt_max',NF90_DOUBLE,sz_dim_id,dmt_max_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_min',NF90_DOUBLE,sz_dim_id,dmt_min_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_naa',NF90_DOUBLE,sz_dim_id,dmt_naa_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_nma',NF90_DOUBLE,sz_dim_id,dmt_nma_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_nmr',NF90_DOUBLE,sz_dim_id,dmt_nmr_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_nwr',NF90_DOUBLE,sz_dim_id,dmt_nwr_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_sma',NF90_DOUBLE,sz_dim_id,dmt_sma_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_smr',NF90_DOUBLE,sz_dim_id,dmt_smr_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_swa',NF90_DOUBLE,sz_dim_id,dmt_swa_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_swr',NF90_DOUBLE,sz_dim_id,dmt_swr_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_vma',NF90_DOUBLE,sz_dim_id,dmt_vma_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_vma_src',NF90_DOUBLE,sz_src_dim_id,dmt_vma_src_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_vmr',NF90_DOUBLE,sz_dim_id,dmt_vmr_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_vwa',NF90_DOUBLE,sz_dim_id,dmt_vwa_id)
      rcd=rcd+nf90_def_var(nc_id,'dmt_vwr',NF90_DOUBLE,sz_dim_id,dmt_vwr_id)
      rcd=rcd+nf90_def_var(nc_id,'dns_aer',NF90_DOUBLE,sz_dim_id,dns_aer_id)
      rcd=rcd+nf90_def_var(nc_id,'gsd_anl',NF90_DOUBLE,sz_dim_id,gsd_anl_id)
      rcd=rcd+nf90_def_var(nc_id,'gsd_anl_src',NF90_DOUBLE,sz_src_dim_id,gsd_anl_src_id)
      rcd=rcd+nf90_def_var(nc_id,'mss_frc_src',NF90_DOUBLE,sz_src_dim_id,mss_frc_src_id)
      rcd=rcd+nf90_def_var(nc_id,'nbr_spc_rsl',NF90_DOUBLE,sz_dim_id,nbr_spc_rsl_id)
      rcd=rcd+nf90_def_var(nc_id,'ovr_src_snk_frc',NF90_DOUBLE,dim_sz_src_sz,ovr_src_snk_frc_id)
      rcd=rcd+nf90_def_var(nc_id,'ovr_src_snk_mss',NF90_DOUBLE,dim_sz_src_sz,ovr_src_snk_mss_id)
! fewer arguments in next line to indicate scalar
      rcd=rcd+nf90_def_var(nc_id,'ovr_src_snk_mss_ttl',NF90_DOUBLE,ovr_src_snk_mss_ttl_id)
      rcd=rcd+nf90_def_var(nc_id,'sfc_spc_rsl',NF90_DOUBLE,sz_dim_id,sfc_spc_rsl_id)
      rcd=rcd+nf90_def_var(nc_id,'sz',NF90_DOUBLE,sz_dim_id,sz_id)
      rcd=rcd+nf90_def_var(nc_id,'sz_grd',NF90_DOUBLE,sz_grd_dim_id,sz_grd_id)
      rcd=rcd+nf90_def_var(nc_id,'sz_src',NF90_DOUBLE,sz_src_dim_id,sz_src_id)
      rcd=rcd+nf90_def_var(nc_id,'vlm_spc_rsl',NF90_DOUBLE,sz_dim_id,vlm_spc_rsl_id)
c     Add english text descriptions
      rcd=rcd+nf90_put_att(nc_id,dmt_max_id,'long_name','Maximum diameter in bin')
      rcd=rcd+nf90_put_att(nc_id,dmt_min_id,'long_name','Minimum diameter in bin')
      rcd=rcd+nf90_put_att(nc_id,dmt_naa_id,'long_name','Number mean diameter analytic')
      rcd=rcd+nf90_put_att(nc_id,dmt_nma_id,'long_name','Number median diameter analytic')
      rcd=rcd+nf90_put_att(nc_id,dmt_nmr_id,'long_name','Number median diameter resolved')
      rcd=rcd+nf90_put_att(nc_id,dmt_nwr_id,'long_name','Number weighted diameter resolved')
      rcd=rcd+nf90_put_att(nc_id,dmt_sma_id,'long_name','Surface median diameter analytic')
      rcd=rcd+nf90_put_att(nc_id,dmt_smr_id,'long_name','Surface area median diameter resolved')
      rcd=rcd+nf90_put_att(nc_id,dmt_swa_id,'long_name','Surface area weighted mean diameter analytic')
      rcd=rcd+nf90_put_att(nc_id,dmt_swr_id,'long_name','Surface area weighted diameter resolved')
      rcd=rcd+nf90_put_att(nc_id,dmt_vma_id,'long_name','Mass median diameter analytic')
      rcd=rcd+nf90_put_att(nc_id,dmt_vma_src_id,'long_name','Mass median diameter of source distribution')
      rcd=rcd+nf90_put_att(nc_id,dmt_vmr_id,'long_name','Mass median diameter resolved')
      rcd=rcd+nf90_put_att(nc_id,dmt_vwa_id,'long_name','Mass weighted mean diameter analytic')
      rcd=rcd+nf90_put_att(nc_id,dmt_vwr_id,'long_name','Mass weighted diameter resolved')
      rcd=rcd+nf90_put_att(nc_id,dns_aer_id,'long_name','Particle density')
      rcd=rcd+nf90_put_att(nc_id,gsd_anl_id,'long_name','Geometric standard deviation')
      rcd=rcd+nf90_put_att(nc_id,gsd_anl_src_id,'long_name','Geometric standard deviation of source distribution')
      rcd=rcd+nf90_put_att(nc_id,mss_frc_src_id,'long_name','Mass fraction of source distribution')
      rcd=rcd+nf90_put_att(nc_id,nbr_spc_rsl_id,'long_name','Specific number concentration resolved')
      rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_frc_id,'long_name','Overlap of src dist. i with sink bin j')
      rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_mss_id,'long_name','Mass overlap of src dist. i with sink bin j')
      rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_mss_ttl_id,'long_name','Total transported mass fraction of dust flux')
      rcd=rcd+nf90_put_att(nc_id,sfc_spc_rsl_id,'long_name','Specific surface area resolved')
      rcd=rcd+nf90_put_att(nc_id,sz_grd_id,'long_name','Size grid interfaces')
      rcd=rcd+nf90_put_att(nc_id,sz_id,'long_name','Nominal size bin center')
      rcd=rcd+nf90_put_att(nc_id,sz_src_id,'long_name','Mass median diameter of source distribution')
      rcd=rcd+nf90_put_att(nc_id,vlm_spc_rsl_id,'long_name','Specific volume resolved')
c     Add units
      rcd=rcd+nf90_put_att(nc_id,dmt_max_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_min_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_naa_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_nma_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_nmr_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_nwr_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_sma_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_smr_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_swa_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_swr_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_vma_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_vma_src_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_vmr_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_vwa_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dmt_vwr_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,dns_aer_id,'units','kilogram meter-3')
      rcd=rcd+nf90_put_att(nc_id,gsd_anl_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,gsd_anl_src_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,mss_frc_src_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,nbr_spc_rsl_id,'units','number kilogram-1')
      rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_frc_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_mss_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_mss_ttl_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,sfc_spc_rsl_id,'units','meter2 kilogram-1')
      rcd=rcd+nf90_put_att(nc_id,sz_grd_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,sz_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,sz_src_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,vlm_spc_rsl_id,'units','meter3 kilogram-1')
c     End define mode
      rcd=rcd+nf90_enddef(nc_id)
c     Write data
      rcd=rcd+nf90_put_var(nc_id,dmt_max_id,dmt_max)
      rcd=rcd+nf90_put_var(nc_id,dmt_min_id,dmt_min)
      rcd=rcd+nf90_put_var(nc_id,dmt_naa_id,dmt_naa)
      rcd=rcd+nf90_put_var(nc_id,dmt_nma_id,dmt_nma)
      rcd=rcd+nf90_put_var(nc_id,dmt_nmr_id,dmt_nmr)
      rcd=rcd+nf90_put_var(nc_id,dmt_nwr_id,dmt_nwr)
      rcd=rcd+nf90_put_var(nc_id,dmt_sma_id,dmt_sma)
      rcd=rcd+nf90_put_var(nc_id,dmt_smr_id,dmt_smr)
      rcd=rcd+nf90_put_var(nc_id,dmt_swa_id,dmt_swa)
      rcd=rcd+nf90_put_var(nc_id,dmt_swr_id,dmt_swr)
      rcd=rcd+nf90_put_var(nc_id,dmt_vma_id,dmt_vma)
      rcd=rcd+nf90_put_var(nc_id,dmt_vma_src_id,dmt_vma_src)
      rcd=rcd+nf90_put_var(nc_id,dmt_vmr_id,dmt_vmr)
      rcd=rcd+nf90_put_var(nc_id,dmt_vwa_id,dmt_vwa)
      rcd=rcd+nf90_put_var(nc_id,dmt_vwr_id,dmt_vwr)
      rcd=rcd+nf90_put_var(nc_id,dns_aer_id,dns_aer)
      rcd=rcd+nf90_put_var(nc_id,gsd_anl_id,gsd_anl)
      rcd=rcd+nf90_put_var(nc_id,gsd_anl_src_id,gsd_anl_src)
      rcd=rcd+nf90_put_var(nc_id,mss_frc_src_id,mss_frc_src)
      rcd=rcd+nf90_put_var(nc_id,nbr_spc_rsl_id,nbr_spc_rsl)
      rcd=rcd+nf90_put_var(nc_id,ovr_src_snk_frc_id,ovr_src_snk_frc)
      rcd=rcd+nf90_put_var(nc_id,ovr_src_snk_mss_id,ovr_src_snk_mss)
      rcd=rcd+nf90_put_var(nc_id,ovr_src_snk_mss_ttl_id,ovr_src_snk_mss_ttl)
      rcd=rcd+nf90_put_var(nc_id,sfc_spc_rsl_id,sfc_spc_rsl)
      rcd=rcd+nf90_put_var(nc_id,sz_grd_id,sz_grd)
      rcd=rcd+nf90_put_var(nc_id,sz_id,sz)
      rcd=rcd+nf90_put_var(nc_id,sz_src_id,sz_src)
      rcd=rcd+nf90_put_var(nc_id,vlm_spc_rsl_id,vlm_spc_rsl)
c     Close output file
      rcd=rcd+nf90_close(nc_id)
      call strnul(fl_out)
      if (nstep == 1) then
         write (6,'(a,a47,1x,a)') prg_nm(1:f_strlen(prg_nm)),
     $        ': Initialized size distribution data archive in',fl_out(1:f_strlen(fl_out))
      endif                     ! endif
      if (rcd /= NF90_NOERR) write (6,'(4a,i4)') prg_nm(1:f_strlen(prg_nm)),
     $     ': ERROR writing ',fl_out(1:f_strlen(fl_out)),
     $     ' in aersz2nc(), rcd = ',rcd
      
      return 
      end subroutine aersz2nc                      ! end aersz2nc()
      
      end module dstszdst
