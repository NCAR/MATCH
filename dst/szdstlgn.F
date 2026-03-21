c     $Header$
      
c     Purpose: szdstlgn.F contains routines to initialize, compute, and manipulate
c     lognormal size distributions
c     diff ~/aer/szdstlgn.F /fs/cgd/data0/zender/ccm_dst/src/dst/szdstlgn.F
      
c     NB: Unless otherwise specified, all functions assume standard, untruncated, lognormal distributions 
      
c     Size distribution nomenclature:
c     Three letter acronym: 
c     First letter [nsvx] = number (nbr), surface area (sfc), volume (vlm), cross-sectional area (xsa)
c     Second letter [amw] = average (i.e., mean) (avg), median (mdn), weighted (wgt)
c     Third letter [ar] = analytic (anl), resolved (rsl)
c     If necessary, a fourth letter option ("o"?) could be used 
c     in the second letter slot to denote true "mode" averages 
c     Not all 24 permutations are unique (much less useful):
c     1. Number mean diameter equals number weighted diameter: naa=nwa, nar=nwr
c     2. Surface area "mean" diameters equal cross-sectional area "mean" diameters 
c     for spheres since factor of 4 difference on top and bottom cancel eachother: 
c     saa=xaa,sma=xma,sar=xar,smr=xmr

c     naa = nbr avg anl 
c     nar = nbr avg rsl
c     nma = nbr mdn anl
c     nmr = nbr mdn rsl
c     nwa = nbr wgt anl
c     nwr = nbr wgt rsl
c     saa = sfc avg anl
c     sar = sfc avg rsl
c     sma = sfc mdn anl
c     smr = sfc mdn rsl
c     swa = sfc wgt anl
c     swr = sfc wgt rsl
c     vaa = vlm avg anl
c     var = vlm avg rsl
c     vma = vlm mdn anl
c     vmr = vlm mdn rsl
c     vwa = vlm wgt anl
c     vwr = vlm wgt rsl
c     xaa = xsa avg anl
c     xar = xsa avg rsl
c     xma = xsa mdn anl
c     xmr = xsa mdn rsl
c     xwa = xsa wgt anl
c     xwr = xsa wgt rsl

      subroutine lgn_evl(abc_nbr,abc_ctr,mdn,gsd, ! I
     $     lgn_dst)             ! O
c     Purpose: Given the median value and geometric standard deviation
c     for a single lognormal size distribution, compute and return the value of 
c     the distribution for each of the input values.
c     This function works for any lognormal size distribution, thus when input
c     values are diameters, this function returns n=dN/dD

      use precision

      implicit none
c     Input
      integer abc_nbr           ! [nbr] Dimension size
      real(r8) mdn                  ! Median value of distribution
      real(r8) abc_ctr(abc_nbr)     ! Abscissae
      real(r8) gsd                  ! [frc] Geometric standard deviation
c     Output
      real(r8) lgn_dst(abc_nbr)     ! Lognormal distribution at abc_ctr
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
      double precision pi       ! [frc] 3
      double precision lngsdsqrttwopi_rcp ! Factor in lognormal distribution
      double precision tmp      ! Factor in lognormal distribution
      double precision xpn      ! Factor in lognormal distribution
c     Main Code
      pi=4.0*atan(1.0d0)        ! [frc] 3
      ln_gsd=log(gsd)
      lngsdsqrttwopi_rcp=1.0/(ln_gsd*sqrt(2.0*pi))
      do idx=1,abc_nbr
         tmp=log(abc_ctr(idx)/mdn)/ln_gsd
         xpn=exp(-0.5*tmp*tmp)
         lgn_dst(idx)=lngsdsqrttwopi_rcp*xpn/abc_ctr(idx)
      enddo                     ! end loop over abc
      return
      end                       ! end lgn_evl()
      
      subroutine naa2nma(aer_nbr,dmt_avg,gsd, ! I
     $     dmt_mdn)             ! O
c     Purpose: Given the mean values and geometric standard deviations
c     for a number of lognormal size distributions, compute and return the median
c     values of each distribution 

      use precision

      implicit none
c     Input
      integer aer_nbr           ! [nbr] Dimension size
      real(r8) dmt_avg(aer_nbr)     ! Mean value
      real(r8) gsd(aer_nbr)         ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_mdn(aer_nbr)     ! Median value
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
c     Main Code
      do idx=1,aer_nbr
         ln_gsd=log(gsd(idx))
         dmt_mdn(idx)=dmt_avg(idx)*exp(-ln_gsd*ln_gsd/2.0)
      enddo                     ! end loop over aer
      return
      end                       ! end naa2nma()
      
      subroutine nma2naa(aer_nbr,dmt_mdn,gsd, ! I
     $     dmt_avg)             ! O
c     Purpose: Given the median values and geometric standard deviations
c     for a number of lognormal size distributions, compute and return the mean
c     value of each distribution 
c     Definition of gsd conforms to SeP97

      use precision

      implicit none
c     Input
      integer aer_nbr           ! [nbr] Dimension size
      real(r8) dmt_mdn(aer_nbr)     ! Median value
      real(r8) gsd(aer_nbr)         ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_avg(aer_nbr)     ! Mean value
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
c     Main Code
      do idx=1,aer_nbr
         ln_gsd=log(gsd(idx))
         dmt_avg(idx)=dmt_mdn(idx)*exp(ln_gsd*ln_gsd/2.0)
      enddo                     ! end loop over aer
      return
      end                       ! end nma2naa()
      
      subroutine nma2swa(aer_nbr,dmt_nma,gsd, ! I
     $     dmt_swa)             ! O
c     Purpose: Given the number median particle sizes and geometric standard deviations
c     for a number of lognormal size distributions, compute and return the 
c     area-weighted mean size of each distribution
c     Definition of gsd conforms to SeP97

      use precision

      implicit none
c     Input
      integer aer_nbr           ! [nbr] Dimension size
      real(r8) dmt_nma(aer_nbr)     ! [m] Number median particle size
      real(r8) gsd(aer_nbr)         ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_swa(aer_nbr)     ! [m] Effective (area weighted) mean particle size
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
c     Main Code
      do idx=1,aer_nbr
         ln_gsd=log(gsd(idx))
         dmt_swa(idx)=dmt_nma(idx)*exp(5.0*ln_gsd*ln_gsd/2.0) ! [m]
      enddo                     ! end loop over aer
      return
      end                       ! end nma2swa()
      
      subroutine sma2nma(aer_nbr,dmt_sma,gsd, ! I
     $     dmt_nma)             ! O
c     Purpose: Given the surface median particle sizes and geometric standard deviations
c     for a number of lognormal size distributions, compute and return the number median size of each distribution 
c     Definition of gsd conforms to SeP97

      use precision

      implicit none
c     Input
      integer aer_nbr           ! [nbr] Dimension size
      real(r8) dmt_sma(aer_nbr)     ! [m] Surface median particle size
      real(r8) gsd(aer_nbr)         ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_nma(aer_nbr)     ! [m] Number median particle size
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
c     Main Code
      do idx=1,aer_nbr
         ln_gsd=log(gsd(idx))
         dmt_nma(idx)=dmt_sma(idx)*exp(-2.0*ln_gsd*ln_gsd) ! [m]
      enddo                     ! end loop over aer
      return
      end                       ! end sma2nma()
      
      subroutine vma2nma(aer_nbr,dmt_vma,gsd, ! I
     $     dmt_nma)             ! O
c     Purpose: Given the mass median particle sizes and geometric standard deviations
c     for a number of lognormal size distributions, compute and return the number median size of each distribution 
c     Definition of gsd conforms to SeP97

      use precision

      implicit none
c     Input
      integer aer_nbr           ! [nbr] Dimension size
      real(r8) dmt_vma(aer_nbr)     ! [m] Mass median particle size
      real(r8) gsd(aer_nbr)         ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_nma(aer_nbr)     ! [m] Number median particle size
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
c     Main Code
      do idx=1,aer_nbr
         ln_gsd=log(gsd(idx))
         dmt_nma(idx)=dmt_vma(idx)*exp(-3.0*ln_gsd*ln_gsd) ! [m]
      enddo                     ! end loop over aer
      return
      end                       ! end vma2nma()
      
      subroutine nma2sma(aer_nbr,dmt_nma,gsd, ! I
     $     dmt_sma)             ! O
c     Purpose: Given the number median particle sizes and geometric standard deviations
c     for a number of lognormal size distributions, compute and return the surface median size of each distribution
c     Definition of gsd conforms to SeP97

      use precision

      implicit none
c     Input
      integer aer_nbr           ! [nbr] Dimension size
      real(r8) dmt_nma(aer_nbr)     ! [m] Number median particle size
      real(r8) gsd(aer_nbr)         ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_sma(aer_nbr)     ! [m] Surface median particle size
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
c     Main Code
      do idx=1,aer_nbr
         ln_gsd=log(gsd(idx))
         dmt_sma(idx)=dmt_nma(idx)*exp(2.0*ln_gsd*ln_gsd) ! [m]
      enddo                     ! end loop over aer
      return
      end                       ! end nma2sma()
      
      subroutine saa2nma(aer_nbr,dmt_das,gsd, ! I
     $     dmt_nma)             ! O
c     Purpose: Given the surface average diameter and geometric standard deviations
c     for a number of lognormal size distributions, compute and return the number median size of each distribution
c     Definition of gsd conforms to SeP97

      use precision

      implicit none
c     Input
      integer aer_nbr           ! [nbr] Dimension size
      real(r8) dmt_das(aer_nbr)     ! [m] Surface average particle size
      real(r8) gsd(aer_nbr)         ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_nma(aer_nbr)     ! [m] Number median particle size
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
c     Main Code
      do idx=1,aer_nbr
         ln_gsd=log(gsd(idx))
         dmt_nma(idx)=dmt_das(idx)*exp(-ln_gsd*ln_gsd) ! [m]
      enddo                     ! end loop over aer
      return
      end                       ! end saa2nma()
      
      subroutine vaa2nma(aer_nbr,dmt_dam,gsd, ! I
     $     dmt_nma)             ! O
c     Purpose: Given the mass average diameter and geometric standard deviations
c     for a number of lognormal size distributions, compute and return the number median size of each distribution
c     Definition of gsd conforms to SeP97

      use precision

      implicit none
c     Input
      integer aer_nbr           ! [nbr] Dimension size
      real(r8) dmt_dam(aer_nbr)     ! [m] Mass average particle size
      real(r8) gsd(aer_nbr)         ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_nma(aer_nbr)     ! [m] Number median particle size
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
c     Main Code
      do idx=1,aer_nbr
         ln_gsd=log(gsd(idx))
         dmt_nma(idx)=dmt_dam(idx)*exp(-3.0*ln_gsd*ln_gsd/2.0) ! [m]
      enddo                     ! end loop over aer
      return
      end                       ! end vaa2nma()
      
      subroutine nma2vma(aer_nbr,dmt_nma,gsd, ! I
     $     dmt_vma)             ! O
c     Purpose: Given the number median particle sizes and geometric standard deviations
c     for a number of lognormal size distributions, compute and return the mass median (same as volume-median) size of each distribution 
c     Definition of gsd conforms to SeP97

      use precision

      implicit none
c     Input
      integer aer_nbr           ! [nbr] Dimension size
      real(r8) dmt_nma(aer_nbr)     ! [m] Number median particle size
      real(r8) gsd(aer_nbr)         ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_vma(aer_nbr)     ! [m] Mass median particle size
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
c     Main Code
      do idx=1,aer_nbr
         ln_gsd=log(gsd(idx))
         dmt_vma(idx)=dmt_nma(idx)*exp(3.0*ln_gsd*ln_gsd) ! [m]
      enddo                     ! end loop over aer
      return
      end                       ! end nma2vma()
      
      subroutine nma2vwa(aer_nbr,dmt_nma,gsd, ! I
     $     dmt_vwa)             ! O
c     Purpose: Given the number median particle sizes and geometric standard deviations
c     for a number of lognormal size distributions, compute and return the mass-weighted (same as volume-weighted) mean size of each distribution 
c     Definition of gsd conforms to SeP97

      use precision

      implicit none
c     Input
      integer aer_nbr           ! [nbr] Dimension size
      real(r8) dmt_nma(aer_nbr)     ! [m] Number median particle size
      real(r8) gsd(aer_nbr)         ! [frc] Geometric standard deviation
c     Output
      real(r8) dmt_vwa(aer_nbr)     ! [m] Mass-weighted mean particle size
c     Local
      integer idx               ! [idx] Counting index
      real(r8) ln_gsd               ! ln(gsd)
c     Main Code
      do idx=1,aer_nbr
         ln_gsd=log(gsd(idx))
         dmt_vwa(idx)=dmt_nma(idx)*exp(7.0*ln_gsd*ln_gsd/2.0) ! [m]
      enddo                     ! end loop over aer
      return
      end                       ! end nma2vwa()
      
      subroutine ovr_src_snk_frc_get(src_nbr,mdn_src,gsd_src, ! I
     $     snk_nbr,dmt_min_snk,dmt_max_snk, ! I
     $     ovr_src_snk_frc)         ! O
c     Purpose: Given one set (the "source") of lognormal distributions, 
c     and one set of bin boundaries (the "sink"), compute and return 
c     the overlap factors between the source distributions and the sink bins.
c     The first group, src, contains src_nbr distributions
c     Each src distribution consists of a specified median diameter and geometric standard deviation
c     The second group, snk, contains snk_nbr bins 
c     Each sink bin consists of a specified minimum and maximum diameter
c     The output is a matrix ovr_src_snk_frc(src_nbr,snk_nbr)
c     Element ovr_src_snk_frc(i,j) is the fraction of size distribution i in group src 
c     that overlaps sink bin j
c     If the i'th source is completely bracketed by sink bins, then
c     \sum_{j=1}^{j=snk_nbr} ovr_src_snk_frc(i,j) = 1
c     If the j'th sink is completely bounded by the sources, then
c     \sum_{i=1}^{i=src_nbr} ovr_src_snk_frc(i,j) = 1
c     Note that input are generic median diameters of lognormal distributions
c     Thus if routine is called with number median diameters it will compute
c     the overlap factors for number concentration, and if it is called with mass
c     median diameters then it will compute the overlap factors for mass.
c     When ovr_src_snk_frc_get() is called with mass distributions parameters as arguments,
c     then ovr_src_snk_frc(src_idx,*) should be multiplied by mss_frc(src_idx) to obtain
c     the absolute mass fraction mapping from source bins to sink bins.

      use precision

      implicit none
c     Input
      integer snk_nbr           ! [nbr] Dimension size
      integer src_nbr           ! [nbr] Dimension size
      real(r8) dmt_max_snk(snk_nbr) ! [m] Maximum diameter in bin
      real(r8) dmt_min_snk(snk_nbr) ! [m] Minimum diameter in bin
      real(r8) gsd_src(src_nbr)     ! [frc] Geometric standard deviation
      real(r8) mdn_src(src_nbr)     ! [m] Mass median particle size
c     Output
      real(r8) ovr_src_snk_frc(src_nbr,snk_nbr) ! [frc] Fractional overlap of src with snk
c     Local
      integer src_idx           ! [idx] Counting index for src
      integer snk_idx           ! [idx] Counting index for snk
      real(r8) ln_gsd               ! ln(gsd)
      real(r8) sqrt2lngsdi          ! [frc] Factor in erf argument
      real(r8) lndmaxjovrdmdni      ! [frc] Factor in erf argument
      real(r8) lndminjovrdmdni      ! [frc] Factor in erf argument
c     External
      real(r8) erf                  ! Error function (libspecfun)
      external erf              ! Error function (libspecfun)
c     Main Code
c     Sanity check
c     19990913: erf() in SGI /usr/lib64/mips4/libftn.so is bogus
      if (abs(0.8427-erf(1.0_r8))/0.8427 > 0.001) then
         write (6,'(a,f12.10)') 'erf(1.0) = ',erf(1.0_r8)
         stop 'dst: ovr_src_snk_frc_get() reports Error function error'
      endif                     ! endif
      if (erf(0.0_r8) /= 0.0) then
         write (6,'(a,f12.10)') 'erf(0.0) = ',erf(0.0_r8)
         stop 'dst: ovr_src_snk_frc_get() reports Error function error'
      endif                     ! endif
      do src_idx=1,src_nbr
         sqrt2lngsdi=sqrt(2.0)*log(gsd_src(src_idx)) ! [frc]
         do snk_idx=1,snk_nbr
            lndmaxjovrdmdni=log(dmt_max_snk(snk_idx)/mdn_src(src_idx)) ! [frc]
            lndminjovrdmdni=log(dmt_min_snk(snk_idx)/mdn_src(src_idx)) ! [frc]
            ovr_src_snk_frc(src_idx,snk_idx)= ! [frc]
     $           0.5*(erf(lndmaxjovrdmdni/sqrt2lngsdi)-
     $           erf(lndminjovrdmdni/sqrt2lngsdi))
         enddo                  ! end loop over snk
      enddo                     ! end loop over src

      return
      end                       ! end ovr_src_snk_frc_get()
