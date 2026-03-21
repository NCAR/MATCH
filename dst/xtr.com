c     $Header$ -*-fortran-*- 

c     Purpose: Common block xtr stores extrapolation flags and strings

c     Usage: 
c     The values of the vairables in these common blocks are set with routine xtr_cmn_ini() in libcsz

c#include <xtr.com> /* Extrapolation constants */

c     Extrapolation types are as follows:
      integer xtr_msk_nbr
      parameter(xtr_msk_nbr=11)

      common / xtr /
     $     xtr_flg_LHS,         ! Initialized in xtr_ini()
     $     xtr_flg_RHS,         ! Initialized in xtr_ini()
     $     xtr_sng_LHS,         ! Initialized in xtr_ini()
     $     xtr_sng_RHS          ! Initialized in xtr_ini()
      logical xtr_flg_LHS(0:xtr_msk_nbr-1) ! Initialized in xtr_ini()
      logical xtr_flg_RHS(0:xtr_msk_nbr-1) ! Initialized in xtr_ini()
      character*80 xtr_sng_LHS  ! Initialized in xtr_ini()
      character*80 xtr_sng_RHS  ! Initialized in xtr_ini()

      integer xtr_prt_msk       ! 
      integer xtr_fll_msk       ! Perhaps fll should imply prt
      integer xtr_fll_nil_msk   ! nil ngh lnr are mutually exclusive
      integer xtr_fll_ngh_msk   ! nil ngh lnr are mutually exclusive
      integer xtr_fll_lnr_msk   ! nil ngh lnr are mutually exclusive
      integer xtr_prt_nil_msk   ! nil ngh lnr are mutually exclusive
      integer xtr_prt_ngh_msk   ! nil ngh lnr are mutually exclusive
      integer xtr_prt_lnr_msk   ! nil ngh lnr are mutually exclusive
      integer xtr_prt_frc_msk   ! nil ngh lnr are mutually exclusive, frc implies prt 
      integer xtr_prt_wgt_msk   ! wgt implies frc
      integer xtr_vrb_msk       !

c     Keys must be unique (like C enums)
      parameter(xtr_prt_msk=0)  ! Partial extrapolation is allowed 
      parameter(xtr_fll_msk=1)  ! Full extrapolation is allowed (implies xtr_prt_msk)
      parameter(xtr_fll_nil_msk=2) ! Set extrapolated value to 0.0
      parameter(xtr_fll_ngh_msk=3) ! Set extrapolated value to value of nearest valid neighbor
      parameter(xtr_fll_lnr_msk=4) ! Perform linear extrapolation using two nearest valid neighbors
      parameter(xtr_prt_nil_msk=5) ! Set extrapolated value to 0.0
      parameter(xtr_prt_ngh_msk=6) ! Set extrapolated value to value of nearest valid neighbor
      parameter(xtr_prt_lnr_msk=7) ! Perform linear extrapolation using two nearest valid neighbors
      parameter(xtr_prt_frc_msk=8) ! Use average value of overlap region
      parameter(xtr_prt_wgt_msk=9) ! Set extrapolated value to average value of overlap region weighted by size of overlap region plus 0.0 weighted by size of non overlap region (implies xtr_frc_msk)
      parameter(xtr_vrb_msk=10)  ! Print verbose warnings whenever extrapolation is performed
      
      integer xtr_err
      integer xtr_fll_nil       
      integer xtr_fll_ngh   
      integer xtr_fll_lnr       
      integer xtr_prt_nil       
      integer xtr_prt_ngh   
      integer xtr_prt_lnr       
      integer xtr_prt_frc       
      integer xtr_prt_wgt   
      integer xtr_vrb
      
c     Metaflags
      parameter(xtr_err=0)      ! No extrapolation allowed under any circumstances
      parameter(xtr_vrb=2**xtr_vrb_msk) ! Print circumstances of all extrapolations

c     Full extrapolation types. Selecting one of these sets xtr_fll_msk=.true.
      parameter(xtr_fll_nil=2**xtr_fll_nil_msk) ! Set extrapolated value to 0.0
      parameter(xtr_fll_ngh=2**xtr_fll_ngh_msk) ! Set extrapolated value to value of nearest valid neighbor
      parameter(xtr_fll_lnr=2**xtr_fll_lnr_msk) ! Linearly extrapolate values using values two nearest valid neighbors

c     Partial extrapolation types. Selecting one of these sets xtr_prt_msk=.true.
      parameter(xtr_prt_nil=2**xtr_prt_nil_msk) ! Set extrapolated value to 0.0
      parameter(xtr_prt_lnr=2**xtr_prt_lnr_msk) ! Linearly extrapolate value
      parameter(xtr_prt_ngh=2**xtr_prt_ngh_msk) ! Set extrapolated value to value of nearest valid neighbor
      parameter(xtr_prt_frc=2**xtr_prt_frc_msk) ! Compute average in overlap region. Use xtr_prt_wgt_msk to determine how to weight final value.
      parameter(xtr_prt_wgt=2**xtr_prt_wgt_msk) ! Set extrapolated value to fraction in valid region weighted by size of valid region plus 0.0 weighted by size of non-overlap region



