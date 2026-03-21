c     $Header$

c     Purpose: library of extrapolation routines

c     Usage: 
c     cd /home/zender/csz; make -W xtr.F OPTS=D libcsz; cd -
c     cd /home/zender/csz; make -W xtr.F libcsz; cd -
c     cd /home/zender/csz; make OPTS=D libcsz; cd -

      subroutine xtr_ini(xtr_typ_LHS,xtr_typ_RHS)
c     Sets extrapolation flags and strings for given extrapolation types
      implicit none
c     Commons
#include <dbg.com> /* Debugging constants */
#include <xtr.com> /* Extrapolation constants */
c     Input
      integer xtr_typ_LHS
      integer xtr_typ_RHS
c     Local
      character*9 xtr_msk_sng(0:xtr_msk_nbr-1) 
      integer idx
      integer itmp
      integer rmn               ! Remainder
      integer mlt               ! Multiple (greatest factor s.t. mlt*msk<=rmn)
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)

c     Main Code
      if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering xtr_ini()'

      if (xtr_typ_LHS /= xtr_typ_RHS) then
         write (6,'(a,i4,a,i4)')
     $        'xtr_ini(): WARNING xtr_typ_LHS = ',xtr_typ_LHS,'  /= xtr_typ_RHS = ',xtr_typ_RHS
      endif

c     Initialize strings
      do idx=0,xtr_msk_nbr-1
         call strini(xtr_msk_sng(idx))
      enddo                     ! end loop over sng
      call f_strcpy(xtr_msk_sng(xtr_prt_msk),'prt ')
      call f_strcpy(xtr_msk_sng(xtr_fll_msk),'fll ')
      call f_strcpy(xtr_msk_sng(xtr_fll_nil_msk),'fll_nil ')
      call f_strcpy(xtr_msk_sng(xtr_fll_ngh_msk),'fll_ngh ')
      call f_strcpy(xtr_msk_sng(xtr_fll_lnr_msk),'fll_lnr ')
      call f_strcpy(xtr_msk_sng(xtr_prt_nil_msk),'prt_nil ')
      call f_strcpy(xtr_msk_sng(xtr_prt_ngh_msk),'prt_ngh ')
      call f_strcpy(xtr_msk_sng(xtr_prt_lnr_msk),'prt_lnr ')
      call f_strcpy(xtr_msk_sng(xtr_prt_frc_msk),'prt_frc ')
      call f_strcpy(xtr_msk_sng(xtr_prt_wgt_msk),'prt_wgt ')
      call f_strcpy(xtr_msk_sng(xtr_vrb_msk),'vrb ')

      call strini(xtr_sng_LHS)
      call strini(xtr_sng_RHS)

c     Proceed in reverse order by mask size
      rmn=xtr_typ_LHS           ! Initialize remainder to original flag
      do idx=xtr_msk_nbr-1,0,-1
         xtr_flg_LHS(idx)=.false.
         itmp=2**idx
         mlt=rmn/itmp
         rmn=mod(rmn,itmp)      ! Next iteration remainder is what is leftover
         if (mlt == 1) xtr_flg_LHS(idx)=.true.
      enddo                     ! end loop over msk

c     Implement and check logical consequences of flags
      if (xtr_flg_LHS(xtr_fll_nil_msk).or.xtr_flg_LHS(xtr_fll_ngh_msk).or.xtr_flg_LHS(xtr_fll_lnr_msk)) then
         xtr_flg_LHS(xtr_fll_msk)=.true.
      else
         xtr_flg_LHS(xtr_fll_msk)=.false.
      endif                     ! endif
      if (xtr_flg_LHS(xtr_prt_nil_msk).or.xtr_flg_LHS(xtr_prt_ngh_msk).or.
     $     xtr_flg_LHS(xtr_prt_lnr_msk).or.xtr_flg_LHS(xtr_prt_frc_msk).or.
     $     xtr_flg_LHS(xtr_prt_wgt_msk)) then
         xtr_flg_LHS(xtr_prt_msk)=.true.
      else
         xtr_flg_LHS(xtr_prt_msk)=.false.
      endif                     ! endif
      if (xtr_flg_LHS(xtr_prt_wgt_msk)) xtr_flg_LHS(xtr_prt_frc_msk)=.true.

c     Assign strings
      do idx=0,xtr_msk_nbr-1
         if (xtr_flg_LHS(idx)) call f_strcat(xtr_sng_LHS,xtr_msk_sng(idx)) 
      enddo                     ! end loop over msk

c     Proceed in reverse order by mask size
      rmn=xtr_typ_RHS           ! Initialize remainder to original flag
      do idx=xtr_msk_nbr-1,0,-1
         xtr_flg_RHS(idx)=.false.
         itmp=2**idx
         mlt=rmn/itmp
         rmn=mod(rmn,itmp)      ! Next iteration remainder is what is leftover
         if (mlt == 1) xtr_flg_RHS(idx)=.true.
      enddo                     ! end loop over msk

c     Implement and check logical consequences of flags
      if (xtr_flg_RHS(xtr_fll_nil_msk).or.xtr_flg_RHS(xtr_fll_ngh_msk).or.xtr_flg_RHS(xtr_fll_lnr_msk)) then
         xtr_flg_RHS(xtr_fll_msk)=.true.
      else
         xtr_flg_RHS(xtr_fll_msk)=.false.
      endif                     ! endif
      if (xtr_flg_RHS(xtr_prt_nil_msk).or.xtr_flg_RHS(xtr_prt_ngh_msk).or.
     $     xtr_flg_RHS(xtr_prt_lnr_msk).or.xtr_flg_RHS(xtr_prt_frc_msk).or.
     $     xtr_flg_RHS(xtr_prt_wgt_msk)) then
         xtr_flg_RHS(xtr_prt_msk)=.true.
      else
         xtr_flg_RHS(xtr_prt_msk)=.false.
      endif                     ! endif
      if (xtr_flg_RHS(xtr_prt_wgt_msk)) xtr_flg_RHS(xtr_prt_frc_msk)=.true.

c     Assign strings
      do idx=0,xtr_msk_nbr-1
         if (xtr_flg_RHS(idx)) call f_strcat(xtr_sng_RHS,xtr_msk_sng(idx)) 
      enddo                     ! end loop over msk

      if (dbg_lvl >= dbg_io) then
         write (6,'(a,i4,a,a)') 'xtr_typ_LHS = ',xtr_typ_LHS,' = ',xtr_sng_LHS(1:f_strlen(xtr_sng_LHS))
         write (6,'(a,i4,a,a)') 'xtr_typ_RHS = ',xtr_typ_RHS,' = ',xtr_sng_RHS(1:f_strlen(xtr_sng_RHS))
      endif                     ! endif dbg
      if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting xtr_ini()'

      return
      end                       ! end xtr_ini()


