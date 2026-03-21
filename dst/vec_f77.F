c     $Header$

c     Purpose: library of Fortran vector manipulation routines

c     Usage: 
c     cd $HOME/f; make -W vec.F OPTS=D libcsz; cd -
c     cd $HOME/f; make -W vec.F libcsz; cd -
c     cd $HOME/f; make OPTS=D libcsz; cd -

      subroutine ntp_vec(in_nbr,crd_in,dat_in,out_nbr,crd_out,dat_out,xtr_typ_LHS,xtr_typ_RHS)
c     Purpose: Project a vector onto another vector.
c     Input and output coordinate vectors must be monotonic.

c     crd_in defines the coordinate grid associated with the values dat_in.
c     crd_out defines a new (user) coordinate grid onto which dat_in is linearly interpolated. 

      use precision

      implicit none
c     Commons
#include <dbg.com> /* Debugging constants */
#include <xtr.com> /* Extrapolation constants */
c     Input
      integer in_nbr
      integer out_nbr
      integer xtr_typ_LHS
      integer xtr_typ_RHS
      real(r8) crd_in(in_nbr)       ! 
      real(r8) crd_out(out_nbr)     ! 
      real(r8) dat_in(in_nbr)       ! 
c     Output
      real(r8) dat_out(out_nbr)     ! 
c     Local
      integer brk_lft_idx
      integer brk_rgt_idx
      integer out_idx
      logical in_dcr
      logical in_ncr
      logical out_dcr
      logical out_ncr
c     External
      integer f_strlen            ! f_strlen() is stored in libcsz

c     Main Code
      if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering ntp_vec()'

c     Set extrapolation flags and strings
      call xtr_ini(xtr_typ_LHS,xtr_typ_RHS)

c     Vet input
      if (in_nbr <= 1) stop 'nbr_dat_in <= 1 in ntp_vec()'
      if (out_nbr < 1) stop 'nbr_dat_out < 1 in ntp_vec()'

c     Find direction of monotonicity
      if (crd_in(2)-crd_in(1) > 0.0) then 
         in_ncr=.true. 
      else 
         in_ncr=.false.
      endif                     ! endif
      in_dcr=.not.in_ncr
      out_ncr=.true. 
      if (out_nbr > 1) then 
         if (crd_out(2)-crd_out(1) < 0.0) out_ncr=.false.
      endif                     ! endif
      out_dcr=.not.out_ncr
      
c     Convert both input and output arrays to monotonic increasing arrays
      if (in_dcr) then
         if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Reversing input grid in ntp_vec()'
         call rvr_vec(crd_in,in_nbr)
         call rvr_vec(dat_in,in_nbr)
      endif                     ! endif in_dcr
      if (out_dcr) then
         if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Reversing output grid in ntp_vec()'
         call rvr_vec(crd_out,out_nbr)
      endif                     ! endif in_dcr

c     Initialize bracketing index
      brk_lft_idx=1
c     Loop over desired output coordinates
      do out_idx=1,out_nbr
c     19990428: Following line causes index out of range error because second condition is tested even though first condition is violated. This may be ignored but is annoying
         do while((brk_lft_idx <= in_nbr).and.(crd_in(brk_lft_idx) < crd_out(out_idx)))
            brk_lft_idx=brk_lft_idx+1
         end do                 ! end while
         brk_lft_idx=brk_lft_idx-1
c     Handle identity interpolation separately to preserve symmetry in extrapolation code 
         if (crd_in(brk_lft_idx+1) == crd_out(out_idx)) then
            dat_out(out_idx)=dat_in(brk_lft_idx+1)
            if (brk_lft_idx == 0) brk_lft_idx=1 ! Reset brk_lft_idx to 1 so next while loop works
c     NB: In C or F90, use "continue" or "cycle", respectively, instead of "goto" 
            goto 100            ! Jump to next iteration
         endif                  ! endif
         if (brk_lft_idx == 0) then
c     LHS Extrapolation required
c     Degenerate case: crd_out(out_idx) < crd_in(1)
            brk_lft_idx=1       ! Reset brk_lft_idx to 1 so next while loop works
            if (xtr_flg_LHS(xtr_vrb_msk)) then
               write (6,'(1p,a,i5,a,i5,a,e10.3,a,i5,a,e10.3,a,i5,a,e10.3,0p)') 
     $              'WARNING: ntp_vec() output value dat_out(',
     $              out_idx,') at coordinate crd_out(',
     $              out_idx,') = ',
     $              crd_out(out_idx),' requires LHS extrapolation beyond leftmost valid coordinate at crd_in(',
     $              brk_lft_idx,') = ',
     $              crd_in(brk_lft_idx),'. Nearest valid datum is dat_in(',
     $              brk_lft_idx,') = ',
     $              dat_in(brk_lft_idx)
            endif               ! endif
c     Extrapolation options are presented in decreasing order of preference
            if (.not.xtr_flg_LHS(xtr_fll_msk)) then
               stop 'ERROR: Full LHS extrapolation required but not permitted in ntp_vec()'
            else if (xtr_flg_LHS(xtr_fll_nil_msk)) then
               dat_out(out_idx)=0.0
            else if (xtr_flg_LHS(xtr_fll_ngh_msk)) then
               dat_out(out_idx)=dat_in(1)
            else if (xtr_flg_LHS(xtr_fll_lnr_msk)) then
               dat_out(out_idx)=dat_in(1)-
     $              (crd_in(1)-crd_out(out_idx))*
     $              (dat_in(2)-dat_in(1))/(crd_in(2)-crd_in(1))
            else
               stop 'Unknown xtr_typ_LHS in ntp_vec()'
            endif               ! endif xtr_typ_LHS
            if (xtr_flg_LHS(xtr_vrb_msk)) write (6,'(a,a,i4,a,1p,e10.3,0p)') 
     $           xtr_sng_LHS(1:f_strlen(xtr_sng_LHS)),' yields dat_out(',out_idx,') = ',dat_out(out_idx)
         else if (brk_lft_idx < in_nbr) then
c     Normal case: crd_out is interpolable
            brk_rgt_idx=brk_lft_idx+1
c     Linearly interpolate
            dat_out(out_idx)=
     $           dat_in(brk_lft_idx)+
     $           (crd_out(out_idx)-crd_in(brk_lft_idx))*
     $           (dat_in(brk_rgt_idx)-dat_in(brk_lft_idx))/
     $           (crd_in(brk_rgt_idx)-crd_in(brk_lft_idx))
         else if (brk_lft_idx == in_nbr) then
c     RHS Extrapolation required
c     Degenerate case: brk_lft_idx is last element of crd_in 
            brk_rgt_idx=brk_lft_idx 
            if (xtr_flg_RHS(xtr_vrb_msk)) then
               write (6,'(1p,a,i5,a,i5,a,e10.3,a,i5,a,e10.3,a,i5,a,e10.3,0p)') 
     $              'WARNING: ntp_vec() output value dat_out(',
     $              out_idx,') at coordinate crd_out(',
     $              out_idx,') = ',
     $              crd_out(out_idx),' requires RHS extrapolation beyond rightmost valid coordinate at crd_in(',
     $              brk_rgt_idx,') = ',
     $              crd_in(brk_rgt_idx),'. Nearest valid datum is dat_in(',
     $              brk_rgt_idx,') = ',
     $              dat_in(brk_rgt_idx)
            endif               ! endif
c     Extrapolation options are presented in decreasing order of preference
            if (.not.xtr_flg_RHS(xtr_fll_msk)) then
               stop 'ERROR: Full RHS extrapolation required but not permitted in ntp_vec()'
            else if (xtr_flg_RHS(xtr_fll_nil_msk)) then
               dat_out(out_idx)=0.0
            else if (xtr_flg_RHS(xtr_fll_ngh_msk)) then
               dat_out(out_idx)=dat_in(in_nbr)
            else if (xtr_flg_RHS(xtr_fll_lnr_msk)) then
               dat_out(out_idx)=dat_in(in_nbr)+
     $              (crd_out(out_idx)-crd_in(in_nbr))*
     $              (dat_in(in_nbr)-dat_in(in_nbr-1))/
     $              (crd_in(in_nbr)-crd_in(in_nbr-1))
            else
               stop 'Unknown xtr_typ_RHS in ntp_vec()'
            endif               ! endif xtr_typ_RHS
            if (xtr_flg_RHS(xtr_vrb_msk)) write (6,'(a,a,i4,a,1p,e10.3,0p)') 
     $           xtr_sng_RHS(1:f_strlen(xtr_sng_RHS)),' yields dat_out(',out_idx,') = ',dat_out(out_idx)
         else
            stop 'ERROR: Unforeseen value of brk_lft_idx in ntp_vec()'
         endif
 100     continue               ! Branch from identity interpolation
      enddo                     ! end loop over output coordinates
      
c     Convert input and output arrays to original direction of monotonicity
      if (in_dcr) then
         if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Un-reversing input grid in ntp_vec()'
         call rvr_vec(crd_in,in_nbr)
         call rvr_vec(dat_in,in_nbr)
      endif                     ! endif in_dcr
      if (out_dcr) then
         if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Un-reversing output grid in ntp_vec()'
         call rvr_vec(crd_out,out_nbr)
         call rvr_vec(dat_out,out_nbr)
      endif                     ! endif in_dcr

      if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting ntp_vec()'
      return
      end                       ! end ntp_vec()
      
      subroutine rbn_vec(in_nbr,grd_in,dat_in, 
     $     out_nbr,grd_out,dat_out,
     $     xtr_typ_LHS,xtr_typ_RHS)
c     Purpose: Rebin a vector from one grid onto another grid
c     Input and output coordinate vectors must both be monotonic 
c     Since they are each monotonic, they are each specified with a single array containing the interface coordinate values

c     crd_min_in and crd_max_in define the coordinate grid associated with the values dat_in
c     crd_min_out and crd_max_out define a new (user) coordinate grid onto which dat_in is interpolated

c     The output (interpolated) values are set so that the output data has the same integral properties as the input data
c     In plain English this means the interpolated values are the average values in some neighborhood of the specified grid point
c     Use this option to interpolate noisy input data, e.g., solar specta or NO2 absorption cross sections, to more regular grids

      use precision

      implicit none
c     Parameters
c     Commons
#include <dbg.com> /* Debugging constants */
#include <xtr.com> /* Extrapolation constants */
c     Input
      integer in_nbr            ! [nbr] max(in_nbr)=49934 for Kur95 1 cm-1 spectrum
      integer out_nbr           ! [nbr] max(out_nbr)=8192 for TGC98 H2OH2O data
      integer xtr_typ_LHS
      integer xtr_typ_RHS
      real(r8) grd_in(in_nbr+1)     ! Input grid
      real(r8) grd_out(out_nbr+1)   ! Output grid
      real(r8) dat_in(in_nbr)       ! Input data
c     Output
      real(r8) dat_out(out_nbr)     ! Input data rebinned to output grid
c     Local
      integer brk_lft_idx(out_nbr+1)
      integer brk_rgt_idx(out_nbr+1)
      integer grd_in_nbr
      integer grd_out_idx
      integer grd_out_nbr
      integer idx
      integer in_idx
      integer out_idx
      integer out_idxp1
      integer trp_idx
      integer trp_srt_idx
      integer trp_nbr
      logical mnt
      logical in_dcr
      logical in_ncr
      logical out_dcr
      logical out_ncr
      real(r8) crd_in(in_nbr)       ! Input coordinates (midpoints of grid)
      real(r8) crd_out(out_nbr)     ! Output coordinates (midpoints of grid)
      real(r8) dat_ntp_in(in_nbr+1) ! Linearly interpolated data on grd_in grid
      real(r8) dat_ntp_out(out_nbr+1) ! Linearly interpolated data on grd_out grid
      real(r8) trp_area(in_nbr) ! Area of each trapezoid defined by input data
      real(r8) crd_dlt_lft
      real(r8) crd_dlt_rgt
      real(r8) crd_dlt_ttl
      real(r8) ovl_avg              ! Overlap average
      real(r8) ovl_frc              ! Overlap fraction
      real(r8) trp_area_lft
      real(r8) trp_area_rgt
      real(r8) trp_area_ttl
c     External
      logical mnt_chk
      integer f_strlen            ! f_strlen() is stored in libcsz

c     Main Code
      if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering rbn_vec()'

c     Vet input
      if (in_nbr < 1) stop 'in_nbr < 1 in rbn_vec()'
      if (out_nbr < 1) stop 'out_nbr < 1 in rbn_vec()'

c     Initialize default values
      grd_in_nbr=in_nbr+1
      grd_out_nbr=out_nbr+1

c     Set extrapolation flags and strings
      call xtr_ini(xtr_typ_LHS,xtr_typ_RHS)
 
c     Check for monotonicity
      mnt=mnt_chk(grd_in,in_nbr+1)
      if (.not.mnt) stop 'grd_in not monotonic in rbn_vec()'
      mnt=mnt_chk(grd_out,out_nbr+1)
      if (.not.mnt) stop 'grd_out not monotonic in rbn_vec()'

c     Find direction of monotonicity
      if (grd_in(2)-grd_in(1) > 0.0) then 
         in_ncr=.true. 
      else 
         in_ncr=.false.
      endif
      in_dcr=.not.in_ncr
      if (grd_out(2)-grd_out(1) > 0.0) then 
         out_ncr=.true. 
      else 
         out_ncr=.false.
      endif
      out_dcr=.not.out_ncr

c     Convert both input and output arrays to monotonic increasing arrays
      if (in_dcr) then
         if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Reversing input grid in rbn_vec()'
         call rvr_vec(grd_in,grd_in_nbr)
         call rvr_vec(dat_in,in_nbr)
      endif                     ! endif in_dcr
      if (out_dcr) then
         if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Reversing output grid in rbn_vec()'
         call rvr_vec(grd_out,grd_out_nbr)
      endif                     ! endif in_dcr

c     Input and output coordinates both increase monotonically
      if (xtr_typ_LHS == xtr_err.and.xtr_typ_RHS == xtr_err) then
c     Vet input more stringently
         if (grd_out(grd_out_nbr) < grd_in(1)) then
            write (6,'(1p,a,i4,a,e10.3,a,e10.3,0p)') 
     $           'ERROR: rbn_vec() has grd_out(grd_out_nbr=',grd_out_nbr,
     $           ') = ',grd_out(grd_out_nbr),
     $           ' < grd_in(1) = ',grd_in(1)
            stop
         endif
         if (grd_out(1) > grd_in(grd_in_nbr)) then
            write (6,'(1p,a,e10.3,a,i4,a,e10.3,0p)') 
     $           'ERROR: rbn_vec() has grd_out(1) = ',grd_out(1),
     $           ' > grd_in(grd_in_nbr = ',grd_in_nbr,
     $           ') = ',grd_in(grd_in_nbr)
            stop
         endif
      endif                     ! endif xtr_err

c     Initialize coordinate arrays
      do in_idx=1,in_nbr
         crd_in(in_idx)=0.5*(grd_in(in_idx)+grd_in(in_idx+1))
      enddo                     ! end loop over input coordinates
      do out_idx=1,out_nbr
         crd_out(out_idx)=0.5*(grd_out(out_idx)+grd_out(out_idx+1))
      enddo                     ! end loop over output coordinates

c     Assemble LHS bracketing indices for each output gridpoint
      do grd_out_idx=1,grd_out_nbr
         brk_lft_idx(grd_out_idx)=1
c     Initialize bracketing index with previous value if possible
         if (grd_out_idx > 1) then
            if (brk_lft_idx(grd_out_idx-1) >= 1.and.
     $           brk_lft_idx(grd_out_idx-1) <= grd_in_nbr) 
     $           brk_lft_idx(grd_out_idx)=brk_lft_idx(grd_out_idx-1)
         endif                  ! endif
         do while((grd_in(brk_lft_idx(grd_out_idx)) < grd_out(grd_out_idx)).and.
     $        (brk_lft_idx(grd_out_idx) <= grd_in_nbr))
            brk_lft_idx(grd_out_idx)=brk_lft_idx(grd_out_idx)+1
         end do                 ! end while
         if (grd_in(brk_lft_idx(grd_out_idx)) /= grd_out(grd_out_idx)) then
c     Decrement bracketing index unless while loop was broken by an equality
            brk_lft_idx(grd_out_idx)=brk_lft_idx(grd_out_idx)-1
         endif                  ! endif
      enddo                     ! end loop over output coordinates

c     Assemble RHS bracketing indices for each output gridpoint
      do grd_out_idx=grd_out_nbr,1,-1
         brk_rgt_idx(grd_out_idx)=grd_in_nbr
c     Initialize bracketing index with previous value if possible
         if (grd_out_idx < grd_out_nbr) then
            if (brk_rgt_idx(grd_out_idx+1) <= grd_in_nbr.and.
     $           brk_rgt_idx(grd_out_idx+1) >= 1) 
     $           brk_rgt_idx(grd_out_idx)=brk_rgt_idx(grd_out_idx+1)
         endif                  ! endif
         do while((grd_in(brk_rgt_idx(grd_out_idx)) > grd_out(grd_out_idx)).and.
     $        (brk_rgt_idx(grd_out_idx) >= 1))
            brk_rgt_idx(grd_out_idx)=brk_rgt_idx(grd_out_idx)-1
         end do                 ! end while
         if (grd_in(brk_rgt_idx(grd_out_idx)) /= grd_out(grd_out_idx)) then
c     Increment bracketing index unless while loop was broken by an equality
            brk_rgt_idx(grd_out_idx)=brk_rgt_idx(grd_out_idx)+1
         endif                  ! endif
      enddo                     ! end loop over output coordinates

c     Vet bracketing indices
      do grd_out_idx=1,grd_out_nbr
         if (brk_lft_idx(grd_out_idx) > brk_rgt_idx(grd_out_idx)) then
            stop 'brk_lft_idx(grd_out_idx) > brk_rgt_idx(grd_out_idx) in rbn_vec()'
         endif                  ! endif
      enddo                     ! end loop over output coordinates

c     Precompute all linearly interpolated values for use below
c     Input values interpolated to input grid (Recall values are input on coordinates, not grid)
      call ntp_vec(in_nbr,crd_in,dat_in,grd_in_nbr,grd_in,dat_ntp_in,xtr_fll_lnr,xtr_fll_lnr)
c     Input values interpolated to output grid
      call ntp_vec(in_nbr,crd_in,dat_in,grd_out_nbr,grd_out,dat_ntp_out,xtr_fll_lnr,xtr_fll_lnr)
c     Precompute area of all trapezoids defined by input grid for use below
      call trp_area_vec(in_nbr,grd_in,dat_ntp_in,trp_area)

c     Reset extrapolation flags and strings
      call xtr_ini(xtr_typ_LHS,xtr_typ_RHS)

c     Assign an output value to each output region
      do out_idx=1,out_nbr
         out_idxp1=out_idx+1
c     Traverse truth table by examining all possible values of brk_lft_idx(out_idx)
         if (brk_lft_idx(out_idx) < 0) then
            stop 'ERROR: brk_lft_idx(out_idx) < 0 in rbn_vec()'
         else if (brk_lft_idx(out_idx) == 0) then
c     LHS extrapolation required
            if (xtr_flg_LHS(xtr_vrb_msk)) then
               write (6,'(1p,a,i5,a,i5,a,e10.3,a,i5,a,e10.3,a,e10.3,a,e10.3,a,e10.3,0p)') 
     $              'WARNING: rbn_vec() output bin value centered at crd_out(',
     $              out_idx,') between interface coordinates grd_out(',
     $              out_idx,') = ',
     $              grd_out(out_idx),' and grd_out(',
     $              out_idxp1,') = ',
     $              grd_out(out_idxp1),
     $              ' requires LHS extrapolation beyond leftmost valid interface coordinates at grd_in(1) = ',
     $              grd_in(1),' and grd_in(2) = ',
     $              grd_in(2),'. Nearest valid datum is dat_in(1) = ',
     $              dat_in(1)
            endif               ! endif xtr_flg_vrb
            if (brk_lft_idx(out_idxp1) == 0) then
c     Fully degenerate case: current output bin has no overlap with input bins
               if (xtr_flg_LHS(xtr_vrb_msk)) then
                  write (6,'(1p,2(a,e10.3),0p)')
     $                 'Full LHS extrpolation required. Distances of bin endpoints from valid data are ',
     $                 grd_in(1)-grd_out(out_idx),' and ',grd_in(1)-grd_out(out_idxp1)
               endif            ! endif xtr_flg_vrb
c     Extrapolation options are presented in decreasing order of preference
               if (.not.xtr_flg_LHS(xtr_fll_msk)) then
                  stop 'ERROR: Full LHS extrapolation required but not permitted in rbn_vec()'
               else if (xtr_flg_LHS(xtr_fll_nil_msk)) then
                  dat_out(out_idx)=0.0
               else if (xtr_flg_LHS(xtr_fll_ngh_msk)) then
                  dat_out(out_idx)=dat_in(1)
               else if (xtr_flg_LHS(xtr_fll_lnr_msk)) then
                  stop 'xtr_fll_lnr not implemented yet in rbn_vec()'
               else
                  stop 'Unknown xtr_typ_LHS in fll branch of rbn_vec()'
               endif            ! endif xtr_typ_LHS
            else if (brk_lft_idx(out_idxp1) > 0) then ! brk_lft_idx(out_idxp1) == 0
c     Half degenerate case: current output bin has partial overlap with input bins
               ovl_frc=(grd_out(out_idxp1)-grd_in(1))/(grd_out(out_idxp1)-grd_out(out_idx))
               if (xtr_flg_LHS(xtr_vrb_msk)) then
                  write (6,'(a,f9.3,a)') 'Partial LHS extrapolation required. Overlap fraction is ',ovl_frc
               endif            ! endif xtr_flg_vrb
c     Extrapolation options are presented in decreasing order of preference
               if (.not.xtr_flg_LHS(xtr_prt_msk)) then
                  stop 'ERROR: Partial LHS extrapolation required but not permitted in rbn_vec()'
               else if (xtr_flg_LHS(xtr_prt_frc_msk)) then
c     Compute average of overlap input bins
                  crd_dlt_rgt=grd_out(out_idxp1)-grd_in(brk_lft_idx(out_idxp1))
                  trp_area_rgt=0.5*crd_dlt_rgt*(dat_ntp_out(out_idxp1)+dat_ntp_in(brk_lft_idx(out_idxp1)))
                  crd_dlt_ttl=crd_dlt_rgt
                  trp_area_ttl=trp_area_rgt
                  trp_nbr=brk_lft_idx(out_idxp1)-1
                  do trp_idx=1,trp_nbr
                     in_idx=trp_idx
                     trp_area_ttl=trp_area_ttl+trp_area(in_idx)
                     crd_dlt_ttl=crd_dlt_ttl+grd_in(in_idx+1)-grd_in(in_idx)
                  end do        ! end loop over inner trapezoids
                  if (crd_dlt_ttl <= 0.0) then 
                     stop 'ERROR: crd_dlt_ttl <= 0.0 in xtr_LHS branch of rbn_vec()'
                  endif         ! endif crd_dlt_ttl <= 0.0
                  ovl_avg=trp_area_ttl/crd_dlt_ttl 
                  if (xtr_flg_LHS(xtr_prt_wgt_msk)) then
                     dat_out(out_idx)=ovl_frc*ovl_avg
                  else
                     dat_out(out_idx)=ovl_avg
                  endif         ! endif wgt
               else if (xtr_flg_LHS(xtr_prt_lnr_msk)) then
                  stop 'xtr_prt_lnr not implemented yet in rbn_vec()'
               else if (xtr_flg_LHS(xtr_prt_ngh_msk)) then
                  dat_out(out_idx)=dat_in(1)
               else if (xtr_flg_LHS(xtr_prt_nil_msk)) then
                  dat_out(out_idx)=0.0
               else
                  stop 'Unknown xtr_typ_LHS in prt branch of rbn_vec()'
               endif            ! endif xtr_typ_LHS
            else                ! brk_lft_idx(out_idxp1) == 0
               write (6,'(a)') 'Unforeseen brk_lft_idx in xtr_typ_LHS in rbn_vec()'
            endif               ! brk_lft_idx(out_idxp1) == 0
            if (xtr_flg_LHS(xtr_vrb_msk)) write (6,'(a,a,i4,a,1p,e10.3,0p)') 
     $           xtr_sng_LHS(1:f_strlen(xtr_sng_LHS)),' yields dat_out(',out_idx,') = ',dat_out(out_idx)
         else if ((brk_lft_idx(out_idx) <= grd_in_nbr).and.
     $           (brk_rgt_idx(out_idxp1) <= grd_in_nbr)) then ! brk_lft_idx(out_idx) == 0
c     Normal case: 99% of execution time for large arrays should be in this block
c     Current and next output gridpoints are both interpolable
c     Find out how many input crds lay between current and next output gridpoints
            trp_nbr=brk_lft_idx(out_idxp1)-brk_rgt_idx(out_idx)
            if (trp_nbr < 0) then
c     Case A: Current and next bracketing indices are the same
c     This means the output crd resolution is finer than the input crd resolution
c     Use trapezoidal are current point (this is the same as using just a LHS trapezoid)
               dat_out(out_idx)=0.5*(dat_ntp_out(out_idx)+dat_ntp_out(out_idxp1))
            else if (trp_nbr <= in_nbr) then
c     Case B: Current and previous bracketing indices differ by one
c     This means the output crd resolution roughly equals the input crd resolution
c     This case can be handled by the following code since the number of inner trapezoids
c     is 0 and the loop will not be executed. 
c     The input/output area will just be the sum of the LHS and RHS trapezoids
               
c     Case C: Current and previous bracketing indices differ by many.
c     This means the output crd resolution is coarser than the input crd resolution.
c     First, figure out area under the RHS and LHS bookend trapezoids:
               crd_dlt_lft=grd_in(brk_rgt_idx(out_idx))-grd_out(out_idx)
               trp_area_lft=0.5*crd_dlt_lft*(dat_ntp_out(out_idx)+dat_ntp_in(brk_rgt_idx(out_idx)))
               crd_dlt_rgt=grd_out(out_idxp1)-grd_in(brk_lft_idx(out_idxp1))
               trp_area_rgt=0.5*crd_dlt_rgt*(dat_ntp_out(out_idxp1)+dat_ntp_in(brk_lft_idx(out_idxp1)))
c     NB: We use trapezoids to integrate. But we are unable to use the classical 
c     Trapezoidal rule because the abscissas are, in general, unevenly spaced.
c     Therefore we use the trapezoidal rule generalized to unevenly spaced abscissas.
c     If the input grid equals the output grid, there is one inner trapezoid.
               
c     Trapezoids are indexed with the Fortran (1-based) convention where
c     trapezoid 1 is bounded by the abscissae crd_in(1) and crd_in(2) and
c     trapezoid N-1 is bounded by the abscissae crd_in(N-1) and crd_in(N).
               crd_dlt_ttl=crd_dlt_lft+crd_dlt_rgt
               trp_area_ttl=trp_area_lft+trp_area_rgt
               if (brk_lft_idx(out_idx) /= brk_rgt_idx(out_idx)) then
                  trp_srt_idx=brk_lft_idx(out_idx)
               else
                  trp_srt_idx=brk_lft_idx(out_idx)-1
               endif            ! endif
               do trp_idx=1,trp_nbr
                  in_idx=trp_srt_idx+trp_idx
                  trp_area_ttl=trp_area_ttl+trp_area(in_idx)
                  crd_dlt_ttl=crd_dlt_ttl+grd_in(in_idx+1)-grd_in(in_idx)
               end do           ! end loop over inner trapezoids
c     This relationship results from requiring output rectangle to be equal in area to 
c     total of input trapezoids (between current and next gridpoints).
               if (crd_dlt_ttl <= 0.0) then 
                  write (6,'(a)') 'ERROR: crd_dlt_ttl <= 0.0 in main branch of rbn_vec()'
               else
                  dat_out(out_idx)=trp_area_ttl/crd_dlt_ttl 
               endif            ! endif crd_dlt_ttl <= 0.0
            else if (trp_nbr > in_nbr) then
               write (6,'(2(a,i5,a,i5))') 
     $              'brk_lft_idx(',out_idx,
     $              ') = ',brk_lft_idx(out_idx),
     $              ', brk_rgt_idx(',out_idx,
     $              ') = ',brk_rgt_idx(out_idx)
               write (6,'(2(a,i5,a,i5))') 
     $              'brk_lft_idx(',out_idxp1,
     $              ') = ',brk_lft_idx(out_idxp1),
     $              ', brk_rgt_idx(',out_idxp1,
     $              ') = ',brk_rgt_idx(out_idxp1)
               write (6,'(a,i5)') 'trp_nbr = ',trp_nbr
               stop 'Unforeseen trp_nbr in rbn_vec()' 
            endif               ! endif trp_nbr == 0
         else if (brk_rgt_idx(out_idxp1) > grd_in_nbr) then ! brk_lft_idx(out_idx) == 0
c     RHS extrapolation required
            if (xtr_flg_RHS(xtr_vrb_msk)) then
               write (6,'(1p,a,i5,a,i5,a,e10.3,a,i5,a,e10.3,a,i5,a,e10.3,a,i5,a,e10.3,a,i5,a,e10.3,0p)') 
     $              'WARNING: rbn_vec() output bin value centered at crd_out(',
     $              out_idx,') between interface coordinates grd_out(',
     $              out_idx,') = ',
     $              grd_out(out_idx),' and grd_out(',
     $              out_idxp1,') = ',
     $              grd_out(out_idxp1),
     $              ' requires RHS extrapolation beyond rightmost valid interface coordinates at grd_in(',
     $              grd_in_nbr-1,') = ',
     $              grd_in(grd_in_nbr-1),' and grd_in(',
     $              grd_in_nbr,') = ',
     $              grd_in(grd_in_nbr),'. Nearest valid datum is dat_in(',
     $              in_nbr,') = ',
     $              dat_in(in_nbr)
            endif               ! endif xtr_flg_vrb
            if (brk_rgt_idx(out_idx) == grd_in_nbr+1) then
c     Fully degenerate case: current output bin has no overlap with input bins
               if (xtr_flg_RHS(xtr_vrb_msk)) then
                  write (6,'(1p,2(a,e10.3),0p)')
     $                 'Full RHS extrpolation required. Distances of bin endpoints from valid data are ',
     $                 grd_out(out_idx)-grd_in(grd_in_nbr),' and ',grd_out(out_idxp1)-grd_in(grd_in_nbr)
               endif            ! endif xtr_flg_vrb
c     Extrapolation options are presented in decreasing order of preference
               if (.not.xtr_flg_RHS(xtr_fll_msk)) then
                  stop 'ERROR: Full RHS extrapolation required but not permitted in rbn_vec()'
               else if (xtr_flg_RHS(xtr_fll_nil_msk)) then
                  dat_out(out_idx)=0.0
               else if (xtr_flg_RHS(xtr_fll_ngh_msk)) then
                  dat_out(out_idx)=dat_in(in_nbr)
               else if (xtr_flg_RHS(xtr_fll_lnr_msk)) then
                  stop 'xtr_fll_lnr not implemented yet in rbn_vec()'
               else
                  stop 'Unknown xtr_typ_RHS in rbn_vec()'
               endif            ! endif xtr_typ_RHS
            else                ! brk_rgt_idx(out_idx) == grd_in_nbr+1
c     Half degenerate case: current output bin has partial overlap with input bins
               ovl_frc=(grd_in(grd_in_nbr)-grd_out(out_idx))/(grd_out(out_idxp1)-grd_out(out_idx))
               if (xtr_flg_RHS(xtr_vrb_msk)) then
                  write (6,'(a,f9.3,a)') 'Partial RHS extrapolation required. Overlap fraction is ',ovl_frc
               endif            ! endif xtr_flg_vrb
c     Extrapolation options are presented in decreasing order of preference
               if (.not.xtr_flg_RHS(xtr_prt_msk)) then
                  stop 'ERROR: Partial RHS extrapolation required but not permitted in rbn_vec()'
               else if (xtr_flg_RHS(xtr_prt_frc_msk)) then
c     Compute average of overlap input bins
                  crd_dlt_lft=grd_in(brk_rgt_idx(out_idx))-grd_out(out_idx)
                  trp_area_lft=0.5*crd_dlt_lft*(dat_ntp_out(out_idx)+dat_ntp_in(brk_rgt_idx(out_idx)))
                  crd_dlt_ttl=crd_dlt_lft
                  trp_area_ttl=trp_area_lft
                  trp_nbr=grd_in_nbr-brk_rgt_idx(out_idx)
                  do trp_idx=1,trp_nbr
                     in_idx=brk_lft_idx(out_idx)+trp_idx
                     trp_area_ttl=trp_area_ttl+trp_area(in_idx)
                     crd_dlt_ttl=crd_dlt_ttl+grd_in(in_idx+1)-grd_in(in_idx)
                  end do        ! end loop over inner trapezoids
                  if (crd_dlt_ttl <= 0.0) then 
                     stop 'ERROR: crd_dlt_ttl <= 0.0 in xtr_RHS branch of rbn_vec()'
                  endif         ! endif crd_dlt_ttl <= 0.0
                  ovl_avg=trp_area_ttl/crd_dlt_ttl 
                  if (xtr_flg_RHS(xtr_prt_wgt_msk)) then
                     dat_out(out_idx)=ovl_frc*ovl_avg
                  else
                     dat_out(out_idx)=ovl_avg
                  endif         ! endif wgt
               else if (xtr_flg_RHS(xtr_prt_lnr_msk)) then
                  stop 'xtr_prt_lnr not implemented yet in rbn_vec()'
               else if (xtr_flg_RHS(xtr_prt_ngh_msk)) then
                  dat_out(out_idx)=dat_in(in_nbr)
               else if (xtr_flg_RHS(xtr_prt_nil_msk)) then
                  dat_out(out_idx)=0.0
               else
                  stop 'Unknown xtr_typ_RHS in prt branch of rbn_vec()'
               endif            ! endif xtr_typ_RHS
            endif               ! brk_rgt_idx(out_idx) == grd_in_nbr+1
            if (xtr_flg_RHS(xtr_vrb_msk)) write (6,'(a,a,i4,a,1p,e10.3,0p)') 
     $           xtr_sng_RHS(1:f_strlen(xtr_sng_RHS)),' yields dat_out(',out_idx,') = ',dat_out(out_idx)
         else
            write (6,'(a)') 'ERROR: Incomplete truth table rbn_vec():'
            write (6,'(2(a,i5,a,i5))') 
     $           'brk_lft_idx(',out_idx,
     $           ') = ',brk_lft_idx(out_idx),
     $           ', brk_rgt_idx(',out_idx,
     $           ') = ',brk_rgt_idx(out_idx)
            stop 'Unforeseen bracketing values in rbn_vec()' 
         endif                  ! endelse RHS extrapolate
      enddo                     ! end loop over output coordinates

      if (dbg_lvl >= dbg_io) then
c     Print input crd
         write (6,'(a)') 'Input crd_in to rbn_vec():'
         write (6,'(a,i4,a,a)') 'xtr_typ_LHS = ',xtr_typ_LHS,' = ',xtr_sng_LHS(1:f_strlen(xtr_sng_LHS))
         write (6,'(a,i4,a,a)') 'xtr_typ_RHS = ',xtr_typ_RHS,' = ',xtr_sng_RHS(1:f_strlen(xtr_sng_RHS))
         write (6,'(5(a,1x))') 'idx','crd_in','crd_in_min','crd_in_max','dat_in'
         do idx=1,in_nbr
            write (6,'(i4,1x,1p,4(e10.3,1x),0p)') 
     $           idx,crd_in(idx),grd_in(idx),grd_in(idx+1),dat_in(idx)
         end do                 ! end loop over crd_out
      endif                     ! endif dbg

      if (dbg_lvl >= dbg_io) then
c     Print input grd
         write (6,'(a)') 'Input grd_in rbn_vec():'
         write (6,'(3(a,1x))') 'idx','grd_in','dat_ntp_in'
         do idx=1,grd_in_nbr
            write (6,'(i4,1x,1p,2(e10.3,1x),0p)') 
     $           idx,grd_in(idx),dat_ntp_in(idx)
         end do                 ! end loop over crd_out
      endif                     ! endif dbg

      if (dbg_lvl >= dbg_io) then
c     Print bracketing indices
         write (6,'(/,a)') 'Grid output from rbn_vec():'
         write (6,'(5(a,1x))') 'idx','grd_out','brk_lft','brk_rgt','dat_ntp_out'
         do idx=1,grd_out_nbr
            write (6,'(i4,1x,1p,1(e10.3,1x,i4,1x,i4,1x),e10.3,0p)') 
     $           idx,
     $           grd_out(idx),brk_lft_idx(idx),brk_rgt_idx(idx),
     $           dat_ntp_out(idx)
         end do                 ! end loop over crd_out
      endif                     ! endif dbg

      if (dbg_lvl >= dbg_io) then
c     Print results of rebinning
         write (6,'(/,a)') 'Coordinate output from rbn_vec():'
         write (6,'(5(a,1x))') 'idx','crd_out','crd_out_min','crd_out_max','dat_out'
         do idx=1,out_nbr
            write (6,'(i4,1x,1p,4(e10.3,1x),0p)') 
     $           idx,crd_out(idx),grd_out(idx),grd_out(idx+1),dat_out(idx)
         end do                 ! end loop over crd_out
      endif                     ! endif dbg

      if (dbg_lvl >= dbg_io) then
c     Print trapezoid info
         write (6,'(/,a)') 'Trapezoid output from rbn_vec():'
         write (6,'(5(a,1x))') 'idx','grd_out','grd_outp1','trp_nbr','trp_idx'
         do out_idx=1,out_nbr
            trp_nbr=brk_lft_idx(out_idx+1)-brk_rgt_idx(out_idx)
            write (6,'(i4,1x,1p,2(e10.3,1x),i2,a2,$)') 
     $           out_idx,grd_out(out_idx),grd_out(out_idx+1),trp_nbr,': '
            if (brk_lft_idx(out_idx) /= brk_rgt_idx(out_idx)) then
               trp_srt_idx=brk_lft_idx(out_idx)
            else
               trp_srt_idx=brk_lft_idx(out_idx)-1
            endif               ! endif
            do trp_idx=1,trp_nbr-1
               in_idx=trp_srt_idx+trp_idx
               write (6,'(i2,a1,$)') in_idx,','
            end do              ! end loop over trp
            in_idx=trp_srt_idx+trp_nbr
            if (trp_nbr > 0) then 
               write (6,'(i2)') trp_srt_idx+trp_nbr 
            else 
               write (6,'(a4)') 'none'
            endif               ! endif
         end do                 ! end loop over crd_out
      endif                     ! endif dbg
      
c     Convert input and output arrays to original direction of monotonicity
      if (in_dcr) then
         if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Un-reversing input grid in rbn_vec()'
         call rvr_vec(grd_in,grd_in_nbr)
         call rvr_vec(dat_in,in_nbr)
      endif                     ! endif in_dcr
      if (out_dcr) then
         if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Un-reversing output grid in rbn_vec()'
         call rvr_vec(grd_out,grd_out_nbr)
         call rvr_vec(dat_out,out_nbr)
      endif                     ! endif in_dcr

      if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting rbn_vec()'
      return
      end                       ! end rbn_vec()
      
      subroutine rvr_vec(dat_out,out_nbr)
c     Purpose: Reverse an array in place

      use precision

      implicit none
c     Parameters
c     Input
      integer out_nbr           ! [nbr] max(out_nbr)=205001 for NO2.F max(bnd_nbr_JPL,bnd_nbr_NCAR,bnd_nbr_NOAA)
c     Input/Output
      real(r8) dat_out(out_nbr)     ! Input grid
c     Local
      integer idx
      real(r8) dat_swp(out_nbr)     ! Swap array
c     Main Code
c     Copy input array into local swap array
      do idx=1,out_nbr
         dat_swp(idx)=dat_out(idx)
      end do                    ! end loop over out_idx
c     Reverse local swap array into output array
      do idx=1,out_nbr
         dat_out(idx)=dat_swp(out_nbr-idx+1)
      end do                    ! end loop over out_idx
      return
      end                       ! end rvr_vec()

      subroutine trp_area_vec(in_nbr,grd_in,dat_in,trp_area)
c     Purpose: Compute area of each trapezoid in input array of trapezoids
c     N trapezoids are centered within N+1 gridpoints
c     Trapezoids are indexed with the Fortran (1-based) convention where
c     Trapezoid 1 is bounded by the abscissae grd_in(1) and grd_in(2)
c     Trapezoid 2 is bounded by the abscissae grd_in(2) and grd_in(3)
c     Trapezoid k is bounded by the abscissae grd_in(k) and grd_in(k+1)
c     Trapezoid N-1 is bounded by the abscissae grd_in(N-1) and grd_in(N)
c     Trapezoid N is bounded by the abscissae grd_in(N) and grd_in(N+1)

      use precision

      implicit none
c     Input
      integer in_nbr
      real(r8) grd_in(in_nbr+1)     ! Input grid
      real(r8) dat_in(in_nbr+1)     ! Input data
c     Output
      real(r8) trp_area(in_nbr)     ! Area of each trapezoid defined by input data
c     Local
      integer in_idx
c     Main Code
      do in_idx=1,in_nbr
         trp_area(in_idx)=0.5*
     $        (grd_in(in_idx+1)-grd_in(in_idx))*
     $        (dat_in(in_idx+1)+dat_in(in_idx))
      end do                    ! end loop over trapezoids
      return
      end                       ! end trp_area_vec()

      function ntp_vec_one(crd_nbr,crd_in,dat_in,crd_out)
c     Purpose: Set up and perform a one point interpolation call to ntp_vec()
c     This function is fairly useful, e.g., at inverting functions
c     The function takes care of the messy overhead needed to call ntp_vec()
c     Prototype:

      use precision

      implicit none

      real(r8) ntp_vec_one       ! Scalar interpolation (libcsz_f77)

c     Commons
#include <dbg.com> /* Debugging constants */
#include <xtr.com> /* Extrapolation constants */
c     Input
      integer crd_nbr
      integer out_nbr
      real(r8) crd_out(1)
      real(r8) crd_in(crd_nbr)
      real(r8) dat_in(crd_nbr)
c     Output
      real(r8) dat_out(1)
c     Local
      integer xtr_typ_LHS
      integer xtr_typ_RHS
c     Main Code
c     Initialize
      out_nbr=1
      xtr_typ_LHS=xtr_fll_lnr
      xtr_typ_RHS=xtr_fll_lnr
      call ntp_vec(crd_nbr,crd_in,dat_in,out_nbr,crd_out,dat_out,xtr_typ_LHS,xtr_typ_RHS)
      ntp_vec_one=dat_out(1)
      return
      end                       ! end ntp_vec_one()

      integer function vec_val2idx(crd,crd_nbr,val)
c     Purpose: Locate index of array member closest to specified value
c     Prototype:
c     integer vec_val2idx       ! Locate index of array member closest to specified value (libcsz_f77)

      use precision

      implicit none
c     Commons
#include <dbg.com> /* Debugging constants */
c     Input
      integer crd_nbr           ! Size of array
      real(r8) crd(crd_nbr)         ! Array
      real(r8) val                  ! Value to search for
c     Local
      integer idx               ! [idx] Counting index
      integer idx_val           ! Index ov value
      real(r8) dst_new              ! Distance from current point to value
      real(r8) dst_old              ! Closest distance yet to value
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
c     Main code
c     Vet input:
      if (crd_nbr < 1) write (6,'(2a,i4)') prg_nm(1:f_strlen(prg_nm)),': ERROR vec_val2idx() reports crd_nbr = ',crd_nbr
      idx_val=1
      dst_old=abs(crd(1)-val)
      do idx=1,crd_nbr
         dst_new=abs(crd(idx)-val)
         if (dst_new < dst_old) then
            idx_val=idx
            dst_old=dst_new
         endif                  ! endif
      enddo                     ! end loop over crd
      vec_val2idx=idx_val
      return
      end                       ! end vec_val2idx()

      subroutine vec_set(vec,vec_nbr,val)
c     Purpose: Set a vector vec of vec_nbr real(r8)'s to the value val
c     NB: Routine has same interface as CCM:control/resetr()

      use precision

      implicit none
c     Input
      integer vec_nbr           ! Size of input vector
      real(r8) val                  ! Value to fill vector with
c     Input/Output
      real(r8) vec(vec_nbr)         ! Input vector
c     Local
      integer idx               ! [idx] Counting index
      do idx=1,vec_nbr
         vec(idx)=val
      end do                    ! end loop over vec
      return
      end                       ! end vec_set()

#if (defined BXM) && (defined SGI)
c     Routines which are in default Cray and/or CCM libraries 
      subroutine whenflt(chk_nbr,crd,idx_srd,val,vld_idx,vld_nbr)
c     Purpose: Return array of indices vld_idx into crd which satisfy crd(vld_idx) < val
c     Only every idx_srd'th value in the input array will be checked
c     Use idx_srd=1 to check every element
c     chk_nbr is the number of elements which are compared
c     chk_nbr does not necessarily equal the size of crd or vld_idx
c     When idx_srd != 1, only a hyperslab of crd is checked
c     Function is rewrite of CCM:srchutil/whenflt()
c     This routine is very flexible, but somewhat confusing

      use precision

      implicit none
c     Input
      integer idx_srd           ! Stride to employ
      integer chk_nbr           ! [nbr] Number of elements to check
      real(r8) crd(*)               ! Values of candidates
      real(r8) val                  ! Target value to use in comparison
c     Output
      integer vld_nbr           ! [nbr] Number of valid points
      integer vld_idx(*)        ! Array of valid indices
c     Local
      integer idx               ! [idx] Counting index
      integer idx_crr           ! Current index into input array
c     Main Code
c     Initialize counter and initial index
      vld_nbr=0
      idx_crr=1
      if (idx_srd < 0) idx_crr=(-idx_srd)*(chk_nbr-1)+1
      do idx=1,chk_nbr
         if (crd(idx_crr) < val) then
            vld_nbr=vld_nbr+1
            vld_idx(vld_nbr)=idx
         end if                 ! endif criteria satisfied
         idx_crr=idx_crr+idx_srd
      end do                    ! end loop over idx
      return
      end                       ! end whenflt()

      subroutine whenfgt(chk_nbr,crd,idx_srd,val,vld_idx,vld_nbr)
c     Purpose: Return array of indices vld_idx into crd which satisfy crd(vld_idx) > val
c     Only every idx_srd'th value in the input array will be checked
c     Use idx_srd=1 to check every element
c     chk_nbr is the number of elements which are compared
c     chk_nbr does not necessarily equal the size of crd or vld_idx
c     When idx_srd != 1, only a hyperslab of crd is checked
c     Function is rewrite of CCM:srchutil/whenfgt()
c     This routine is very flexible, but somewhat confusing

      use precision

      implicit none
c     Input
      integer idx_srd           ! Stride to employ
      integer chk_nbr           ! [nbr] Number of elements to check
      real(r8) crd(*)               ! Values of candidates
      real(r8) val                  ! Target value to use in comparison
c     Output
      integer vld_nbr           ! [nbr] Number of valid points
      integer vld_idx(*)        ! Array of valid indices
c     Local
      integer idx               ! [idx] Counting index
      integer idx_crr           ! Current index into input array
c     Main Code
c     Initialize counter and initial index
      vld_nbr=0
      idx_crr=1
      if (idx_srd < 0) idx_crr=(-idx_srd)*(chk_nbr-1)+1
      do idx=1,chk_nbr
         if (crd(idx_crr) > val) then
            vld_nbr=vld_nbr+1
            vld_idx(vld_nbr)=idx
         end if                 ! endif criteria satisfied
         idx_crr=idx_crr+idx_srd
      end do                    ! end loop over idx
      return
      end                       ! end whenfgt()
#endif /* CRAY or CCM */
      
