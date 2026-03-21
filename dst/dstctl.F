c     $Header$
      
c     Purpose: dstctl.F contains control routines for the dust model

c     Requires dst.h for DST_MSS_BDG and many other tokens
#include <dst.h> /* Dust preprocessor tokens */

      subroutine dst_bnr()
c     Purpose: Print banner with dust parameterization information
c     dst_bnr() is called by CCM:dst/dst_szdst_ini(), MATCH:main()
      implicit none
#include <dstgrd.h> /* Dust grid sizes */
#include <dbg.com> /* Debugging constants */
c     Local
      character CVS_Date*28
      character CVS_Revision*16
      character CVS_Id*80       ! CVS ID
      character CVS_Name*30
      character cmd_ln*200      ! Command line
      character prg_ID*200      ! Program ID
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
c     Initialize defaults
      CVS_Date='$Date$'
      CVS_Id='$Id$'
      CVS_Revision='$Revision$'
      CVS_Name='$Name$'
c      dbg_lvl=dbg_off           ! dbg_lvl allocated in dbg.com, Option D
      call cmd_ln_sng(cmd_ln)
      call prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
c     Describe the dust parameterization being used
      write (6,'(a)') '------------------------------------------------------------'
      write (6,'(a)') 'Mineral dust parameterization is active:'
      write (6,'(a)') prg_ID(1:f_strlen(prg_ID))
      write (6,'(a)') '------------------------------------------------------------'
c      write (6,'(a)') 'Mineral Dust Version ' // CVS_Name(8:30) // 'dated ' // CVS_Date(8:26) // ' GMT'
      write (6,'(i2,a,i2,a,i2)') 
     $     dst_nbr,' mineral dust bins are advected ' //
     $     'species ',dst_idx_srt,' through ',dst_idx_end
      write (6,'(a)') 'Mobilization parameterization:'
      write (6,'(2x,a)') 
     $     'Dynamic source scheme is active'
      write (6,'(a)') 'Radiative mode:'
      write (6,'(2x,a)') 
#ifdef DST_RAD
     $     'Radiative effects of dust are computed'
#else /* not DST_RAD */
     $     'Radiative effects of dust are not computed'
#endif /* not DST_RAD */
      write (6,'(2x,a)') 
#ifdef DST_FDB
     $     'Dust radiative effects are interactive, not diagnostic'
#else /* not DST_FDB */
     $     'Dust radiative effects are diagnostic, not interactive'
#endif /* not DST_FDB */
      write (6,'(2x,a)') 'Diagnostic optical depths will be archived'
      write (6,'(a)') 'Land surface mode:'
      write (6,'(2x,a)') 
#ifdef DST_LSM
     $     'Interactive LSM surface soil moisture is active'
#else /* not DST_LSM */
     $     'Interactive LSM surface soil moisture is turned off'
#endif /* not DST_LSM */
      write (6,'(a)') 'Chemistry mode:'
      write (6,'(2x,a)') 
#ifdef DST_CHM
     $     'Heterogeneous chemistry is active (but incomplete)'
#else /* not DST_CHM */
     $     'Heterogeneous chemistry is turned off'
#endif /* not DST_CHM */
      write (6,'(a)') 'Debugging modes:'
      write (6,'(2x,a)') 
#ifdef DST_DBG
     $     'Unvectorized run-time debugging blocks are active'
#else /* not DST_DBG */
     $     'Unvectorized run-time debugging blocks are turned off'
#endif /* not DST_DBG */
      write (6,'(2x,a)') 
#ifdef DST_MSS_BDG
     $     'Diagnostic global mass budget is active'
#else /* not DST_MSS_BDG */
     $     'Diagnostic global mass budget is turned off'
#endif /* not DST_MSS_BDG */
      return
      end                       ! end dst_bnr()

