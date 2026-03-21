c     $Header$ -*-fortran-*-

c     Purpose: Common block dstodx stores information needed to compute optical depth
c     dstodx is initialized in dst_odx_cmn_ini() which is called by dst_msc_cmn_ini()

c     Usage:
c#include <dstodx.h> /* Optical depth information */

c     Optical depth info
      real(r8) ext_cff_mss_dst_dgn(dst_nbr) ! [m2 kg-1] Mass extinction coefficient of dust for diagnostic optical depth
      common / dstodx /         ! dstodx is initialized in dst_odx_cmn_ini()
     $     ext_cff_mss_dst_dgn  ! [m2 kg-1] Mass extinction coefficient of dust for diagnostic optical depth
