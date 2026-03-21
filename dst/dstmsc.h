c     $Header$ -*-fortran-*-

c     Purpose: Miscellaneous common blocks
c     These common blocks are initialized in dst_msc_cmn_ini() which is called by CCM:inti() 

c     Usage: 
c     #include <dstmsc.h> /* Miscellaneous variables for dust parameterization */
c     params.h needed for DST_NBR

c     Debug info
      common / dst_dbg_cmn /    ! dst_dbg_cmn is initialized in dst_msc_cmn_ini()
     $     lat_dbg,             ! [idx] Latitude  index for diagnostics
     $     lev_dbg,             ! [idx] Level     index for diagnostics
     $     lon_dbg,             ! [idx] Longitude index for diagnostics
     $     time_dbg             ! [idx] Time      index for diagnostics
      integer lat_dbg           ! [idx] Latitude  index for diagnostics
      integer lev_dbg           ! [idx] Level     index for diagnostics
      integer lon_dbg           ! [idx] Longitude index for diagnostics
      integer time_dbg          ! [idx] Time      index for diagnostics

c     Physics info
      real(r8) flx_mss_sgn          ! [kg m-2 s-1] Minimum significant flux of dust
      real(r8) flx_mss_mxm          ! [kg m-2 s-1] Maximum plausible flux of dust
      real(r8) q_dst_mxm            ! [kg kg-1] Maximum plausible mass mixing ratio of dust
      real(r8) q_dst_sgn            ! [kg kg-1] Minimum significant mass mixing ratio of dust
      real(r8) z_scv(dst_nbr)       ! (kg dust/kg rain)/(kg dust/kg air) Scavenging ratio

      common / dst_phys /       ! dst_phys is initialized in dst_msc_cmn_ini()
     $     flx_mss_sgn,         ! [kg m-2 s-1] Minimum significant flux of dust
     $     flx_mss_mxm,         ! [kg m-2 s-1] Maximum plausible flux of dust
     $     q_dst_mxm,           ! [kg kg-1] Maximum plausible mass mixing ratio of dust
     $     q_dst_sgn,           ! [kg kg-1] Minimum significant mass mixing ratio of dust
     $     z_scv                ! (kg dust/kg rain)/(kg dust/kg air) Scavenging ratio
      

