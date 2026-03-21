!c     common block for storing sea salt opt depth
      real(r8) tsslt2d (plond,plat) ! sea salt optical depth
      real(r8) sslt_profiles( plond, plev, plat )   !CTB
      common /ssltcom/ tsslt2d, sslt_profiles  ! new way !CTB
