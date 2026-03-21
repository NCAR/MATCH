#include <params.h>

c Set kind of real variables to have at least 12 digits of precision.
c     integer, parameter :: REALKIND = selected_real_kind( p = 12 )
c      integer, parameter :: REALKIND = selected_real_kind( p = RPREC )

c     Basic grid point resolution parameters

      integer, parameter ::
     $  plon = PLON                       ! number of longitudes
     $, plat = PLAT                       ! number of latitudes
     $, plev = PLEV                       ! number of vertical levels
     $, pcnst = PCNST                     ! number of advected constituents
     $, pnats = PNATS                     ! number of non-advected trace species
     $, plevp = plev + 1
     $, plevd = 2*plev
     $, nxpt = 1                          ! no. of pts outside active domain of interpolant
     $, jintmx = 1                        ! number of extra latitudes in polar region
     $, plond = plon + 1 + 2*nxpt
     $, platd = plat + 2*nxpt + 2*jintmx
     $, i1 = nxpt + 1                     ! model starting longitude (3-d)
     $, j1 = jintmx + nxpt + 1            ! model starting latitude (3-d)
     $, j1m = j1 - 1                      ! model starting offset (3-d)
     $, padv = 2                          ! number of miscellaneous advected fields
     $, mxdynflds = 42                    ! maximum number of dynamics input fields
     $, mxoutflds = 7*pcnst + pnats + 100 ! maximum number of history output fields

      logical, parameter ::
     $  masterproc = .true.               ! for CCM compatibility.  true => shared memory code
