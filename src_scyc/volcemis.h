c $Id$
c SO2 emissions from non-eruptive volcanoes.

#include <params.h>

      integer plon, plat
      parameter( plon=PLON, plat=PLAT )

      real(r8)
     $  so2(plon,plat)

      common / volcemis_R /
     $  so2

      integer
     $  hgtidx(plon,plat)
     $, ncid      ! ID for netCDF file

      common / volcemis_I /
     $  hgtidx
     $, ncid

