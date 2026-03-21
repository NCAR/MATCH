c $Id$

      real(r8)
     $  lat        ! grid latitudes in radians
     $, lon        ! grid longitudes in radians
     $, latDeg     ! grid latitudes in degrees
     $, lonDeg     ! grid longitudes in degrees

      common /grid_R/ lat(plat), lon(plon), latDeg(plat), lonDeg(plon)
