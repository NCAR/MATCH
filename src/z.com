c $Id$
c Data for geopotential height calculations

      common /z0/ rair, gravit  

      real(r8)
     $  rair       ! Gas constant for dry air
     $, gravit     ! acceleration due to gravity (m/s^2)
