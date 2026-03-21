c $Id$
c Data for potential temperature calculations

      common /pottem0/ cappa, p0

      real(r8)
     $  cappa      ! R / Cp for dry air
     $, p0         ! reference pressure (1.e5 Pa)
