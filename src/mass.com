c $Id$
c Data for mass integrals.

      common /mass/ rgrav, cwava, gw, etamid

      real(r8)
     $  rgrav        ! 1./gravit
     $, cwava        ! 1./(gravit*plon)
     $, gw(plat)     ! Gauss weights (these sum to 2)
     $, etamid(plev) ! hybrid coordinate a layer midpoints
