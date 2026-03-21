c $Id$
c Data for dry deposition calculation.

      real(r8)
     $  rair                ! Gas constant for dry air (J/K/kg)
     $, gravit              ! Gravitational acceleration
     $, phi(plat)           ! grid latitudes (radians)

      common /drydep0/ rair, gravit, phi
