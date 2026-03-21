c $Id$
c Data for dry mass adjustment.

      common /drymass/ mdryatm

      real(r8)
     $  mdryatm      ! atmospheric dry mass (kg/m^2)
