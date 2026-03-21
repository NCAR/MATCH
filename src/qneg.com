c $Id$
c Minimum mass mixing ratio for constituents.

      common /qneg0/ qmin(pcnst)

      real(r8)
     $  qmin     ! Global minimum constituent concentration
