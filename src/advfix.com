c $Id$
c Data for advection mass fixer.

      logical
     $  limfix       ! true => limit the fixer

      common /advfix_L/ limfix

      real(r8)
     $  etamid(plev) ! hybrid coordinate a layer midpoints
     $, alpha(pcnst) ! Fixer coeffiecients

      common /advfix_R/ etamid, alpha
