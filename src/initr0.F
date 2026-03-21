c $Id$
      subroutine initr0( q, as, etamid )

c Initialize idealized tracer concentrations.

      use precision

      implicit none
c---------------------------------------------------------------------
#include <pmgrid.h>
c---------------------------------------------------------------------

c Input argument:
      real(r8) q(plond,plev)           ! specific humidity
     $, etamid(plev)  ! hyam + hybm

c Output argument:
      real(r8)
     $  as(plond,plev,pcnst,platd) ! advected species

c Local variables:

      integer
     $  i, j, k, m   ! indices

      real(r8) weight, apres

c---------------------------------------------------------------------

c Init to zero
      as(1:plon,:,:,1:plat) = 0.0

#ifndef AER
#ifdef TESTTRAC
c Ozone-like tracer.
      do k = 1, plev
         apres = 1000.*etamid(k) ! approximate pres in mb
         if (apres.lt.50.) then
            as(1:plon,k,1,1:plat) = 1.
         else if (apres.gt.100) then
            as(1:plon,k,1,1:plat) = 0.
         else
            weight = min(1._r8,max(0._r8,(apres-50)/(100.-50.)))
            as(1:plon,k,1,1:plat) = 1.-weight
         endif

c Upper troposphere (about 250 mb)
         if (apres.gt.200.and.apres.lt.300.) then
            as(1:plon,k,2,1:plat) = 1.
         endif

c Middle troposphere (about 500 mb)
         if (apres.gt.450..and.apres.lt.550.) then
            as(1:plon,k,3,1:plat) = 1.
         endif

      end do

c Boundary Layer
      as(1:plon,plev-1:plev,4,1:plat) = 1.

c Unit tracer
      as(1:plon,:,5,1:plat) = 1.

#endif /* TESTTRAC */
#endif

      return
      end
