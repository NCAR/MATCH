c $Id$
c Hybrid level definitions: p = a*p0 + b*ps
c        interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
c        midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps

      real(r8)
     $  hyai    ! ps0 component of hybrid coordinate - interfaces
     $, hyam    ! ps0 component of hybrid coordinate - midpoints
     $, hybi    ! ps component of hybrid coordinate - interfaces
     $, hybm    ! ps component of hybrid coordinate - midpoints
     $, ps0     ! base state sfc pressure for level definitions
     $, dhya    ! dhya(k) = hyai(k+1) - hyai(k)
     $, dhyb    ! dhyb(k) = hybi(k+1) - hybi(k)
     $, ptop    ! pressure of top interface

      common /plevs/
     $     hyai(plevp) ,hyam(plev) ,hybi(plevp) ,hybm(plev) ,ps0
     $,    dhya(plev), dhyb(plev), ptop
