c $Id$
c Data for etadot calculation.

      integer
     $  nprlev        ! index of lowest pure pressure interface.
     $, nfltpts(plat) ! number of points in filter

      common /etadot0/
     $  nprlev, nfltpts

      logical
     $  filter        ! true => apply filter to calculated etadot values

      common /etadot1/
     $  filter

      real(r8)
     $  etamid(plev)  ! hybrid coords of layer midpoints
     $, hybi(plevp)
     $, hybd(plev)    ! diffs of hybrid coords across layers
     $, dlam          ! longitude grid spacing (radians)
     $, arad          ! radius of earth (m)
     $, phi(platd)    ! Latitudes (radians) on the extended grid.

      common /etadot2/
     $  etamid, hybi, hybd, dlam, arad, phi
