c $Id$
c Data for divergence calculation.

      common /div0/
     $  wall    ,w1s     ,w2s     ,w1n     ,w2n     ,ddlam   ,ddmu
     $, cw

      real(r8)
     $  wall(plat+1) ! position of walls for cells containing Gaussian latitudes
     $, w1s(plat)    ! weight for interpolation to south cell wall
     $, w2s(plat)    ! weight for interpolation to south cell wall
     $, w1n(plat)    ! weight for interpolation to north cell wall
     $, w2n(plat)    ! weight for interpolation to north cell wall
     $, ddlam(plat)  ! 1./( 2.*dlam*rearth*cos( phi(j) ) )
     $, ddmu(plat)   ! 1./( rearth*dmu(j) )
     $, cw(plat+1)   ! cos( wall(j) )
