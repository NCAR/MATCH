c $Id$
      common / cld0 / phi, coslat, cappa, gravit, rair, latvap, cpair
     $,               rhow, tfh2o, k700

      real(r8)
     $  phi(plat)    ! Gaussian latitudes (radians)
     $, coslat(plat) ! cosine of latitude
     $, cappa        ! R/Cp
     $, gravit       ! Gravitational acceleration
     $, rair         ! Gas constant for dry air
     $, latvap       ! latent heat of vaporization for liq. water
     $, cpair        ! Specific heat of dry air
     $, rhow         ! water density (g/cm^3)
     $, tfh2o        ! freezing point of water at STP (K)

      integer
     $  k700         ! Level closest to 700mb over ocean.
