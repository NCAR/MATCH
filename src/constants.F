      module constants

      use precision

      implicit none
      save
      private
      public ::
     $  init_const      ! initialize constants module

      real(r8), parameter, public ::
     $  cpwv   = 1.81e3           ! Specific heat of water vapor
     $, cpair  = 1004.64          ! Specific heat of dry air
     $, rair   = 287.04           ! Gas constant for air
     $, rh2o   = 461.             ! Gas constant for water vapor
     $, zvir   = rh2o/rair - 1.   ! rh2o/rair - 1.
     $, gravit = 9.80616          ! Gravitational acceleration (m/s^2)
     $, rgravit= 1./gravit        ! Gravitational acceleration (m/s^2)
     $, rearth = 6.37122e+6       ! radius of earth (m)
     $, ra     = 1./rearth        ! 1./rearth
     $, epsilo = 0.622            ! ratio of mol. wt. water to dry air
     $, latvap = 2.5104e06        ! latent heat of vaporization for liq. water (J/kg)
     $, latice = 3.336e5          ! latent heat of sublimation for ice
     $, cappa  = rair/cpair       ! R/Cp (air)
     $, rhoh2o = 1.e3             ! density of liquid water (STP) (kg/m^3)
     $, tfh2o  = 273.16           ! freezing point of water at STP (K)

      real(r8), public ::
     $  pi                        ! 3.14159...

c##############################################################################
      contains
c##############################################################################

      subroutine init_const()

      implicit none

      pi = 4.d0*atan( 1.d0 )

      end subroutine init_const

      end module constants
