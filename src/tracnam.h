! $Id$
! Names for advected and non-advected trace species, and
! surface fluxes and process tendencies for the advected species.
! This block is initialized in inirun.

       character*8 tracnam(pcnst)   ! specie names
       character*8 sflxnam(pcnst)   ! surface flux names
!      character
!     $  tracnam(pcnst)*8   ! specie names
       character*8  tanam(pcnst)     ! total advection tendency names
       character*8 vdnam(pcnst)    ! vertical diffusion tendency names
       character*8 dcnam(pcnst)     ! convection tendency names
       character*8 xflxnam(pcnst+2) ! fluxes on x cell walls, tracers + air + Q
       character*8 yflxnam(pcnst+2) ! fluxes on y cell walls, tracers + air + Q
       character*8 zflxnam(pcnst+2) ! fluxes on z cell walls, tracers + air + Q
#if defined STRAT_CHEM || defined TROP_CHEM
       character*8 natsnam(pnats)   ! names of non-advected trace species
#else
       character*8 natsnam(1)
#endif

      common / tracnm_C / tracnam, sflxnam, tanam, vdnam, dcnam, xflxnam, yflxnam, zflxnam, natsnam
