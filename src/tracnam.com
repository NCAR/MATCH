c $Id$
c Names for advected and non-advected trace species, and
c surface fluxes and process tendencies for the advected species.
c This block is initialized in inirun.

      character
     $  tracnam(pcnst)*8   ! specie names
     $, sflxnam(pcnst)*8   ! surface flux names
     $, tanam(pcnst)*8     ! total advection tendency names
     $, vdnam(pcnst)*8     ! vertical diffusion tendency names
     $, dcnam(pcnst)*8     ! convection tendency names
     $, xflxnam(pcnst+2)*8 ! fluxes on x cell walls, tracers + air + Q
     $, yflxnam(pcnst+2)*8 ! fluxes on y cell walls, tracers + air + Q
     $, zflxnam(pcnst+2)*8 ! fluxes on z cell walls, tracers + air + Q
#if defined STRAT_CHEM || defined TROP_CHEM
     $, natsnam(pnats)*8   ! names of non-advected trace species
#else
     $, natsnam(1)*8
#endif

      common / tracnm_C /
     $  tracnam, sflxnam, tanam, vdnam, dcnam
     $, xflxnam, yflxnam, zflxnam
     $, natsnam
