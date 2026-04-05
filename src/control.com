c $Id$
c Control parameters.

      logical
     $  incore        ! t => incore history buffer and nas buffer
     $, divdiff       ! t => diagnose vertical diffusion
     $, arvdiff       ! t => archived vertical diffusion
     $, diconvccm     ! t => diagnose CCM convection
     $, arconvccm     ! t => archived CCM convection
     $, dicldphys     ! t => diagnose cloud physics
     $, arcldphys     ! t => archived cloud physics
     $, stratchem     ! t => stratospheric chemistry
     $, tropchem      ! t => tropospheric chemistry
     $, chembgt       ! t => perform budget computations for trop chem
     $, vdiffuse      ! t => vertical diffusion scheme enabled
     $, convect       ! t => convection scheme enabled
     $, chemistry     ! t => chemistry scheme enabled
     $, physics       ! t => vdiffuse or convect or chemistry or dicldphys
     $, fixmas        ! t => apply SLT mass fixer
     $, calcedot      ! t => calculate etadot
     $, advqth        ! t => advect Q and THETA (potential temperature)
     $, rdomega       ! t => read omega from dynamics input files
     $, satvap        ! t => enable saturated vapor pressure calculation
     $, tvbnd         ! t => emissions from time varying boundary dataset
     $, dowetdep      ! t => apply wet deposition parameterization
     $, addqrad       ! t => add radiative heating tendency to the temperature
                      !      field supplied to the convection routine
     $, limqin        ! t => limit input Q s.t. it doesn't exceed saturation
     $, globmass      ! t => write global tracer mass at each level into netCDF
                      !      file at each timestep
     $, dodrydep      ! t => apply dry deposition parameterization
     $, readQIC       ! t => read initial conditions for the advected Q field
                      !      rather than initializing from meteorology.
     $, useAltConv    ! t => use alternative convective parameterization
     $, doAODasm      ! t => invoke AOD assimilation code
     $, wind_ps_adj   ! t => adjust the winds to enforce consistency with input
                      !      surface pressure field.
      logical
     $  begstep     ! t => current timestep is beginning timestep for run
     $, laststep    ! t => current timestep is final timestep for run
     $, endofrun    ! t => last timestep or shutdown requested
     $, wrhstts     ! t => write a history file time sample this timestep
     $, fullhst     ! t => history file full after writing time sample
     $, closehst    ! t => close history file (flush output buffers)
     $, wrestart    ! t => write a restart file this timestep

      common /control_L/ incore
     $, divdiff, arvdiff, diconvccm, arconvccm, dicldphys, arcldphys
     $, stratchem, tropchem, chembgt, vdiffuse, convect
     $, chemistry, fixmas, calcedot, physics, advqth
     $, rdomega, satvap, tvbnd, dowetdep, addqrad, limqin, globmass
     $, dodrydep, readQIC, useAltConv, doAODasm, wind_ps_adj
     $, begstep, laststep, endofrun
     $, wrhstts, fullhst, closehst, wrestart

c     Control time sequence.
c     *** Note on terminology ***  timestep N advances the integration from
c     time index N-1 to time index N.

      integer
     $  rstflg      ! 0 for initial run, 1 for restart
     $, itim0       ! beginning time index for run
     $, curstep     ! current timestep for run
     $, nestep      ! ending timestep for run
     $, delt        ! timestep interval in seconds
     $, icday       ! day of initial conditions (time index 0)
     $, icdaysec    ! seconds relative to icday
     $, icdate      ! date of initial conditions (time index 0) in yyyymmdd format
     $, icdatesec   ! seconds relative to icdate
     $, mpdate      ! date at mid-point of current timestep
     $, mpdatesec   ! seconds relative to mpdate
     $, curday      ! day at end of current timestep
     $, curdaysec   ! seconds relative to curday
     $, curdate     ! date at end of current timestep
     $, curdatesec  ! seconds relative to curdate
     $, basedate    ! date at time = 0.0
     $, basedatesec ! seconds relative to basedate

      common /control_I/ rstflg, itim0, nestep, curstep, delt
     $, icday, icdaysec, icdate, icdatesec, mpdate, mpdatesec
     $, curday, curdaysec, curdate, curdatesec, basedate, basedatesec

c     Updating advected Q and potential temperature.
      real(r8)
     $  qrelax      ! Relaxation factor for "nudging" predicted Q field
c                   ! towards input Q field.  qrelax=1 => Q is reset to the
c                   ! input values each timestep.  qrelax=0 => Q evolves
c                   ! without being affected by the input Q field except for
c                   ! being initialized by it.  The default value is ~ .1.
     $, vwc_scale   ! [frc] Soil moisture rescaling slope for dust emission
c                   ! vwc_sfc = vwc_scale*vwc_sfc + vwc_offset
c                   ! Default 1.0 (identity). Set per-month when tuning
c                   ! against MODIS/VIIRS column AOD.
     $, vwc_offset  ! [m3 m-3] Soil moisture rescaling offset for dust emission
c                   ! Default 0.0 (identity).
     $, sslt_scale  ! [frc] Sea salt emission flux multiplier
c                   ! Default 1.0 (no change). Set per-month when tuning
c                   ! against MODIS/VIIRS column AOD over clean-ocean regions.

      common /control2/ qrelax, vwc_scale, vwc_offset, sslt_scale




