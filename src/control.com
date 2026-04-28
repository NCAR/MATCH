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

c     Sea salt regional tuning: compile-time limit on number of bands.
      integer max_sslt_bands
      parameter (max_sslt_bands = 10)

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
     $, sslt_scale  ! [frc] Sea salt emission flux global flux multiplier
c                   ! Default 1.0 (no change). Composes multiplicatively with
c                   ! any per-band scaling in sslt_bands.
     $, dst_rgn_scale(7) ! [frc] Per-region dust emission multiplier.
c                   ! Regions (lat/lon bounds in dstmbl.F):
c                   !   1 = Sahara/Sahel   (10-30N,  20W-30E)
c                   !   2 = Arabia         (15-30N,  30-55E)
c                   !   3 = Central Asia   (30-50N,  55-90E)
c                   !   4 = Gobi/Taklamakan(35-50N,  90-115E)
c                   !   5 = Australia      (35-15S, 115-145E)
c                   !   6 = SW N.America   (25-40N, 115W-100W)
c                   !   7 = Patagonia      (55-35S,  75W-60W)
c                   ! Default 1.0 (no per-region change).
     $, sslt_bands(3, max_sslt_bands) ! Sea salt regional tuning bands.
c                   ! Each active column (ib = 1..n_sslt_bands) is a triple
c                   !   (lat_min_deg, lat_max_deg, scale).
c                   ! Applied in src/getsslt.F: a column with latitude in
c                   ! [lat_min, lat_max] gets its sslt multiplied by scale
c                   ! (and by sslt_scale globally). Columns outside every
c                   ! band get scale 1.0. Bands are validated at startup in
c                   ! main.F — latitudes must lie in [-90, 90], lat_min must
c                   ! be strictly less than lat_max, and bands must not
c                   ! overlap (touching at the boundary is allowed).
c                   ! Default 0.0 for all entries (inactive until
c                   ! n_sslt_bands > 0 and entries populated).
     $, so2_rgn_scale(4)  ! [frc] Per-region SOx (SO2+SO4) emission multiplier.
c                   ! Applied in src_scyc/sulemis.F90 to all SOx-derived sflx
c                   ! components (98% SO2 + 2% SO4, both 0m and >100m).
c                   ! Volcanic SO2 is unaffected. Regions:
c                   !   1 = E. Asia       (20-50N, 100-145E)
c                   !   2 = S. Asia       ( 5-35N,  65-95E)
c                   !   3 = Europe        (35-65N,  10W-50E)
c                   !   4 = N. America    (25-60N, 125W-65W)
c                   ! Default 1.0 (no per-region change).
     $, dms_rgn_scale(6)  ! [frc] Per-region DMS emission multiplier.
c                   ! Applied in src_scyc/sulemis.F90 to the DMS sflx slot.
c                   ! Regions (ocean basins):
c                   !   1 = N. Pacific    (  0-60N, 120E-100W)
c                   !   2 = S. Pacific    (60S- 0 , 145E- 70W)
c                   !   3 = N. Atlantic   (  0-60N, 100W- 20E)
c                   !   4 = S. Atlantic   (60S- 0 ,  70W- 20E)
c                   !   5 = Indian        (30S-30N,  30E-120E)
c                   !   6 = Southern      (90S-60S, all)
c                   ! Default 1.0 (no per-region change).
     $, oc_rgn_scale(9)   ! [frc] Per-region organic carbon emission multiplier.
c                   ! Applied in src/caer.F90 to the OC components
c                   ! (BBOCSF + FFOCSF + NOCSF) before composing sflx.
c                   ! Regions (tested in order; first match wins):
c                   !   1 = Amazon        (15S-10N,  80W- 50W)
c                   !   2 = S. Africa     (35S- 5N,  10E- 45E)
c                   !   3 = SE Asia/Indo  (10S-20N,  95E-150E)
c                   !   4 = Australia     (40S-10S, 110E-155E)
c                   !   5 = E. Asia       (20N-50N, 100E-145E)
c                   !   6 = S. Asia       ( 5N-35N,  65E- 95E)
c                   !   7 = Europe        (35N-65N,  10W- 50E)
c                   !   8 = N. America    (25N-60N, 125W- 65W)
c                   !   9 = Boreal NH     (50N-75N, all)
c                   ! Default 1.0 (no per-region change).
     $, bc_rgn_scale(9)   ! [frc] Per-region black carbon emission multiplier.
c                   ! Applied in src/caer.F90 to the BC components
c                   ! (BBBCSF + FFBCSF) before composing sflx. Regions match
c                   ! oc_rgn_scale (see above). Default 1.0.
      integer
     $  n_sslt_bands ! Number of active entries in sslt_bands, 0..max_sslt_bands.
c                    ! Default 0 disables regional scaling.

      common /control2/ qrelax, vwc_scale, vwc_offset, sslt_scale
     $, dst_rgn_scale, sslt_bands
     $, so2_rgn_scale, dms_rgn_scale, oc_rgn_scale, bc_rgn_scale
      common /control3/ n_sslt_bands




