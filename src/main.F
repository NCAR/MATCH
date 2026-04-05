c $Id$
c $Name$
#undef DEBUG

#ifdef DST
c++csz
#include <dst.h> /* Dust preprocessor tokens */
c--csz
#endif
      program main

c Model for Atmospheric Transport and CHemistry.

      use acsf, only: acsfint
      use advect, only: adv_dr
      use dyninp, only: chkdyn, interp1, interp2
      use histout, only: wrhd, disphst, nexthist
      use cldwat, only: inimc, inimland
#ifdef SCYC
c++scyc
      use scyc, only: doScyc
      use acbnd, only: acbndint
      use dmsbnd, only: dmsbndint
      use soxbnd, only: soxbndint, soxsfconst, soxsfyr
c--scyc
#endif
#ifdef CAER
c++caer
      use caer, only: docaer
      use caerbnd, only: caerbndint
c--caer
#endif
#ifdef ASSIM
c++assim
      use aod_assim_3d, only: Assimilate_AOD
c--assim
#endif
      use precision
#ifdef DST
      use dstszdst
      use dsttvbds
#endif

      use calendar

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <control.com>
c-----------------------------------------------------------------------
c++csz
#ifdef DST
c     Initialize all dust model block data
#include <dstgrd.h> /* Dust grid sizes */
      external dstlsmbd         ! Initialize LSM block data
      external dstszdstbd       ! Initialize aerosol block data
#endif /* not DST */
c--csz

      integer
     $  i
     $, nt          ! time index for advected species array (current)
     $, ntp1        ! time index for advected species array (forecast)
     $, ncdate      ! date at end of current timestep in yymmdd format
     $, ncsec       ! seconds relative to ncdate
     $, mpdat       ! date at midpoint of current timestep in yymmdd format
     $, mpsec       ! seconds relative to mpdat
     $, nstep       ! timestep counter
     $, nadv        ! number of miscellaneous advected fields

      real(r8)
     $  dtime       ! timestep interval in seconds
     $, calday      ! Julian calandar day for current timestep

      real(r8)
     $  uv(plond,plev,2,platd)      ! u,v winds
     $, w(plon,plevp,plat)          ! vertical wind (etadot) (1/s)
     $, wr(plon,plevp,plat)         ! raw vertical wind (etadot) (1/s)
     $, mca(plon,plevp,plat)        ! raw mass convergence above (etadot) (Pa/s)
     $, mcar(plon,plevp,plat)       ! mass convergence above (etadot) (Pa/s)
     $, omga(plond,plev,plat)       ! vertical velocity (omega) (Pa/s)
     $, ps(plond,plat,2)            ! surface pressure
     $, psdry(plond,plat,2)         ! surface pressure of dry atm
     $, oro(plond,plat)             ! orography
     $, phis(plond,plat)            ! surface geopotential (m**2/s**2)
     $, as(plond,plev,pcnst,platd,2)! advected trace species
     $, hw1(pcnst)                  ! constituent mass before advection.
     $, trmk(plev,pcnst)            ! global average mass (kg/m2)
     $, mixk(plev,pcnst)            ! global average m.r. (kg/kg)

      real(r8), dimension(:,:,:,:,:), allocatable ::
     $   x                          ! misc advected variables

      real(r8)
     $  t(plond,plev,plat)     ! temperature
     $, q(plond,plev,plat,2)   ! specific humidity

c     Move this stack storage to static to avoid stack size limits.
      save uv, w, omga, as, t, q

      real(r8)
     $  ts(plond,plat)       ! surface temperature
     $, taux(plond,plat)     ! stress in x direction
     $, tauy(plond,plat)     ! stress in y direction
     $, hflx(plond,plat)     ! surface sensible heat flux (w/m2)
     $, qflx(plond,plat)     ! surface water vapor flux (kg/m2/s)
     $, fsds(plond,plat)     ! downward solar flux at surface (W/m2)
     $, fsus(plond,plat)     ! upward solar flux at surface (W/m2)
     $, flds(plond,plat)     ! downward longwave flux at surface (W/m2)
     $, vwcsfc(plond,plat)   ! volumetric soil water content at surface (frac)
     $, snowh(plond,plat)    ! water equivalent snow depth (m)

c     Vertical Diffusion.
      real(r8), dimension(:,:,:), allocatable ::
     $  cgs       ! counter-gradient coefficient
     $, kvh       ! vertical diffusion coefficient

c     Moist convection.
      real(r8), dimension(:,:,:), allocatable ::
     $  zmu       ! mu2 from conv_ccm, kg/m2/s
     $, zmd       ! md2 from conv_ccm, kg/m2/s
     $, zeu       ! eu2 from conv_ccm, 1/s
     $, hketa     ! convective mass flux, Hack
     $, hkbeta    ! overshoot parameter, Hack
     $, cmfdqr    ! dq/dt due to convective rainout 
     $, conicw    ! convective in-cloud water mixing ratio
     $, qrad      ! radiative heating tendency (K/s)

c     Prognostic cloud water.
      real(r8), dimension(:,:,:), allocatable ::
     $  cldfr     ! cloud fraction
     $, concld    ! convective cloud fraction
     $, cme       ! rate of cond-evap within the cloud
     $, prain     ! rate of release of stratiform precip (1/s)
     $, evapr     ! rate of evaporation of precipitation (1/s)
     $, cwat      ! cloud water (kg/kg)

#ifdef SCYC
c++scyc
      integer mpyear
c--scyc
#endif

c Externals:
!     integer newdate
      logical
     $  adj_winds       ! adjust winds to enforce mass conservation
     $, use_slt         ! t => use SLT advection

c++assim
      logical AODasm
c--assim

      character vtag*64  ! version tag expanded from cvs keyword



c-----------------------------------------------------------------------


c     print out version tag from CVS
      vtag = '$Name$'
      write(*,*)'MATCH SARB version ',vtag

c     Initialize time level indices.
      nt   = 1
      ntp1 = 2

c     Initialize date and field variables.

      call inirun
#ifdef DEBUG
      write(6,*)'DEBUG main 0. after inirun '
#endif

      dtime = float( delt )

c     Output initial mass diagnostics
      if ( globmass ) then
         call gamk2( as(i1,1,1,j1,nt), psdry(:,:,nt), trmk, mixk )
         call wrmass2nc( trmk, mixk )
      end if
c++csz
#ifdef DST_MSS_BDG
c     Calling bdg_update() before main loop writes an nstep=0 record containing initial values to dst_mss_bdg.nc
c     This places dst_mss_bdg.nc one timestep out of sync with MATCH history tapes since MATCH never writes nstep=0 initial values
cc     Date at end of timestep
c      call getcurdate(ncdate,ncsec)
cc     Output global average mass budget for this timestep
c      call bdg_update(itim0,ncdate,ncsec)
#endif /* not DST_MSS_BDG */
c--csz

c     ***********************
c     * BEGIN TIME STEPPING *
c     ***********************


      do 1000 nstep = itim0+1, nestep
#ifdef DEBUG
         write(6,*)'DEBUG main 0.3 starting time step'
#endif

c        Set flags to control the current timestep.
         call ctlstep( nstep )

c        Set beginning of timestep values in misc advected species array.
         if ( advqth ) then
c           Initialize misc advected variables with Q and THETA.
c           *** N.B. *** It is assumed here that q and t have not been
c           modified in the physics scan because after the first time step
c           the values in q and t are the ones from interpolation to the
c           end of the previous timestep.
c           x(1,1,1,.)=q, x(1,1,2,.)= potential temperature from t
            call setqth( ps(:,:,nt), q(:,:,:,nt), t, x(i1,1,1,j1,1)
     $,                  x(i1,1,1,j1,2) )
         end if
c        Set date at midpoint of timestep.
         call getmpdate( mpdat, mpsec )
         calday = caldayr( mpdat, mpsec )

c        Interpolate surface flux data.
         call acsfint( calday )

#ifdef SCYC
c++scyc, interpolate boundary data
         if ( doScyc() ) then
c           soxsfconst true => use one year of soxsf data ( soxsfyr) 
            if ( .not. soxsfconst ) soxsfyr = mpdate/10000
c-drb            mpyear = mpdate/10000
            if ( soxsfyr .lt. 1000 ) soxsfyr = soxsfyr + 1900
            call soxbndint( soxsfyr, calday )
            call dmsbndint( calday )
            call acbndint( calday )
         end if
c--scyc
#endif
c++csz
#ifdef DST
c     Time interpolate seasonal cycle of surface dust properties
         call dst_tvbds_ntp(calday) 
#endif /* not DST */
c--csz
#ifdef CAER
c++caer
         if ( docaer() ) then
            call caerbndint( calday )
         end if
c--caer
#endif

c        Interpolate data to midpoint of timestep.
         call chkdyn( mpdat, mpsec, phis, oro )
         call interp2( mpdat, mpsec, uv(i1,1,1,j1), w, cgs, kvh
     $,                ts, taux, tauy, hflx, qflx
c    $,                fsds, fsus, flds, vwcsfc, snowh, qrad, omga )
     $,                fsds, fsus, flds, vwcsfc, snowh, omga )

c        Date at end of timestep.
         call getcurdate( ncdate, ncsec )

c        Interpolate data to end of the timestep.
         call chkdyn( ncdate, ncsec, phis, oro )
         call interp1( ncdate, ncsec, ps(:,:,ntp1), t, q(:,:,:,ntp1)
     $,                zmu, zmd, zeu, hketa, hkbeta, cmfdqr, conicw
     $,                cldfr, concld, cme, prain, evapr, cwat )
#ifdef DEBUG
         call pgavg(' main: after interp1, ps,ntp1 ',ps(:,:,ntp1))
#endif

         if ( limqin ) then
c           Do not allow input specific humidity to exceed saturation.
            call limsh( ps(:,:,ntp1), t, q(:,:,:,ntp1) )
         end if

#ifdef DEBUG
         call pgavg(' main: before adjdry, ps,ntp1 ',ps(:,:,ntp1))
#endif
c        Adjust surface pressure to maintain constant dry mass.
         call adjdry( ps(:,:,ntp1), q(:,:,:,ntp1), psdry(:,:,ntp1)
     $,               .false. )
c        Overwrite ps with the dry ps.
         ps(1:plon,:,ntp1) = psdry(1:plon,:,ntp1)

#ifdef DEBUG
         call pgavg(' main: after adjdry, psdry and ps,ntp1 ',ps(:,:,ntp1))
         write(6,*)'DEBUG main 0.8'
#endif         

c        Write header records if a time sample is to be written to the 
c        history file.
         if ( wrhstts ) call wrhd( nstep )

c        ******************
c        * Advect Tracers *
c        ******************o

#ifdef DEBUG
         write(6,*)'DEBUG main 1.0 before adv_dr'
         call pgavg(' main: before adv_dr, psdry,ntp1 ',psdry(:,:,ntp1))
         call pgavg(' main: before adv_dr, psdry,nt ',psdry(:,:,nt))
#endif
         call adv_dr( nstep, dtime, uv, w
     $,               psdry(:,:,nt),   as(1,1,1,1,nt)
     $,               psdry(:,:,ntp1), as(1,1,1,1,ntp1)
     $,               hw1, nadv, x, x(1,1,1,1,2), mcar, mca, wr )
#ifdef DEBUG
         call pgavg(' main: after adv_dr, psdry,ntp1 ',psdry(:,:,ntp1))
         call pgavg(' main: after adv_dr, psdry,nt ',psdry(:,:,nt))
#endif


         if ( adj_winds() ) then  ! psdry is updated
c           Overwrite ps with the dry ps.
            ps(1:plon,:,ntp1) = psdry(1:plon,:,ntp1)
         end if
#ifdef DEBUG
         call pgavg(' main: after adv_dr, ps,ntp1 ',ps(:,:,ntp1))
#endif

         if ( fixmas ) then
c           Evaluate coefficients for advection mass fixer.
            call advfixcf3( pcnst, psdry(:,:,ntp1), as(i1,1,1,j1,nt)
     $,                     as(i1,1,1,j1,ntp1), hw1 )
         end if
#ifdef DEBUG
         call pgavg(' main: after advfixcf3, psdry,nt ',psdry(:,:,ntp1))
#endif

         call fixneg (as(i1,1,1,j1,ntp1))

c        ****************************************************************
c        * Physics/Chemistry parameterizations and write history buffer *
c        * to disk.                                                     *
c        ****************************************************************


         call physcan( dtime, calday, ncdate, ncsec, nstep
     $,                ps(:,:,ntp1), oro, phis, x(i1,1,1,j1,1), t
     $,                q(:,:,:,ntp1)
     $,                uv(i1,1,1,j1), ts, taux, tauy, hflx, qflx
c     $,                fsds, fsus, flds, vwcsfc, snowh, qrad
     $,                fsds, fsus, flds, vwcsfc, snowh
     $,                w, omga, cgs, kvh
     $,                zmu, zmd, zeu, hketa, hkbeta, cmfdqr, conicw
     $,                cldfr, concld, cme, prain, evapr, cwat
     $,                as(i1,1,1,j1,nt), as(i1,1,1,j1,ntp1)
     $,                x(i1,1,1,j1,2)
     $,                wr, mcar, mca 
     $)



c        Output mass diagnostics
         if ( globmass ) then
            call gamk2( as(i1,1,1,j1,ntp1), psdry(:,:,ntp1), trmk, mixk)
            call wrmass2nc( trmk, mixk )
         end if
#ifdef ASSIM
c++assim


         if (AODasm() ) then
c-drb            call Assimilate_AOD( ncdate, ncsec,ps(:,:,nt), x(i1,1,1,j1,2), t
c-drb   ps,t (x(i1,1,2,j1,2) changed to be consistent with physcan:physlic:optdep calls
            call Assimilate_AOD( ncdate, ncsec,ps(:,:,ntp1), x(i1,1,1,j1,2), x(i1,1,2,j1,2)
     $,                          as(i1,1,1,j1,ntp1), delt )
!     $,                          as(i1,1,1,j1,ntp1) )
         end if
c--assim
#endif


c++csz
#ifdef DST_MSS_BDG
c     Output global average mass budget for this timestep
         call bdg_update(nstep,ncdate,ncsec)
#endif /* not DST_MSS_BDG */
c--csz


c        Dispose the history file and open the next one.
         if ( closehst ) then
            call disphst( laststep, .false. )
            if ( fullhst ) call nexthist( laststep )
         end if



c        Write restart files and dispose.
         if ( wrestart ) then
            call wrrst( nstep, closehst, as(1,1,1,j1,ntp1), x(1,1,1,1,2)
     $,                 psdry(:,:,ntp1)  )
         end if

c        Output to summary file.
         write(*,*) 'done timestep ', nstep ,'  date = ', ncdate
     $,                '  seconds = ', ncsec

c        ************************
c        * Check for end of run *
c        ************************


         if ( endofrun ) then
c           Close mass diagnostics file
            if ( globmass ) then
               call clmass2nc( )
            end if
c++csz
#ifdef DST
            call dst_tvbds_close()

#endif /* not DST */
#ifdef DST_MSS_BDG
            call bdg_close()

#endif /* not DST_MSS_BDG */
c--csz

            stop
         endif

c        Toggle time level indices.

         nt   = ntp1
         ntp1 = mod( nt, 2 ) + 1


 1000 continue                  ! End loop over timesteps.

c##############################################################################
      contains
c##############################################################################

      subroutine inirun
   
c Initialize configuration, parameterizations, dynamics and trace specie arrays.

      use constants
      use acsf, only: iniacsf
      use advect, only: adv_type, ini_adv
      use dyninp, only: inidyn
      use histout, only: inihist
      use conv_hack, only: mfinti
      use conv_pjr, only: conv_ini
      use cloud 
      use vdiff_mod, only: vdinti
      use surface, only: inisflx

#ifdef SCYC
c++scyc
      use scyc, only: scycini, doscyc
      use acbnd, only: acbndini
      use dmsbnd, only: dmsbndini
      use soxbnd, only: soxbndini, soxsfconst, soxsfyr
      use sulchem, only: inisulchem
      use volcemission
c--scyc
#endif
#ifdef CAER
c++caer
      use caer, only: caerini, docaer
      use caerbnd, only: caerbndini
c--scyc
#endif

#ifdef ASSIM
c++assim
      use sat_obs, only: sat_obs_ini
c--assim
#endif

      use calendar

      implicit none
c---------------------------------------------------------------------
#include <control.com>
c---------------------------------------------------------------------
#include <tracnam.h>
c---------------------------------------------------------------------

c Local variables:

      character
     $  trnum*3         ! for converting integers to character strings
     $, case*16         ! case name from the environment variable CASE
     $, logname*16      ! user's name from the environment variable LOGNAME

      integer
     $  i, j, k, m   ! indices
     $, iret, istat  ! return values
     $, inerror   ! keep track of input errors
     $, lundyn    ! unit number for dynamics input
     $, luntra    ! unit number for tracer initial conditions input
     $, lunnas    ! unit number for out of core storage of non-advected species
     $, lunhst    ! unit number for output history file
     $, lunhbf    ! unit number for out of core history buffer
     $, lunrst0   ! unit number for restart info file
     $, lunrst1   ! unit number for restart data file
     $, lunsflx   ! unit number for surface flux data
     $, lunoro    ! unit number for oro field when not read from dynamics
     $, ncdate    ! current date in yymmdd format
     $, ncsec     ! current seconds relative to ncdate
     $, itimrst   ! time index of restart data
     $, rstdate   ! date of restart data
     $, rstdatesec! seconds relative to rstdate
     $, rstday    ! day of restart data
     $, rstdaysec ! seconds relative to rstday
     $, ntshst    ! number of time samples in initial history file

      real(r8)
     $  dlam         ! longitude spacing
     $, lam(plond)   ! longitudes in extended SLT grid (radians)
     $, phi(platd)   ! Gauss latitudes in extended SLT grid (radians)
     $, latwts(plat) ! Latitude weights.
     $, calday       ! calandar day for current timestep

      real(r8)
     $  hyai(plevp)   ! p0 component of hybrid coordinate - interfaces
     $, hyam(plev)    ! p0 component of hybrid coordinate - midpoints
     $, hybi(plevp)   ! ps component of hybrid coordinate - interfaces
     $, hybm(plev)    ! ps component of hybrid coordinate - midpoints
     $, hybd(plev)    ! difference of the hybi's
     $, etamid(plev)  ! hyam + hybm
     $, etaint(plevp) ! hyai + hybi
     $, hypi(plevp)   ! p0*etaint
     $, p0            ! reference pressure in Pa.
      parameter( p0 = 1.e5 )

      integer
     $  nprlev        ! index of lowest pure pressure interface

c     Control variables not stored in control.com (they are passed to the
c     individual modules).

      logical
     $  rdt    ! read T from the dynamics files
     $, rdoro  ! read ORO from the dynamics files
     $, rdphis ! read PHIS from dynamics files

      data rdt, rdphis / 2 * .false. /
      data rdoro / .true. /

c     Arrays to pass any user specified names or scale factors for the dynamics
c     input fields to inidyn.
      character
     $  dynnm(mxdynflds)*32 ! names of dynamics fields:
      real(r8)
     $  dynsc(mxdynflds)    ! scale factor to convert units in input files to
c                           ! the internal units.

c *******************
c * Namelist input: *
c *******************

c     namelist variables stored in /control0/
c      integer
c     $  rstflg      ! 0 for initial run, 1 for "restart"
c     $, nestep      ! ending timestep (pos) or day (neg)
c     $, delt        ! timestep interval in seconds
c     $, icdate      ! date of initial conditions in yyyymmdd format (time index 0)
c     $, icdatesec   ! seconds relative to icdate
c     $, icday       ! day at initial time (time index 0)
c     $, icdaysec    ! seconds relative to icday
c      logical
c     $  incore      ! t => incore history buffer and nas buffer
c     $, fixmas      ! t => apply global mass fixer
c     $, tvbnd       ! t => emissions from time varying boundary dataset
c     $, dowetdep    ! t => apply wet deposition parameterization
c     $, addqrad     ! t => add radiative heating tendency to the temperature
c                           field supplied to the convection routine
c     $, limqin      ! t => limit input Q s.t. it doesn't exceed saturation
c     $, globmass    ! t => write global tracer mass at each level into netCDF
c                    !      file at each timestep
c     $, dodrydep    ! t => apply dry deposition parameterization
c     $, readQIC     ! t => read initial conditions for the advected Q field
c                    !      rather than initializing from meteorology.
c     $, useAltConv  ! t => use alternative convective parameterization
c     $, wind_ps_adj ! t => adjust the winds to enforce consistency with input
c                    !      surface pressure field.
c     $, calcedot    ! t => calculate etadot
c     $, rdomega     ! t => read omega from dynamics input files
c      real(r8)
c     $  qrelax      ! Relaxation factor for "nudging" predicted Q field
c                    ! towards input Q field.  qrelax=0 => MATCH computes
c                    ! complete hydrologic cycle.  qrelax=1 => MATCH resets
c                    ! Q to input values at beginning of each timestep.

      integer
     $  trdate     ! date of as and nas initial conditions in yymmdd format
     $, trsec      ! seconds relative to trdate
     $, dynoffset  ! offset (seconds) added to (date,sec) for dynamics data
     $, nhtfrq     ! freq of hist file writes in timesteps (pos) or hours (neg)
     $, ndens      ! write density for output history tapes
     $, mfilt      ! max number of time samples written to hist tape
     $, irt        ! retention period in days on mass store for history files.
c                  ! If 0 (default) then output files not transfered.
     $, rstfrq     ! freq of restart file writes in timesteps (pos) or day (neg)
     $, rstrrt     ! retention period in days on mass store for restart files.

      logical
     $  limfix     ! true => limit fixer
     $, filtedot   ! true => apply filter to calculated etadot
     $, esice      ! true => saturation vapor pressures account for ice phase
     $, zmtrigon   ! true => turn on trigger mods for Zhang convection
     $, tracfrht   ! true => initialize tracers from a history file
     $, async      ! true => output files transferred to ms asynchronously
     $, linkhst    ! true => link completed hist files to linkdir/linkfn
     $, rmout      ! true => remove hist file after linking and/or transferring
c                  !         to mass store
     $, rmdyn      ! true => remove old dynamics files
     $, rmtra      ! true => remove tracer initial conditions file
     $, rmrst      ! true => remove restart data files

c     Names and scale factors for dynamics input data.
      character(len=32) ::
     $  unm=' '         !  1) zonal wind component                             
     $, vnm=' '         !  2) meridional wind component                        
     $, psnm=' '        !  3) surface pressure                                 
     $, oronm=' '       !  4) orography                                        
     $, phisnm=' '      !  5) surface geopotential                             
     $, etadotnm=' '    !  6) etadot (vertical wind component, d(eta)/dt )     
     $, omeganm=' '     !  7) omega  (vertical wind component, d(pressure)/dt )
     $, tnm=' '         !  8) temperature                                      
     $, qnm=' '         !  9) specific humidity                                
     $, tsnm=' '        ! 10) surface temperature (K)  Note: setting tsnm=tnm tells
c                             MATCH to use the bottom level of the T field for TS.
     $, tauxnm=' '      ! 11) surface zonal stress (n/m2)                      
     $, tauynm=' '      ! 12) surface meridional stress (n/m2)                 
     $, shflxnm=' '     ! 13) surface heat flux (w/m2)                         
     $, qflxnm=' '      ! 14) surface moisture flux (kg/m2/s)                  
     $, kvhnm=' '       ! 15) vertical diffusion coefficient                   
     $, cgsnm=' '       ! 16) counter-gradient coefficient                     
     $, zmunm=' '       ! 17) zhang mu2 from conv_ccm, kg/m2/s                 
     $, zmdnm=' '       ! 18) zhang md2 from conv_ccm, kg/m2/s                 
     $, zeunm=' '       ! 19) zhang eu2 from conv_ccm, 1/s                     
     $, hketanm=' '     ! 20) hack convective mass flux (kg/m^2/s)             
     $, hkbetanm=' '    ! 21) hack overshoot parameter (fraction)              
     $, qradnm=' '      ! 22) radiative heating tendency (k/s)                 
     $, cmfdqrnm=' '    ! 23) dq/dt due to convective rainout (kg/kg/s)
     $, conicwnm=' '    ! 24) convective in-cloud water mixing ratio (kg/kg)
     $, cldfrnm=' '     ! 25) total cloud fraction
     $, concldnm=' '    ! 26) convective cloud fraction
     $, cmenm=' '       ! 27) rate of cond-evap within the cloud (kg/kg/s)
     $, prainnm=' '     ! 28) dq/dt due to stratiform precip (kg/kg/s)
     $, evaprnm=' '     ! 29) rate of evaporation of strat precip (kg/kg/s)
     $, cwatnm=' '      ! 30) grid box cloud water (kg/kg)
     $, fsdsnm=' '      ! 31) downward solar flux at surface (W/m2)                  
     $, snowhnm=' '     ! 32) water equivalent snow depth (m)
     $, fsusnm=' '      ! 33) upward solar flux at surface (W/m2)
     $, fldsnm=' '      ! 34) downward longwave flux at surface (W/m2)
     $, vwcsfcnm=' '    ! 35) volumetric soil water content at surface (frac)
      real(r8) ::
     $  usc=1., vsc=1., pssc=1., orosc=1., phissc=1., etadotsc=1.
     $, omegasc=1., tsc=1., qsc=1., tssc=1., tauxsc=1., tauysc=1.
     $, shflxsc=1., qflxsc=1., kvhsc=1., cgssc=1., zmusc=1., zmdsc=1.
     $, zeusc=1., hketasc=1., hkbetasc=1., qradsc=1., cmfdqrsc=1.
     $, conicwsc=1., cldfrsc=1., concldsc=1., cmesc=1., prainsc=1.
     $, evaprsc=1., cwatsc=1., fsdssc=1.
     $, snowhsc=1., fsussc=1., fldssc=1., vwcsfcsc=1.

c     Names used to read initial conditions from a history file and write to
c     output history file.
c     These variables stored in /tracnm/
c      character
c     $  tracnam(pcnst)*8 ! advected species' names
c     $, sflxnam(pcnst)*8 ! surface flux names (for advected species)
c     $, natsnam(pnats)*8 ! non-advected species' names
c     $, tanam(pcnst)*8   ! total advection tendency names
c     $, vdnam(pcnst)*8   ! vertical diffusion tendency names
c     $, dcnam(pcnst)*8   ! moist convection tendency names
c     $, xflxnam(pcnst+2)*8 ! fluxes on x cell walls, tracers + air + Q
c     $, yflxnam(pcnst+2)*8 ! fluxes on y cell walls, tracers + air + Q
c     $, zflxnam(pcnst+2)*8 ! fluxes on z cell walls, tracers + air + Q

      real(r8)
     $  qmin(pcnst)     ! minimum tracer mixing ratio values.

      character*8
     $  outinst(mxoutflds)  ! names of instantaneous output fields
     $, outtimav(mxoutflds) ! names of time averaged output fields

c     netCDF surface flux data
      character
     $  lpthncsf*168       ! local pathname for netCDF surface flux data
     $, ncsfvar(pcnst)*32  ! ncsfvar(i) is the variable name in the netCDF
c                          ! file to use for tracnam(i).  Leave value blank
c                          ! if there is no data for tracnam(i).
      character
     $  lpthdyn*168     ! local path for dynamics data
     $, rpthdyn*168     ! remote path for dynamics data
     $, lpthtra*168     ! local path for tracer initial conditions
     $, rpthtra*168     ! remote path for tracer initial conditions
     $, rrootd*16       ! Default remote root directory for output history and
c                       ! restart files.  Defaults to the environment variable
c                       ! LOGNAME.
     $, lpthhst*168     ! local path for history file
     $, rpthhst*168     ! remote path for history file
     $, lpthrst0*128    ! local pathname for restart info file
     $, rpthrst0*128    ! remote pathname for restart info file
     $, lpthrst1*128    ! local pathname for restart data file
     $, rpthrst1*128    ! remote pathname for restart data file
     $, title*80        ! case title
     $, msvol*8         ! mass store virtual volume name
     $, wpasswd*8       ! MS write password
     $, lndirhst*128    ! link directory for output hist tapes
     $, caltype*16     ! caltype type.  Currently '365' or 'gregorian'

#ifdef SCYC
c++scyc
      logical volcEmis
      character*4 ncnam(5)
! soxsfconst and soxsfyr declared in soxbnd.F90
!      logical soxsfconst ! true => surface flux constant. Set year with soxsfyr
!      integer soxsfyr    ! used if soxsfconst, defaults to most recent year on SOx file
c--scyc
#endif
#ifdef ASSIM
c++assim
      character*120 assimDir ! Directory containing satellite data
c      logical doAODasm       ! t => invoke AOD assimilation code
c--assim
#endif

      namelist /nlist/
     $  incore, fixmas, limfix, filtedot, esice, zmtrigon, limqin
     $, qrelax
     $, dcnam, xflxnam, yflxnam, zflxnam
     $, rstflg, nestep, delt, icdate, icdatesec, trdate, trsec
     $, dynoffset, caltype, icday, icdaysec
     $, tracfrht,  qmin 
     $, lpthdyn, rpthdyn, lpthtra, rpthtra, rrootd, lpthhst, rpthhst
     $, title, nhtfrq, ndens, mfilt, irt, async
     $, msvol, wpasswd, linkhst, lndirhst, rmout, rmdyn, rmtra
     $, rstfrq, lpthrst0, rpthrst0, lpthrst1, rpthrst1, rmrst, rstrrt
     $, tvbnd, dowetdep, rdoro, rdomega, addqrad, globmass
     $, dodrydep, readQIC, useAltConv, wind_ps_adj, calcedot
     $, vwc_scale, vwc_offset, sslt_scale
     $, lpthncsf, ncsfvar
     $, unm, vnm, psnm, oronm, phisnm, etadotnm, omeganm, tnm, qnm, tsnm
     $, tauxnm, tauynm, shflxnm, qflxnm, kvhnm, cgsnm, zmunm, zmdnm
     $, zeunm, hketanm, hkbetanm, qradnm, cmfdqrnm, conicwnm
     $, cldfrnm, concldnm, cmenm, prainnm, evaprnm, cwatnm
     $, fsdsnm, snowhnm, fsusnm, fldsnm
     $, vwcsfcnm
     $, usc, vsc, pssc, orosc, phissc, etadotsc, omegasc, tsc, qsc, tssc
     $, tauxsc, tauysc, shflxsc, qflxsc, kvhsc, cgssc, zmusc, zmdsc
     $, zeusc, hketasc, hkbetasc, qradsc, cmfdqrsc, conicwsc
     $, cldfrsc, concldsc, cmesc, prainsc, evaprsc, cwatsc
     $, fsdssc, snowhsc, fsussc, fldssc
     $, vwcsfcsc,  tanam, vdnam
#ifdef SCYC
c++scyc
     $, volcEmis, soxsfconst, soxsfyr
c--scyc
#endif
#ifdef ASSIM
     $,     assimDir,doAODasm
#endif
     $, tracnam
     $, sflxnam
     $, natsnam
     $, outinst
     $, outtimav


      real(r8) temp

c Externals:
      character*16
     $  upcase    ! Convert from lower case to upper case.
!      real(r8)
!     $  caldayr   ! Caltype day (includes fractional part of day)
      logical
     $  adj_winds ! t => adjust winds to enforce mass conservation
c-----------------------------------------------------------------------

      call init_const()

c     Set defaults for namelist variables.  Initialize local variables.

      inerror = 0

      rstflg    = 0
      itim0     = 0
      delt      = 0
      nestep    = 0
      icdate    = 0
      icdatesec = 0
      icday     = 0
      icdaysec  = 0
      trdate    = -1
      trsec     = 0
      lndirhst  = ' '
      lpthncsf  = ' '
      incore    = .true.
      filtedot  = .true.
      esice     = .true.
      zmtrigon  = .true.
      limqin    = .false.
      tracfrht  = .false.
      dynoffset = 0
      caltype  = 'gregorian'
      tracnam = 'not-set'
      sflxnam = 'not-set'
      tanam = 'not-set'
      vdnam = 'not-set'
      dcnam    = 'not-set'
      xflxnam  = 'not-set'
      yflxnam  = 'not-set'
      zflxnam  = 'not-set'

      do m = 1, pcnst
         qmin(m) = 0.0
         write( trnum, '(i3)' ) m+100
c-drb         tracnam(m)  = 'TR'//trnum(2:3)
c-drb         sflxnam(m)  = 'SF'//trnum(2:3)
         tracnam(m)  = 'DEFAULT'
         sflxnam(m)  = 'DEFAULT'
         ncsfvar(m)  = ' '
         tanam(m)    = 'TA'//trnum(2:3)
         vdnam(m)    = 'VD'//trnum(2:3)
         dcnam(m)    = 'DC'//trnum(2:3)
         xflxnam(m)  = 'XFLX'//trnum(2:3)
         yflxnam(m)  = 'YFLX'//trnum(2:3)
         zflxnam(m)  = 'ZFLX'//trnum(2:3)
      end do
      xflxnam(pcnst+1) = 'XFLXAIR'
      yflxnam(pcnst+1) = 'YFLXAIR'
      zflxnam(pcnst+1) = 'ZFLXAIR'
      xflxnam(pcnst+2) = 'XFLXQ'
      yflxnam(pcnst+2) = 'YFLXQ'
      zflxnam(pcnst+2) = 'ZFLXQ'
      do m = 1, mxoutflds
         outinst(m) = '        '
         outtimav(m) = '        '
      end do

      rrootd   = '                '
      ntshst   = 0
      title    = 'MATCH'
      nhtfrq   = 0
      ndens    = 1
      mfilt    = 1
      irt      = 0
      rstrrt   = 0
      async    = .true.
      msvol    = '        '
      wpasswd  = '        '
      linkhst  = .false.
      rmout    = .false.
      rmdyn    = .false.
      rmtra    = .false.
      rmrst    = .false.
      rstfrq   = 0
      lpthdyn  = ' '
      rpthdyn  = ' '
      lpthtra  = ' '
      rpthtra  = ' '
      lpthhst  = ' '
      rpthhst  = ' '
      rpthrst0 = ' '
      lpthrst0 = ' '
      rpthrst1 = ' '
      lpthrst1 = ' '
c     Logical unit numbers for Fortran I/O

      lundyn = 1
      luntra = 2
      lunnas = 3
      lunhbf = 4
      lunhst = 7
      lunrst0 = 10
      lunrst1 = 11
      lunsflx = 12
      lunoro  = 13

c     Initialize control variables.

      vwc_scale  = 1.0
      vwc_offset = 0.0
      sslt_scale = 1.0
      divdiff = .false.
      arvdiff = .false.
      diconvccm = .false.
      arconvccm = .false.
      dicldphys = .false.
      arcldphys = .false.
      vdiffuse = .false.
      convect = .false.
      chemistry = .false.
      physics = .false.
      advqth = .false.
      rdomega = .false.
      satvap = .false.
      tvbnd = .false.
      dowetdep = .false.
      addqrad = .false.
      globmass = .false.
      dodrydep = .false.
      readQIC  = .false.
      useAltConv = .false.
#ifdef SCYC
c++scyc
      volcEmis = .false.
c--scyc
#endif
#ifdef ASSIM
c++assim
      doAODasm = .true.
      assimDir = 'not-set'
c--assim
#endif
      calcedot = .true.
      if ( adv_type() == 'spitfire' ) then
         wind_ps_adj = .true.
         fixmas = .false.
      else
         wind_ps_adj = .false.
         fixmas = .true.
      end if
      limfix    = .true.

#ifdef DI_VDIFF
      divdiff = .true.
      advqth = .true.
      satvap = .true.
#endif
#ifdef AR_VDIFF
      arvdiff = .true.
#endif
#ifdef DI_CONV_CCM
      diconvccm = .true.
      advqth = .true.
#endif
#ifdef AR_CONV_CCM
      arconvccm = .true.
#endif
#ifdef DI_CLOUD_PHYS
      dicldphys = .true.
      advqth = .true.
      tracnam(pcnst) = 'CWAT    '
#endif
#ifdef AR_CLOUD_PHYS
      arcldphys = .true.
#endif

      if ( divdiff .or. arvdiff ) vdiffuse = .true.
      if ( diconvccm .or. arconvccm ) convect = .true.

      if ( vdiffuse .or. convect .or. chemistry ) then
         physics = .true.
         rdt = .true.
         rdphis = .true.
      end if

      qrelax = 0.
      if ( .not. dicldphys ) qrelax = 1.

c     Check for inconsistent options in CPP variable specifications.

      if ( diconvccm .and. .not. divdiff ) then
         write(*,*)'inirun: diagnosed CCM convection uses pblh, tpert'
         write(*,*)'  and qpert calculated by diagnosed vertical'
         write(*,*)'  diffusion scheme.  Must define DI_VDIFF.'
         inerror = inerror + 1
      end if

      if ( dicldphys .and. .not. diconvccm ) then
         write(*,*)'inirun: diagnosed cloud physics requires diagnosed'
         write(*,*)'  convection.  Must define DI_CONV_CCM'
         inerror = inerror + 1
      end if

      if ( dowetdep .and. (.not. dicldphys .and. .not. arcldphys)) then
         write(*,*)'inirun: wet deposition requires cloud physics.'
         write(*,*)'  Must define DI_CLOUD_PHYS or AR_CLOUD_PHYS.'
         inerror = inerror + 1
      end if

      if ( dicldphys .and. pcnst .eq. 1 ) then
         write(*,*)'inirun: diagnosed cloud physics uses tracer pcnst'
         write(*,*)'  for the cloud water.  No space has been allocated'
         write(*,*)'  for user defined tracers.  Must define PCNST > 1'
         inerror = inerror + 1
      end if


c     *************************************
c     * Read namelist from standard input *
c     *************************************
      write (6,*) ' ready to read nlist '
      read( *, nlist )
      write( *, nlist )
c      read (*, tracernames)
c      write (*, tracernames)
c      read (*, sflxnames)
c      write (*, sflxnames)
c      read (*, inames)
c      write (*, inames)
c      read (*, tanames)
c      write (*,tanames)

c     Some checking for invalid namelist input:

      if ( delt .eq. 0 ) then
         write(*,*)'inirun: Must specify timestep size.'
         write(*,*)'        Use the DELT namelist variable.'
         inerror = inerror + 1
      end if

c     Cast "nestep", "nhtfrq" and "rstfrq" in terms of time step
      if ( nestep .lt. 0 ) nestep = -nestep*86400/delt
      if ( rstfrq .lt. 0 ) rstfrq = -rstfrq*86400/delt
      if ( nhtfrq .lt. 0 ) nhtfrq = -nhtfrq*3600/delt

      if ( icdate .eq. 0 ) then
         write(*,*)'inirun: Must specify date of initial conditions.'
         write(*,*)'        Use the ICDATE namelist variable.'
         inerror = inerror + 1
      end if
      if ( trdate .eq. -1 ) then
         trdate = icdate
         trsec = icdatesec
      end if
      if ( nestep .eq. 0 ) then
         write(*,*)'inirun: Must specify non-zero number of timesteps.'
         write(*,*)'        Use the NESTEP namelist variable.'
         inerror = inerror + 1
      end if
      if ( rstfrq .eq. 0 ) rstfrq = mfilt*nhtfrq
      if ( outinst(1) .eq. '        ' .and.
     $     outtimav(1) .eq. '        ' ) then
         write(*,*)'inirun: Must specify output fields using either'
         write(*,*)'        outinst or outtimav (namelist)'
         inerror = inerror + 1
      end if
      if ( lpthdyn(1:1) .eq. ' ' .and. rpthdyn(1:1) .eq. ' ' ) then
         write(*,*)'inirun: Must specify dynamics datasets using the'
         write(*,*)'        lpthdyn and/or rpthdyn namelist variables'
         inerror = inerror + 1
      end if
      if ( tracfrht ) then
         if ( lpthtra(1:1) .eq. ' ' .and. rpthtra(1:1) .eq. ' ' ) then
            write(*,*)'inirun: Must specify tracer initial condition'
            write(*,*)'        datasets using the lpthtra and/or'
            write(*,*)'        rpthtra namelist variables'
            inerror = inerror + 1
         end if
      end if
      if ( .not. dicldphys .and. qrelax .ne. 1 ) then
         write(*,*)'inirun: must run with fully diagnosed physics to'
         write(*,*)'  compute the hydrologic cycle.  Set qrelax=1.0.'
         inerror = inerror + 1
      end if
      if ( wind_ps_adj .and. .not. calcedot ) then
         write(*,*)'inirun: doing wind/ps adjustments to enforce mass'
         write(*,*)'  conservation ==> etadots must be calculated.'
         write(*,*)'  Either set wind_ps_adj=.false. or set calcedot'
         write(*,*)'  to .true. (which is the default value).'
         inerror = inerror + 1
      end if
      if ( wind_ps_adj .and. (adv_type() == 'slt') ) then
         write(*,*)'inirun: inconsistent options specified.  Cannot'
         write(*,*)'  do wind/ps adjustments with SLT advection.'
         write(*,*)'  Either set wind_ps_adj=.false. (default for SLT)'
         write(*,*)'  or use SPITFIRE advection.'
         inerror = inerror + 1
      end if
!      if ( soxsfconst .and. (soxsfyr == 0 ) ) then
!         write(*,*)'inirun: constant sox surface flux requested without'
!         write(*,*)'      specifying which year to use.'
!         write(*,*)'      Either set soxsfconst = .true. (default) or'
!         write(*,*)'      set soxsfyr = 1990,2005,2020,2035,2065 or 2095'
!         inerror = inerror + 1
!      end if

c     Get environment variables for CASE and LOGNAME
      call getenv( 'CASE', case )
      if ( case .eq. ' ' ) then
         write(*,*)'inirun: CASE not found.  Set to debug'
         case = 'debug'
      end if
      logname = ' '
      call getenv( 'LOGNAME', logname )
      write (6,*) ' logname >',logname,'<'
      if ( logname .eq. ' ' ) then
         write(*,*)'inirun: LOGNAME not set'
      else
         if ( rrootd(1:1) .eq. ' ' ) rrootd = upcase( logname )
      end if

      if ( rpthhst(1:1) .eq. ' '  .and.  irt .ne. 0  .and. 
     $     rrootd(1:1) .eq. ' ' ) then
         write(*,*)'inirun: Do not have remote path for disposing'
         write(*,*)'        history files.  Use namelist variables'
         write(*,*)'        RPTHHST or RROOTD'
         inerror = inerror + 1
      end if
      if ( (rpthrst0(1:1) .eq. ' ' .or. rpthrst1(1:1) .eq. ' ')
     $     .and.  irt .ne. 0  .and.  rrootd(1:1) .eq. ' ' ) then
         write(*,*)'inirun: Do not have remote path for disposing'
         write(*,*)'        restrart files.  Use namelist variables'
         write(*,*)'        RPTHRST0/RPTHRST1 or RROOTD'
         inerror = inerror + 1
      end if

c     Exit if input error(s) found.
      if ( inerror .gt. 0 ) stop

      nadv = 0
      if ( advqth ) nadv = padv

c     Allocate storage for field arrays.
      if ( advqth ) then
         allocate( x(plond,plev,padv,platd,2), stat=istat )
         call alloc_err( istat, 'inirun', 'x', plond*plev*padv*platd*2 )
         x = 0.
      end if
      if ( vdiffuse ) then
         allocate( cgs(plond,plevp,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'cgs', plond*plevp*plat )
         allocate( kvh(plond,plevp,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'kvh', plond*plevp*plat )
      end if
      if ( diconvccm .or. arconvccm ) then
         allocate( zmu(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'zmu', plond*plev*plat )
         allocate( zmd(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'zmd', plond*plev*plat )
         allocate( zeu(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'zeu', plond*plev*plat )
         allocate( hketa(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'hketa', plond*plev*plat )
         allocate( hkbeta(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'hkbeta', plond*plev*plat )
         allocate( cmfdqr(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'cmfdqr', plond*plev*plat )
         allocate( conicw(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'conicw', plond*plev*plat )
      end if
      if ( addqrad ) then
         allocate( qrad(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'qrad', plond*plev*plat )
      else
         allocate( qrad(1,1,1), stat=istat )
         call alloc_err( istat, 'inirun', 'qrad', 1 )
      end if
      if ( dicldphys .or. arcldphys ) then
         allocate( cldfr(plond,plevp,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'cldfr', plond*plevp*plat )
         allocate( concld(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'concld', plond*plev*plat )
         allocate( cme(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'cme', plond*plev*plat )
         allocate( prain(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'prain', plond*plev*plat )
         allocate( evapr(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'evapr', plond*plev*plat )
         allocate( cwat(plond,plev,plat), stat=istat )
         call alloc_err( istat, 'inirun', 'cwat', plond*plev*plat )
      end if

c     Initialize caltype calculations.
      call inicalendar( caltype )

c     Initialize dynamics input data names and scale factors.
      dynnm = ' '
      dynnm(1) = unm      
      dynnm(2) = vnm      
      dynnm(3) = psnm     
      dynnm(4) = oronm    
      dynnm(5) = phisnm   
      dynnm(6) = etadotnm 
      dynnm(7) = omeganm  
      dynnm(8) = tnm      
      dynnm(9) = qnm      
      dynnm(10) = tsnm
      dynnm(11) = tauxnm   
      dynnm(12) = tauynm   
      dynnm(13) = shflxnm  
      dynnm(14) = qflxnm   
      dynnm(15) = kvhnm    
      dynnm(16) = cgsnm    
      dynnm(17) = zmunm    
      dynnm(18) = zmdnm    
      dynnm(19) = zeunm    
      dynnm(20) = hketanm 
      dynnm(21) = hkbetanm
      dynnm(22) = qradnm   
      dynnm(23) = cmfdqrnm
      dynnm(24) = conicwnm
      dynnm(25) = cldfrnm
      dynnm(26) = concldnm
      dynnm(27) = cmenm
      dynnm(28) = prainnm
      dynnm(29) = evaprnm
      dynnm(30) = cwatnm
      dynnm(31) = fsdsnm
      dynnm(32) = snowhnm
      dynnm(33) = fsusnm
      dynnm(34) = fldsnm
      dynnm(35) = vwcsfcnm
      dynsc = 1.
      dynsc(1) = usc      
      dynsc(2) = vsc      
      dynsc(3) = pssc     
      dynsc(4) = orosc    
      dynsc(5) = phissc   
      dynsc(6) = etadotsc 
      dynsc(7) = omegasc  
      dynsc(8) = tsc      
      dynsc(9) = qsc      
      dynsc(10) = tssc
      dynsc(11) = tauxsc   
      dynsc(12) = tauysc   
      dynsc(13) = shflxsc  
      dynsc(14) = qflxsc   
      dynsc(15) = kvhsc    
      dynsc(16) = cgssc    
      dynsc(17) = zmusc    
      dynsc(18) = zmdsc    
      dynsc(19) = zeusc    
      dynsc(20) = hketasc 
      dynsc(21) = hkbetasc
      dynsc(22) = qradsc   
      dynsc(23) = cmfdqrsc
      dynsc(24) = conicwsc
      dynsc(25) = cldfrsc
      dynsc(26) = concldsc
      dynsc(27) = cmesc
      dynsc(28) = prainsc
      dynsc(29) = evaprsc
      dynsc(30) = cwatsc
      dynsc(31) = fsdssc
      dynsc(32) = snowhsc
      dynsc(33) = fsussc
      dynsc(34) = fldssc
      dynsc(35) = vwcsfcsc

c     Initialize restart.
      call inirest( rstflg, lunrst0, lunrst1, rstrrt, rstfrq
     $,             nestep, rrootd, case, lpthrst0, rpthrst0
     $,             lpthrst1, rpthrst1, async, rmrst
     $,             diconvccm, dicldphys, advqth
     $,             itimrst, rstdate, rstdatesec, rstday, rstdaysec
     $,             ntshst, lpthdyn, rpthdyn, lpthhst, rpthhst
     $)

c     Initialize control variables.
      call inictl( itimrst, rstdate, rstdatesec, rstday, rstdaysec )
      call getcurdate( ncdate, ncsec )

      write (6,*) ' b inidyn: rdoro ', rdoro
c     Initialize dynamics data and get coordinate data.
      call inidyn( plon, plat, plev, plond, platd, plevp, mxdynflds
     $,            lpthdyn, rpthdyn, ncdate, ncsec, dynoffset, lundyn
     $,            rmdyn
     $,            rdt, rdoro, rdphis, calcedot, rdomega
     $,            arvdiff, divdiff, arconvccm, addqrad
     $,            arcldphys
     $,            dynnm, dynsc
     $,            phi(j1), latwts, hyam, hybm, hyai, hybi, phis, oro
     $)

      write (6,*) ' hyai is ', hyai
      write (6,*) ' hybi is ', hybi
      write (6,*) ' hyam is ', hyam
      write (6,*) ' hybm is ', hybm
      if ( ( adv_type() .eq. 'spitfire' ) .and. ( hybi(1).ne.0. ) ) then
         write (6,*) 'WARNING: top interface set to zero as interrum fix for spitfire ', hybi(1)
         hybi(1) = 0.
      endif
      if ( ( adv_type() .eq. 'slt' ) .and. ( hybi(1).eq.0. ) ) then
         hybi(1) = 0.01*hybm(1)
         write (6,*) 'WARNING: top interface was zero, set to ', hybi(1), 'for slt'
      endif

      
      temp = 0.
      do j = 1,plat
         temp = temp + latwts(j)
      end do
      write (6,66)  temp
   66 format (' sum of gaussian weights is ', f15.13)

c     Check that input data interval contains initial time.  If not, read
c     time samples until the interval is found.
      call chkdyn( ncdate, ncsec, phis, oro )
c     Interpolate dynamics to initial time.
      call interp1( ncdate, ncsec, ps(:,:,nt), t, q(:,:,:,nt)
     $,             zmu, zmd, zeu, hketa, hkbeta, cmfdqr, conicw
     $,             cldfr, concld, cme, prain, evapr, cwat
     $)

      call interp2( ncdate, ncsec, uv(i1,1,1,j1), w, cgs, kvh
     $,             ts, taux, tauy, hflx, qflx
c    $,             fsds, fsus, flds, vwcsfc, snowh, qrad
     $,             fsds, fsus, flds, vwcsfc, snowh
     $,             omga
     $)

c     Initialize horizontal grid coordinates on the extended SLT grid.
      call grdxy( pi, dlam, lam, phi, latwts )

      if ( .not. rdoro ) then
         write(*,*)'inirun: ****WARNING**** the oro field has not been'
         write(*,*)'        initialized.  It is required by the surface'
         write(*,*)'        flux and diagnosed cloud physics modules'
c         call inioro( lunoro, phi(j1), oro )
      end if

c     Initialize vertical coordinates.
      call hycoef( hyai, hyam, hybi, hybm
     $,            hybd, nprlev, etamid, etaint )
      do k = 1, plevp
         hypi(k) = etaint(k)*p0
      end do

c     Initialize pressure levels calculation.
      call iniplevs( hyai, hyam, hybi, hybm )

      if ( satvap .or. limqin ) then
c        Initialize saturation vapor pressure tables.
         call esinti( epsilo, latvap, latice, rh2o, cpair, esice )
      end if

      if ( limqin ) then
c        Do not let input specific humidity exceed saturation.
         call limsh( ps(:,:,nt), t, q(:,:,:,nt) )
      end if

c     Initialize mass integral calculations.
      call inimass( gravit, latwts, etamid )

c     Initialize virtual temperature calculation.
      call inivtem( zvir )

c     potential temperature calculation.
      call inipottem( cappa )

c     Initialize geopotential height calculation.
      call iniz( rair, gravit )

c     Initialize surface flux calculation.
      call inisflx( rair, gravit, rearth, lam(i1), phi(j1), latwts, oro)
      calday = caldayr( ncdate, ncsec )
      call iniacsf( lpthncsf, calday, ncsfvar )

c     Initialize advection scheme.
      call ini_adv( dlam, lam, phi, latwts, etaint, etamid, hybi
     $,             nadv )
      if ( fixmas ) then
c        Initialize advection mass fixer.
         call iniadvfx( limfix, etamid )
      end if

c     Initialize minimum mixing ratio value enforcer.
      call iniqneg( qmin )

      if ( arconvccm .or. diconvccm ) then
c        Initialize constants used in both Zhang and Hack schemes (comcmf.h).
         call mfinti(rair, cpair, gravit, latvap, rhoh2o, hypi )
         if ( useAltConv ) then
c           Initialization for modified ZM convection w/ triggers
            call conv_ini( gravit, cpair, rair, hypi, zmtrigon )
         end if
      end if

c     Initialize cloud physics calculations.
      if ( arcldphys .or. dicldphys ) then
         call inicld( phi(j1), cappa, gravit, rair, latvap, cpair
     $,               rhoh2o/1000., tfh2o, etamid )
      end if
      if ( dicldphys ) then
         call inimc( etamid )
         call inimland( oro, phi(j1) )
      end if

      if ( calcedot .and. .not. wind_ps_adj ) then
c        Initialize divergence calculation.
         call inidiv( pi, phi(j1), latwts, dlam, rearth )

c        Initialize etadot calculation.
         call inietd( rearth, dlam, nprlev, hybi, hybd, etamid
     $,               filtedot, pi, phi )
      end if

      if ( vdiffuse ) then
c        Initialize vertical diffusion.
         call vdinti( cpwv, cpair, gravit, rair, qmin, hyam, hybm )
      end if

#ifdef DRYDEP
c     Initialize Dry Deposition.
      if ( dodrydep ) then
         call inidrydep( rair, gravit, phi(j1) )
      end if
#endif

#ifdef SCYC
c++scyc
      ! add namelist controls for the first 2 args later.
      call scycini( .true., 'SMITH', volcEmis )
      if ( doScyc() ) then

c-drb         ncyear = ncdate/10000
c-drb         if ( ncyear .lt. 1000 ) ncyear = ncyear + 1900
c-drb         call soxbndini( ncyear )

         if ( (soxsfconst).and.(soxsfyr.eq.0) ) soxsfyr = ncdate/10000
         if ( .not. soxsfconst ) soxsfyr = ncdate/10000
         if ( soxsfyr .lt. 1000 ) soxsfyr = soxsfyr + 1900

         call soxbndini( soxsfyr )
         call dmsbndini( calday )
         if ( volcEmis ) call volcemisini()
         call inisulchem()
         ncnam(1) = 'O3'
         ncnam(2) = 'OH'
         ncnam(3) = 'HO2'
         ncnam(4) = 'NO3'
         ncnam(5) = 'H2O2'
         call acbndini( 'oxid_3d.nc', calday, 5, ncnam )
      end if
c--scyc
#endif
c++csz
#ifdef DST
c     Initialize dust names
c     NB: dst_nm_cmn_ini() must be called before MATCH:src/deffld()
      call dst_nm_cmn_ini()
c     Initialize time-invariant physical constants
      call dst_cst_cmn_ini(rair,gravit,rearth)
c     Initialize size grid
      call dst_szdst_ini()
c     Initialize miscellaneous common blocks
      call dst_msc_cmn_ini()
#ifdef DST_MSS_BDG
c     NB: MATCH reads delt as integer, CCM reads dtime as real(r8)
c     Initialize mass budget common block
      call bdg_cmn_ini(gravit,latwts,real(delt,r8))
#endif /* not DST_MSS_BDG */
c     Initialize time-invariant boundary data from netCDF file
      call dst_tibds_ini('dst_bnd.nc')
c     NB: dst_tibds_ini() opens and closes the file, while
c     dst_tvbds_ini() opens the file and leaves it open
c     Initialize time-varying seasonal cycle data from netCDF file
      call dst_tvbds_ini('dst_bnd.nc',calday)
c     Initialize fields derived from external datasets
c     call dst_lsm_ini()
#endif /* not DST */
c--csz
#ifdef CAER
c++caer
      if ( docaer() ) then
         call caerini()
         call caerbndini( calday )
      end if
c--caer
#endif
#ifdef AER
c++sslt
c     initialize the code for calculating a 10m windspeed
      call initu10 (gravit, cpair, cpwv, rair, zvir)
c--sslt
#endif
#ifdef ASSIM
c++assim
      if ( doAODasm ) then
         write(6,*)'Initializing assimilation'
         call sat_obs_ini(doAODasm,assimDir,ncdate, ncsec )
      end if
c--assim
#endif

c     Initialize output history file and buffer.
      call inihist( lpthhst, rpthhst, rrootd, case, lunhst
     $,             ntshst, outinst, outtimav, nestep, nhtfrq
     $,             ndens, mfilt, icday, icdaysec, icdate
     $,             icdatesec, delt, rstflg, title, pi
     $,             phi(j1), latwts, hyam, hybm, hyai
     $,             hybi, irt, async, rmout, wpasswd
     $)

c     Initialize global mass monitor.  (Done here because the control module
c     is queried for time info.)
      if ( globmass ) then
         call inimass2nc( etamid, tracnam )
      end if

      if ( rstflg .eq. 0 ) then

c        Initialize tracers.
         if ( tracfrht ) then
            call initrac( lpthtra, rpthtra, trdate, trsec, luntra
     $,                   rmtra, chemistry, readQIC
     $,                   as(1,1,1,1,nt), x(1,1,1,1,2) )
         else
            call initr0( q(:,:,:,nt), as(i1,1,1,j1,nt), etamid )
         end if

c        Initialize dry mass adjustment calculation.
         call inidry( ps(:,:,nt), q(:,:,:,nt), psdry(:,:,nt) )

c        Overwrite ps with the dry ps.
         ps(1:plon,:,nt) = psdry(1:plon,:,nt)

      else if ( rstflg .eq. 1 ) then

c        Read restart data file.
         call rdrst( as(1,1,1,j1,nt), x(1,1,1,1,2), psdry(:,:,nt) )
         write (6,*) ' aft rdrst q, t', x(i1+70-1,18,:,j1+2-1,:)

         if ( .not. adj_winds() ) then
c           Adjust surface pressure to maintain constant dry mass.
c           (Need the atm dry mass read from restart file.)
            call adjdry( ps(:,:,nt), q(:,:,:,nt), psdry(:,:,nt), .true.)
         end if

c        Overwrite ps with the dry ps.
         ps(1:plon,:,nt) = psdry(1:plon,:,nt)

      end if

      end subroutine inirun

      subroutine fixneg (array)

      use precision
 
      implicit none

c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------


      real(r8) array (plond,plev,pcnst,plat) ! constituent field

c----------------------------- local variables -------------------------
      integer i,j,k,m                                                                                   
      integer ip1,jp1,kp1, im1, jm1, km1
      real(r8) surround

      do m = 1,pcnst
         do k = 1,plev
            km1 = max(1,k-1)
            kp1 = min(plev,k+1)
            do j = 1,plat
               jm1 = max(1,j-1)
               jp1 = min(plat,j+1)
               do i = 1,plon
                  im1 = max(1,i-1)
                  ip1 = min(plon,i+1)
                  surround = max(abs(array(i,k,m,jp1)),
     $                           abs(array(ip1,k,m,j)),
     $                           abs(array(im1,k,m,j)),
     $                           abs(array(i,kp1,m,j)),
     $                           abs(array(i,km1,m,j)),
     $                           abs(array(i,k,m,jm1))
     $                          )
                  if (abs(array(i,k,m,j)).lt.1.e-7*surround) then
                     array(i,k,m,j) = 0.
                  endif
               end do
            end do
         end do
!         write (6,*) ' fixneg range for m ',m, minval(array(1:plon,1:plev,m,1:plat)), 
!     $                                         maxval(array(1:plon,1:plev,m,1:plat))
      end do

      end subroutine fixneg

      end program main
