c $Id$
      module dyninp

      use calendar
      use precision
      use netcdf

c Module for reading meteorological driving data in netCDF format.

      implicit none
      save

c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

      character(len=168) ::
     $  lpath             ! local path
     $, rpath             ! remote pathname
      character(len=32) ::
     $  nam(mxdynflds)    ! names of dynamics fields:
     $, tname             ! work variable

      logical ::
     $  rmdyn      ! .true. => remove old dynamics files
     $, rdt        ! read T
     $, rdoro      ! read ORO
     $, rdphis     ! read PHIS
     $, calcedot   ! f => read ETADOT
     $, rdomega    ! read OMEGA
     $, limqin     ! true => limit input Q s.t. it doesn't exceed sat mr.
     $, arvdiff    ! read CGS and KVH for archived vertical diffusion
     $, divdiff    ! read variables used to diagnose vertical diffusion
                   ! (TAUX,TAUY,SHFLX,LHFLX)
     $, arconvccm  ! read ZMU, ZMD, ZEU, CONBETA, CONETA for archived CCM
                   ! convection
     $, arcldphys  ! read fields for archived cloud physics
     $, addqrad    ! t => add radiative heating tendency to the temperature
                   !      field supplied to the convection routine
     $, readfsds   ! t => read FSDS
     $, readsnowh  ! t => read SNOWH
     $, readfsus   ! t => read FSUS
     $, readflds   ! t => read FLDS
     $, readvwcsfc ! t => read VWCSFC

      integer ::
     $  ncid              ! netCDF file ID
     $, vid(mxdynflds)    ! IDs for dynamics fields
     $, nflds             ! number of fields to be read from input file
     $, namidx(mxdynflds) ! indices into nam and vid arrays for all fields
                          ! to be read from input file
     $, datevid(2)        ! IDs for date and datesec
     $, nlon              ! number of longitudes
     $, nlat              ! number of latitudes
     $, nlev              ! number of model layers
     $, ntim              ! number of time samples in current file
     $, dyndate1          ! date previous to current date (in yyyymmdd format)
     $, dynsec1           ! seconds relative to dyndate1
     $, ndyn1             ! index in data arrays for previous time sample
     $, dyndate2          ! current date (in yyyymmdd format)
     $, dynsec2           ! seconds relative to dyndate2
     $, ndyn2             ! index in data arrays for current time sample
     $, curr_ts           ! time index for most recently read data
     $, dynoffset         ! offset (seconds) added to (date,sec) info

      real(r8) ::
     $  uvin(plon,plev,2,plat,2) ! u,v winds (input, m/s)
     $, psin(plon,plat,2)        ! surface pressure (input, Pa)
     $, tin(plon,plev,plat,2)    ! temperature (input, K)
     $, qin(plon,plev,plat,2)    ! specific humidity (input, kg/kg)
     $, tsin(plon,plat,2)        ! surface temperature (input, K)
     $, tauxin(plon,plat,2)      ! surface zonal stress (input, N/m2)      
     $, tauyin(plon,plat,2)      ! surface meridional stress (input, N/m2) 
     $, hflxin(plon,plat,2)      ! surface heat flux (input, W/m2)         
     $, qflxin(plon,plat,2)      ! surface moisture flux (input, kg/m2/s)  
     $, fsdsin(plon,plat,2)      ! downward solar flux at surface (input, W/m2)  
     $, snowhin(plon,plat,2)     ! water equivalent snow depth (input, m)  
     $, fsusin(plon,plat,2)      ! upward solar flux at surface (input, W/m2)  
     $, fldsin(plon,plat,2)      ! downward longwave flux at surface (input, W/m2)  
     $, vwcsfcin(plon,plat,2)    ! volumetric soil water content at surface (frac)
     $, scfac(mxdynflds)         ! scale factor to convert units in input files
                                 ! to the units required by the interface.

      real(r8), dimension(:,:,:,:), allocatable ::
     $  win         ! vertical wind ("etadot", input, 1/s)
     $, omgain      ! vertical wind ("omega", input, Pa/s)
     $, cgsin       ! counter-gradient coefficient (input)  
     $, kvhin       ! vertical diffusion coefficient (input)
     $, zmuin       ! mu2 from conv_ccm, kg/m2/s            
     $, zmdin       ! md2 from conv_ccm, kg/m2/s            
     $, zeuin       ! eu2 from conv_ccm, 1/s                
     $, betain      ! overshoot parameter (input, fraction) 
     $, etain       ! convective mass flux (input, kg/m^2/s)
     $, qradin      ! radiative heating tendency (K/s)
     $, cmfdqrin    ! dq/dt due to convective rainout 
     $, conicwin    ! convective in-cloud water mixing ratio
     $, cldfrin     ! cloud fraction
     $, concldin    ! convective cloud fraction
     $, cmein       ! rate of cond-evap within the cloud
     $, prainin     ! rate of release of stratiform precip (1/s)
     $, evaprin     ! rate of evaporation of precipitation (1/s)
     $, cwatin      ! cloud water (kg/kg)

      integer, dimension(:), allocatable ::
     $  date      ! dates in input file
     $, datesec   ! seconds relative to date

      private
      public ::
     $  inidyn    ! initialize the dynamics input module
     $, chkdyn    ! check that necessary input data is available
     $, interp1   ! Interpolate the fields whose values are required at the
                  ! begining of a timestep.
     $, interp2   ! Interpolate the fields whose values are required at the
                  ! timestep midpoints.
     $, qdyn      ! Query the dynamics module.

c##############################################################################
      contains
c##############################################################################

      subroutine inidyn( xplon, xplat, xplev, xplond, xplatd, xplevp
     $,                  xmxdynflds
     $,                  lpth, rpth, icdate, icsec, offset, lun, remove
     $,                  xrdt, xrdoro, xrdphis, xcalcedot, xrdomega
     $,                  xarvdiff, xdivdiff, xarconvccm, xaddqrad
     $,                  xarcldphys
     $,                  dynnm, dynsc
     $,                  phi, latwts, hyam, hybm, hyai, hybi, phis, oro
     $)

c Initialize the dynamics input module.  The dynamics fields will be
c linearly interpolated to the initial time (icdate,icsec).
c
c NetCDF version.
c
c N.B.  Floating point data is extracted from the netCDF file utilizing
c       the interface that assumes an 8-byte local datatype.  Hence the
c       r8 parameter must be specifed to provide 8-byte reals.

      implicit none
c-----------------------------------------------------------------------
#include <netcdf.inc>
c-----------------------------------------------------------------------

c Input arguments:
c     These args are dummies here.  They are provided for interfaces
c     written in C which don't have access to the parameters in pmgrid.h.
      integer xplon, xplat, xplev, xplond, xplatd, xplevp, xmxdynflds

      character(len=*) ::
     $  lpth                ! local path for dynamics data
     $, rpth                ! remote path for dynamics data
      character(len=32) ::
     $  dynnm(mxdynflds)    ! names of dynamics fields:
      integer ::
     $  icdate     ! date of initial conditions in yyyymmdd format
     $, icsec      ! seconds relative to date for initial conditions
     $, offset     ! add offset (seconds) to dates read from headers
     $, lun        ! unit number for dynamics input (not used in netCDF version)
      logical ::
     $  remove     ! t => remove old dynamics files
     $, xrdt       ! read T

     $, xrdoro     ! read ORO from the dynamics files
     $, xrdphis    ! read PHIS
     $, xcalcedot  ! f => read ETADOT
     $, xrdomega   ! read OMEGA
     $, xarvdiff   ! read CGS and KVH for archived vertical diffusion
     $, xdivdiff   ! read variables used to diagnose vertical diffusion
c                  ! (TAUX,TAUY,SHFLX,LHFLX)
     $, xarconvccm ! read ZMU, ZMD, ZEU, CONBETA, CONETA for archived CCM
c                    convection
     $, xaddqrad   ! t => read radiative heating tendency
     $, xarcldphys ! read fields for archived cloud physics
      real(r8) ::
     $  dynsc(mxdynflds) ! scale factor to convert units in input files to
c                        ! the internal units.

c Output arguments:
      real(r8) ::
     $  phi(plat)          ! Latitude coordinate values ordered south to north.
     $, latwts(plat)       ! Latitude weights for use in integrals over surface area.
     $, hyam(plev)         ! Hybrid A coefficient at layer midpoints.
     $, hybm(plev)         ! Hybrid B coefficient at layer midpoints.
     $, hyai(plevp)        ! Hybrid A coefficient at layer interfaces.
     $, hybi(plevp)        ! Hybrid B coefficient at layer interfaces.
     $, phis(plond,plat)   ! surface geopotential  units: m^2/s^2
     $, oro(plond,plat)    ! orography, values: 0 over ocean, 1 over land,
c                                               2 over sea ice.

c Local variables:
      integer ::
     $  i, ii, j, k
     $, istat
     $, did(4), vid1
     $, n3d, n3di
      real(r8) ::
     $  tmp

      character*32 tname
      integer ndims

c Externals:
      integer :: glc, lastsl, getfile
!      real(r8):: diffdat
      logical :: dryDep
c-----------------------------------------------------------------------
c     Set variables in /dyninp/.

c     Default names for input data:
      nam = 'not-set'
      nam( 1) = 'U'       ! zonal wind component                             
      nam( 2) = 'V'       ! meridional wind component                        
      nam( 3) = 'PS'      ! surface pressure                                 
      nam( 4) = 'ORO'     ! orography                                        
      nam( 5) = 'PHIS'    ! surface geopotential                             
      nam( 6) = 'ETADOT'  ! etadot (vertical wind component, d(eta)/dt )     
      nam( 7) = 'OMEGA'   ! omega  (vertical wind component, d(pressure)/dt )
      nam( 8) = 'T'       ! temperature                                      
      nam( 9) = 'Q'       ! specific humidity                                
      nam(10) = 'TS'      ! surface temperature (K)
      nam(11) = 'TAUX'    ! surface zonal stress (N/m2)                      
      nam(12) = 'TAUY'    ! surface meridional stress (N/m2)                 
      nam(13) = 'SHFLX'   ! surface heat flux (W/m2)                         
      nam(14) = 'QFLX'    ! surface moisture flux (kg/m2/s)                  
      nam(15) = 'KVH'     ! vertical diffusion coefficient                   
      nam(16) = 'CGS'     ! counter-gradient coefficient                     
      nam(17) = 'ZMU'     ! zhang mu2 from conv_ccm, kg/m2/s                 
      nam(18) = 'ZMD'     ! zhang md2 from conv_ccm, kg/m2/s                 
      nam(19) = 'ZEU'     ! zhang eu2 from conv_ccm, 1/s                     
      nam(20) = 'CONETA'  ! hack convective mass flux (kg/m^2/s)             
      nam(21) = 'CONBETA' ! hack overshoot parameter (fraction)              
      nam(22) = 'QRAD'    ! radiative heating tendency (K/s)                 
      nam(23) = 'CMFDQR'  ! dq/dt due to convective rainout (kg/kg/s)
      nam(24) = 'CONICW'  ! convective in-cloud water mixing ratio (kg/kg)
      nam(25) = 'CLOUD'   ! total cloud fraction
      nam(26) = 'CONCLD'  ! convective cloud fraction
      nam(27) = 'CME'     ! rate of cond-evap within the cloud (kg/kg/s)
      nam(28) = 'PRAIN'   ! dq/dt due to stratiform precip (kg/kg/s)
      nam(29) = 'EVAPR'   ! rate of evaporation of strat precip (kg/kg/s)
      nam(30) = 'CWAT'    ! grid box cloud water (kg/kg)
      nam(31) = 'FSDS'    ! downward solar flux at surface (W/m2)
      nam(32) = 'SNOWH'   ! water equivalent snow depth (m)
      nam(33) = 'FSUS'    ! upward solar flux at surface (W/m2)
      nam(34) = 'FLDS'    ! downward longwave flux at surface (W/m2)
      nam(35) = 'VWC_SFC'  ! volumetric soil water content at surface (frac)

      do i = 1, mxdynflds
         if ( dynnm(i)(1:1) .ne. ' ' ) nam(i) = dynnm(i)
         scfac(i) = dynsc(i)  ! dynsc default set to 1. in inirun()
      end do

c     If local path not specified, take filename from remote path and put
c     it in subdirectory DYN
      if ( lpth(1:1) .eq. ' ' ) then
         lpath = 'DYN/' // rpth(lastsl(rpth)+1:glc(rpth))
      else
         lpath = lpth
      end if
      rpath    = rpth
      rmdyn    = remove
      dynoffset= offset
      rdt      = xrdt
      rdoro    = xrdoro
      rdphis   = xrdphis
      calcedot = xcalcedot
      rdomega  = xrdomega
      arvdiff  = xarvdiff
      divdiff  = xdivdiff
      arconvccm= xarconvccm
      arcldphys= xarcldphys
      addqrad  = xaddqrad
#ifdef DST      
      readfsds = .true.
      readfsus = .true.
      readflds = .true.
      readvwcsfc = .true.
      readsnowh = .true.
#else
      readfsds = .false.
      readfsus = .false.
      readflds = .false.
      readvwcsfc = .false.
      readsnowh = .false.
#endif
c changed default dry dep so that fsds and snowh not required
c      if ( dryDep() ) then
c         readfsds = .true.
c         readsnowh = .true.
c      end if

c     Get dynamics file.
      istat = getfile( lpath, rpath )
      if ( istat .ne. 0 ) then
         write(*,*)'inidyn: getfile returned ', istat
         call endrun
      end if

c     Open file for reading.
      call handle_ncerr( nf90_open( lpath, NF90_NOWRITE, ncid )
     $, 'inidyn: error opening file '//trim(lpath) )

c     Make list of indices into nam and vid arrays for all fields
c     that need to be read from the input file.
      namidx(1) = 1              ! U
      namidx(2) = 2              ! V
      namidx(3) = 3              ! PS
      i = 3
      if ( rdoro ) then
         namidx(i+1) = 4         ! ORO
         i = i + 1
      end if
      if ( rdphis ) then
         namidx(i+1) = 5         ! PHIS
         i = i + 1
      end if
      if ( .not. calcedot ) then
         namidx(i+1) = 6         ! ETADOT
         i = i + 1
      end if
      if ( rdomega ) then
         namidx(i+1) = 7         ! OMEGA
         i = i + 1
      end if
      if ( rdt ) then
         namidx(i+1) = 8         ! T
         i = i + 1
      end if
      namidx(i+1) = 9            ! Q
      i = i + 1
      if ( divdiff ) then
         namidx(i+1) = 10        ! TS
         namidx(i+2) = 11        ! TAUX
         namidx(i+3) = 12        ! TAUY
         namidx(i+4) = 13        ! SHFLX
         namidx(i+5) = 14        ! QFLX
         i = i + 5
      end if
      if ( readfsds ) then
         namidx(i+1) = 31        ! FSDS
         i = i + 1
      end if
      if ( readsnowh ) then
         namidx(i+1) = 32        ! SNOWH
         i = i + 1
      end if
      if ( readfsus ) then
         namidx(i+1) = 33        ! FSUS
         i = i + 1
      end if
      if ( readflds ) then
         namidx(i+1) = 34        ! FLDS
         i = i + 1
      end if
      if ( readvwcsfc ) then
         namidx(i+1) = 35        ! VWCSFC
         i = i + 1
      end if
      if ( arvdiff ) then
         namidx(i+1) = 15        ! KVH
         namidx(i+2) = 16        ! CGS
         i = i + 2
      end if
      if ( arconvccm ) then
         namidx(i+1) = 17        ! ZMU
         namidx(i+2) = 18        ! ZMD
         namidx(i+3) = 19        ! ZEU
         namidx(i+4) = 20        ! CONETA
         namidx(i+5) = 21        ! CONBETA
         i = i + 5
      end if
      if ( addqrad ) then
         namidx(i+1) = 22        ! QRAD
         i = i + 1
      end if
      if ( arcldphys ) then
         namidx(i+1) = 23        ! CMFDQR
         namidx(i+2) = 24        ! CONICW
         namidx(i+3) = 25        ! CLOUD
         namidx(i+4) = 26        ! CONCLD
         namidx(i+5) = 27        ! CME
         namidx(i+6) = 28        ! PRAIN
         namidx(i+7) = 29        ! EVAPR
         namidx(i+8) = 30        ! CWAT
         i = i + 8
      end if
      nflds = i

c     Get variable IDs for fields required to be in the input data.
      do i = 1, nflds
         ii = namidx(i)
         call handle_ncerr( nf90_inq_varid( ncid, nam(ii), vid(ii) )
     $,  'inidyn: '//trim(nam(ii))//' not found in dynamics input file')
      end do
c     Determine dimensions by looking at zonal velocity variable.
!     call handle_ncerr( nf90_inquire_variable( ncid, vid(2), name=tname, dimids=did(1:4), ndims=ndims )
      call handle_ncerr( nf_inq_vardimid( ncid, vid(1), did )
     $, 'inidyn: getting dimids for zonal velocity' )
      call handle_ncerr( nf90_inquire_dimension( ncid, did(1), len=nlon )
     $, 'inidyn: getting nlon' )
      call handle_ncerr( nf90_inquire_dimension( ncid, did(2), len=nlat )
     $, 'inidyn: getting nlat' )
      call handle_ncerr( nf90_inquire_dimension( ncid, did(3), len=nlev )
     $, 'inidyn: getting nlev' )
      call handle_ncerr( nf90_inquire_dimension( ncid, did(4), len=ntim )
     $, 'inidyn: getting ntim' )

c     Check that data resolution agrees with user specified parameters.
      if ( nlon .ne. plon .or. nlat .ne. plat .or. nlev .ne. plev ) then
         write(*,*)'inidyn: resolution of dynamics file does not match model'
         write(*,*)'in dynamics file       : nlon, nlat, nlev= '
     $        , nlon, nlat, nlev
         write(*,*)'model (set in params.h): plon, plat, plev= '
     $        , plon, plat, plev
         call endrun
      end if

c     Get date info.
      call handle_ncerr( nf90_inq_varid( ncid, 'date', datevid(1) )
     $, 'inidyn: date not found in dynamics input file' )
      call handle_ncerr( nf90_inq_varid( ncid, 'datesec', datevid(2) )
     $, 'inidyn: datesec not found in dynamics input file' )
      allocate( date(ntim) )
      allocate( datesec(ntim) )
      call handle_ncerr( nf90_get_var( ncid, datevid(1), date )
     $, 'inidyn: getting date' )
      call handle_ncerr( nf90_get_var( ncid, datevid(2), datesec )
     $, 'inidyn: getting datesec' )

c     Adjust dates by the specified offset.
      do i = 1, ntim
         call addsec2dat( dynoffset, date(i), datesec(i), date(i)
     $,                   datesec(i) )
      end do


c     Find latest date that is earlier than (icdate,icsec).
      curr_ts = lotim( icdate, icsec )
      if ( curr_ts .eq. 0 ) then
         write(*,*)'inidyn: first time sample in dynamics file after IC: '
         write(*,*)'dynamics file:',date(1), datesec(1)
         write(*,*)'requested start time:',icdate,icsec
         call endrun
      end if

c     Compute Gauss latitudes and weights.  On return; phi contains the
c     sine of the latitudes starting closest to the north pole and going
c     toward the south; w contains the corresponding Gauss weights.
      call gauaw( phi, latwts, plat )
c     Reorder and compute latitude values.  Don't need to reorder weights
c     since they are symmetric.
      do j = 1,plat/2
         tmp = asin( phi(j) )
         phi(j) = asin( phi(plat-j+1) )
         phi(plat-j+1) = tmp
      end do

c     Get vertical hybid coordinate coefficients.
      call handle_ncerr( nf90_inq_varid( ncid, 'hyam', vid1 )
     $, 'inidyn: hyam not found in dynamics input file' )
      call handle_ncerr( nf90_get_var( ncid, vid1, hyam )
     $, 'inidyn: getting hyam' )
      call handle_ncerr( nf90_inq_varid( ncid, 'hybm', vid1 )
     $, 'inidyn: hybm not found in dynamics input file' )
      call handle_ncerr( nf90_get_var( ncid, vid1, hybm )
     $, 'inidyn: getting hybm' )
      call handle_ncerr( nf90_inq_varid( ncid, 'hyai', vid1 )
     $, 'inidyn: hyai not found in dynamics input file' )
      call handle_ncerr( nf90_get_var( ncid, vid1, hyai )
     $, 'inidyn: getting hyai' )
      call handle_ncerr( nf90_inq_varid( ncid, 'hybi', vid1 )
     $, 'inidyn: hybi not found in dynamics input file' )
      call handle_ncerr( nf90_get_var( ncid, vid1, hybi )
     $, 'inidyn: getting hybi' )

c     Allocate memory for optional 3D fields.
      n3d = plon*plev*plat
      n3di = plon*plevp*plat
      if ( .not. calcedot ) then
         allocate( win(plon,plevp,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'win', 2*n3di )
      end if
      if ( rdomega ) then
         allocate( omgain(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'omgain', 2*n3d )
      end if
      if ( arvdiff ) then
         allocate( cgsin(plon,plevp,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'cgsin', 2*n3di )
         allocate( kvhin(plon,plevp,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'kvhin', 2*n3di )
      end if
      if ( arconvccm ) then
         allocate( zmuin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'zmuin', 2*n3d )
         allocate( zmdin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'zmdin', 2*n3d )
         allocate( zeuin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'zeuin', 2*n3d )
         allocate( betain(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'betain', 2*n3d )
         allocate( etain(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'etain', 2*n3d )
      end if
      if ( addqrad ) then
         allocate( qradin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'qradin', 2*n3d )
      end if
      if ( arcldphys ) then
         allocate( cmfdqrin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'cmfdqrin', 2*n3d )
         allocate( conicwin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'conicwin', 2*n3d )
         allocate( cldfrin(plon,plevp,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'cldfrin', 2*n3di )
         allocate( concldin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'concldin', 2*n3d )
         allocate( cmein(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'cmein', 2*n3d )
         allocate( prainin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'prainin', 2*n3d )
         allocate( evaprin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'evaprin', 2*n3d )
         allocate( cwatin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'cwatin', 2*n3d )
      end if

c     Read data into time level 1.
      ndyn1    = 1
      dyndate1 = date(curr_ts)
      dynsec1  = datesec(curr_ts)
      call readData( curr_ts, ndyn1, phis, oro )

c     If next time sample is not in current file then get next file.
      if ( curr_ts .eq. ntim ) call nextdyn()
      curr_ts  = curr_ts + 1

c     Check that current time sample is an upper bound to the IC time.
      if(diffdat(icdate,icsec,date(curr_ts),datesec(curr_ts)).lt.0.)then
         write(*,*)'inidyn: IC time not bounded by specified dynamics '
     $,            'input file.  Check namelist settings.'
         call endrun
      end if

c     Read data into time level 2.
      ndyn2    = 2
      dyndate2 = date(curr_ts)
      dynsec2  = datesec(curr_ts)
      call readData( curr_ts, ndyn2, phis, oro )

      end subroutine inidyn

c##############################################################################

      subroutine readData( itim, idyn, phis, oro )

c Read the field data for the specified time sample

      use netcdf
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <netcdf.inc>
c-----------------------------------------------------------------------

c Input arguments:
      integer ::
     $  itim     ! time index for requested record from netCDF file
     $, idyn     ! index in dynamics input array

c Output arguments:
      real(r8) ::
     $  phis(plond,plat)   ! surface geopotential  units: m^2/s^2
     $, oro(plond,plat)    ! orography, values: 0 over ocean, 1 over land,

c Local variables:
      integer ::
     $  m
     $, i
     $, count(4) = (/ plon, plat, plev, 1 /)
     $, counti(4) = (/ plon, plat, plevp, 1 /)
     $, start(4) = (/ 1, 1, 1, 1 /)
     $, count2(3) = (/ plon, plat, 1 /)
     $, start2(3) = (/ 1, 1, 1 /)

      real(r8) ::
     $  tmp2d(plon,plat)         ! temp 2d field
     $, tmp3d(plon,plat,plevp)   ! temp 3d field at midpoints or interfaces
     $, qmin = 2.e-12
c-----------------------------------------------------------------------

      start(4) = itim
      start2(3) = itim

c     ZONAL WIND COMPONENT, U, m/s
      m = 1
      write (6,*) ' m, nam ', m, nam(m)
      call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start, count)
     $,                  'readData: getting '//trim(nam(m)) )
      call mv_uv( 1, scfac(m), tmp3d, uvin(1,1,1,1,idyn) )

c     MERIDIONAL WIND COMPONENT, V, m/s
      m = 2
      write (6,*) ' m, nam ', m, nam(m)
      call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start, count)
     $,                  'readData: getting '//trim(nam(m)) )
      call mv_uv( 2, scfac(m), tmp3d, uvin(1,1,1,1,idyn) )

c     SURFACE PRESSURE, PS, Pa
      m = 3
      write (6,*) ' m, nam ', m, nam(m)
      call handle_ncerr( nf90_get_var(ncid, vid(m), psin(:,:,idyn), 
     $                   start2, count2 ), 'readData: getting '//trim(nam(m)) )
      if ( scfac(m) .ne. 1.0 )
     $     call scale( scfac(m), plon, nlon, nlat, psin(1,1,idyn) )

      if ( rdoro ) then
c        OROGRAPHY, ORO
         m = 4
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var(ncid, vid(m), tmp2d, start2
     $,            count2 ), 'readData: getting '//trim(nam(m)) )
         call mv_2d( scfac(m), tmp2d, oro )
      end if

      if ( rdphis ) then
c        SURFACE GEOPOTENTIAL, PHIS, m2/s2
         m = 5
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var(ncid, vid(m), tmp2d, start2
     $,            count2 ), 'readData: getting '//trim(nam(m)) )
         call mv_2d( scfac(m), tmp2d, phis )
      end if

      if ( .not. calcedot ) then
c        VERTICAL WIND COMPONENT, D(ETA)/DT, ETADOT, 1/s
         m = 6
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            counti ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plevp, scfac(m), tmp3d, win(1,1,1,idyn) )
      end if

      if ( rdomega ) then
c        VERTICAL WIND COMPONENT, D(PRESSURE)/DT, OMEGA, Pa/s
         m = 7
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, omgain(1,1,1,idyn) )
      end if

      if ( rdt ) then
c        TEMPERATURE, T, K
         m = 8
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, tin(1,1,1,idyn) )
      end if

c     SPECIFIC HUMIDITY, Q, kg/kg
      m = 9
      write (6,*) ' m, nam ', m, nam(m)
      call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
      call mv_3d( plev, scfac(m), tmp3d, qin(1,1,1,idyn) )
c     Limit Q to be >= qmin.  Some of the convection/cloud codes blow
c     up when given q=0.
      where ( qin(:,:,:,idyn) < qmin ) qin(:,:,:,idyn) = qmin

      if ( divdiff ) then
c        SURFACE TEMPERATURE, TS, K
         if ( nam(10) .ne. nam(8) ) then
            m = 10
            write (6,*) ' m, nam ', m, nam(m)
            call handle_ncerr( nf90_get_var( ncid, vid(m), tsin(:,:,idyn), 
     $           start2, count2 ),'readData: getting '//trim(nam(m)))
            if ( scfac(m) .ne. 1.0 )
     $          call scale( scfac(m), plon, nlon, nlat, tsin(1,1,idyn) )
         else
            tsin(:,:,idyn)  = tin(:,plev,:,idyn)
         end if

c        SURFACE ZONAL STRESS, TAUX, N/m2
         m = 11
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tauxin(:,:,idyn), 
     $        start2, count2 ),'readData: getting '//trim(nam(m)))
         if ( scfac(m) .ne. 1.0 )
     $        call scale( scfac(m), plon, nlon, nlat, tauxin(1,1,idyn) )

c        SURFACE MERIDIONAL STRESS, TAUY, N/m2
         m = 12
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tauyin(:,:,idyn),
     $        start2,   count2 ),'readData: getting '//trim(nam(m)))
         if ( scfac(m) .ne. 1.0 )
     $        call scale( scfac(m), plon, nlon, nlat, tauyin(1,1,idyn) )

c        SURFACE HEAT FLUX, SHFLX, W/m2
         m = 13
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), hflxin(:,:,idyn),
     $        start2,  count2 ),'readData: getting '//trim(nam(m)))
         if ( scfac(m) .ne. 1.0 )
     $        call scale( scfac(m), plon, nlon, nlat, hflxin(:,:,idyn) )

c        SURFACE MOISTURE FLUX, QFLX, kg/m2/s
         m = 14
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), qflxin(:,:,idyn),
     $        start2,   count2 ),'readData: getting '//trim(nam(m)))
         if ( scfac(m) .ne. 1.0 )
     $        call scale( scfac(m), plon, nlon, nlat, qflxin(1,1,idyn) )
      end if

      if ( readfsds ) then
c        DOWNWARD SOLAR FLUX AT SURFACE, FSDS, W/m2
         m = 31
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), fsdsin(:,:,idyn),
     $        start2,   count2 ),'readData: getting '//trim(nam(m)))
         if ( scfac(m) .ne. 1.0 )
     $        call scale( scfac(m), plon, nlon, nlat, fsdsin(1,1,idyn) )
      else
         write (6,*) ' ZEROING OUT FSDS '
         fsdsin = 0.
      end if

      if ( readsnowh ) then
c        WATER EQUIVALENT SNOW DEPTH, SNOWH, m
         m = 32
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), snowhin(:,:,idyn),
     $        start2,  count2),'readData: getting '//trim(nam(m)))
         if ( scfac(m) .ne. 1.0 )
     $        call scale( scfac(m), plon, nlon, nlat, snowhin(1,1,idyn))
      else
         write (6,*) ' ZEROING OUT SNOWH '
         snowhin = 0.
      end if

      if ( readfsus ) then
c        UPWARD SOLAR FLUX AT SURFACE, FSUS, W/m2
         m = 33
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), fsusin(:,:,idyn),
     $        start2,   count2 ),'readData: getting '//trim(nam(m)))
         if ( scfac(m) .ne. 1.0 )
     $        call scale( scfac(m), plon, nlon, nlat, fsusin(1,1,idyn) )
      else
         write (6,*) ' ZEROING OUT FSUS '
         fsusin = 0.
      end if

      if ( readflds ) then
c        DOWNWARD LONGWAVE FLUX AT SURFACE, FLDS, W/m2
         m = 34
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m),  fldsin(:,:,idyn),
     $        start2, count2 ),'readData: getting '//trim(nam(m)))
         if ( scfac(m) .ne. 1.0 )
     $        call scale( scfac(m), plon, nlon, nlat, fldsin(1,1,idyn) )
      else
         write (6,*) ' ZEROING OUT FLDSIN '
         fldsin = 0.
      end if

      if ( readvwcsfc ) then
c        VOLUMETRIC SOIL WATER CONTENT AT SURFACE, VWCSFC, (frac)
         m = 35
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), vwcsfcin(:,:,idyn), 
     $        start2,   count2),'readData: getting '//trim(nam(m)))
         if ( scfac(m) .ne. 1.0 )
     $        call scale( scfac(m), plon, nlon, nlat,vwcsfcin(1,1,idyn))
         do m = 1,nlat
            do i = 1,plon
               if (vwcsfcin(i,m,idyn).lt.-0.01) then
                  vwcsfcin(i,m,idyn) = 1.
               endif
            end do
         end do
      else
         write (6,*) ' ZEROING OUT VWCSFC '
         vwcsfcin = 0.
      end if

      if ( arvdiff ) then

c        VERTICAL DIFFUSION COEFFICIENT, KVH
         m = 15
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            counti ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plevp, scfac(m), tmp3d, kvhin(1,1,1,idyn) )

c        COUNTER-GRADIENT COEFFICIENT, CGS
         m = 16
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d,
     $        start, counti), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plevp, scfac(m), tmp3d, cgsin(1,1,1,idyn) )

c        Set cgs and kvh at top of atmosphere to 0.
c        Set kvh at bottom of atmosphere to 0.
         kvhin(:,1,:,idyn) = 0.0
         cgsin(:,1,:,idyn) = 0.0
         kvhin(:,plevp,:,idyn) = 0.0
      end if

      if ( arconvccm ) then

c        ZHANG MU2 FROM CONV_CCM, ZMU, kg/m2/s
         m = 17
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, zmuin(1,1,1,idyn) )

c        ZHANG MD2 FROM CONV_CCM, ZMD, kg/m2/s
         m = 18
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d,  start,
     $        count), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, zmdin(1,1,1,idyn) )

c        ZHANG EU2 FROM CONV_CCM, ZEU, 1/s
         m = 19
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, zeuin(1,1,1,idyn) )

c        HACK CONVECTIVE MASS FLUX, CONETA, kg/m^2/s
         m = 20
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, etain(1,1,1,idyn) )

c        HACK OVERSHOOT PARAMETER, CONBETA, fraction
         m = 21
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, betain(1,1,1,idyn) )
      end if

      if ( addqrad ) then

c        RADIATIVE HEATING TENDENCY, QRAD, K/s
         m = 22
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, qradin(1,1,1,idyn) )
      end if

      if ( arcldphys ) then

c        DQ/DT DUE TO CONVECTIVE RAINOUT, CMFDQR, kg/kg/s
         m = 23
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, cmfdqrin(1,1,1,idyn) )

c        CONVECTIVE IN-CLOUD WATER MIXING RATIO, CONICW, kg/kg
         m = 24
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, conicwin(1,1,1,idyn) )

c        TOTAL CLOUD FRACTION, CLOUD, (top 18 levels output from CCM3)
         m = 25
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, cldfrin(1,1,1,idyn) )

c        CONVECTIVE CLOUD FRACTION, CONCLD
         m = 26
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, concldin(1,1,1,idyn) )

c        RATE OF COND-EVAP WITHIN THE CLOUD, CME, kg/kg/s
         m = 27
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, cmein(1,1,1,idyn) )

c        DQ/DT DUE TO STRATIFORM PRECIP, PRAIN, kg/kg/s
         m = 28
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, prainin(1,1,1,idyn) )

c        RATE OF EVAPORATION OF STRATIFORM PRECIP, EVAPR, kg/kg/s
         m = 29
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, evaprin(1,1,1,idyn) )

c        GRID BOX CLOUD WATER, CWAT, kg/kg
         m = 30
         write (6,*) ' m, nam ', m, nam(m)
         call handle_ncerr( nf90_get_var( ncid, vid(m), tmp3d, start
     $,            count ), 'readData: getting '//trim(nam(m)) )
         call mv_3d( plev, scfac(m), tmp3d, cwatin(1,1,1,idyn) )
      end if

      write(*,*)'readData: read time sample for date = ', date(itim)
     $,            '  seconds = ', datesec(itim)

      end subroutine readData

c##############################################################################

      subroutine mv_uv( uvidx, scale, data, dynin )

c Reorder and scale input U and V.

      implicit none

c Input arguments:
      integer ::
     $  uvidx                   ! 1 => U, 2 => V

      real(r8) ::
     $  scale                   ! scale input units to m/s
     $, data(plon,plat,plev)    ! input U or V

c Output argument:
      real(r8) ::
     $  dynin(plon,plev,2,plat) ! internal storage for time samples

c Local variables:
      integer i, j, k
c-----------------------------------------------------------------------

      do j = 1, plat
         do k = 1, plev
            do i = 1, plon
               dynin(i,k,uvidx,j) = scale * data(i,j,k)
            end do
         end do
      end do
      end subroutine mv_uv

c##############################################################################

      subroutine mv_3d( nlevs, scale, data, dynin )

c Reorder and scale 3D data fields.

      implicit none

c Input arguments:
      integer ::
     $  nlevs                   ! number of levels in 3D field
      real(r8) ::
     $  scale                   ! scale input units to m/s
     $, data(plon,plat,nlevs)   ! input data

c Output argument:
      real(r8) ::
     $  dynin(plon,nlevs,plat) ! internal storage for time samples

c Local variables:
      integer i, j, k
c-----------------------------------------------------------------------

      do j = 1, plat
         do k = 1, nlevs
            do i = 1, plon
               dynin(i,k,j) = scale * data(i,j,k)
            end do
         end do
      end do
      end subroutine mv_3d

c##############################################################################

      subroutine mv_2d( scale, data, out )

c Reorder and scale 2D data fields.

      implicit none

c Input arguments:
      real(r8) ::
     $  scale             ! scale input units
     $, data(plon,plat)   ! input field

c Output argument:
      real(r8) ::
     $  out(plond,plat)

c Local variables:
      integer i, j
c-----------------------------------------------------------------------

      out(plon+1:plond,:) = 0.
      do j = 1, plat
         do i = 1, plon
            out(i,j) = scale * data(i,j)
         end do
      end do
      end subroutine mv_2d

c##############################################################################

      subroutine chkdyn( ncdate, ncsec, phis, oro )

c Check that date to interpolate the data to is contained in the interval
c of dynamics data currently in memory.
c If not, then read the next time sample from the dynamics file, acquiring
c the next file if necessary.

      implicit none

c Input arguments:
      integer ::
     $  ncdate  ! Interpolate the input data to (ncdate,ncsec)
     $, ncsec   !

c Output arguments: (only changed if new dynamics data read)
      real(r8) ::
     $  phis(plond,plat)            ! surface geopotential
     $, oro(plond,plat)             ! orography

c Local variables:
      integer ::
     $  istat
     $, lo, hi

c Externals:
!     real(r8):: diffdat
c-----------------------------------------------------------------------

      if ( diffdat( ncdate, ncsec, dyndate2, dynsec2 ) .ge. 0. ) return

c     (ncdate,ncsec) is outside of current dynamics time interval.  Get the
c     time sample(s) for the new time interval.

c     Make the upper bound of the current time interval the lower bound of the
c     new time interval.  In the "usual" case we won't need to read new
c     data into this slot.
      ndyn1 = ndyn2

c     Is the new upper bound time sample in the current dynamics file?
      lo = lotim( ncdate, ncsec )
      hi = lo + 1
      if ( hi .le. ntim ) then   ! all required data is present

c        Is the current upper bound time sample the same as the lower bound
c        time sample for the new time interval?  If not we need to acquire
c        new lower bound data.
         if ( curr_ts .eq. lo ) then
            dyndate1 = dyndate2
            dynsec1  = dynsec2
         else
            dyndate1 = date(lo)
            dynsec1  = datesec(lo)
            call readData( lo, ndyn1, phis, oro )
         end if

      else

c        Need next file.  If last time sample in current file is not the
c        current upper bound then read it into the lower bound position
c        before getting next file.
         if ( curr_ts .eq. ntim ) then
            dyndate1 = dyndate2
            dynsec1  = dynsec2
         else
            dyndate1 = date(ntim)
            dynsec1  = datesec(ntim)
            call readData( ntim, ndyn1, phis, oro )
         end if

         call nextdyn()
         lo = lotim( ncdate, ncsec )
         hi = lo + 1

c        Make sure the new upper bound time sample is in this file.
         if ( hi .gt. ntim ) then
            write(*,*)'chkdyn: current time not bounded by data on next'
     $,               ' dynamics input file.'
            call endrun
         end if

c        Is the current lower bound time sample (i.e., the last time sample
c        from the previous file) the correct one?  If not we need to acquire
c        new lower bound data.
         if ( lo .gt. 0 ) then
            dyndate1 = date(lo)
            dynsec1  = datesec(lo)
            call readData( lo, ndyn1, phis, oro )
         end if
      end if

c     Read data for new upper bound.
      ndyn2    = mod( ndyn1, 2 ) + 1
      dyndate2 = date(hi)
      dynsec2  = datesec(hi)
      call readData( hi, ndyn2, phis, oro )
      curr_ts = hi

      end subroutine chkdyn

c##############################################################################

      subroutine interp1( ncdate, ncsec, ps, t, q
     $,                   zmu, zmd, zeu, hketa, hkbeta, cmfdqr, conicw
     $,                   cldfr, concld, cme, prain, evapr, cwat
     $)

c Interpolate the fields whose values are required at the begining of
c a timestep.

      implicit none

c Input arguments:
      integer ::
     $  ncdate  ! Interpolate the input data to (ncdate,ncsec)
     $, ncsec   !

c Output arguments:
      real(r8) ::
     $  ps(plond,plat)          ! surface pressure (time interpolant)
     $, t(plond,plev,plat)      ! temperature (time interpolant)
     $, q(plond,plev,plat)      ! specific humidity (time interpolant)
     $, zmu(plond,plev,plat)    ! mu2 from conv_ccm, kg/m2/s
     $, zmd(plond,plev,plat)    ! md2 from conv_ccm, kg/m2/s
     $, zeu(plond,plev,plat)    ! eu2 from conv_ccm, 1/s
     $, hketa(plond,plev,plat)  ! Hack convective mass flux
     $, hkbeta(plond,plev,plat) ! Hack overshoot parameter
     $, cmfdqr(plond,plev,plat) ! dq/dt due to convective rainout 
     $, conicw(plond,plev,plat) ! convective in-cloud water mixing ratio (kg/kg)
     $, cldfr(plond,plevp,plat) ! cloud fraction
     $, concld(plond,plev,plat) ! convective cloud fraction
     $, cme(plond,plev,plat)    ! rate of cond-evap within the cloud
     $, prain(plond,plev,plat)  ! rate of release of stratiform precip (1/s)
     $, evapr(plond,plev,plat)  ! rate of evaporation of precipitation (1/s)
     $, cwat(plond,plev,plat)   ! cloud water (kg/kg)

c Local variables:
      integer, parameter :: one = 1
      real(r8), parameter :: zero = 0.0
      integer ::
     $  j     ! index
      real(r8) ::
     $  dynday2 ! (dyndate2,dynsec2) - (dyndate1,dynsec1) in days
     $, intpday ! (ncdate,ncsec) - (dyndate1,dynsec1) in days

c Externals:
!     real(r8):: diffdat
c-----------------------------------------------------------------------

      dynday2 = diffdat( dyndate1, dynsec1, dyndate2, dynsec2 )
      intpday = diffdat( dyndate1, dynsec1, ncdate, ncsec )

! OpenMP directive
!$OMP PARALLEL DO
      do j = 1, plat

         call intp2d( plon, plond, one, zero, dynday2, intpday
     $,               psin(1,j,ndyn1), psin(1,j,ndyn2), ps(1,j) )

         call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,               qin(1,1,j,ndyn1), qin(1,1,j,ndyn2), q(1,1,j) )

         if ( rdt ) then
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  tin(1,1,j,ndyn1), tin(1,1,j,ndyn2), t(1,1,j) )
         end if

         if ( arconvccm ) then
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  zmuin(1,1,j,ndyn1), zmuin(1,1,j,ndyn2)
     $,                  zmu(1,1,j) )
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  zmdin(1,1,j,ndyn1), zmdin(1,1,j,ndyn2)
     $,                  zmd(1,1,j) )
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  zeuin(1,1,j,ndyn1), zeuin(1,1,j,ndyn2)
     $,                  zeu(1,1,j) )
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  etain(1,1,j,ndyn1), etain(1,1,j,ndyn2)
     $,                  hketa(1,1,j) )
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  betain(1,1,j,ndyn1), betain(1,1,j,ndyn2)
     $,                  hkbeta(1,1,j) )
         end if

         if ( arcldphys ) then
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  cmfdqrin(1,1,j,ndyn1), cmfdqrin(1,1,j,ndyn2)
     $,                  cmfdqr(1,1,j) )
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  conicwin(1,1,j,ndyn1), conicwin(1,1,j,ndyn2)
     $,                  conicw(1,1,j) )
            call intp2d( plon, plond, plevp, zero, dynday2, intpday
     $,                  cldfrin(1,1,j,ndyn1), cldfrin(1,1,j,ndyn2)
     $,                  cldfr(1,1,j) )
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  concldin(1,1,j,ndyn1), concldin(1,1,j,ndyn2)
     $,                  concld(1,1,j) )
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  cmein(1,1,j,ndyn1), cmein(1,1,j,ndyn2)
     $,                  cme(1,1,j) )
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  prainin(1,1,j,ndyn1), prainin(1,1,j,ndyn2)
     $,                  prain(1,1,j) )
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  evaprin(1,1,j,ndyn1), evaprin(1,1,j,ndyn2)
     $,                  evapr(1,1,j) )
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  cwatin(1,1,j,ndyn1), cwatin(1,1,j,ndyn2)
     $,                  cwat(1,1,j) )
         end if
      end do

      end subroutine interp1

c##############################################################################

      subroutine interp2( ncdate, ncsec, uv, w, cgs, kvh
     $,                   ts, taux, tauy, hflx, qflx, fsds
c    $,                   fsus, flds, vwcsfc, snowh, qrad
     $,                   fsus, flds, vwcsfc, snowh
     $,                   omga
     $)

c Interpolate the fields whose values are required at the midpoint of
c a timestep.

      implicit none

c Input arguments:
      integer ::
     $  ncdate  ! Interpolate the input data to (ncdate,ncsec)
     $, ncsec   !

c Output arguments:
      real(r8) ::
c    $  uv(plond,plev,2,platd) ! u,v winds (time interpolant)
     $  uv(plond,plev,2,plat)  ! u,v winds (time interpolant)
     $, w(plon,plevp,plat)     ! vertical wind (time interpolant)
     $, cgs(plond,plevp,plat)  ! interpolated counter-gradient coefficient
     $, kvh(plond,plevp,plat)  ! interpolated vertical diffusion coefficient
     $, ts(plond,plat)         ! surface temperature
     $, taux(plond,plat)       ! stress in x direction
     $, tauy(plond,plat)       ! stress in y direction
     $, hflx(plond,plat)       ! surface sensible heat flux (w/m2)
     $, qflx(plond,plat)       ! surface water vapor flux (kg/m2/s)
     $, fsds(plond,plat)       ! downward solar flux at surface (W/m2)
     $, fsus(plond,plat)       ! upward solar flux at surface (W/m2)
     $, flds(plond,plat)       ! downward longwave flux at surface (W/m2)
     $, vwcsfc(plond,plat)     ! volumetric soil water content at surface (frac)
     $, snowh(plond,plat)      ! water equivalent snow depth (m)
c    $, qrad(plond,plev,plat)  ! radiative heating tendency (K/s)
     $, omga(plond,plev,plat)  ! vertical velocity (omega) (Pa/s)

c Local variables:
      integer, parameter :: one = 1
      real(r8), parameter :: zero = 0.0
      integer
     $  j     ! index
      real(r8) ::
     $  dynday2 ! (dyndate2,dynsec2) - (dyndate1,dynsec1) in days
     $, intpday ! (ncdate,ncsec) - (dyndate1,dynsec1) in days

c Externals:
!     real(r8):: diffdat
      integer system
c-----------------------------------------------------------------------

      dynday2 = diffdat( dyndate1, dynsec1, dyndate2, dynsec2 )
      intpday = diffdat( dyndate1, dynsec1, ncdate, ncsec )

! OpenMP directive
!$OMP PARALLEL DO
      do j = 1, plat
         call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,               uvin(1,1,1,j,ndyn1), uvin(1,1,1,j,ndyn2)
     $,               uv(1,1,1,j) )
         call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,               uvin(1,1,2,j,ndyn1), uvin(1,1,2,j,ndyn2)
     $,               uv(1,1,2,j) )

         if ( .not. calcedot ) then
            call intp2d( plon, plon, plevp, zero, dynday2, intpday
     $,                  win(1,1,j,ndyn1), win(1,1,j,ndyn2)
     $,                  w(1,1,j) )
         end if

         if ( rdomega ) then
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  omgain(1,1,j,ndyn1), omgain(1,1,j,ndyn2)
     $,                  omga(1,1,j) )
         end if

         if ( divdiff ) then
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                  tsin(1,j,ndyn1), tsin(1,j,ndyn2), ts(1,j) )
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                  tauxin(1,j,ndyn1),tauxin(1,j,ndyn2),taux(1,j) )
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                  tauyin(1,j,ndyn1),tauyin(1,j,ndyn2),tauy(1,j) )
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                  hflxin(1,j,ndyn1),hflxin(1,j,ndyn2),hflx(1,j))
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                  qflxin(1,j,ndyn1),qflxin(1,j,ndyn2),qflx(1,j))
         end if

         if ( readfsds ) then
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                  fsdsin(1,j,ndyn1),fsdsin(1,j,ndyn2),fsds(1,j))
         end if

         if ( readsnowh ) then
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                snowhin(1,j,ndyn1),snowhin(1,j,ndyn2),snowh(1,j))
         end if

         if ( readfsus ) then
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                  fsusin(1,j,ndyn1),fsusin(1,j,ndyn2),fsus(1,j))
         else
            fsus(:,j) = 0.
         end if

         if ( readflds ) then
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                  fldsin(1,j,ndyn1),fldsin(1,j,ndyn2),flds(1,j))
         else
            flds(:,j) = 0.
         end if

         if ( readvwcsfc ) then
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                  vwcsfcin(1,j,ndyn1),vwcsfcin(1,j,ndyn2)
     $,                  vwcsfc(1,j))
         end if

         if ( arvdiff ) then
            call intp2d( plon, plond, plevp, zero, dynday2, intpday
     $,                  cgsin(1,1,j,ndyn1), cgsin(1,1,j,ndyn2)
     $,                  cgs(1,1,j) )
            call intp2d( plon, plond, plevp, zero, dynday2, intpday
     $,                  kvhin(1,1,j,ndyn1), kvhin(1,1,j,ndyn2)
     $,                  kvh(1,1,j) )
         end if

c         if ( addqrad ) then
c            call intp2d( plon, plond, plev, zero, dynday2, intpday
c     $,                  qradin(1,1,j,ndyn1), qradin(1,1,j,ndyn2)
c     $,                  qrad(1,1,j) )
c        end if

      end do

      end subroutine interp2

c##############################################################################

      subroutine nextdyn()

c Close [and remove] the current dynamics file.  Open the next dynamics file.
c Acquire file from mass store if local version doesn't exist.

      implicit none
c-----------------------------------------------------------------------
#include <netcdf.inc>
c-----------------------------------------------------------------------

c Local variables:
      integer ::
     $  i, istat
     $, recdid, nrec
      character(len=128) ::
     $  ctmp

c Externals:
      integer glc, incstr, getfile, system
c-----------------------------------------------------------------------

c     Close current dynamics file.
      call handle_ncerr( nf90_close( ncid )
     $, 'nextdyn: error closing file '//trim(lpath) )

c     Remove if requested.
      if ( rmdyn ) then
         ctmp = 'rm -f ' // trim(lpath)
         write(*,*)'nextdyn: system issuing command - '
         write(*,*) trim( ctmp )
         istat = system( ctmp )
      endif

c     Build filenames for next dynamics file.  Assume that netCDF files
c     have a .nc extension.
      if ( rpath(1:1) .ne. ' ' ) then
         istat = incstr( rpath(:glc(rpath)-3), 1 )
         if ( istat .ne. 0 ) then
            write(*,*)'nextdyn: incstr returned ', istat
            write(*,*)'         while trying to increment ',rpath
            stop
         end if
      end if
      istat = incstr( lpath(:glc(lpath)-3), 1 )
      if ( istat .ne. 0 ) then
         write(*,*)'nextdyn: incstr returned ', istat
         write(*,*)'         while trying to increment ',lpath
         stop
      end if

c     Get dynamics file.
      istat = getfile( lpath, rpath )
      if ( istat .ne. 0 ) then
         write(*,*)'nextdyn: getfile returned ', istat
         call endrun
      end if

c     Open file for reading.
      call handle_ncerr( nf90_open( lpath, NF90_NOWRITE, ncid )
     $, 'nextdyn: error opening file '//trim(lpath) )

c     Update date info.
c      call handle_ncerr( nf90_inquire( ncid, unlimiteddimid=recdid )
c     $, 'nextdyn: inquiring record dimension ID' )
      call handle_ncerr( nf90_inq_dimid( ncid, 'time', recdid ) )
      call handle_ncerr( nf90_inquire_dimension( ncid, recdid, len=nrec )
     $, 'nextdyn: inquiring record dimension length' )
      if ( nrec .ne. ntim ) then
         deallocate( date, datesec )
         ntim = nrec
         allocate( date(ntim), datesec(ntim) )
      end if
      call handle_ncerr( NF90_INQ_VARID( ncid, 'date', datevid(1) )
     $, 'nextdyn: date not found in dynamics input file' )
      call handle_ncerr( NF90_INQ_VARID( ncid, 'datesec', datevid(2) )
     $, 'nextdyn: datesec not found in dynamics input file' )
      call handle_ncerr( NF90_GET_VAR( ncid, datevid(1), date )
     $, 'nextdyn: getting date' )
      call handle_ncerr( NF90_GET_VAR( ncid, datevid(2), datesec )
     $, 'nextdyn: getting datesec' )

c     Adjust dates by the specified offset.
      do i = 1, ntim
         call addsec2dat( dynoffset, date(i), datesec(i), date(i)
     $,                   datesec(i) )
      end do

c     Reset current time sample index.
      curr_ts = 0

      end subroutine nextdyn

c##############################################################################

      subroutine qdyn( lpth, rpth, tim_idx )

c Query the dynamics module.

      implicit none

c Output arguments:
      character(len=*) ::
     $  lpth       ! local path for dynamics data
     $, rpth       ! remote path for dynamics data
      integer ::
     $  tim_idx    ! time index for most recently read data
c-----------------------------------------------------------------------

      lpth = lpath
      rpth = rpath
      tim_idx = curr_ts
      return
      end subroutine qdyn

c##############################################################################

      subroutine intp2d( nd1, nd2, nlev, t1, t2, tint, f1, f2, fint )

c Linearly interpolate between f1(t1) and f2(t2) to fint(tint).

      implicit none

c Input arguments:
      integer ::
     $  nd1
     $, nd2
     $, nlev
      real(r8) ::
     $  t1               ! time level of f1
     $, t2               ! time level of f2
     $, tint             ! interpolant time
     $, f1(nd1,nlev)     ! field at time t1
     $, f2(nd1,nlev)     ! field at time t2

c Output arguments:
      real(r8) ::
     $  fint(nd2,nlev)  ! field at time tint

c Local variables:
      integer :: i, j
      real(r8) :: factor
c------------------------------------------------------------------------------

      factor = ( tint - t1 )/( t2 - t1)

      do j = 1, nlev
         do i = 1, nd1
            fint(i,j) = f1(i,j) + ( f2(i,j) - f1(i,j) )*factor
         end do
      end do

      end subroutine intp2d

c##############################################################################

      integer function lotim( cdate, csec )

c Return the index of the dynamics time sample that is the lower
c bound of the interval that contains the input date.  If
c (cdate,csec) is earlier than the first time sample then 0 is
c returned.  If (cdate,csec) is later than the last time sample then
c that index is returned.  If (cdate,csec) is equal to the date of a
c dynamics time sample then that index is returned.

      implicit none

c Input arguments:
      integer ::
     $  cdate    ! date in yyyymmdd
     $, csec     ! seconds relative to date

c Local variables:
      integer ::
     $  i

c Externals:
!     real(r8):: diffdat
c-----------------------------------------------------------------------

c     Find latest date that is earlier than or equal to (date,sec).
      do i = 1, ntim
         if ( diffdat(cdate, csec, date(i), datesec(i)) .gt. 0. ) then
            lotim = i-1
            exit
         end if
         if ( i .eq. ntim ) lotim = ntim
      end do

      end function lotim

c##############################################################################

      end module dyninp
