c $Id$
      module dyninp

c Module for reading meteorological driving data in CCM history format.

      implicit none
      save
      private
      public ::
     $  inidyn    ! initialize the dynamics input module
     $, chkdyn    ! check that necessary input data is available
     $, interp1   ! Interpolate the fields whose values are required at the
                  ! begining of a timestep.
     $, interp2   ! Interpolate the fields whose values are required at the
                  ! timestep midpoints.
     $, qdyn      ! Query the dynamics module.

c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

      integer, parameter ::
     $  mxfldsIn = 200    ! max number of fields in input history files

      character
     $  lpath*168         ! local path
     $, rpath*168         ! remote pathname
     $, nam(mxdynflds)*32 ! names of dynamics fields:
     $, mcflds(2,mxfldsIn)*8 ! field info from character header records

      logical
     $  cosbdyn    ! .true. => file is COS blocked
     $, rmdyn      ! .true. => remove old dynamics files
     $, rdt        ! read T
     $, rdoro      ! read ORO
     $, rdphis     ! read PHIS
     $, calcedot   ! f => read ETADOT
     $, rdomega    ! read OMEGA
     $, limqin     ! true => limit input Q s.t. it does not exceed sat mr.
     $, arvdiff    ! read CGS and KVH for archived vertical diffusion
     $, divdiff    ! read variables used to diagnose vertical diffusion
c                  ! (TAUX,TAUY,SHFLX,LHFLX)
     $, arconvccm  ! read ZMU, ZMD, ZEU, CONBETA, CONETA for archived CCM
c                    convection
     $, arcldphys  ! read fields for archived cloud physics
     $, addqrad    ! t => add radiative heating tendency to the temperature
c                         field supplied to the convection routine
     $, readfsds   ! t => read FSDS
     $, readsnowh  ! t => read SNOWH
     $, readfsus   ! t => read FSUS
     $, readflds   ! t => read FLDS
     $, readvwcsfc ! t => read VWCSFC

      integer
     $  lundyn     ! Fortran unit number attached to file
     $, lenhdi     ! length of integer header record (number of integers)
     $, lenhdc     ! length of character header record (number of character*8)
     $, lenhdr     ! length of real header record (number of reals)
     $, maxsiz     ! length of data record (number of reals)
     $, nfld       ! number of fields in data records
     $, nlon       ! number of longitudes
     $, nlonw      ! number of longitudes written to history tape
     $, nlat       ! number of latitudes
     $, nlev       ! number of model layers
     $, dyndate1   ! date previous to current date (in yyyymmdd format)
     $, dynsec1    ! seconds relative to dyndate1
     $, ndyn1      ! index in data arrays for previous time sample
     $, dyndate2   ! current date (in yyyymmdd format)
     $, dynsec2    ! seconds relative to dyndate2
     $, ndyn2      ! index in data arrays for current time sample
     $, dynoffset  ! offset (seconds) added to (date,sec) info
     $, mfilh      ! time sample number of the one most recently read.
c                  ! (maintained internally... not read from integer header)
     $, ndens      ! packing density for data records
     $, mflds(3,mxfldsIn) ! field info from integer header records

      real(r8)
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
     $, scfac(mxdynflds)         ! scale factor to convert units in input files to
                                 ! the internal units.

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

c Initialize module data, dynamics arrays, and a few other variables from
c the dynamics file headers.  The dynamics fields will be linearly
c interpolated to the initial time (icdate,icsec).

      implicit none

c Input arguments:
c     These args are dummies here.  They are provided for interfaces
c     written in C which do not have access to the parameters in pmgrid.h.
      integer xplon, xplat, xplev, xplond, xplatd, xplevp, xmxdynflds

      character
     $  lpth*(*)            ! local path for dynamics data
     $, rpth*(*)            ! remote path for dynamics data
     $, dynnm(mxdynflds)*32 ! names of dynamics fields:
      integer
     $  icdate     ! date of initial conditions in yyyymmdd format
     $, icsec      ! seconds relative to date for initial conditions
     $, offset     ! add offset (seconds) to dates read from headers
     $, lun        ! unit number for dynamics input
      logical
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
      real(r8)
     $  dynsc(mxdynflds) ! scale factor to convert units in input files to
c                        ! the internal units.

c Output arguments:
      real(r8)
     $  phi(plat)          ! Latitude coordinate values ordered south to north.
     $, latwts(plat)       ! Latitude weights for use in integrals over surface area.
     $, hyam(plev)         ! Hybrid A coefficient at layer midpoints.
     $, hybm(plev)         ! Hybrid B coefficient at layer midpoints.
     $, hyai(plevp)        ! Hybrid A coefficient at layer interfaces.
     $, hybi(plevp)        ! Hybrid B coefficient at layer interfaces.
     $, phis(plond,plat)   ! surface geopotential  units: m^2/s^2
     $, oro(plond,plat)    ! orography, values: 0 over ocean, 1 over land,
                           !                    2 over sea ice.

c Local variables:
      integer
     $  i, j, k
     $, istat, istat2 ! return status
     $, mpsig   ! position in real header of first word of vertical coord list
     $, mpflds  ! position in integer header of field info list
     $, mpcfld  ! position in character header of field info list
      integer
     $  ntrm    ! spectral truncation parameter M
     $, ntrn    ! spectral truncation parameter N
     $, ntrk    ! spectral truncation parameter K

      real(r8)
     $  tmp

c Externals:
      integer glc, lastsl, attach, rdhd
      real(DBLKIND) diffdat
      logical dryDep
c-----------------------------------------------------------------------
c     Set modulue variables.

c     Default names for input data:
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

c     Default scale factors for input data:
      do i = 1, mxdynflds
         scfac(i) = 1.0
      end do

c     If local path not specified, take filename from remote path and put
c     it in subdirectory DYN
      if ( lpth(1:1) .eq. ' ' ) then
         lpath = 'DYN/' // rpth(lastsl(rpth)+1:glc(rpth))
      else
         lpath = lpth
      end if
      rpath    = rpth
      lundyn   = lun
      rmdyn    = remove
      ndyn1    = 1
      ndyn2    = 2
      dynoffset= offset
      mfilh    = 0
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
      readfsds = .false.
      readsnowh = .false.
c changed default dry dep so that fsds and snowh not required
c      if ( dryDep() ) then
c         readfsds = .true.
c         readsnowh = .true.
c      end if
      readfsus = .false.
      readflds = .false.
      readvwcsfc = .true.

      do i = 1, mxdynflds
         if ( dynnm(i)(1:1) .ne. ' ' ) nam(i) = dynnm(i)
         scfac(i) = dynsc(i)
      end do

c     Get dynamics history tape.
      istat = attach( lpath, rpath, lundyn, .true., cosbdyn )
      if ( istat .ne. 0 ) then
         write(*,*)'inidyn: attach returned ', istat
         write(*,*)' lpath=',lpath
         write(*,*)' rpath=',rpath
         stop
      end if
      if (cosbdyn .and. RPREC < 7) then
         write (6,*) ' MATCH is not yet configured to read cray files in 4 byte float mode '
         stop
       endif

c     Read fixed part of first header record to get size and position info.
      call inirdhd( lundyn, cosbdyn, lenhdi, lenhdc, lenhdr, maxsiz
     $,             nlon, nlonw, nlat, nlev, ntrm, ntrn, ntrk, nfld
     $,             dyndate1, dynsec1, mpsig, mpflds, mpcfld, ndens )
      call addsec2dat( dynoffset, dyndate1, dynsec1, dyndate1, dynsec1 )

c     Check that data resolution agrees with user specified parameters.
      if ( nlon .ne. plon .or. nlat .ne. plat .or. nlev .ne. plev ) then
         write(*,*)'inidyn: nlon, nlat, nlev= ', nlon, nlat, nlev
         write(*,*)'inidyn: plon, plat, plev= ', plon, plat, plev
         stop
      end if

c     Make sure that initial time is not earlier than time of first time
c     sample on dynamics tape.
      if ( diffdat(  icdate, icsec, dyndate1, dynsec1 ) .gt. 0. ) then
         write(*,*)'inidyn: first time sample= ', dyndate1, dynsec1
         stop
      end if

c     Make sure fixed length arrays are large enough.
      if ( nfld .gt. mxfldsIn ) then
         write(*,*)'inidyn: MXFLDSIN must be increased to ', nfld
         stop
      end if

c     Compute Gauss latitudes and weights.  On return; phi contains the
c     sine of the latitudes starting closest to the north pole and going
c     toward the south; w contains the corresponding Gauss weights.
      call gauaw( phi, latwts, plat )
c     Reorder and compute latitude values.  Do not need to reorder weights
c     since they are symmetric.
      do j = 1,plat/2
         tmp = asin( phi(j) )
         phi(j) = asin( phi(plat-j+1) )
         phi(plat-j+1) = tmp
      end do

c     Initial reading of complete header records.
      istat = rdhd( lundyn, cosbdyn, lenhdi, lenhdc, lenhdr, nfld
     $,             dyndate1, dynsec1, mflds, mcflds
     $,             nlev, hyam, hybm, hyai, hybi )
      if ( istat .ne. 0 ) then
         write(*,*)'inidyn: rdhd returned ', istat
         stop
      end if
      call addsec2dat( dynoffset, dyndate1, dynsec1, dyndate1, dynsec1 )

c     Allocate memory for optional 3D fields.
      if ( .not. calcedot ) then
         allocate( win(plon,plevp,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'win', plon*plevp*plat*2 )
      end if
      if ( rdomega ) then
         allocate( omgain(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'omgain', plon*plev*plat*2 )
      end if
      if ( arvdiff ) then
         allocate( cgsin(plon,plevp,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'cgsin', plon*plevp*plat*2 )
         allocate( kvhin(plon,plevp,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'kvhin', plon*plevp*plat*2 )
      end if
      if ( arconvccm ) then
         allocate( zmuin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'zmuin', plon*plev*plat*2 )
         allocate( zmdin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'zmdin', plon*plev*plat*2 )
         allocate( zeuin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'zeuin', plon*plev*plat*2 )
         allocate( betain(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'betain', plon*plev*plat*2 )
         allocate( etain(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'etain', plon*plev*plat*2 )
      end if
      if ( addqrad ) then
         allocate( qradin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'qradin', plon*plev*plat*2 )
      end if
      if ( arcldphys ) then
         allocate( cmfdqrin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'cmfdqrin', plon*plev*plat*2 )
         allocate( conicwin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'conicwin', plon*plev*plat*2 )
         allocate( cldfrin(plon,plevp,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'cldfrin', plon*plevp*plat*2 )
         allocate( concldin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'concldin', plon*plev*plat*2 )
         allocate( cmein(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'cmein', plon*plev*plat*2 )
         allocate( prainin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'prainin', plon*plev*plat*2 )
         allocate( evaprin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'evaprin', plon*plev*plat*2 )
         allocate( cwatin(plon,plev,plat,2), stat=istat )
         call alloc_err( istat, 'inidyn', 'cwatin', plon*plev*plat*2 )
      end if

c     Read data into time level 1.
      call getdyn( ndyn1, phis, oro )

c     Read header records for 2nd time sample.
      istat = rdhd( lundyn, cosbdyn, lenhdi, lenhdc, lenhdr, nfld
     $,             dyndate2, dynsec2, mflds, mcflds
     $,             nlev, hyam, hybm, hyai, hybi )
      if ( istat .eq. 1 ) then
         write(*,*)'inidyn: error returned by rdhd'
         stop
      else if ( istat .eq. -1 ) then
         call nextdyn()
         istat2 = rdhd( lundyn, cosbdyn, lenhdi, lenhdc, lenhdr, nfld
     $,                 dyndate2, dynsec2, mflds, mcflds
     $,                 nlev, hyam, hybm, hyai, hybi )
         if ( istat2 .ne. 0 ) then
            write(*,*)'inidyn: error returned by rdhd for newly '
     $,               'acquired dynamics file'
            stop
         end if
      end if
      call addsec2dat( dynoffset, dyndate2, dynsec2, dyndate2, dynsec2 )

c     Read data into time level 2.
      call getdyn( ndyn2, phis, oro )

      end subroutine inidyn

c##############################################################################

      subroutine dynslic( ti, phis, oro )

c Read the data record for the current latitude and extract required fields.
c *** N.B. *** This code will break if it is built for single precision real(r8)s 
c              and the data is cray data because the variable drecin(maxsiz)
c              won't be large enough to hold the full data record which will
c              be 8*maxsize bytes.

      implicit none

c Input argument:
      integer
     $  ti                      ! time index for module arrays

c Output arguments:
      real(r8)
     $  phis(plond,plat)        ! surface geopotential
     $, oro(plond,plat)         ! orography

c Local variables:
      integer, parameter :: one = 1
      real(r8) drecin(maxsiz)  ! storage for a data record from the hist tape
      real(r8) :: qmin = 2.e-12
      integer
     $  lat        ! latitude number of data record
     $, i, k, m    ! indices
c-----------------------------------------------------------------------

      call rddat( lundyn, cosbdyn, drecin, maxsiz, ndens )

      lat = int( drecin(1) + .01 ) ! type 43 history tape

c     Extract appropriate slice and transfer to input arrays.

c     ZONAL WIND COMPONENT
      m = 1
      call mkslic( trim(nam(m)), uvin(1,1,1,lat,ti)
     $,            nlev, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
      if ( scfac(m) .ne. 1.0 )
     $ call scale( scfac(m), plon, nlon, nlev, uvin(1,1,1,lat,ti) )

c     MERIDIONAL WIND COMPONENT
      m = 2
      call mkslic( trim(nam(m)), uvin(1,1,2,lat,ti)
     $,            nlev, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
      if ( scfac(m) .ne. 1.0 )
     $ call scale( scfac(m), plon, nlon, nlev, uvin(1,1,2,lat,ti) )

c     SURFACE PRESSURE
      m = 3
      call mkslic( trim(nam(m)), psin(1,lat,ti)
     $,            one, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
      if ( scfac(m) .ne. 1.0 )
     $ call scale( scfac(m), plon, nlon, one, psin(1,lat,ti) )

      if ( rdoro ) then
c        OROGRAPHY
         m = 4
         call mkslic( trim(nam(m)), oro(1,lat)
     $,               one, mflds, mcflds, nfld
     $,               drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, oro(1,lat) )
      end if

      if ( rdphis ) then
c        SURFACE GEOPOTENTIAL
         m = 5
         call mkslic( trim(nam(m)), phis(1,lat)
     $,               one, mflds, mcflds, nfld
     $,               drecin, maxsiz, nlon, nlonw, plond, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, phis(1,lat) )
      end if

      if ( .not. calcedot ) then
c        ETADOT (VERTICAL WIND COMPONENT, D(ETA)/DT )
         m = 6
         call mkslic( trim(nam(m)), win(1,2,lat,ti)
     $,               nlev, mflds, mcflds, nfld
     $,               drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, win(1,2,lat,ti) )
c        Set etadot at top of atmosphere to 0.
         win(:,1,lat,ti) = 0.0
      end if

      if ( rdomega ) then
c        OMEGA  (VERTICAL WIND COMPONENT, D(PRESSURE)/DT )
         m = 7
         call mkslic( trim(nam(m)), omgain(1,1,lat,ti)
     $,               nlev, mflds, mcflds, nfld
     $,               drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, omgain(1,1,lat,ti) )
      end if

      if ( rdt ) then
c        TEMPERATURE
         m = 8
         call mkslic( trim(nam(m)), tin(1,1,lat,ti)
     $,               nlev, mflds, mcflds, nfld
     $,               drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, tin(1,1,lat,ti) )
      end if

c     SPECIFIC HUMIDITY
      m = 9
      call mkslic( trim(nam(m)), qin(1,1,lat,ti)
     $,            nlev, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
      if ( scfac(m) .ne. 1.0 )
     $     call scale( scfac(m), plon, nlon, nlev, qin(1,1,lat,ti) )
c     Limit Q to be >= qmin.  Some of the convection/cloud codes blow
c     up when given q=0.
      where ( qin(:,:,lat,ti) < qmin ) qin(:,:,lat,ti) = qmin

      if ( divdiff ) then

c        SURFACE TEMPERATURE (K)
         m = 10
         if ( nam(m) .ne. nam(8) ) then
            call mkslic( trim(nam(m)), tsin(1,lat,ti)
     $,                  one, mflds, mcflds, nfld
     $,                  drecin, maxsiz, nlon, nlonw, plon, nlev )
            if ( scfac(m) .ne. 1.0 )
     $       call scale( scfac(m), plon, nlon, one, tsin(1,lat,ti) )
         else
            tsin(:,lat,ti)  = tin(:,plev,lat,ti)
         end if

c        SURFACE ZONAL STRESS (N/m2)
         m = 11
         call mkslic( trim(nam(m)), tauxin(1,lat,ti)
     $,            one, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, tauxin(1,lat,ti) )

c        SURFACE MERIDIONAL STRESS (N/m2)
         m = 12
         call mkslic( trim(nam(m)), tauyin(1,lat,ti)
     $,            one, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, tauyin(1,lat,ti) )

c        SURFACE HEAT FLUX (W/m2)
         m = 13
         call mkslic( trim(nam(m)), hflxin(1,lat,ti)
     $,            one, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, hflxin(1,lat,ti) )
      end if

      if ( divdiff ) then
c        SURFACE MOISTURE FLUX (kg/m2/s)
         m = 14
         call mkslic( trim(nam(m)), qflxin(1,lat,ti)
     $,               one, mflds, mcflds, nfld
     $,               drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, qflxin(1,lat,ti) )
      end if

      if ( readfsds ) then
c        DOWNWARD SOLAR FLUX AT SURFACE (W/m2)
         m = 31
         call mkslic( trim(nam(m)), fsdsin(1,lat,ti)
     $,            one, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, fsdsin(1,lat,ti) )
      end if

      if ( readsnowh ) then
c        WATER EQUIVALENT SNOW DEPTH
         m = 32
         call mkslic( trim(nam(m)), snowhin(1,lat,ti)
     $,            one, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, snowhin(1,lat,ti) )
      end if

      if ( readfsus ) then
c        UPWARD SOLAR FLUX AT SURFACE (W/m2)
         m = 33
         call mkslic( trim(nam(m)), fsusin(1,lat,ti)
     $,            one, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, fsusin(1,lat,ti) )
      end if

      if ( readflds ) then
c        DOWNWARD LONGWAVE FLUX AT SURFACE (W/m2)
         m = 34
         call mkslic( trim(nam(m)), fldsin(1,lat,ti)
     $,            one, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, fldsin(1,lat,ti) )
      end if

      if ( readvwcsfc ) then
c        VOLUMETRIC SOIL WATER CONTENT AT SURFACE (frac)
         m = 35
         call mkslic( trim(nam(m)), vwcsfcin(1,lat,ti)
     $,            one, mflds, mcflds, nfld
     $,            drecin, maxsiz, nlon, nlonw, plon, nlev )
         do i = 1,plon
            if (vwcsfcin(i,lat,ti).lt.-9998) then
               vwcsfcin(i,lat,ti) = 1.e30
            endif
         end do
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, one, vwcsfcin(1,lat,ti) )
      end if

      if ( arvdiff ) then

c        VERTICAL DIFFUSION COEFFICIENT
         m = 15
         call mkslic( trim(nam(m)), kvhin(1,2,lat,ti)
     $,               nlev, mflds, mcflds, nfld
     $,               drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, kvhin(1,2,lat,ti) )

c        COUNTER-GRADIENT COEFFICIENT
         m = 16
         call mkslic( trim(nam(m)), cgsin(1,2,lat,ti)
     $,               nlev, mflds, mcflds, nfld
     $,               drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, cgsin(1,2,lat,ti) )

c        Set cgs and kvh at top of atmosphere to 0.
c        Set kvh at bottom of atmosphere to 0.
         kvhin(:,1,lat,ti) = 0.0
         kvhin(:,plevp,lat,ti) = 0.0
         cgsin(:,1,lat,ti) = 0.0
      end if

      if ( arconvccm ) then

c        ZHANG MU2 FROM CONV_CCM, kg/m2/s
         m = 17
         call mkslic( trim(nam(m)), zmuin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, zmuin(1,1,lat,ti) )

c        ZHANG MD2 FROM CONV_CCM, kg/m2/s
         m = 18
         call mkslic( trim(nam(m)), zmdin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, zmdin(1,1,lat,ti) )

c        ZHANG EU2 FROM CONV_CCM, 1/s
         m = 19
         call mkslic( trim(nam(m)), zeuin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, zeuin(1,1,lat,ti) )

c        HACK CONVECTIVE MASS FLUX (kg/m^2/s)
         m = 20
         call mkslic( trim(nam(m)), etain(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, etain(1,1,lat,ti) )

c        HACK OVERSHOOT PARAMETER (fraction)
         m = 21
         call mkslic( trim(nam(m)), betain(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, betain(1,1,lat,ti) )
      end if

      if ( addqrad ) then

c        RADIATIVE HEATING TENDENCY (K/s)                 
         m = 22
         call mkslic( trim(nam(m)), qradin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, qradin(1,1,lat,ti) )
      end if

      if ( arcldphys ) then
c        DQ/DT DUE TO CONVECTIVE RAINOUT (kg/kg/s)
         m = 23
         call mkslic( trim(nam(m)), cmfdqrin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, cmfdqrin(1,1,lat,ti) )

c        CONVECTIVE IN-CLOUD WATER MIXING RATIO (kg/kg)
         m = 24
         call mkslic( trim(nam(m)), conicwin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, conicwin(1,1,lat,ti) )

c        TOTAL CLOUD FRACTION (top 18 levels output from CCM3)
         m = 25
         call mkslic( trim(nam(m)), cldfrin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, cldfrin(1,1,lat,ti) )
c        Set level plevp to zero (it's not used)
         cldfrin(:,plevp,lat,ti) = 0.0

c        CONVECTIVE CLOUD FRACTION
         m = 26
         call mkslic( trim(nam(m)), concldin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, concldin(1,1,lat,ti) )

c        RATE OF COND-EVAP WITHIN THE CLOUD
         m = 27
         call mkslic( trim(nam(m)), cmein(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, cmein(1,1,lat,ti) )

c        DQ/DT DUE TO STRATIFORM PRECIP
         m = 28
         call mkslic( trim(nam(m)), prainin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, prainin(1,1,lat,ti) )

c        RATE OF EVAPORATION OF STRATIFORM PRECIP
         m = 29
         call mkslic( trim(nam(m)), evaprin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, evaprin(1,1,lat,ti) )

c        GRID BOX CLOUD WATER
         m = 30
         call mkslic( trim(nam(m)), cwatin(1,1,lat,ti)
     $,               nlev, mflds, mcflds
     $,               nfld, drecin, maxsiz, nlon, nlonw, plon, nlev )
         if ( scfac(m) .ne. 1.0 )
     $    call scale( scfac(m), plon, nlon, nlev, cwatin(1,1,lat,ti) )
      end if

      end subroutine dynslic

c##############################################################################

      subroutine chkdyn( ncdate, ncsec, phis, oro )

c Check that date to interpolate the data to is contained in the interval
c of dynamics data currently in memory.
c If not, then read the next time sample from the dynamics file, acquiring
c the next file if necessary.

      implicit none

c Input arguments:
      integer
     $  ncdate  ! Interpolate the input data to (ncdate,ncsec)
     $, ncsec   !

c Output arguments: (only changed if new dynamics data read)
      real(r8)
     $  phis(plond,plat)            ! surface geopotential
     $, oro(plond,plat)             ! orography

c Local variables:
      integer
     $  istat
     $, istat2
      logical
     $  gotnext  ! true => next dynamics file has been acquired
      real(r8)
     $  hyam(plev),hybm(plev),hyai(plevp),hybi(plevp)

c Externals:
      integer rdhd
      real(DBLKIND) diffdat
c-----------------------------------------------------------------------

      gotnext = .false.

   10 continue

c     If (ncdate,ncsec) is outside of dynamics time interval, then get the
c     next time sample.

      if ( diffdat( dyndate2, dynsec2, ncdate, ncsec ) .gt. 0. ) then

c        Read header of next time sample.  If EOF encountered, acquire
c        next dynamics file.

         dyndate1 = dyndate2
         dynsec1 = dynsec2
         ndyn1 = ndyn2
         ndyn2 = mod( ndyn1, 2 ) + 1
         istat = rdhd( lundyn, cosbdyn, lenhdi, lenhdc, lenhdr, nfld
     $,                dyndate2, dynsec2, mflds, mcflds
     $,                nlev, hyam, hybm, hyai, hybi )
         if ( istat .eq. 1 ) then
            write(*,*)'chkdyn: error returned by rdhd'
            stop
         else if ( istat .eq. -1 .and. .not.gotnext ) then
            call nextdyn()
            gotnext = .true.
            istat2 = rdhd( lundyn, cosbdyn, lenhdi, lenhdc, lenhdr, nfld
     $,                    dyndate2, dynsec2, mflds, mcflds
     $,                    nlev, hyam, hybm, hyai, hybi )
            if ( istat2 .ne. 0 ) then
               write(*,*)'chkdyn: error returned by rdhd for newly '
     $,                  'acquired dynamics file'
               stop
            end if
         else if ( istat .eq. -1 .and. gotnext ) then
            write(*,*)'chkdyn: required time sample not found on next '
     $,               'dynamics file'
            stop
         end if
         call addsec2dat( dynoffset, dyndate2, dynsec2, dyndate2
     $,                   dynsec2 )

         write(*,*)'chkdyn: read time sample for date = ', dyndate2
     $,            '  seconds = ', dynsec2

c        Get "new" dynamics fields.
         call getdyn( ndyn2, phis, oro )

c        Jump back to make sure that (ncdate,ncsec) is contained in the
c        dynamics time interval.
         go to 10

      endif

      end subroutine chkdyn

c##############################################################################

      subroutine getdyn( ti, phis, oro )

c Multitasked loop over call to get the required dynamics fields from
c a single data record.  Assume the history tape is positioned past the
c header records and in front of the first data record when getdyn is
c called.

      implicit none

c Input argument:
      integer
     $  ti                   ! time index for module arrays

c Output arguments:
      real(r8)
     $  phis(plond,plat)     ! surface geopotential
     $, oro(plond,plat)      ! orography

c Local variables:
      integer j
c-----------------------------------------------------------------------

      do j = 1, plat
         call dynslic( ti, phis, oro )
      end do

c     Update time sample counter.
      mfilh = mfilh + 1

      end subroutine getdyn

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
      real(DBLKIND) :: diffdat
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
     $,                   fsus, flds, vwcsfc, snowh, qrad
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
     $  uv(plond,plev,2,platd) ! u,v winds (time interpolant)
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
     $, qrad(plond,plev,plat)  ! radiative heating tendency (K/s)
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
      real(DBLKIND) :: diffdat
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
         end if

         if ( readflds ) then
            call intp2d( plon, plond, one, zero, dynday2, intpday
     $,                  fldsin(1,j,ndyn1),fldsin(1,j,ndyn2),flds(1,j))
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

         if ( addqrad ) then
            call intp2d( plon, plond, plev, zero, dynday2, intpday
     $,                  qradin(1,1,j,ndyn1), qradin(1,1,j,ndyn2)
     $,                  qrad(1,1,j) )
         end if

      end do

      end subroutine interp2

c##############################################################################

      subroutine nextdyn()

c Close [and remove] the current dynamics file.  Open the next dynamics file.
c Acquire file from mass store if local version doesn't exist.

      implicit none

c Local variables:
      integer
     $  istat  ! return status
     $, i      ! index
      character
     $  ctmp*128    ! temp string to build commands in

c Externals:
      integer glc, fsystem, incstr, attach

c Library functions:
      integer
     $  crayclose  ! close a COS blocked file on a non-Cray machine
c-----------------------------------------------------------------------

c     Close current dynamics file.

#ifdef CRAY
      close( lundyn )
#else
      if ( cosbdyn ) then
         istat = crayclose( lundyn )
      else
         close( lundyn )
      end if
#endif

c     Remove if requested.
      if ( rmdyn ) then
         do i = 1, len(ctmp)
            ctmp(i:i) = ' '
         end do
         ctmp = 'rm -f ' // lpath(:glc(lpath))
         write(*,*)'nextdyn: fsystem issuing command - '
         write(*,*) ctmp(:glc(ctmp))
         istat = fsystem( ctmp )
      endif

c     Build filenames for next dynamics file.
      if ( rpath(1:1) .ne. ' ' ) then
         istat = incstr( rpath, 1 )
         if ( istat .ne. 0 ) then
            write(*,*)'nextdyn: incstr returned ', istat
            write(*,*)'         while trying to increment ',rpath
            stop
         end if
      end if
      istat = incstr( lpath, 1 )
      if ( istat .ne. 0 ) then
         write(*,*)'nextdyn: incstr returned ', istat
         write(*,*)'         while trying to increment ',lpath
         stop
      end if

c     Attach next dynamics file.
      istat = attach( lpath, rpath, lundyn, .true., cosbdyn )
      if ( istat .ne. 0 ) then
         write(*,*)'nextdyn: attach returned ', istat
         stop
      end if

c     Reset time sample counter.

      mfilh = 0

      end subroutine nextdyn

c##############################################################################

      subroutine qdyn( lpth, rpth, ntim )

c Query the dynamics module.

      implicit none

c Output arguments:
      character
     $  lpth*(*)   ! local path for dynamics data
     $, rpth*(*)   ! remote path for dynamics data
      integer
     $  ntim        ! time sample number
c-----------------------------------------------------------------------

      lpth = lpath
      rpth = rpath
      ntim = mfilh

      end subroutine qdyn

c##############################################################################

      subroutine intp2d( nd1, nd2, nlev, t1, t2, tint, f1, f2, fint )

c Linearly interpolate between f1(t1) and f2(t2) to fint(tint).

      implicit none

c Input arguments:
      integer
     $  nd1
     $, nd2
     $, nlev
      real(r8)
     $  t1               ! time level of f1
     $, t2               ! time level of f2
     $, tint             ! interpolant time
     $, f1(nd1,nlev)     ! field at time t1
     $, f2(nd1,nlev)     ! field at time t2

c Output arguments:
      real(r8)
     $  fint(nd2,nlev)  ! field at time tint

c Local variables:
      integer i, j
      real(r8) factor
c------------------------------------------------------------------------------

      factor = ( tint - t1 )/( t2 - t1)

      do j = 1, nlev
         do i = 1, nd1
            fint(i,j) = f1(i,j) + ( f2(i,j) - f1(i,j) )*factor
         end do
      end do

      end subroutine intp2d

c##############################################################################

      end module dyninp
