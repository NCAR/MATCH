c $Id$
#ifdef DST
c++csz
#include <dst.h> /* Dust preprocessor tokens */
c--csz
#endif
      module histout

c Module for writing output fields in CCM history format.

      use precision
      use rec_io

      implicit none
      save
      private
      public ::
     $  inihist   ! initialize the history output module
     $, outfld    ! update history for fields written to output file
     $, qhstc     ! query histout module for character data
     $, qhsti     ! query histout module for integer data
     $, wrhd      ! write CCM header records to output file
     $, rhstrst   ! read histout restart data
     $, whstrst   ! write histout restart data
     $, disphst   ! copy current history output file to remote storage
     $, nexthist  ! prepare to write the next history output file
     $, appndhst  ! open history file and position it at the end

c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

      integer, parameter ::
     $  mxdiagflds=500     ! max number of fields in table containing info for the
                           ! diagnostic fields (set in deffld() )
     $, lenhdr=3*(2*PLEV+1)+2*PLAT   ! length of real header record (number of reals)
     $, mpflds=38          ! start index in integer header record for field info
     $, mpcfld=90          ! start index in character header record for field info

      character
     $  lpath*168    ! local path for history files
     $, rpath*168    ! remote path for history files
      character
     $  wpasswd*8        ! MS write password
     $, lnhstc*80        ! mss pathname of current history file
     $, lnhstf*80        ! mss pathname of first history file for this case
     $, mcstit*80        ! case title
     $, flddat(3,mxdiagflds)*8 ! name,levels,units for diagnostic fields
     $, diagnams(mxdiagflds)*8 ! names for diagnostic fields

      logical
     $  async      ! true => output hist files transferred to ms asynchronously
     $, rmout      ! true => remove hist file after linking and/or transferring
c                  !         to mass store

      integer
     $  nfld        ! total number of user requested output fields
     $, ninstfld    ! number of user requested instantaneous output fields
     $, ntimavfld   ! number of user requested time averaged output fields
     $, ndiagfld    ! num of fields in table for diagnostic fields
     $, fld_idx     ! current index into flddat array
     $, lenhdi      ! length of integer header record (number of integers)
     $, lenhdc      ! length of character header record (number of character*8)
      integer
     $  lunhst        ! Fortran unit number attached to file
     $, maxsiz        ! length of output data records (may be packed)
     $, unpksiz       ! length of unpacked output data record
     $, hbuflen       ! length of history buffer (one latitude slice)
     $, nestep        ! ending timestep for model run
     $, nhtfrq        ! freq of hist tape writes in timesteps
     $, ndens         ! write density for output history tapes
     $, mfilt         ! max number of time samples written to hist tape
     $, ifcnt         ! counter for fields as they are initialized in inifld.
c                     ! used to make sure number of fields specified in fldlst
c                     ! is equal to the parameter nfld set above.
     $, tafcnt        ! counter for time averaged fields
     $, irt           ! retention period in days for mass store files.
     $, ndcurf        ! day of first time sample on current history file
     $, nscurf        ! seconds relative to ndcurf
     $, ncdatf        ! date of first time sample on current history file
     $, ncsecf        ! seconds relative to ncdatf

      integer, dimension(:), allocatable ::
     $  hdi         ! integer header record
     $, hbufp       ! offsets into a history buffer latitude slice for each field

      character(len=8), dimension(:), allocatable ::
     $  hdc         ! character header record

      real(r8) ::  hdr(lenhdr) ! real header record

      real(r8) :: fnorm       ! normalizer for time averaged fields

      real(r8), dimension(:,:), allocatable ::
     $  hbuf3d      ! history buffer


c##############################################################################
      contains
c##############################################################################

      subroutine inihist( lpth, rpth, rrootd, case, iunhst
     $,                   ntshst, outinst, outtimav, xnestep, xnhtfrq
     $,                   xndens, xmfilt, stday, stsec, icdate
     $,                   icsec, dtime, rstflg, title, pi
     $,                   latrad, gw, hyam, hybm, hyai
     $,                   hybi, xirt, xasync, xrmout, xwpasswd
     $)

c Initialize private data.  Open first history file.  Initialize history buffer.
c For a restart run, the private data that are updated as the run progresses,
c and the history buffer, will be initialized by calling rhstrst() from the restart
c module.

      use calendar

      implicit none

c-----------------------------------------------------------------------
#include <tracnam.h>
c-----------------------------------------------------------------------


c Input arguments:
      character
     $  lpth*(*)      ! local path for history file
     $, rpth*(*)      ! remote path for history file
     $, rrootd*(*)    ! Default remote root directory
     $, case*(*)      ! case identifier
     $, outinst(mxoutflds)*8  ! instantaneous output fields names
     $, outtimav(mxoutflds)*8 ! time averaged output fields names
     $, title*(*)     ! case title
     $, xwpasswd*(*)  ! MS write password
      integer
     $  iunhst     ! Fortran unit for output history file
     $, ntshst     ! number of time samples in the starting hist file
     $, xnestep    ! ending timestep for model run
     $, xnhtfrq    ! freq of hist tape writes in timesteps
     $, xndens     ! write density for output history tapes
     $, xmfilt     ! max number of time samples written to hist tape
     $, stday      ! starting day of initial run
     $, stsec      ! seconds relative to stday
     $, icdate     ! date of initial conditions in yyyymmdd format
     $, icsec      ! seconds relative to icdate
     $, dtime      ! timestep interval in seconds
     $, rstflg     ! 0 for initial run, 1 for "restart"
     $, xirt       ! retention period in days for mass store files.  If set to
                   ! 0 (default) then output files not transfered to mass store
      logical
     $  xasync     ! true => output hist files transferred to ms asynchronously
     $, xrmout     ! true => remove hist file after transferring to mass store
      real(r8) ::
     $  pi                ! 3.14159...
     $, latrad(plat)      ! latitudes in radians
     $, gw(plat)          ! Gauss weights
     $, hyam(plev)        ! Hybrid A coefficient at layer midpoints.
     $, hybm(plev)        ! Hybrid B coefficient at layer midpoints.
     $, hyai(plevp)       ! Hybrid A coefficient at layer interfaces.
     $, hybi(plevp)       ! Hybrid B coefficient at layer interfaces.

c Local Variables:
      integer
     $  ioff, aoff, boff  ! offsets
     $, i, j, k           ! indicies
     $, junk
     $, istat             ! return status
     $, ihd1(27)   ! for reading part of the integer header record
     $, ntrm       ! spectral truncation parameter M
     $, ntrn       ! spectral truncation parameter N
     $, ntrk       ! spectral truncation parameter K
      logical
     $  cosbhst   ! .true. => history file is COS blocked.  This info is
c                 ! not used since this file is only written in the format
c                 ! of the local machine.
      real(r8) ::
     $  hya(2*plev+1)  ! hybrid A coefficients in real header order
     $, hyb(2*plev+1)  ! hybrid B coefficients in real header order

c Externals:
      integer attach, glc, lastsl
c-----------------------------------------------------------------------

c     Define fields available for output to history file.
      fld_idx = 0
      flddat = ' '
      call deffld()


c     Initialize /histout/ variables directly from input arguments.
      if ( lpth(1:1) .eq. ' ' ) then
         lpath = 'h0001'
      else
         lpath = lpth
      end if
      if ( rpth(1:1) .eq. ' ' ) then
         rpath = '/' // rrootd(:glc(rrootd)) // '/match/' //
     $           case(:glc(case)) // '/hist/' //
     $           lpath(lastsl(lpath)+1:glc(lpath))
      else
         rpath = rpth
      end if
      lunhst   = iunhst
      nestep   = xnestep
      nhtfrq   = xnhtfrq
      ndens    = xndens
      mfilt    = xmfilt
      ifcnt    = 0
      tafcnt   = 0
      hbuflen  = 0
      irt      = xirt
      async    = xasync
      rmout    = xrmout
      wpasswd  = xwpasswd

c     Initialize field specific values in header records.

c     Count diagnostic fields and fill array that holds names for use in
c     the findfld routine.
      ndiagfld = 0
      do i = 1, mxdiagflds
         if ( flddat(1,i) .ne. '        ' ) then
            ndiagfld = ndiagfld + 1
            diagnams(i) = flddat(1,i)
         end if
      end do


c++drb   User can request ALL time average or instanteous fields 
c     If first field in outinst or outtimav is 'ALL' then modify 
      if (( outinst(1).eq.'ALL').or.(outtimav(1).eq.'ALL')) then
         if (ndiagfld+2*pcnst.gt.mxoutflds) then
            write(6,*)'requested ALL diag fields (currently ',ndiagfld,')'
     $           ,'plus tracnam and sflxnam ',2*pcnst
     $           ,'exceeds mxoutflds ',mxoutflds,' defined in pmgrid.h'
     $           ,'calling abort'
            call abort()
         end if
      end if

c     list to include all available variables for instantaneous
      if ( outinst(1).eq.'ALL') then    
         do i = 1,ndiagfld              ! Add all diag fields
            outinst(i) = diagnams(i)
         end do
                                        ! Add the tracer fields (named)
         j = ndiagfld+1                 ! j = next index in outinst to write to     
         do i = 1,pcnst-1               ! CWAT is last pcnst, is in diagflds
            if (tracnam(i).ne.'DEFAULT') then
               outinst(j) = tracnam(i)
               j = j+1
            end if
         end do
                                        ! Now add the surface flux fields
         do i = 1,pcnst
            if (sflxnam(i).ne.'DEFAULT') then
               outinst(j) = sflxnam(i)
               j = j + 1
            end if
         end do
         write(6,*)'OUTINST = ALL gives field list:',outinst
      end if

c     list to include all available variables for time average
      if ( outtimav(1).eq.'ALL') then
         do i = 1,ndiagfld              ! Add all diag fields
            outtimav(i) = diagnams(i)
         end do
                                        ! Add the tracer fields (named)
         j = ndiagfld+1                 ! j = next index in outtimav to write to
         do i = 1,pcnst-1   ! CWAT is last pcnst, is in diagflds
            if (tracnam(i).ne.'DEFAULT') then
               outtimav(j) = tracnam(i)
               j = j+1
            end if
         end do
         do i = 1,pcnst                 ! Now add the surface flux fields
            if (sflxnam(i).ne.'DEFAULT') then
               outtimav(j) = sflxnam(i)
               j = j+1
            end if
         end do
         write(6,*)'OUTTIMAV = ALL gives field list:',outtimav
      end if
c--drb


c     Check that the requested fields are known, and that there are no
c     duplicates.  Count the requested fields.
      call fldparse( outinst, outtimav )

      nfld = ninstfld + ntimavfld  ! number of fields written to output file
      if ( nfld .eq. 0 ) then
         write(*,*)'inihist: no output fields requested.'
         write(*,*)'         Use the namelist variables outinst'
     $,            ' and/or outtimav.'
         stop
      end if

c     Allocate memory for header records.
      lenhdi = 37 + 3*nfld
      lenhdc = 89 + 2*nfld
      allocate( hdi(lenhdi), stat=istat )
      call alloc_err( istat, 'inihist', 'hdi', lenhdi )
      allocate( hdc(lenhdc), stat=istat )
      call alloc_err( istat, 'inihist', 'hdc', lenhdc )
      allocate( hbufp(nfld), stat=istat )
      call alloc_err( istat, 'inihist', 'hbufp', nfld )
      hbufp = 0   ! initialize pointer to zero

c     Add field into to integer and char header records.  Compute offsets into
c     data records and history buffer latitude slices.
      call fldlst( ninstfld, outinst, 0 )
      call fldlst( ntimavfld, outtimav, 1 )

c     Integer header.
      hdi(1)  = lenhdi    ! length of integer header
      hdi(2)  = 43        ! code for ccm2 history tape
      hdi(3)  = ntshst    ! counter for time samples in this file
      hdi(4)  = mfilt     ! max number of time samples in a file
      hdi(5)  = 3         ! number of header records
      hdi(6)  = maxsiz    ! length of data record
      hdi(7)  = unpksiz   ! length of unpacked data record
      hdi(8)  = 0         ! unused
      hdi(9)  = plon      ! number of longitudes
      hdi(10) = plon      ! number of longitude values written to data record
      hdi(11) = plat      ! number of latitudes
      hdi(12) = plev      ! number of vertical levels

c     Spectral truncation parameters are set based on horizontal resolution.
      call getspecres( plon, plat, ntrm, ntrn, ntrk )
      hdi(13) = ntrm      ! spectral truncation parameter M
      hdi(14) = ntrn      ! spectral truncation parameter N
      hdi(15) = ntrk      ! spectral truncation parameter K
      hdi(16) = nfld      ! number of fields in each time sample
      hdi(18) = 0         ! timestep number at start of run
      hdi(20) = 0         ! base day 
      hdi(21) = 0         ! seconds relative to base day
      if ( stday .eq. 0 .and. stsec .eq. 0 ) then
         hdi(24) = icdate       ! base date
         hdi(25) = icsec        ! seconds relative to base date
      else
         call addsec2dat( -(stday*86400+stsec), icdate, icsec, hdi(24)
     $,                   hdi(25) )
      end if
      hdi(28) = dtime     ! model timestep in seconds
      hdi(29) = nhtfrq    ! frequency of time sample writes in timesteps
      hdi(30) = rstflg    ! flag for type of run
      hdi(31) = lenhdc    ! length of character header
      hdi(32) = lenhdr    ! length of real header
      hdi(33) = 1         ! position in real header of vert coord list
      hdi(34) = hdi(33)+6*plev+3 ! position in real header of lats (degrees)
      hdi(35) = hdi(34)+plat ! position in real header of gauss weights
      hdi(36) = mpflds    ! position in int header of mflds
      hdi(37) = mpcfld    ! position in char header of mcflds
      if ( rstflg .eq. 0 ) then
         hdi(17) = 0      ! current iteration number
         hdi(19) = 0      ! timesteps since last time sample written
         hdi(22) = stday  ! current day corresponding to current timestep
         hdi(23) = stsec  ! seconds relative to current day
         hdi(26) = icdate ! current date corresponding to current timestep
         hdi(27) = icsec  ! current seconds relative to current date
      end if

c     Character header.
      do i = 1, 89
         hdc(i) = '        '
      end do

      hdc(1) = case    ! case identifier
      mcstit = title   ! case title
      lnhstc = rpath(:glc(rpath))
      do i = 1, 10
         hdc(1+i) = mcstit(8*(i-1)+1:8*i)
         hdc(11+i) = lnhstc(8*(i-1)+1:8*i)
      end do
      if ( rstflg .eq. 0 ) then
         lnhstf = lnhstc
         do i = 1, 10
            hdc(24+i) = lnhstf(8*(i-1)+1:8*i)
         end do
      end if

c     Real header.
c     Shuffle the hybrid coeffs.
      do k = 1, plev
         hya(2*k-1) = hyai(k)
         hyb(2*k-1) = hybi(k)
         hya(2*k)   = hyam(k)
         hyb(2*k)   = hybm(k)
      end do
      hya(2*plev+1) = hyai(plevp)
      hyb(2*plev+1) = hybi(plevp)

      ioff = hdi(33) - 1
      aoff = 2*plev + 1
      boff = 4*plev + 2
      do i = 1, 2*plev+1
         hdr(ioff     +i) = hya(i) + hyb(i)
         hdr(ioff+aoff+i) = hya(i)
         hdr(ioff+boff+i) = hyb(i)
      end do

      ioff = hdi(34) - 1
      do i = 1, plat
         hdr(ioff+i) = latrad(i) * 180./pi   ! latitudes (S to N, degrees)
      end do

      ioff = hdi(35) - 1
      do i = 1, plat
         hdr(ioff+i) = gw(i)                 ! Gauss weights
      end do

      if ( ntshst .eq. 0 ) then

c        Open new history file.
         open( lunhst, file=lpath(:glc(lpath)), form='unformatted'
     $,        status='unknown' )

      else

c        Acquire most recent history file for this run.
         istat = attach( lpath, rpath, lunhst, .true., cosbhst )
         if ( istat .ne. 0 ) then
            write(*,*)'inihst: attach returned ', istat
            stop
         end if

c        Position the history file.  Read part of integer header for the 
c        first time sample to get the day/date info.  Skip the rest.

         read( lunhst, end=100 ) ihd1
         ndcurf = ihd1(22)
         nscurf = ihd1(23)
         ncdatf = ihd1(26)
         ncsecf = ihd1(27)
         do j = 1, 2 + plat           ! read 2 header and plat data records
            read( lunhst, end=100 ) junk
         end do

         if ( ntshst .gt. 1 ) then
            do i = 1, ntshst-1
               do j = 1, 3 + plat     ! read 3 header and plat data records
                  read( lunhst, end=100 ) junk
               end do
            end do
         end if

      end if

c     Initialize history buffer for time averaged fields.
      if ( hbuflen .ne. 0 ) then
c        Allocate memory for history buffer and initialize.
         allocate( hbuf3d(hbuflen,plat), stat=istat )
         call alloc_err( istat, 'inihist', 'hbuf3d', hbuflen*plat )
         if ( rstflg .eq. 0 ) then
            hbuf3d = 0.0
         end if
      end if

      return

c     Error exit

  100 continue
      write(*,*)'inihist: EOF encountered trying to position the'
      write(*,*)'         history file after ', ntshst,' time samples'
      stop
      
      end subroutine inihist

c##############################################################################

      subroutine addfld( nm, lvl, units )

c Add field information to master field list.

      implicit none

c Input arguments.
      character(len=*), intent(in) ::
     $  nm            ! Field name
     $, lvl           ! Field level keyword
     $, units         ! Field dimensional units
c----------------------------------------------------------------------

      fld_idx = fld_idx + 1
      if ( fld_idx .gt. mxdiagflds ) then
         write(*,*)'addfld: mxdiagflds must be set larger than '
     $,            mxdiagflds, ' in histout.F'
         stop
      end if
      flddat(1,fld_idx) = nm
      flddat(2,fld_idx) = lvl
      flddat(3,fld_idx) = units

      end subroutine addfld

c##############################################################################

      subroutine deffld()

c Define fields that may be written to the output history file.

#ifdef SCYC
c++scyc
      use scyc, only: doscyc, add_volc_emis
c--scyc
#endif
#ifdef CAER
c++caer
      use caer, only: docaer, caer_ncomp, caer_idx1
c--caer
#endif

      implicit none

c++csz
#ifdef DST
#include <dstgrd.h> /* Dust grid sizes */
#include <dstnm.h> /* Nomenclature for outfld() */
c     Local
      integer fld_nbr_non_dst   ! Number of standard (non-dust) fields
      integer m                 ! Counting index
#endif /* not DST */
c--csz
c++caer
c----------------------------------------------------------------------
#include <tracnam.h>
c----------------------------------------------------------------------
      integer mm, nphob
c--caer
c----------------------------------------------------------------------

      call addfld( 'PS',       'single',  'PA' )          ! surface pressure
      call addfld( 'U',        'multmid', 'M/S' )         ! zonal wind
      call addfld( 'V',        'multmid', 'M/S' )         ! meridional wind
      call addfld( 'ETADOT',   'multint', '1/S' )         ! vertical velocity in eta coord
      call addfld( 'ETDOTR',   'multint', '1/S' )         ! vertical velocity in eta coord before correction
      call addfld( 'OMEGA',    'multmid', 'pa/s' )        ! vertical velocity in pressure coord
      call addfld( 'MCA  ',    'multmid', 'pa/s' )        ! Mass convergence above
      call addfld( 'MCAR ',    'multmid', 'pa/s' )        ! Mass convergence above before correction
      call addfld( 'T',        'multmid', 'K' )           ! temperature
      call addfld( 'TS',       'single',  'K' )           ! surface temperature
      call addfld( 'Q',        'multmid', 'KG/KG' )       ! advected specific humidity
      call addfld( 'RELHUM',   'multmid', 'fraction' )    ! relative humidity
      call addfld( 'QIN',      'multmid', 'kg/kg' )       ! specific humidity from input data
      call addfld( 'TCWV',     'single',  'kg/m2' )       ! total column water vapor
      call addfld( 'O3',       'multmid',  'kg/kg' )      ! 3D ozone
      call addfld( 'Z',        'multmid', 'm' )           ! geopotential height above surface at layer midpoints
      call addfld( 'ZI',       'multint', 'm' )           ! geopotential height above surface at interfaces
      call addfld( 'ORO',      'single',  'FLAG' )        ! Orography
      call addfld( 'PHIS',     'single',  'M2/S2' )       ! surface geopotential
      call addfld( 'TAUX',     'single',  'N/m2' )        ! surface zonal stress (N/m2)
      call addfld( 'TAUY',     'single',  'N/m2' )        ! surface meridional stress (N/m2)
      call addfld( 'SHFLX',    'single',  'W/m2' )        ! surface sensible heat flux (W/m2)
      call addfld( 'QFLX',     'single',  'kg/m2/s' )     ! surface moisture flux (kg/m2/s)
      call addfld( 'FSDS',     'single',  'W/m2' )        ! downward solar flux at surface (W/m2)
      call addfld( 'FSUS',     'single',  'W/m2' )        ! upward solar flux at surface (W/m2)
      call addfld( 'FLDS',     'single',  'W/m2' )        ! downward longwave flux at surface (W/m2)
      call addfld( 'VWCSFC',   'single',  'fraction' )    ! volumetric soil water content at surface (frac)
      call addfld( 'SNOWH',    'single',  'm' )           ! water equivalent accumulated snow depth (m)
      call addfld( 'KVH',      'multint', 'M2/S' )        ! vertical diffusion coeff for heat
      call addfld( 'CGS',      'multint', 'S/M2' )        ! vertical diffusion counter gradient coeff
      call addfld( 'PBLH',     'single',  'M' )           ! PBL height
      call addfld( 'TPERT',    'single',  'K' )           ! tpert from pbl
      call addfld( 'QPERT',    'single',  'KG/KG' )       ! qpert from pbl
      call addfld( 'ZMU',      'multmid', 'KG/M2/S' )     ! upward mass flux (Zhang)
      call addfld( 'ZMD',      'multmid', 'KG/M2/S' )     ! downward mass flux (Zhang)
      call addfld( 'ZEU',      'multmid', '1/S' )         ! entrainment in updraft (Zhang)
      call addfld( 'ZDU',      'multmid', '1/S' )         ! detrainment in updraft (Zhang)
      call addfld( 'CONBETA',  'multmid', 'FRACTION' )    ! overshoot parameter (Hack)
      call addfld( 'CONETA',   'multmid', 'KG/M2/S' )     ! convective mass flux (Hack)
      call addfld( 'CMFDQR',   'multmid', 'KG/KG/S' )     ! convective rainout
      call addfld( 'CMFMC',    'multmid', 'KG/KG/S' )     ! moist convection cloud mass flux
      call addfld( 'CONICW',   'multmid', 'KG/KG' )       ! convective in-cloud water m.r.
      call addfld( 'PRECC',    'single',  'M/S' )         ! convective precip
      call addfld( 'PRECT',    'single',  'm/s' )         ! total precip
      call addfld( 'ZCMFDQR',  'multmid', 'KG/KG/S' )     ! convective rainout (Zhang)
      call addfld( 'ZCMFMC',   'multmid', 'KG/KG/S' )     ! moist convection cloud mass flux (Zhang)
      call addfld( 'ZPRECC',   'single',  'M/S' )         ! convective precip (Zhang)
      call addfld( 'CLOUD',    'multmid', 'FRAC' )        ! Total Cloud Fraction
      call addfld( 'CONCLD',   'multmid', 'FRAC' )        ! Cloud Fraction (convective)
      call addfld( 'CLDTOT',   'single',  'FRAC' )        ! Vertically-integrated  random overlap, total cloud amount
      call addfld( 'CLDLOW',   'single',  'FRAC' )        ! Vertically-integrated  random overlap, low cloud amount
      call addfld( 'CLDMED',   'single',  'FRAC' )        ! Vertically-integrated  random overlap, mid-level cloud amount
      call addfld( 'CLDHGH',   'single',  'FRAC' )        ! Vertically-integrated  random overlap, high cloud amount
      call addfld( 'CME',      'multmid', 'KG/KG/S' )     ! rate of cond-evap within cloud
      call addfld( 'EVAPR',    'multmid', 'KG/KG/S' )     ! rate of evap of falling precip
      call addfld( 'PRAIN',    'multmid', 'KG/KG/S' )     ! rate of conversion of condensate to precipitation
      call addfld( 'PRECL',    'single',  'M/S' )         ! stratiform precip
      call addfld( 'CWAT',     'multmid', 'KG/KG' )       ! Grid box cloud Water
      call addfld( 'XTHETA',   'multmid', 'K' )           ! advected THETA
      call addfld( 'DQR',      'multmid', 'kg/kg/s' )     ! total Q tendency due to rain
      call addfld( 'TAQ',      'multmid', 'kg/kg/s' )     ! total advection tendency of Q
      call addfld( 'TATHETA',  'multmid', 'kg/kg/s' )     ! total advection tendency of potential temperature
      call addfld( 'VDQ',      'multmid', 'kg/kg/s' )     ! vertical diffusion tendency of Q
      call addfld( 'DCQ',      'multmid', 'kg/kg/s' )     ! moist convection tendency of Q
!      call addfld( 'PDMFUEN',  'multmid', 'KG/M2/S' )     ! entrainment rate for updraft (Tiedtke)
!      call addfld( 'PDMFUDE',  'multmid', 'KG/M2/S' )     ! detrainment rate for updraft (Tiedtke)
!      call addfld( 'PDMFDEN',  'multmid', 'KG/M2/S' )     ! entrainment rate for downdraft (Tiedtke)
!      call addfld( 'PDMFDDE',  'multmid', 'KG/M2/S' )     ! detrainment rate for downdraft (Tiedtke)
!      call addfld( 'PMFU',     'multmid', 'KG/M2/S' )     ! updraft mass flux (Tiedtke)
!      call addfld( 'PMFD',     'multmid', 'KG/M2/S' )     ! downdraft mass flux (Tiedtke)
!      call addfld( 'PCBOT',    'single',  'PA' )          ! pressure at bottom of cloud
!      call addfld( 'PCTOP',    'single',  'PA' )          ! pressure at top of cloud
!      call addfld( 'DEEPC',    'single',  'FRACTION' )    ! Deep Convection Diagnosed?  1= Deep, 0= not deep
!      call addfld( 'SHALC',    'single',  'FRACTION' )    ! Shallow Convection Diagnosed? 1=Shallow, 0=not shallow
!      call addfld( 'MIDC',     'single',  'FRACTION' )    ! Mid-level Convection Diagnosed?  1= Mid, 0= not mid
!      call addfld( 'CFREQ',    'single',  'FRACTION' )    ! Convection at all diagnosed (either of above 3)?

      call addfld( 'U10',      'single',  'm/s' )         ! 10 m. wind speed
C      call addfld( 'TAUSSLT',  'single',  'Opt.Dpth' )    ! opt depth for sea salt
      call addfld( 'SSLTOD',  'single',  'Opt.Dpth' )    ! opt depth for sea salt
      call addfld( 'SSLT',    'multmid', 'kg/kg' )        ! sea salt mixing ratio
      call addfld( 'AEROD',    'single',  ' ' )
      call addfld( 'O3INT',       'single',  'Dobson units' )   ! integrated ozone in Dobson  units

#ifdef SCYC
c++scyc
      if ( doscyc() ) then
         if (add_volc_emis()) call addfld( 'SO2VOLC',  'multmid', 'kg/kg/s' )
!        call addfld( 'SO2DVEL',  'single',  'm/s' )
         call addfld( 'SO2DRY',   'single',  'kg/m2/s' )
         call addfld( 'SO2SRCG',  'multmid', 'kg/kg/s' )
         call addfld( 'SO2SRCG2', 'multmid', 'kg/kg/s' )
         call addfld( 'SO2WET',   'multmid', 'kg/kg/s' )
!        call addfld( 'SO4DVEL',  'single',  'm/s' )
         call addfld( 'SO4DRY',   'single',  'kg/m2/s' )
         call addfld( 'SO4SRC',   'multmid', 'kg/kg/s' ) ! total production of SO4
         call addfld( 'SO4SRCG',  'multmid', 'kg/kg/s' ) ! SO4 source due to gas phase chemistry
         call addfld( 'SO4WET',   'multmid', 'kg/kg/s' )
!        call addfld( 'H2O2DVEL', 'single',  'm/s' )
         call addfld( 'H2O2DRY',  'single',  'kg/m2/s' )
         call addfld( 'H2O2SRC',  'multmid', 'kg/kg/s' )
         call addfld( 'H2O2SNKA', 'multmid', 'kg/kg/s' )
         call addfld( 'H2O2SNKG', 'multmid', 'kg/kg/s' )
         call addfld( 'H2O2WET',  'multmid', 'kg/kg/s' )
         call addfld( 'PH',       'multmid', ' ' ) ! pH of cloud water
         call addfld( 'DMSSNK',   'multmid', 'kg/kg/s' )
         call addfld( 'SO4OD',    'single',  ' ' )
#ifdef ASSIM
         call addfld( 'SO4ASM',   'multmid', 'kg/kg/s' ) !assimilation tendency
                  ! VLCOD optical depth due to Pinatubo interpolated from volcaod.nc file
         call addfld( 'VLCOD',  'single', 'optical depth' )
                  ! VLCODASM optical depth due to Pinatubo removed from assimilation
         call addfld( 'VLCODASM',  'single', 'optical depth' )         
         call addfld( 'DBAODOUT',  'single', 'optical depth' )         
         call addfld( 'DBAODIN',  'single', 'optical depth' )         
         call addfld( 'DBAODDIF',  'single', 'optical depth' )         
         call addfld( 'ASMRATIO',  'multmid', 'ratio' )    
         call addfld( 'ASMAOD',  'single', 'optical depth' )         

#endif  
      end if
c--scyc
#endif
c++csz
#ifdef DST
      fld_nbr_non_dst=fld_idx   ! Number of standard (non-dust) fields
c     NB: Fields MUST be padded to 8 characters when passed to addfld()
c     Not doing so will cause MATCH to fail, at least on SGI
c     Currently there are 16 dust fields (independent of dst_nbr)
c     Thus we add 16 to mxdiagflds in histout.com
c     NB: pmgrid.h:mxoutflds must also be large enough
      call addfld('DSTMPC  ','single','KG/M2   ')
      call addfld('DSTODXC ','single','FRACTION')
      call addfld('DSTQ    ','multmid','KG/KG   ')
      call addfld('DSTSFDRY','single','KG/M2/S ')
      call addfld('DSTSFGRV','single','KG/M2/S ')
      call addfld('DSTSFMBL','single','KG/M2/S ')
      call addfld('DSTSFPCP','single','KG/M2/S ')
      call addfld('DSTSFTRB','single','KG/M2/S ')
      call addfld('DSTSSEVP','multmid','KG/KGS  ')
      call addfld('DSTSSDRY','multmid','KG/KGS  ')
      call addfld('DSTSSPCP','multmid','KG/KGS  ')
      call addfld('FRC_WET ','single','FRC     ')
      call addfld('TPT_GND ','single','K       ')
      call addfld('VWC_SFC ','single','M3/M3   ')
      call addfld('WND_FRC ','single','M/S     ')
      call addfld('WND_FRCS','single','M/S     ')
      call addfld('WND_FRCT','single','M/S     ')
      call addfld('WND_RFR ','single','M/S     ')
      call addfld('WND_RFRT','single','M/S     ')
!     call addfld('PRECT2  ','single','M/S     ')
c     Following field names were initialized in dst_nm_cmn_ini()
c     One field of each type is created for each dust size bin
c     Thus, if you add/remove fields here, remember to adjust mxdiagflds in histout.com
      do m=1,dst_nbr
         call addfld(odxc_dst_nm(m),'single','FRACTION')
         call addfld(mpc_dst_nm(m),'single','KG/M2   ')
         call addfld(flx_mss_dry_sfc_nm(m),'single','KG/M2/S ')
         call addfld(flx_mss_grv_sfc_nm(m),'single','KG/M2/S ')
         call addfld(flx_mss_mbl_sfc_nm(m),'single','KG/M2/S ')
         call addfld(flx_mss_pcp_sfc_nm(m),'single','KG/M2/S ')
         call addfld(flx_mss_trb_sfc_nm(m),'single','KG/M2/S ')
#ifdef ASSIM
         call addfld( 'DST01ASM',  'multmid', 'kg/kg/s' )
         call addfld( 'DST02ASM',  'multmid', 'kg/kg/s' )
         call addfld( 'DST03ASM',  'multmid', 'kg/kg/s' )
         call addfld( 'DST04ASM',  'multmid', 'kg/kg/s' )
#endif
      end do
      write (6,'(a,i3,a,i3)') 'dst: mxoutflds = ',mxoutflds,', mxdiagflds = ',mxdiagflds
      write (6,'(a,i3,a,i3,a)') 'dst: deffld() defines ',fld_idx,' fields, ',fld_idx-fld_nbr_non_dst,' of which are dust-related'
#endif /* not DST */
c--csz
#ifdef CAER
c++caer
      if ( docaer() ) then
         nphob = caer_ncomp()/2
         do m = 1, nphob
            mm = caer_idx1() + m - 1 ! hydrophobic species
            call addfld( trim(tracnam(mm))//'SRC', 'multmid', 'kg/kg/s' )
            call addfld( trim(tracnam(mm))//'DRY', 'single',  'kg/m2/s' )
            call addfld( trim(tracnam(mm))//'OD',  'single',  ' ' )
            call addfld( trim(tracnam(mm))//'ASM',  'multmid', 'kg/kg/s' )
         end do
         do m = 1, nphob
            mm = caer_idx1() + nphob + m - 1 ! hydrophilic species
            call addfld( trim(tracnam(mm))//'DRY', 'single',  'kg/m2/s' )
            call addfld( trim(tracnam(mm))//'WET', 'multmid', 'kg/kg/s' )
            call addfld( trim(tracnam(mm))//'OD',  'single',  ' ' )
            call addfld( trim(tracnam(mm))//'ASM',  'multmid', 'kg/kg/s' )
         end do
      end if
c--caer
#endif
      end subroutine deffld

c##############################################################################

      subroutine fldlst( mxflds, outnam, timerep )

c For each requested output field find out what type of field it is, then add
c appropriate info to the integer and character header records and compute
c the offsets into the output data record and history buffer latitude slice.

      implicit none
c----------------------------------------------------------------------
#include <tracnam.h>
c----------------------------------------------------------------------

c Input arguments:
      integer
     $  mxflds
     $, timerep    ! 0 => instant, 1 => averaged
      character
     $  outnam(mxflds)*8

c Local variables:
      integer
     $  single  ! flag for single level field
     $, multint ! flag for multi level field at layer interfaces
     $, multmid ! flag for multi level field at layer midpoints
      parameter(
     $  single=0, multint=1, multmid=2 )
      integer
     $  i       ! index
     $, indx
     $, levrep
c----------------------------------------------------------------------

      do 10 i = 1, mxflds

         call lookup( outnam(i), pcnst, tracnam, indx )
         if ( indx .gt. 0 ) then
            call inifld( outnam(i), multmid , timerep, 'kg/kg' )
            go to 10
         end if

         call lookup( outnam(i), pcnst, sflxnam, indx )
         if ( indx .gt. 0 ) then
            call inifld( outnam(i), single, timerep, 'kg/m2/s' )
            go to 10
         end if

         call lookup( outnam(i), pnats, natsnam, indx )
         if ( indx .gt. 0 ) then
            call inifld( outnam(i), multmid , timerep, 'kg/kg' )
            go to 10
         end if

         call lookup( outnam(i), pcnst, tanam, indx )
         if ( indx .gt. 0 ) then
            call inifld( outnam(i), multmid , timerep, 'kg/kg/s' )
            go to 10
         end if

         call lookup( outnam(i), pcnst, vdnam, indx )
         if ( indx .gt. 0 ) then
            call inifld( outnam(i), multmid , timerep, 'kg/kg/s' )
            go to 10
         end if

         call lookup( outnam(i), pcnst, dcnam, indx )
         if ( indx .gt. 0 ) then
            call inifld( outnam(i), multmid , timerep, 'kg/kg/s' )
            go to 10
         end if

         call lookup( outnam(i), pcnst+2, xflxnam, indx )
         if ( indx .gt. 0 ) then
            call inifld( outnam(i), multmid , timerep, 'kg/m2/s' )
            go to 10
         end if

         call lookup( outnam(i), pcnst+2, yflxnam, indx )
         if ( indx .gt. 0 ) then
            call inifld( outnam(i), multmid , timerep, 'kg/m2/s' )
            go to 10
         end if

         call lookup( outnam(i), pcnst+2, zflxnam, indx )
         if ( indx .gt. 0 ) then
            call inifld( outnam(i), multint , timerep, 'kg/m2/s' )
            go to 10
         end if

         call lookup( outnam(i), ndiagfld, diagnams, indx )
         if ( indx .gt. 0 ) then
            if ( flddat(2,indx) .eq. 'single' ) then
               levrep = single
            else if ( flddat(2,indx) .eq. 'multint' ) then
               levrep = multint
            else if ( flddat(2,indx) .eq. 'multmid' ) then
               levrep = multmid
            end if
            call inifld( outnam(i), levrep, timerep, flddat(3,indx) )
            go to 10
         end if

   10 continue

      end subroutine fldlst

c##############################################################################

      subroutine fldparse( outinst, outtimav )

c Parse the output field lists to determine the number of output fields.
c Check lists for inconsistencies, such as duplicated requests or fields
c that are not contained in the tracnam or diagnosed field lists.

      implicit none
c----------------------------------------------------------------------
#include <tracnam.h>
c----------------------------------------------------------------------

c Input/Output arguments:
      character
     $  outinst(mxoutflds)*8  ! instantaneous output fields
     $, outtimav(mxoutflds)*8 ! time averaged output fields
      
c Local variables:
      integer
     $  i,j
      logical
     $  found
c-----------------------------------------------------------------------

      ninstfld = 0
      do 10 i = 1, mxoutflds
         if ( outinst(i) .eq. '        ' ) go to 20

         call findfld( mxoutflds, outinst, i, pcnst, tracnam, ninstfld
     $,                found )
         if ( found ) go to 10
         call findfld( mxoutflds, outinst, i, pcnst, sflxnam, ninstfld
     $,                found )
         if ( found ) go to 10
         call findfld( mxoutflds, outinst, i, pnats, natsnam, ninstfld
     $,                found )
         if ( found ) go to 10
         call findfld( mxoutflds, outinst, i, pcnst, tanam, ninstfld
     $,                found )
         if ( found ) go to 10
         call findfld( mxoutflds, outinst, i, pcnst, vdnam, ninstfld
     $,                found )
         if ( found ) go to 10
         call findfld( mxoutflds, outinst, i, pcnst, dcnam, ninstfld
     $,                found )
         if ( found ) go to 10
#ifdef CALC_ADV_FLUX
         call findfld( mxoutflds, outinst, i, pcnst+2, xflxnam
     $,                ninstfld, found )
         if ( found ) go to 10
         call findfld( mxoutflds, outinst, i, pcnst+2, yflxnam
     $,                ninstfld, found )
         if ( found ) go to 10
         call findfld( mxoutflds, outinst, i, pcnst+2, zflxnam
     $,                ninstfld, found )
         if ( found ) go to 10
#endif
         call findfld( mxoutflds, outinst, i, ndiagfld, diagnams
     $,                ninstfld, found )
         if ( found ) go to 10

         write(*,*)'fldparse: requested instantaneous field ',outinst(i)
     $,            ' not found... stopping.'
         stop
   10 continue
   20 continue

      ntimavfld = 0
      do 30 i = 1, mxoutflds
         if ( outtimav(i) .eq. '        ' ) go to 40

         call findfld( mxoutflds, outtimav, i, pcnst, tracnam, ntimavfld
     $,                found )
         if ( found ) go to 30
         call findfld( mxoutflds, outtimav, i, pcnst, sflxnam, ntimavfld
     $,                found )
         if ( found ) go to 30
         call findfld( mxoutflds, outtimav, i, pnats, natsnam, ntimavfld
     $,                found )
         if ( found ) go to 30
         call findfld( mxoutflds, outtimav, i, pcnst, tanam, ntimavfld
     $,                found )
         if ( found ) go to 30
         call findfld( mxoutflds, outtimav, i, pcnst, vdnam, ntimavfld
     $,                found )
         if ( found ) go to 30
         call findfld( mxoutflds, outtimav, i, pcnst, dcnam, ntimavfld
     $,                found )
         if ( found ) go to 30
#ifdef CALC_ADV_FLUX
         call findfld( mxoutflds, outtimav, i, pcnst+2, xflxnam
     $,                ntimavfld, found )
         if ( found ) go to 30
         call findfld( mxoutflds, outtimav, i, pcnst+2, yflxnam
     $,                ntimavfld, found )
         if ( found ) go to 30
         call findfld( mxoutflds, outtimav, i, pcnst+2, zflxnam
     $,                ntimavfld, found )
         if ( found ) go to 30
#endif
         call findfld( mxoutflds, outtimav, i, ndiagfld, diagnams
     $,                ntimavfld, found )
         if ( found ) go to 30

         write(*,*)'fldparse: requested time average field ',outtimav(i)
     $,            ' not found... stopping.'
         stop
   30 continue
   40 continue

c     Check for duplication between time average and instantaneous lists.

      if ( ninstfld .gt. 0 .and. ntimavfld .gt. 0 ) then
         do i = 1, ntimavfld
            do j = 1, ninstfld
               if ( outtimav(i) .eq. outinst(j) ) then
                  write(*,*)'fldparse: ',outtimav(i),' requested as '
     $,              'both instantaneous and time averaged... stopping.'
                  stop
               end if
            end do
         end do
      end if

      end subroutine fldparse

c##############################################################################

      subroutine findfld( lenreq, reqlst, ireq, lenlook, looklst, nfld
     $,                   found )

c The field at the current position in the request list (reqlst(ireq))
c is being searched for in looklst.  If it is found, then search through
c the previously updated part of the reqlst to make sure the field hasn't
c already been specified.  If it hasn't, then update the reqlst with the 
c new field.  Updating the reqlst is done to overwrite any duplicate
c specifications.  That way when the list is used to construct the history
c tape header we don't need to go through this again.

      implicit none

c Input arguments:
      integer
     $  lenreq
     $, ireq     ! current position in reqlist
     $, lenlook
      character
     $  looklst(lenlook)*8  ! a list of fields that the model can output

c Input/Output arguments:
      character
     $  reqlst(lenreq)*8   ! user requested fields
      integer
     $  nfld                ! current number of valid fields in reqlist

c Output arguments:
      logical
     $  found   ! t => valid user requested field
c               ! f => requested field not found in input looklst

c Local variables:
      integer
     $  j, k
c-----------------------------------------------------------------------

      found = .false.

      if ( lenlook .gt. 0 ) then
         do j = 1, lenlook
            if ( reqlst(ireq) .eq. looklst(j) ) then
               found = .true.
c              Make sure field hasn't already been requested.
               if ( nfld .gt. 0 ) then
                  do k = 1, nfld
                     if ( reqlst(ireq) .eq. reqlst(k) ) go to 10
                  end do
                  nfld = nfld + 1
c                 Update list to remove duplicate field requests
                  reqlst(nfld) = reqlst(ireq)
               else
                  nfld = 1
               end if
               go to 10
            end if
         end do
   10    continue
      end if

      
      end subroutine findfld

c##############################################################################

      subroutine lookup( name, lenlook, looklst, indx )

c Look for the input name in the array looklst.  Return array index.

      implicit none

c Input arguments:
      integer
     $  lenlook
      character
     $  name*8
     $, looklst(lenlook)*8

c Output arguments:
      integer
     $  indx   ! looklst(indx) .eq. name,  return 0 if name not found.

c Local variables:
      integer
     $  j
c-----------------------------------------------------------------------

      indx = 0

      if ( lenlook .gt. 0 ) then
         do j = 1, lenlook
            if ( name .eq. looklst(j) ) then
               indx = j
               return
            end if
         end do
      end if

      end subroutine lookup

c##############################################################################

      subroutine fldinfo( name, found, mflds, hbufst )

c Get field info from integer header record and start index into history buffer
c needed to write history tape and update counters.

      implicit none

c Input arguments:
      character
     $  name*(*)          ! requested field

c Output arguments:
      logical
     $  found     ! t => requested field found
      integer
     $  mflds(3)  ! integer header info for requested field
     $, hbufst    ! start index for requested field in history buffer slice

c Local variables:
      integer
     $  j
c-----------------------------------------------------------------------

      found = .false.

      do j = 1, nfld
         if ( hdc(mpcfld+2*(j-1)) .eq. name ) then
            found = .true.
            mflds(1) = hdi(mpflds  +3*(j-1))
            mflds(2) = hdi(mpflds+1+3*(j-1))
            mflds(3) = hdi(mpflds+2+3*(j-1))
            hbufst = hbufp(j)
            return
         end if
      end do

      end subroutine fldinfo

c##############################################################################

      subroutine appndhst( )

c Re-open history file and position it at the end.  We do the positioning
c ourselves because opening a file in "APPEND" mode is not part of f77 std.

      implicit none

c Local variables:
      integer
     $  i, j, junk

c Externals:
      integer glc
c-----------------------------------------------------------------------

      open( lunhst, file=lpath(:glc(lpath)), status='OLD'
     $,     form='UNFORMATTED' )

c     Position the history file.
      do i = 1, hdi(3)
         do j = 1, 3 + plat     ! read 3 header and plat data records
            read( lunhst, end=100 ) junk
         end do
      end do

      return

c     Error exit
  100 continue
      write(*,*)'appndhst: EOF encountered trying to position the'
      write(*,*)'          history file after ',hdi(3),' time samples'
      stop

      end subroutine appndhst

c##############################################################################

      subroutine disphst( laststep, noremove )

c Copy the history file to remote storage device.  Remove the local copy
c if requested.

      implicit none

c Input arguments:
      logical
     $  laststep     ! true => current timestep is final timestep for run
     $, noremove     ! true => don't remove history file after disposing.
c                    ! (used to override the user specified value)

c Local variables:
      logical
     $  remove      ! true => remove history file after disposing
      character
     $  com*80      ! comment field for mss

c Externals:
      integer glc
c-----------------------------------------------------------------------

      remove = rmout
      if ( laststep .or. noremove ) remove = .false.
      call ddicom( ndcurf, nscurf, hdi(22), hdi(23), ncdatf, ncsecf
     $,           hdi(26), hdi(27), com )
      call dispose( lunhst, lpath, rpath, irt
     $,             async, com, wpasswd, remove )

      end subroutine disphst

c##############################################################################

      subroutine inifld( fname, levflg, accflg, units )

c Initialize field information in the integer and character header records.
c Initialize field offsets for both the output data records and the internal
c history buffer.

      implicit none

c Input arguments:
      character
     $  fname*(*)   !  name of field
     $, units*(*)   !  units of field
      integer
     $  levflg      ! number of levels flag
c                   ! 0 = single level, 1 = multi at layer interfaces
c                   ! 2 = multi at layer midpoints
     $, accflg      ! accumulation flag
c                   ! 0 = instantaneous, 1 = averaged

c Local variables:
      integer
     $  ioff    ! offset into integer header record
     $, coff    ! offset into character header record
     $, ilen    ! field length
     $, iplen   ! packed field length
     $, ilenpr  ! field length of previous field
     $, iplenpr ! packed field length of previous field
     $, hbufprev! history buffer offset of previous time averaged field
      save ilenpr, iplenpr, hbufprev
c-----------------------------------------------------------------------

      ifcnt = ifcnt + 1

      if ( ifcnt .gt. nfld ) then
         write(*,*)'inifld: parameter nfld not set correctly in '
     $,            'inihist'
         stop
      end if

      if ( accflg .eq. 1 ) tafcnt = tafcnt + 1

      ioff = mpflds - 1 + 3*(ifcnt-1)
      coff = mpcfld - 1 + 2*(ifcnt-1)
      hdc(coff+1) = fname
      hdc(coff+2) = units
      hdi(ioff+1) = accflg*10 + levflg
      hdi(ioff+3) = ndens

c     Calculate lengths of field (unpacked and packed)
      ilen = plon
      iplen = plon
      if ( ndens .gt. 1 ) iplen = 3 + (plon-1)/ndens
      if ( levflg .ge. 1 ) then
         ilen = ilen*plev
         iplen = iplen*plev
      end if

c     Calc start indices and update length for output buffer.
      if ( ifcnt .eq. 1 ) then
         hdi(ioff+2) = 3
         maxsiz = 2 + iplen
         unpksiz = 2 + ilen
      else
         hdi(ioff+2) = hdi(ioff-1) + iplenpr  ! hdi(ioff-1) is start index for previous field
         maxsiz = maxsiz + iplen
         unpksiz = unpksiz + ilen
      end if
c     Save packed length of field to set offset for next field.
      iplenpr = iplen

c     Calc start indices and update length for history buffer.
      if ( accflg .eq. 1 ) then
         if ( tafcnt .eq. 1 ) then
            hbufp(ifcnt) = 1
            hbuflen = ilen
         else
            hbufp(ifcnt) = hbufprev + ilenpr
            hbuflen = hbuflen + ilen
         end if
c        Save unpacked length of field and offset to set offset for next field.
         ilenpr = ilen
         hbufprev = hbufp(ifcnt)
      end if


      end subroutine inifld

c##############################################################################

      subroutine nexthist( endofrun )

c Determine pathnames for the next history file. Update CCM header records.
c Open next history file if we are not at the end of the run.

      implicit none

c Input argument:
      logical
     $  endofrun    ! true => last timestep or shutdown requested

c Local variables:
      integer
     $  istat
     $, i

c Externals:
      integer glc, incstr
c-----------------------------------------------------------------------

c     Build filenames for next history file.
      istat = incstr( rpath, 1 )
      if ( istat .ne. 0 ) then
         write(*,*)'nexthist: incstr returned ', istat
         write(*,*)'         while trying to increment ',rpath
         stop
      end if
      istat = incstr( lpath, 1 )
      if ( istat .ne. 0 ) then
         write(*,*)'nexthist: incstr returned ', istat
         write(*,*)'         while trying to increment ',lpath
         stop
      end if

c     Reset counter for number of time samples in a file.
      hdi(3) = 0

c     Update remote filename of current history file.
      lnhstc = rpath
      do i = 1, 10
         hdc(11+i) = lnhstc(8*(i-1)+1:8*i)
      end do

c     Open next history file.
      if ( .not. endofrun ) then
         open( lunhst, file=lpath(:glc(lpath)), form='unformatted'
     $,        status='unknown' )
      end if

      end subroutine nexthist

c##############################################################################

      subroutine outfld( fname, field, idim, lat, obuf )

c For instantaneous fields, write to output buffer if a time sample is to
c be written this timestep.  Otherwise return.
c For time averaged fields, update the history buffer unless it is time
c to write a time sample.  In that case, do the averaging and write to
c the output buffer.

      implicit none

c Input arguments:
      character
     $  fname*(*)     ! name of output field
      integer
     $  idim          ! first dimension of field array
     $, lat           ! latitude index (S to N)
      real(r8)
     $  field(idim,*) ! field data

c Output argument:
      real(r8)
     $  obuf(*)       ! output buffer

c Local variables:
      integer
     $  i,j,k    ! indices
     $, ip       ! pointer into hbuf3d(*,lat)
     $, nlev     ! number of levels in field
     $, mflds(3) ! integer header info for requested field
      logical
     $  found    ! t => requested field found

c Externals:
      logical
     $  outputTS ! true => time sample to be written this timestep
c-----------------------------------------------------------------------

      call fldinfo( fname, found, mflds, ip )

      if ( found ) then

         if ( mflds(1)/10 .eq. 0 ) then              ! instantaneous

            if ( outputTS() ) then                      ! write output buffer

               call wrobuf( mflds, idim, field, obuf )

            else
               return
            end if

         else                                        ! time average

            nlev = 1
            if ( mod(mflds(1),10) .gt. 0 ) then
               nlev = plev
            end if


            if ( outputTS() ) then                      ! write output buffer

               j = 0
               do k = 1, nlev
                  do i = 1, plon
                     hbuf3d(ip+j,lat) = (hbuf3d(ip+j,lat) + field(i,k))
     $                    *fnorm
                     j = j + 1
                  end do
               end do

               call wrobuf( mflds, plon, hbuf3d(ip,lat), obuf )

               do i = 1, plon*nlev
                  hbuf3d(ip+i-1,lat) = 0.0
               end do
!              if (lat.eq.plat) write (6,*) 'accumed, renormalized and wrote accumulator to obuf '

            else                        ! accumulate field in history buffer

!              if (lat.eq.plat) write (6,*) ' accumulating field in hbuf3d '
               j = 0
               do k = 1, nlev
                  do i = 1, plon
                     hbuf3d(ip+j,lat) = hbuf3d(ip+j,lat) + field(i,k)
                     j = j + 1
                  end do
               end do

            end if
         end if
      end if

      end subroutine outfld

c##############################################################################

      integer function qhstc( name, cval )

c Query history file module for character data.

      implicit none

c Input argument:
      character
     $  name*(*)  ! variable name of requested character string

c Output argument:
      character
     $  cval*(*)  ! requested character string

c Return values:
c  0 => success
c -1 => variable name not found
c-----------------------------------------------------------------------

      qhstc = 0

      if ( name .eq. 'LPATH' ) then
         cval = lpath
         return
      else if ( name .eq. 'RPATH' ) then
         cval = rpath
         return
      end if

      qhstc = -1

      end function qhstc

c##############################################################################

      integer function qhsti( name, ival )

c Query history file module for integer data.

      implicit none

c Input argument:
      character
     $  name*(*)  ! variable name of requested integer data

c Output argument:
      integer
     $  ival      ! requested value

c Return values:
c  0 => success
c -1 => variable name not found
c-----------------------------------------------------------------------

      qhsti = 0

      if ( name .eq. 'NTSHST' ) then
         ival = hdi(3)
         return
      else if ( name .eq. 'MXTSHST' ) then
         ival = mfilt
         return
      else if ( name .eq. 'NHTFRQ' ) then
         ival = nhtfrq
         return
      else if ( name .eq. 'HBUFLEN' ) then
         ival = hbuflen
         return
      else if ( name .eq. 'MAXSIZ' ) then
         ival = maxsiz
         return
      else if ( name .eq. 'LUNHST' ) then
         ival = lunhst
         return
      else if ( name .eq. 'NDBASE' ) then
         ival = hdi(20)
         return
      else if ( name .eq. 'NSBASE' ) then
         ival = hdi(21)
         return
      else if ( name .eq. 'NDCUR' ) then
         ival = hdi(22)
         return
      else if ( name .eq. 'NSCUR' ) then
         ival = hdi(23)
         return
      else if ( name .eq. 'NBDATE' ) then
         ival = hdi(24)
         return
      else if ( name .eq. 'NBSEC' ) then
         ival = hdi(25)
         return
      else if ( name .eq. 'NCDATE' ) then
         ival = hdi(26)
         return
      else if ( name .eq. 'NCSEC' ) then
         ival = hdi(27)
         return
      end if

      qhsti = -1

      end function qhsti

c##############################################################################

      subroutine rhstrst( lun )

c Read restart data for history file.

c Header info that is updated during the run: (except for mfilht which
c is written on the info file)
c
c  hdi(17)  ! current iteration number
c  hdi(19)  ! timesteps since last time sample written
c  hdi(22)  ! current day corresponding to current timestep
c  hdi(23)  ! seconds relative to current day
c  hdi(26)  ! current date corresponding to current timestep
c  hdi(27)  ! current seconds relative to current date
c  lnhstf   ! mss pathname of first history file for this case

c  History buffer data and accumulation info.

      implicit none

c Input argument:
      integer
     $  lun     ! logical unit for fortran read

c Local variables:
      integer
     $  i, j       ! indices
c-----------------------------------------------------------------------

      read( lun ) hdi(17),hdi(19),hdi(22),hdi(23),hdi(26),hdi(27)
     $,           lnhstf
      do i = 1, 10
         hdc(24+i) = lnhstf(8*(i-1)+1:8*i)
      end do

      if ( hbuflen .gt. 0 ) then
         do j = 1, plat
            call rdarr( lun, hbuflen, hbuf3d(1,j) )
         end do         
      end if

      end subroutine rhstrst

c##############################################################################

      subroutine whstrst( lun )

c Write restart data for history file.

c Header info that is updated each time a time sample is written:
c (except for mfilht which is written on the info file)
c  hdi(17)  ! current iteration number
c  hdi(19)  ! timesteps since last time sample written
c  hdi(22)  ! current day corresponding to current timestep
c  hdi(23)  ! seconds relative to current day
c  hdi(26)  ! current date corresponding to current timestep
c  hdi(27)  ! current seconds relative to current date

c  lnhstf   ! mss pathname of first history file for this case

      implicit none

c Input argument:
      integer
     $  lun     ! logical unit for fortran write

c Local variables:
      integer
     $  j       ! latitude index
c-----------------------------------------------------------------------

      write( lun ) hdi(17),hdi(19),hdi(22),hdi(23),hdi(26),hdi(27)
     $,            lnhstf

      if ( hbuflen .gt. 0 ) then
         do j = 1, plat
!            write (6,*) ' whstrst: lat, min, max ', j,
!     $           minval(hbuf3d(1:hbuflen,j)), maxval(hbuf3d(1:hbuflen,j))
            call wrarr( lun, hbuflen, hbuf3d(1,j) )
         end do         
      end if

      end subroutine whstrst

c##############################################################################

      subroutine wrhd( nstep )

c Update and write header records.

      use calendar

      implicit none

c Input arguments:
      integer
     $  nstep        ! Current timestep number.

c Local variables:
      integer
     $  itim    ! time index at the end of timestep nstep
     $, nitslf  ! number of timesteps since last time sample
     $, totsec  ! total seconds since last time sample
     $, iday    ! days since last time sample
     $, isec    ! seconds since last time sample relative to eday
     $, ndcur   ! day of current time sample
     $, nscur   ! seconds relative to ndcur
     $, ncdate  ! date of current time sample
     $, ncsec   ! seconds relative to ncdate

c Externals:
!     integer newdate
c-----------------------------------------------------------------------

      itim = nstep

c     Update index of this time sample in the current history file.
      hdi(3) = hdi(3) + 1

c     Update number of timesteps since last time sample written to history file
      nitslf = itim - hdi(17)
      hdi(19) = nitslf
      fnorm = 1. / amax0( nitslf, 1 )

c     Update time index of the current history buffer.
      hdi(17) = itim

c     Time since last time sample written to history file.
      totsec = nitslf * hdi(28)
      iday = totsec/86400
      isec = mod( totsec, 86400 )

c     Update day
      ndcur = hdi(22) + iday
      nscur = hdi(23) + isec
      if ( nscur .ge. 86400 ) then
         ndcur = ndcur + nscur/86400
         nscur = mod( nscur, 86400 )
      end if
      hdi(22) = ndcur
      hdi(23) = nscur

c     Update date
      ncsec = hdi(27) + isec
      if ( ncsec .ge. 86400 ) then
         iday = iday + ncsec/86400
         ncsec = mod( ncsec, 86400 )
      end if
      ncdate = newdate( hdi(26), iday )
      hdi(26) = ncdate
      hdi(27) = ncsec

      write(*,*)'wrhd: writing history file header for time sample ',
     $     hdi(3), ' at date = ', ncdate, '  seconds = ', ncsec

c     Save day and date from first time sample on a file for the
c     comment string.
      if ( hdi(3) .eq. 1 ) then
         ndcurf = ndcur
         nscurf = nscur
         ncdatf = ncdate
         ncsecf = ncsec
      endif

c     Write out header records.  wrintrec used for integer and wrch8rec
c     used for character headers
c     because they are dynamically allocated, so their declarations don't
c     contain the correct size.
      call wrintrec( lunhst, lenhdi, hdi )
      call wrch8rec( lunhst, lenhdc, hdc )
      write( lunhst ) hdr

      end subroutine wrhd

c##############################################################################

      subroutine wrobuf( mflds, dim1, field, obuf )

c Write a latitude slice of data to the output buffer.  Packing is done
c here if it was requested.

       implicit none

c Input arguments:
      integer
     $  mflds(3)       ! integer header info for field being output
     $, dim1           ! 1st dimension of field array
      real(r8)
     $  field(dim1,*)  ! field values

c Output arguments:
      real(r8)
     $  obuf(*)        ! output buffer

c Local variables:
      integer
     $  i, k      ! indices
     $, iq        ! pointer into output buffer
     $, nlev      ! number of levels
     $, nden      ! packing density for current field
     $, iqlen     ! packed length for one latitude circle of data
      real(r8)
     $  amn, amx  ! dummy args for packaf
c-----------------------------------------------------------------------

      iq = mflds(2)
      nlev = 1
      if ( mod( mflds(1), 10 ) .gt. 0 ) nlev = plev
      nden = mflds(3)
      iqlen = 3 + (plon-1)/nden
c     Copy field to output data buffer, packing if requested.
      if ( nden .gt. 1 ) then
         do k = 1, nlev
            call packaf( field(1,k), plon, obuf(iq), nden, amn, amx )
            iq = iq + iqlen
         end do
      else
         do k = 1, nlev
            do i = 1, plon
               obuf(iq+i-1) = field(i,k)
            end do
            iq = iq + plon
         end do
      endif

      end subroutine wrobuf

c##############################################################################

      subroutine getspecres( plon, plat, ntrm, ntrn, ntrk )
c spectral resolution based on horizontal resolution
      implicit none
      integer plon, plat, ntrm, ntrn, ntrk
      if ( plon .eq. 16 .and. plat .eq. 8 ) then  ! T5
         ntrm = 5
         ntrn = 5
         ntrk = 5
      else if ( plon .eq. 64 .and. plat .eq. 32 ) then ! T21
         ntrm = 21
         ntrn = 21
         ntrk = 21
      else if ( plon .eq. 128 .and. plat .eq. 64 ) then  ! T42
         ntrm = 42
         ntrn = 42
         ntrk = 42
      else if ( plon .eq. 64 .and. plat .eq. 64 ) then ! T42x21
         ntrm = 21
         ntrn = 42
         ntrk = 42
      else if ( plon .eq. 192 .and. plat .eq. 94 ) then ! T62
         ntrm = 62
         ntrn = 62
         ntrk = 62
      else if ( plon .eq. 320 .and. plat .eq. 160 ) then ! T106
         ntrm = 106
         ntrn = 106
         ntrk = 106
      else if ( plon .eq. 384 .and. plat .eq. 190 ) then ! T126
         ntrm = 126
         ntrn = 126
         ntrk = 126
      else if ( plon .eq. 512 .and. plat .eq. 256 ) then ! T170
         ntrm = 170
         ntrn = 170
         ntrk = 170
      else
         ntrm = 0
         ntrn = 0
         ntrk = 0
         write(*,*)'getspecres: cannot determine spectral resolution'
         write(*,*)'   to write to output history  file.'
         write(*,*)'   Setting (ntrm,ntrn,ntrk) = (0,0,0)'
      end if

      end subroutine getspecres

c##############################################################################

      end module histout
