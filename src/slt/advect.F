! $Name$
! $Date$

      module advect

c Semi-Lagrangian Transport (SLT) algorithm.

      use precision
      use constants, only: ra

      implicit none
      save
      private
      public ::
     $  adv_type ! returns 'slt'
     $, ini_adv  ! initialize parameterization
     $, adv_dr   ! SLT advection driver
     $, radvrst  ! read restart data
     $, wadvrst  ! Write restart data.

c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

      character(len=*), parameter ::
     $  type = 'slt'

      integer, parameter::
     $  ppdy = 4             ! length of interpolation grid stencil
     $, pmap = 20000         ! Dimension of artificial evenly spaced vertical grid arrays.

      logical, parameter::
     $  plimdr = .true.      ! flag to limit derivatives

      real(r8) ::
     $ dlam                  ! Length of increment in longitude grid.
     $, lam(plond)            ! Longitude values in the extended grid.
     $, sinlam(plon)          ! sin(lam) on non-extended grid
     $, coslam(plon)          ! cos(lam) on non-extended grid
     $, phi(platd)            ! Latitude values in the extended grid.
     $, dphi(platd)           ! Interval between latitudes in the extended grid.
     $, etamid(plev)          ! Hybrid coordinate at layer midpoints.
     $, etaint(plevp)         ! Hybrid coordinate at layer interfaces.
     $, detam(plev)           ! del-eta (between layer midpoints)
     $, detai(plevp)          ! del-eta (between layer interfaces)
     $, lbasdy(4,2,platd)     ! deriv basis functions in y-direction
     $, lbasdz(4,2,plev)      ! deriv basis functions for full levels in z
     $, lbassd(4,2,plevp)     ! deriv basis functinos for half levels in z
     $, lbasiy(4,2,platd)     ! interp basis functions in y-direction

      real(r8), dimension(:,:,:), allocatable ::
     $  lammp                 ! x-coordinate of trajectory midpoint (plon,plev,plat)
     $, phimp                 ! y-coordinate of trajectory midpoint
     $, sigmp                 ! z-coordinate of trajectory midpoint

      integer ::
     $  kdpmpf(pmap)          ! artificial full vertical grid indices
     $, kdpmph(pmap)          ! artificial half vertical grid indices

      real(r8), dimension(:,:,:,:), allocatable ::
     $  xxl  ! left x-deriv of misc transported variables
     $, xxr  ! right x-deriv of misc transported variables

c##############################################################################
      contains
c##############################################################################

      function adv_type()

! Return advection type.

      implicit none
      character(len=len(type)) :: adv_type
      adv_type = type
      end function adv_type

c##############################################################################

      subroutine ini_adv( xdlam, xlam, xphi, gw
     $,                   xetaint, xetamid, hybi, nadv )

c Initialize SLT advection module.

      implicit none

c Input arguments:
      real(r8)
     $ xdlam           ! longitude spacing
     $, xlam(plond)     ! longitudes in extended SLT grid (radians)
     $, xphi(platd)     ! Gauss latitudes in extended SLT grid (radians)
     $, gw(plat)        ! gauss weights for latitude grid (not used)
     $, xetaint(plevp)  ! hybrid coords at layer interfaces
     $, xetamid(plev)   ! hybrid coords at layer midpoints
     $, hybi(plevp)     ! sigma component of hybrid coordinate at interfaces (not used)
      integer
     $  nadv            ! number of misc transported species

c Local variables:
      integer
     $  i, ix, j, k             ! indices
     $, istat
      real(DBLKIND)
     $  etamln(plev)        ! log(etamid)
     $, etailn(plevp)       ! log(etaint)
     $, detamln(plev)       ! dlog(etamid)
     $, detailn(plevp)      ! dlog(etaint)
c-----------------------------------------------------------------------

      dlam = xdlam
      do i = 1, plond
         lam(i) = xlam(i)
      end do
      do j = 1, platd
         phi(j) = xphi(j)
      end do

      do k = 1, plev
         etamid(k) = xetamid(k)
      end do
      do k = 1, plevp
         etaint(k) = xetaint(k)
      end do

      write (6,*) ' eta top, bot and diff ', 
     $     etaint(plevp), etaint(1), etaint(plevp)-etaint(1) 

c     Sine and cosine of longitude.

      do i = 1, plon
         ix = nxpt + i
         sinlam(i) = sin( lam(ix) )
         coslam(i) = cos( lam(ix) )
      end do

c     Basis functions for computing Lagrangian cubic derivatives
c     on unequally spaced latitude and vertical grids.

      call basdy(phi,lbasdy)
      call basdz(plev,etamid,lbasdz)
      call basdz(plevp,etaint,lbassd)

c     Basis functions for computing weights for Lagrangian cubic
c     interpolation on unequally spaced latitude grids.

      call basiy(phi,lbasiy)

c     Compute interval lengths in latitudinal grid

      do j = 1,platd-1
         dphi(j) = phi(j+1) - phi(j)
      end do

c     Compute interval lengths in vertical grids.


      do k = 1,plev
         etamln(k) = log(etamid(k))
      end do
      do k = 1,plevp
         etailn(k) = log(max(etaint(k),etamid(1)*0.01))
      end do
      do k = 1,plev-1
         detam  (k) = etamid(k+1) - etamid(k)
         detamln(k) = etamln(k+1) - etamln(k)
      end do
      detamln(plev) = 0.0
      do k = 1,plev
         detai  (k) = etaint(k+1) - etaint(k)
         detailn(k) = etailn(k+1) - etailn(k)
      end do
      detailn(plevp) = 0.0

c     Build artificial evenly spaced vertical grid for use in determining
c     vertical position of departure point.
c     Build one grid for full model levels and one for half levels.

      call vrtmap(plev ,etamln,detamln,kdpmpf)
      call vrtmap(plevp,etailn,detailn,kdpmph)

c     Allocate memory for derivatives of misc transported variables.
      allocate( lammp(plon,plev,plat), stat=istat )
      call alloc_err( istat, 'ini_adv', 'lammp', plon*plev*plat)
      allocate( phimp(plon,plev,plat), stat=istat )
      call alloc_err( istat, 'ini_adv', 'phimp', plon*plev*plat)
      allocate( sigmp(plon,plev,plat), stat=istat )
      call alloc_err( istat, 'ini_adv', 'sigmp', plon*plev*plat)

c     Initialize mid-point trajectory estimates
c       (Initial guess will be the arrival point)
      do j = 1, plat
         do k = 1, plev
            do i = 1, plon
               lammp(i,k,j) = real(i-1,r8)*dlam
               phimp(i,k,j) = phi(j + nxpt + jintmx)
               sigmp(i,k,j) = etamid(k)
            end do
         end do
      end do

c     Allocate memory for derivatives of misc transported variables.
      if ( nadv .gt. 0 ) then
         allocate( xxl(plond,plev,nadv,platd), stat=istat )
         call alloc_err( istat, 'ini_adv', 'xxl', plond*plev*nadv*platd)
         allocate( xxr(plond,plev,nadv,platd), stat=istat )
         call alloc_err( istat, 'ini_adv', 'xxr', plond*plev*nadv*platd)
      end if

      end subroutine ini_adv

c##############################################################################

      subroutine adv_dr( nstep, dtime, ufld, wfld
     $,                  ps, as1, ps2, as2
     $,                  hw1, nadv, x1, x2 
     $,                  mcar, mca, wr)

c Latitude loop over semi-lagrangian transport forecasts.

      implicit none

c Input arguments:
      integer 
     $  nstep        ! current model timestep
     $, nadv         ! number of advected variables (non-species)
      real(r8)
     $  dtime                       ! time step in seconds
     $, ufld(plond,plev,2,platd)    ! u,v wind components
     $, ps(plond,plat)              ! dry surface pressure at begining of timestep
     $, as1(plond,plev,pcnst,platd) ! advected species before slt
     $, ps2(plond,plat)             ! dry surface pressure at end of timestep.  (not used)
     $, x1(plond,plev,padv,platd)   ! advected variables before slt
     $, wr(plon,plevp,plat) ! raw vertical wind (etadot) (1/s)
     $, mca(plon,plevp,plat)        ! raw mass convergence above (etadot) (Pa/s)
     $, mcar(plon,plevp,plat)       ! mass convergence above (etadot) (Pa/s)


c Input/Output arguments:
      real(r8)
     $  wfld(plon,plevp,plat)       ! vertical wind ("eta-dot").  Is either
                                    ! input from dynamics files or calculated.

c Output arguments:
      real(r8)
     $  as2(plond,plev,pcnst,platd) ! advected species after slt
     $, x2(plond,plev,padv,platd)   ! advected variables after slt
     $, hw1(pcnst)                  ! global mass before slt

c Local variables:
      integer
     $  itermx            ! max # of iters for trajectory calc
     $, itermn            ! min # of iters for trajectory calc
      parameter( itermx = 4, itermn = 1 )

      integer
     $  iter     ! number of iterations for departure point calculation
     $, j        ! latitude index
     $, jcen     ! lat index (extended grid) of forecast
     $, m        ! constituent index

      real(r8)
     $  uxl(plond,plev,2,platd)     ! left  x-deriv of u/v
     $, uxr(plond,plev,2,platd)     ! right  x-deriv of u/v
     $, qxl(plond,plev,pcnst,platd) ! left  x-deriv of constituents
     $, qxr(plond,plev,pcnst,platd) ! right  x-deriv of constituents
      real(r8)
     $  pmid(plond,plev)            ! pressure at model levels
     $, pint(plond,plevp)           ! pressure at interfaces
     $, pdel(plond,plev)            ! pressure difference between interfaces
     $, d(plon,plev)                ! divergence

      real(r8)
     $  hw1j(pcnst,plat)            ! latitude contribution to global mass

c Externals:
      logical
     $  calc_etadot     ! estimate etadot internally
c------------------------------------------------------------------------------
      wr = 0.
      mca = 0.
      mcar = 0.
      hw1j = 0.0

      if (nstep .le. 1) then
         iter = itermx
      else
         iter = itermn
      end if

c     Extend arrays.
      call sltini( sinlam, coslam, ufld, as1, uxl, uxr, qxl, qxr
     $,            nadv, x1 )

! OpenMP directive
!$OMP PARALLEL DO PRIVATE(JCEN, PINT, PMID, PDEL, D)
      do j = 1, plat
         jcen = j + nxpt + jintmx

c        Compute global mass integrals "hw1"
         call plevs0( ps(1,j), pint, pmid, pdel )
         call gamdry2( j, pcnst, as1(i1,1,1,jcen), pdel, hw1j(1,j) )

         if ( calc_etadot() ) then

c           Approximate the divergence field.
            call dcal( j, ufld(1,1,1,jcen-1), ufld(1,1,2,jcen-1), d )

c           Calculate etadot.
            call etdcal( j, ufld(i1,1,1,jcen), ufld(i1,1,2,jcen), d, ps
     $,                  pdel, pmid, wfld(1,1,j) )
         end if

c        Semi-Lagrangian Transport driver.

         call sltb1( jcen, j, dtime, iter, ufld, uxl, uxr
     $,              wfld(1,1,j), as1, qxl, qxr, lam, phi, dphi
     $,              etamid, etaint, detam, detai, lbasdy, lbasdz
     $,              lbassd, lbasiy, kdpmpf, kdpmph, lammp(1,1,j)
     $,              phimp(1,1,j), sigmp(1,1,j), as2(i1,1,1,jcen)
     $,              nadv, x1, x2(i1,1,1,jcen) )

      end do

c Compute global integral hw1 on one cpu for identical answers
c single-threaded vs. multitasked

      do m = 1, pcnst
         hw1(m) = 0.
      end do
      do m = 1, pcnst
         do j = 1, plat
            hw1(m) = hw1(m) + hw1j(m,j)
         end do
      end do

      end subroutine adv_dr

c##############################################################################

      subroutine radvrst( lun )

c Read restart data for slt module.

      implicit none

c Input argument:
      integer
     $  lun    ! logical unit for fortran i/o.
c-----------------------------------------------------------------------

c     Departure points.

      read( lun ) lammp, phimp, sigmp

      end subroutine radvrst

c##############################################################################

      subroutine wadvrst( lun )

c Write restart data from slt module.

      implicit none

c Input argument:

      integer
     $  lun    ! logical unit for fortran write.
c-----------------------------------------------------------------------

c     Departure points.
      write( lun ) lammp, phimp, sigmp

      end subroutine wadvrst

c##############################################################################

      subroutine sltb1(jcen    ,jgc     ,dt      ,
     $                 iterdp  ,ub      ,uxl     ,uxr     ,wb      ,
     $                 fb      ,fxl     ,fxr     ,lam     ,phib    ,
     $                 dphib   ,sig     ,sigh    ,dsig    ,dsigh   ,
     $                 lbasdy  ,lbasdz  ,lbassd  ,lbasiy  ,kdpmpf  ,
     $                 kdpmph  ,lammp   ,phimp   ,sigmp   ,fbout   ,
     $                 nadv, x1, x2 )
C-----------------------------------------------------------------------
C
C Drive the slt algorithm on a given latitude slice in the extended
C data arrays using information from the entire latitudinal extent
C of the arrays.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none

      real(r8) phigs                      ! cut-off latitude: about 70 degrees
      parameter (phigs = 1.221730)
C
C Input arguments
C
      integer
     $     jcen,                            ! index of lat slice(extend)
     $     jgc                              ! index of lat slice (model)
      real(r8) dt                               ! time step (seconds)
      integer
     $     iterdp                           ! iteration count
     $,    nadv                             ! number of misc advected quantities

      real(r8) ub (plond,plev,        2,platd), ! u in extended grid
     $     uxl(plond,plev,        2,platd), ! left  x-deriv of ub/vb
     $     uxr(plond,plev,        2,platd), ! right x-deriv of ub/vb
     $     wb(plon,plevp),                  ! eta-dot
     $     fb (plond,plev,  pcnst,platd),   ! constituent fields
     $     fxl(plond,plev,  pcnst,platd),   ! left  x-deriv of fb
     $     fxr(plond,plev,  pcnst,platd),   ! right x-deriv of fb
     $     x1 (plond,plev,padv,platd),      ! misc advected quantities
     $     lam  (plond),                    ! long. coord of model grid
     $     phib (platd),                    ! lat.  coord of model grid
     $     dphib(platd),                    ! increment between lats.
     $     sig  (plev),                     ! vertical full levels
     $     sigh (plevp),                    ! vertical half levels
     $     dsig (plev),                     ! increment between full levs
     $     dsigh(plevp),                    ! increment between half levs
     $     lbasdy(4,2,platd),               ! lat deriv weights
     $     lbasdz(4,2,plev),                ! vert full level deriv wts
     $     lbassd(4,2,plevp),               ! vert half level deriv wts
     $     lbasiy(4,2,platd)                ! lat interp wts(lagrng)

      integer
     $     kdpmpf(pmap),                    ! artificial vert grid indices
     $     kdpmph(pmap)                     ! artificial vert grid indices
C
C Input/Output arguments
C
      real(r8) lammp(plon,plev),                ! long coord of mid-point
     $     phimp(plon,plev),                ! lat  coord of mid-point
     $     sigmp(plon,plev)                 ! vert coord of mid-point
C
C Output arguments
C
      real(r8) fbout(plond,plev,pcnst)          ! advected constituents
     $, x2(plond,plev,padv)                 ! misc advected quantities
C
C-----------------------------------------------------------------------
C
C  jcen    Latitude index in extended grid corresponding to lat slice
C          being forecasted.
C  jgc     Latitude index in model    grid corresponding to lat slice 
C          being forecasted.
C  dt      Time interval that parameterizes the parcel trajectory.
C  iterdp  Number of iterations used for departure point calculation.
C  ub      u-velocity component.
C  uxl     uxl ( , , 1, ) = x-derivatives of u at the left  (west) edge
C          of given interval
C          uxl ( , , 2, ) = x-derivatives of v at the left  (west) edge
C          of given interval
C  uxr     uxr ( , , 1, ) = x-derivatives of u at the right (east) edge
C          of given interval
C          uxr ( , , 2, ) = x-derivatives of v at the right (east) edge
C          of given interval
C  wb      z-velocity component (eta-dot).
C  fb      Scalar components to be transported.
C  fxl     x-derivatives at the left  edge of each interval containing 
C          the departure point.
C  fxr     x-derivatives at the right edge of each interval containing 
C          the departure point.
C  lam     Longitude values for the extended grid.
C  phib    Latitude  values for the extended grid.
C  dphib   Interval between latitudes in the extended grid.
C  sig     Hybrid eta values at the "full-index" levels.
C  sigh    Half-index eta-levels including sigh(i,1) = eta(1/2) = 0.0
C          and sigh(i,plev+1) = eta(plev+1/2) = 1.  Note that in general
C          sigh(i,k) .lt. sig(i,k)  where sig(i,k) is the hybrid value
C          at the k_th full-index level.
C  dsig    Interval lengths in full-index hybrid level grid.
C  dsigh   Interval lengths in half-index hybrid level grid.
C  lbasdy  Weights for Lagrange cubic derivative estimates on the
C          unequally spaced latitude grid.
C  lbasdz  Weights for Lagrange cubic derivative estimates on the
C          unequally spaced vertical grid (full levels).
C  lbassd  Weights for Lagrange cubic derivative estimates on the
C          unequally spaced vertical grid (half levels).
C  lbasiy  Weights for Lagrange cubic interpolation on the unequally
C          spaced latitude grid.
C  kdpmpf  indices of artificial grid mapped into the full level grid
C  kdpmph  indices of artificial grid mapped into the half level grid
C  lammp   Longitude coordinates of the trajectory mid-points of the
C          parcels that correspond to the global grid points contained
C          in the latitude slice being forecasted.  On entry lammp
C          is an initial guess.
C  phimp   Latitude coordinates of the trajectory mid-points of the
C          parcels that correspond to the global grid points contained
C          in the latitude slice being forecasted.  On entry phimp
C          is an initial guess.
C  sigmp   Hybrid value at the trajectory midpoint for each gridpoint
C          in a vertical slice from the global grid.  On entry sigmp is
C          an initial guess.
C  fbout   Extended array only one latitude of which, however, is filled
C          with forecasted (transported) values.  This routine must be
C          called multiple times to fill the entire array.  This is
C          done to facilitate multi-tasking.
C  diag    Time tendency due to horizontal advection for each scalar
C          component at the global grid points that correspond to the
C          latitude slice being forecasted.  As in "fbout", only one
C          latitude slice at a time is filled.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     istat,
     $     m,                          ! constituent index
     $     idp(plon,plev),             ! zonal      dep point index
     $     jdp(plon,plev),             ! meridional dep point index
     $     kdp(plon,plev)              ! vertical   dep point index
      real(r8) fhr(plon,plev,pcnst),       ! horizontal interpolants
     $     lamdp(plon,plev),           ! zonal      departure pt. coord.
     $     phidp(plon,plev),           ! meridional departure pt. coord.
     $     sigdp(plon,plev),           ! vertical   departure pt. coord.
     $     fhst(plon,plev,pcnst),      ! derivative at top of interval
     $     fhsb(plon,plev,pcnst),      ! derivative at bot of interval
     $     wst(plon,plevp),            ! w derivative at top of interval
     $     wsb(plon,plevp),            ! w derivative at bot of interval
     $     fint(plon,plev,ppdy,pcnst), ! work space
     $     fyb(plon,plev,pcnst),       ! work space
     $     fyt(plon,plev,pcnst)        ! work space

      logical
     $     locgeo                      ! flag indicating coordinate sys

c     Analogs of fhr, fhst, fhsb, fint, fyb, fyt for misc adv. quantities.
      real(r8), allocatable, dimension(:,:,:) ::
     $  xhr, xhst, xhsb, xyb, xyt
      real(r8), allocatable, dimension(:,:,:,:) ::
     $  xint
C-----------------------------------------------------------------------

c     Allocate memory for misc advected quantities.
      if ( nadv .gt. 0 ) then
         allocate( xhr(plon,plev,nadv), stat=istat )
         call alloc_err( istat, 'sltb1', 'xhr', plon*plev*nadv )
         allocate( xhst(plon,plev,nadv), stat=istat )
         call alloc_err( istat, 'sltb1', 'xhst', plon*plev*nadv )
         allocate( xhsb(plon,plev,nadv), stat=istat )
         call alloc_err( istat, 'sltb1', 'xhsb', plon*plev*nadv )
         allocate( xyb(plon,plev,nadv), stat=istat )
         call alloc_err( istat, 'sltb1', 'xyb', plon*plev*nadv )
         allocate( xyt(plon,plev,nadv), stat=istat )
         call alloc_err( istat, 'sltb1', 'xyt', plon*plev*nadv )
         allocate( xint(plon,plev,ppdy,nadv), stat=istat )
         call alloc_err( istat, 'sltb1', 'xint', plon*plev*ppdy*nadv )
      else
         allocate( xhr(1,1,1), stat=istat )
         call alloc_err( istat, 'sltb1', 'xhr', 1*1*1 )
         allocate( xhst(1,1,1), stat=istat )
         call alloc_err( istat, 'sltb1', 'xhst', 1*1*1 )
         allocate( xhsb(1,1,1), stat=istat )
         call alloc_err( istat, 'sltb1', 'xhsb', 1*1*1 )
         allocate( xyb(1,1,1), stat=istat )
         call alloc_err( istat, 'sltb1', 'xyb', 1*1*1 )
         allocate( xyt(1,1,1), stat=istat )
         call alloc_err( istat, 'sltb1', 'xyt', 1*1*1 )
         allocate( xint(1,1,1,1), stat=istat )
         call alloc_err( istat, 'sltb1', 'xint', 1*1**1 )
      end if

      xhr = 0.
      xhst = 0.
      xhsb = 0.
C
C Horizontal interpolation
C
C Compute departure points and corresponding indices.
C Poleward of latitude phigs (radians), perform the computation in
C local geodesic coordinates.
C Equatorward of latitude phigs, perform the computation in global
C spherical coordinates
C
      locgeo = abs(phib(jcen)).ge.phigs
C
      call sphdep(jcen    ,jgc     ,dt      ,iterdp  ,
     $            locgeo  ,ub      ,uxl     ,uxr     ,lam     ,
     $            phib    ,lbasiy  ,lammp   ,phimp   ,lamdp   ,
     $            phidp   ,idp     ,jdp     )
C
C Interpolate scalar fields to the departure points.
C
      call hrintp(pcnst   ,pcnst   ,fb      ,fxl     ,fxr     ,
     $            lam     ,phib    ,dphib   ,lbasdy  ,lamdp   ,
     $            phidp   ,idp     ,jdp     ,plimdr  ,fint    ,
     $            fyb     ,fyt     ,fhr
     $,           nadv, x1, xint, xyb, xyt, xhr )
C
C Vertical interpolation.
C Compute vertical derivatives of vertical wind
C
      call cubzdr(plon    ,plevp   ,wb      ,lbassd  ,wst     ,
     $            wsb     )
C
C Compute departure points and corresponding indices.
C
      call vrtdep(dt      ,iterdp  ,wb      ,wst     ,
     $            wsb     ,sig     ,sigh    ,dsigh   ,kdpmpf  ,
     $            kdpmph  ,sigmp   ,sigdp   ,kdp     )
C
C Vertical derivatives of scalar fields.
C Loop over constituents.
C
      fhst = 0.
      do m = 1,pcnst
         call cubzdr(plon    ,plev    ,fhr(1,1,m),lbasdz  ,fhst(1,1,m),
     $               fhsb(1,1,m)      )
      end do
      if( plimdr )then
         call limdz(pcnst,fhr,dsig,fhst,fhsb)
      end if

      if ( nadv .gt. 0 ) then
         do m = 1, nadv
            call cubzdr( plon, plev, xhr(1,1,m), lbasdz, xhst(1,1,m)
     $,                  xhsb(1,1,m) )
         end do
         if( plimdr )then
            call limdz(nadv,xhr,dsig,xhst,xhsb)
         end if
      end if
C
C Vertical interpolation of scalar fields.
C
      call herzin(plev    ,pcnst   ,fhr     ,fhst    ,fhsb    ,
     $            sig     ,dsig    ,sigdp   ,kdp     ,fbout
     $,           nadv, xhr, xhst, xhsb, x2 )
C

c     Free memory for misc advected quantities.
      if ( nadv .gt. 0 ) then
         deallocate( xhr, xhst, xhsb, xyb, xyt, xint, stat=istat )
         if ( istat .ne. 0 ) then
            write(*,*)'ERROR trying to deallocate memory in sltb1.'
            call endrun
         end if
      end if

      end subroutine sltb1

c##############################################################################

      subroutine bandij(dlam    ,phib    ,lamp    ,phip    ,iband   ,
     $                  jband   )
C-----------------------------------------------------------------------
C
C Calculate longitude and latitude indices that identify the
C intervals in the extended grid that contain the departure points.
C Upon entry, all dep. points should be within jintmx intervals of the
C Northern- and Southern-most model latitudes. Note: the algorithm
C relies on certain relationships of the intervals in the Gaussian grid.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C
C-----------------------------------------------------------------------
      implicit none
      integer pgls              ! length of latitude slice
      parameter (pgls = plon*plev) 
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) dlam,                ! longitude increment
     $     phib(platd),         ! latitude  coordinates of model grid
     $     lamp(pgls),          ! longitude coordinates of dep. points
     $     phip(pgls)           ! latitude  coordinates of dep. points
C
C Output arguments
C
      integer
     $     iband(pgls),         ! longitude index of dep. points
     $     jband(pgls)          ! latitude  index of dep. points
C
C-----------------------------------------------------------------------
C
C  dlam    Length of increment in equally spaced longitude grid (radians)
C  phib    Latitude values for the extended grid.
C  lamp    Longitude coordinates of the points.  It is assumed that
C                     0.0 .le. lamp(i) .lt. 2*pi .
C  phip    Latitude coordinates of the points.
C  iband   Longitude index of the points.  This index points into
C          the extended arrays, e.g.,
C                   lam(iband(i)) .le. lamp(i) .lt. lam(iband(i)+1) .
C  jband   Latitude index of the points.  This index points into
C          the extended arrays, e.g.,
C                   phib(jband(i)) .le. phip(i) .lt. phib(jband(i)+1) .
C
C---------------------------Local workspace-----------------------------
C
      integer
     $     i                    ! index
      real(r8) dphibr,              ! reciprocal of an approximate del phi
     $     phibs                ! latitude of southern-most latitude
C
      dphibr = 1./( phib(platd/2+1) - phib(platd/2) )
      phibs  = phib(1)
C
C Longitude indices.
C
      do i = 1,pgls
         iband(i) = i1 + int( lamp(i)/dlam )
      end do
C
C Latitude indices.
C
      do i = 1,pgls
         jband(i) = int ( (phip(i) - phibs)*dphibr + 1. )
         if( phip(i) .ge. phib(jband(i)+1) ) jband(i) = jband(i) + 1
      end do
C
      end subroutine bandij

c##############################################################################

      subroutine basdy(phi,lbasdy)
C-----------------------------------------------------------------------
C
C Compute weights for the calculation of derivative estimates at the two
C center points of the four point stencil for each interval in the
C unequally spaced latitude grid. Estimates are from differentiating
C a Lagrange cubic polynomial through the four point stencil.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C-----------------------------------------------------------------------
      integer
     $     jfirst,              ! first index to be computed
     $     jlast                ! last  index to be computed
      parameter(jfirst = nxpt + 1,
     $          jlast  = platd - nxpt - 1)
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) phi(platd)           ! latitude coordinates of model grid
C
C Output arguments
C
      real(r8) lbasdy(4,2,platd)    ! derivative estimate weights
C
C-----------------------------------------------------------------------
C
C  phi     Latitude values in the extended grid.
C  lbasdy  Weights for derivative estimates based on Lagrange cubic
C          polynomial on the unequally spaced latitude grid.
C          If grid interval j (in extended grid) is surrounded by
C          a 4 point stencil, then the derivative at the "bottom"
C          of the interval uses the weights lbasdy(1,1,j),
C          lbasdy(2,1,j), lbasdy(3,1,j), and lbasdy(4,1,j).
C          The derivative at the "top" of the interval
C          uses lbasdy(1,2,j), lbasdy(2,2,j), lbasdy(3,2,j),
C          and lbasdy(4,2,j).
C
C---------------------------Local variables-----------------------------
C
      integer
     $     jj                   ! index
C
C-----------------------------------------------------------------------
C
      do jj = jfirst,jlast
         call lcdbas( phi(jj-1), lbasdy(1,1,jj), lbasdy(1,2,jj) )
      end do
C
      end subroutine basdy

c##############################################################################

      subroutine basdz(pkdim,sig,lbasdz)
C-----------------------------------------------------------------------
C
C Compute weights for the calculation of derivative estimates at two
C center points of the four point stencil for each interval in the
C unequally spaced vertical grid (as defined by the array sig).
C Estimates are from differentiating a Lagrange cubic polynomial
C through the four point stencil.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      integer
     $     pkdim                ! vertical dimension
      real(r8) sig(pkdim)           ! sigma levels
C
C Output arguments
C
      real(r8) lbasdz(4,2,pkdim)    ! vertical interpolation weights
C
C-----------------------------------------------------------------------
C
C  pkdim   Number of grid points in vertical grid.
C  sig     Sigma values in the vertical grid.
C  lbasdz  Weights for derivative estimates based on Lagrange cubic
C          polynomial on the unequally spaced vertical grid.
C          If grid interval j is surrounded by a 4 point stencil,
C          then the derivative at the "top" of the interval (smaller
C          sigma value) uses the weights lbasdz(1,1,j),lbasdz(2,1,j),
C          lbasdz(3,1,j), and lbasdz(4,1,j).  The derivative at the
C          "bottom" of the interval uses lbasdz(1,2,j), lbasdz(2,2,j),
C          lbasdz(3,2,j), and lbasdz(4,2,j).  (Recall the vertical
C          level indices increase from the top of the atmosphere
C          towards the bottom.)
C
C---------------------------Local variables-----------------------------
C
      integer
     $     kk                   ! index
C
C-----------------------------------------------------------------------
C
      do kk = 2,pkdim-2
         call lcdbas( sig(kk-1), lbasdz(1,1,kk), lbasdz(1,2,kk) )
      end do
C
      end subroutine basdz

c##############################################################################

      subroutine basiy(phi,lbasiy)
C-----------------------------------------------------------------------
C
C Compute weights used in Lagrange cubic polynomial interpolation in 
C the central interval of a four point stencil.  Done for each interval
C in the unequally spaced latitude grid.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C-----------------------------------------------------------------------
      integer
     $     jfirst,              ! index of first latitude interval
     $     jlast                ! index of last  latitude interval
      parameter(jfirst = nxpt + 1,
     $          jlast  = platd - nxpt - 1)
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) phi(platd)           ! grid values in extended grid
C
C Output arguments
C
      real(r8) lbasiy(4,2,platd)    ! Weights for Lagrange cubic interp
C
C---------------------------Local variables-----------------------------
C
      integer
     $     jj                   ! index
C
C-----------------------------------------------------------------------
C
      do jj = jfirst,jlast
         call lcbas( phi(jj-1),lbasiy(1,1,jj),lbasiy(1,2,jj) )
      end do
C
      end subroutine basiy

c##############################################################################

      subroutine copyac(a       ,b       ,kno)
C-----------------------------------------------------------------------
C
C Copy data from array a to array b
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath, April 1992
C Standardized:      L. Bath, June 1992
C
C-----------------------------------------------------------------------
      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) a(*)       ! Input array name
      integer kno     ! Number of words to copy
C
C Output arguments
C
      real(r8) b(*)       ! Output array name
C
C---------------------------Local variables-----------------------------
C
      integer i       ! Loop index
C
C-----------------------------------------------------------------------
C
      do i=1,kno
         b(i)=a(i)
      end do
C
      end subroutine copyac

c##############################################################################

      subroutine cubxdr(pidim   ,ibeg    ,len     ,dx      ,f       ,
     $                  fxl     ,fxr     )
C-----------------------------------------------------------------------
C
C Compute Lagrangian cubic derivative estimates for data on an equally
C spaced grid.  Suppose grid interval i is centered in a 4 point
C stencil consisting of grid points i-1, i, i+1, and i+2.  Then the
C derivative at the left edge of the interval (i.e., grid point i)
C is stored in fxl(i), and the derivative at the right edge of the
C interval (i.e., grid point i+1) is stored in fxr(i).  Note that
C fxl(i) is not necessarily equal to fxr(i-1) even though both of
C these values are estimates of the derivative at grid point i.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      integer
     $     pidim,               ! dimension
     $     ibeg,                ! starting index to perform computation
     $     len                  ! length over which to perform comp.
      real(r8) dx,                  ! grid interval
     $     f(pidim)             ! input field values
C
C Output arguments
C
      real(r8) fxl(pidim),          ! left  derivative of interval i in "f"
     $     fxr(pidim)           ! right derivative of interval i in "f"
C-----------------------------------------------------------------------
C
C  pidim   Length of f, fxl, and fxr.
C  ibeg    First interval of grid for which derivatives are computed.
C  len     Number of grid intervals for which derivatives are computed.
C          (There are pidim - 1 intervals between the pidim gridpoints
C          represented in f, fxl, and fxr.)
C  dx      Value of grid spacing.
C  f       Values on equally spaced grid for which derivatives are
C          computed.
C  fxl     fxl(i) is the derivative at the left  edge of interval i.
C  fxr     fxr(i) is the derivative at the right edge of interval i.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                   ! index
     $     iend                 ! index denoting end of computation
      real(r8) rdx6                 ! normalization weight
C
C-----------------------------------------------------------------------
C
      iend = ibeg + len - 1
      rdx6 = 1._r8/(6.*dx)
C
      do i = ibeg,iend
         fxl(i) = ( -2.*f(i-1) - 3.*f(i) + 6.*f(i+1) -    f(i+2) )*rdx6
         fxr(i) = (     f(i-1) - 6.*f(i) + 3.*f(i+1) + 2.*f(i+2) )*rdx6
      end do
C
      end subroutine cubxdr

c##############################################################################

      subroutine cubydr(pf      ,fint    ,wdy     ,jdp     ,fyb     ,
     $                  fyt
     $,                 nadv, xint, xyb, xyt )
C-----------------------------------------------------------------------
C
C Compute Lagrangian cubic derivative estimates at both ends of the
C intervals in the y coordinate (unequally spaced) containing the
C departure points for the latitude slice being forecasted.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
      integer  pgls
      parameter(pgls = plon*plev) ! number of pts in vertical slice
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer pf                ! number of constituent fields
     $,    nadv                 ! number of misc advected quantities
      real(r8) fint(pgls,ppdy,pf),  ! constituent x- interpolants
     $     xint(pgls,ppdy,padv),! x interpolants
     $     wdy(4,2,platd)       ! latitude interpolation weights
      integer jdp(pgls)         ! indices of latitude intervals
C
C Output arguments
C
      real(r8) fyb(pgls,pf),        ! Derivative at south end of interval
     $     fyt(pgls,pf)         ! Derivative at north end of interval
     $,    xyb (pgls,padv)      ! y-derivatives at bot of int.
     $,    xyt (pgls,padv)      ! y-derivatives at top of int.
C-----------------------------------------------------------------------
C
C  pf      Number of fields being interpolated.
C  fint    (fint(i,k,j,m),j=1,ppdy) contains the x interpolants at each
C          latitude needed for the y derivative estimates at the
C          endpoints of the interval that contains the departure point
C          for grid point (i,k).  The last index of fint allows for
C          interpolation of multiple fields.  fint is generated by a
C          call to herxin.
C  wdy     Weights for Lagrange cubic derivative estimates on the
C          unequally spaced latitude grid.  If grid interval j (in
C          extended array) is surrounded by a 4 point stencil, then
C          the derivative at the "bottom" of the interval uses the
C          weights wdy(1,1,j),wdy(2,1,j), wdy(3,1,j), and wdy(4,1,j).
C          The derivative at the "top" of the interval uses wdy(1,2,j),
C          wdy(2,2,j), wdy(3,2,j), and wdy(4,2,j).
C  jdp     jdp(i,k) is the index of the y-interval that contains the
C          departure point corresponding to global grid point (i,k) in
C          the latitude slice being forecasted.
C          Suppose yb contains the y-coordinates of the extended array
C          and ydp(i,k) is the y-coordinate of the departure point
C          corresponding to grid point (i,k).  Then,
C                yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
C  fyb     fyb(i,k,.) is the derivative at the bottom of the y interval
C          that contains the departure point of global grid point (i,k).
C  fyt     fyt(i,k,.) is the derivative at the top of the y interval
C          that contains the departure point of global grid point (i,k).
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                   ! index
     $     m                    ! index
      real(r8) wdytem(pgls,4,2)     ! work array to gather weights
C
C-----------------------------------------------------------------------
C
C For sake of efficiency,
C   if (pf .eq. 1) then
C      run first piece of code
C   else
C      run second piece of code
C   end if
C
C-----------------------------------------------------------------------
C
      if( pf+nadv .eq. 1) then
         do i = 1,pgls
            fyb(i,1) = wdy(1,1,jdp(i))*fint(i,1,1) +
     $                 wdy(2,1,jdp(i))*fint(i,2,1) +
     $                 wdy(3,1,jdp(i))*fint(i,3,1) +
     $                 wdy(4,1,jdp(i))*fint(i,4,1)
C
            fyt(i,1) = wdy(1,2,jdp(i))*fint(i,1,1) +
     $                 wdy(2,2,jdp(i))*fint(i,2,1) +
     $                 wdy(3,2,jdp(i))*fint(i,3,1) +
     $                 wdy(4,2,jdp(i))*fint(i,4,1)
         end do
      else
C
C Load temp arrays with weights for Lagrangian cubic derivative
C estimates on unequally spaced grid.
C
         do i = 1,pgls
            wdytem(i,1,1) = wdy(1,1,jdp(i))
            wdytem(i,2,1) = wdy(2,1,jdp(i))
            wdytem(i,3,1) = wdy(3,1,jdp(i))
            wdytem(i,4,1) = wdy(4,1,jdp(i))
            wdytem(i,1,2) = wdy(1,2,jdp(i))
            wdytem(i,2,2) = wdy(2,2,jdp(i))
            wdytem(i,3,2) = wdy(3,2,jdp(i))
            wdytem(i,4,2) = wdy(4,2,jdp(i))
         end do
C
C Loop over fields.
C
         do m = 1,pf
            do i = 1,pgls
               fyb(i,m) = wdytem(i,1,1)*fint(i,1,m) +
     $                    wdytem(i,2,1)*fint(i,2,m) +
     $                    wdytem(i,3,1)*fint(i,3,m) +
     $                    wdytem(i,4,1)*fint(i,4,m)
C
               fyt(i,m) = wdytem(i,1,2)*fint(i,1,m) +
     $                    wdytem(i,2,2)*fint(i,2,m) +
     $                    wdytem(i,3,2)*fint(i,3,m) +
     $                    wdytem(i,4,2)*fint(i,4,m)
            end do
         end do

         if ( nadv .gt. 0 ) then
            do m = 1,nadv
               do i = 1,pgls
                  xyb(i,m) = wdytem(i,1,1)*xint(i,1,m) +
     $                       wdytem(i,2,1)*xint(i,2,m) +
     $                       wdytem(i,3,1)*xint(i,3,m) +
     $                       wdytem(i,4,1)*xint(i,4,m)
C
                  xyt(i,m) = wdytem(i,1,2)*xint(i,1,m) +
     $                       wdytem(i,2,2)*xint(i,2,m) +
     $                       wdytem(i,3,2)*xint(i,3,m) +
     $                       wdytem(i,4,2)*xint(i,4,m)
               end do
            end do
         end if

      end if
C
      end subroutine cubydr

c##############################################################################

      subroutine cubzdr(pidim   ,pkdim   ,f       ,lbasdz  ,dfz1    ,
     $                  dfz2    )
C-----------------------------------------------------------------------
C
C Vertical derivative estimates for a vertical slice using Lagrangian
C cubic formulas.  Derivatives are set to zero at the top and bottom.
C At the "inner nodes" of the top and bottom intervals, a "one sided"
C estimate is used.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments 
C
      integer
     $     pidim,               ! longitudinal dimension
     $     pkdim                ! vertical     dimension
      real(r8) f(pidim,pkdim),      ! constituent field
     $     lbasdz(4,2,pkdim)    ! vertical interpolation weights
C
C Output arguments
C
      real(r8) dfz1(pidim,pkdim),   ! derivative at top of interval
     $     dfz2(pidim,pkdim)    ! derivative at bot of interval
C-----------------------------------------------------------------------
C
C  pidim   Horizontal dimension of arrays. (pidim = plon)
C  pkdim   Vertical dimension of arrays.
C  f       Vertical slice of data for which derivative estimates are
C          made
C  lbasdz  Lagrangian cubic basis functions for evaluating the
C          derivatives on the unequally spaced vertical grid.
C  dfz1    dfz1 contains derivative estimates at the "top" edges of the
C          intervals in the f array.
C  dfz2    dfz2 contains derivative estimates at the "bottom" edges of
C          the intervals in the f array.
C
C---------------------------Local variables-----------------------------
C
      integer i,k               ! indices
C
C-----------------------------------------------------------------------
C
      do k=2,pkdim-2
         do i=1,pidim
C
C Lagrangian derivative estimates (cubic) for the two center nodes in a
C four node stencil.
C
            dfz1(i,k) = lbasdz(1,1,k)*f(i,k-1) +
     $                  lbasdz(2,1,k)*f(i,k)   +
     $                  lbasdz(3,1,k)*f(i,k+1) +
     $                  lbasdz(4,1,k)*f(i,k+2)
C
            dfz2(i,k) = lbasdz(1,2,k)*f(i,k-1) +
     $                  lbasdz(2,2,k)*f(i,k)   +
     $                  lbasdz(3,2,k)*f(i,k+1) +
     $                  lbasdz(4,2,k)*f(i,k+2)
         end do
      end do
C
C Constrain derivatives to zero at top and bottom of vertical grid.
C At the interior nodes of the intervals at the top and bottom of the
C vertical grid, use the derivative estimate at that same node for the
C adjacent interval.  (This is a "one-sided" estimate for that node.)
C
      do i=1,pidim
         dfz1(i,1) = 0.0
         dfz2(i,1) = dfz1(i,2)
         dfz1(i,pkdim-1) = dfz2(i,pkdim-2)
         dfz2(i,pkdim-1) = 0.0
      end do
C
      end subroutine cubzdr

c##############################################################################

      subroutine extx(pkcnst,pkdim,fb)
C-----------------------------------------------------------------------
C
C Copy data to the longitude extensions of the extended array
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
      integer i2pi              ! start of eastern long. extension
      parameter( i2pi = nxpt+plon+1 ) 
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     pkcnst,              ! dimension construct for 3-D arrays
     $     pkdim                ! vertical dimension
C
C Output arguments
C
      real(r8) fb(plond,pkdim*pkcnst,platd) ! array of constituents
C
C---------------------------Local variables-----------------------------
C
      integer i,                ! longitude index
     $        j,                ! latitude  index
     $        k                 ! vertical  index
C
C-----------------------------------------------------------------------
C
C Fill west edge points.
C
      if(nxpt .ge. 1) then
         do j=1,platd
            do k=1,pkdim
               do i=1,nxpt
                  fb(i,k,j) = fb(i+plon,k,j)
               end do
            end do
         end do
      end if
C
C Fill east edge points
C
      do j=1,platd
         do k=1,pkdim
            do i=i2pi,plond
               fb(i,k,j) = fb(i-plon,k,j)
            end do
         end do
      end do
C
      end subroutine extx

c##############################################################################

      subroutine extys(pkcnst,pkdim,fb)
C-----------------------------------------------------------------------
C
C Fill latitude extensions of a scalar extended array.
C This is done in 2 steps:
C
C   1) interpolate to the pole points; use the mean field value on the
C      Gaussian latitude closest to the pole.
C   2) add latitude lines beyond the poles.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
      integer
     $     plon2,        ! plon/2
     $     istart,              ! index to start computation
     $     istop,               ! index to stop  computation
     $     js,                  ! index of southernmost model lat
     $     jn                   ! index of northernmost model lat
      parameter(plon2 = plon/2,
     $         istart = nxpt+1,
     $         istop  = nxpt+plon,
     $         js = 1    + nxpt + jintmx,
     $         jn = plat + nxpt + jintmx)
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     pkcnst,              ! dimensioning construct for 3-D arrays
     $     pkdim                ! vertical dimension
C
C Output arguments
C
      real(r8) fb(plond,pkdim*pkcnst,platd) ! constituent fields
C-----------------------------------------------------------------------
C
C  pkcnst  dimensioning construct for 3-D arrays
C  pkdim   Vertical dimension
C  fb      Output is same as on entry except with the pole latitude
C          and extensions beyond it filled.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                   ! index
     $     j,                   ! index
     $     k                    ! index
      real(r8) zave,                ! accumulator for zonal averaging
     $     flplon               ! float(plon)
C
C-----------------------------------------------------------------------
C
      flplon = real(plon,r8)
C
C Fill north pole line.
C
      do k=1,pkdim
         zave = 0.0
         do i=istart,istop
            zave = zave + fb(i,k,jn  )
         end do
         zave = zave/flplon
         do i=istart,istop
            fb(i,k,jn+1) = zave
         end do
      end do
C
C Fill northern lines beyond pole line.
C
      if( jn+2 .le. platd )then
         do j=jn+2,platd
            do k=1,pkdim
cdir$ ivdep
               do i=istart,istart+plon2-1
                  fb(      i,k,j) = fb(plon2+i,k,2*jn+2-j)
                  fb(plon2+i,k,j) = fb(      i,k,2*jn+2-j)
               end do
            end do
         end do
      end if
C
C Fill south pole line.
C
      do k=1,pkdim
         zave = 0.0
         do i = istart,istop
            zave = zave + fb(i,k,js  )
         end do
         zave = zave/flplon
         do i=istart,istop
            fb(i,k,js-1) = zave
         end do
      end do
C
C Fill southern lines beyond pole line.
C
      if( js-2 .ge. 1 )then
         do j=1,js-2
            do k=1,pkdim
cdir$ ivdep
               do i=istart,istart+plon2-1
                  fb(      i,k,j) = fb(plon2+i,k,2*js-2-j)
                  fb(plon2+i,k,j) = fb(      i,k,2*js-2-j)
               end do
            end do
         end do
      end if
C
      end subroutine extys

c##############################################################################

      subroutine extyv(pkcnst  ,pkdim   ,coslam  ,sinlam  ,vb      )
C-----------------------------------------------------------------------
C
C Fill latitude extensions of a vector component extended array.
C This is done in 2 steps:
C
C   1) interpolate to the pole points;
C      use coefficients for zonal wave number 1 on the Gaussian
C      latitude closest to the pole.
C   2) add latitude lines beyond the poles.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C-----------------------------------------------------------------------
C
      integer
     $     plon2,               ! plon/2
     $     istart,              ! index to start computation
     $     istop,               ! index to stop  computation
     $     js,                  ! index of southernmost model lat
     $     jn                   ! index of northernmost model lat
      parameter(plon2 = plon/2,
     $         istart = nxpt+1,
     $         istop  = nxpt+plon,
     $         js = 1    + nxpt + jintmx,
     $         jn = plat + nxpt + jintmx)
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     pkcnst,              ! dimensioning construct for 3-D arrays
     $     pkdim                ! vertical dimension
      real(r8) coslam(plon),        ! cos lamda
     $     sinlam(plon)         ! sin lamda
C
C Output arguments
C
      real(r8) vb(plond,pkdim*pkcnst,platd) ! one component of vector field
C-----------------------------------------------------------------------
C
C  pkcnst  dimensioning construct for 3-D arrays
C  pkdim   Vertical dimension
C  coslam  Cos of longitude at x-grid points (global grid).
C  sinlam  Sin of longitude at x-grid points (global grid).
C  vb      Output is same as on entry except with the pole latitude
C          and extensions beyond it filled.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                   ! index
     $     ig,                  ! index
     $     j,                   ! index
     $     k                    ! index
      real(r8) zavec,               ! accumulator for wavenumber 1
     $     zaves,               ! accumulator for wavenumber 1
     $     fplon2               ! float(plon2)
C
C-----------------------------------------------------------------------
C
      fplon2 = real(plon2,r8)
C
C Fill north pole line.
C
      do k = 1,pkdim
         zavec = 0.0
         zaves = 0.0
         ig = 0
         do i = istart,istop
            ig    = ig + 1
            zavec = zavec + vb(i,k,jn  )*coslam(ig)
            zaves = zaves + vb(i,k,jn  )*sinlam(ig)
         end do
         zavec = zavec/fplon2
         zaves = zaves/fplon2
         ig    = 0
         do i = istart,istop
            ig           = ig + 1
            vb(i,k,jn+1) = zavec*coslam(ig) + zaves*sinlam(ig)
         end do
      end do
C
C Fill northern lines beyond pole line.
C
      if( jn+2 .le. platd )then
         do j = jn+2,platd
            do k = 1,pkdim
cdir$ ivdep
               do i = istart,istart+plon2-1
                  vb(      i,k,j) = -vb(plon2+i,k,2*jn+2-j)
                  vb(plon2+i,k,j) = -vb(      i,k,2*jn+2-j)
               end do
            end do
         end do
      end if
C 
C Fill south pole line.
C
      do k = 1,pkdim
         zavec = 0.0
         zaves = 0.0
         ig = 0
         do i = istart,istop
            ig    = ig + 1
            zavec = zavec + vb(i,k,js  )*coslam(ig)
            zaves = zaves + vb(i,k,js  )*sinlam(ig)
         end do
         zavec = zavec/fplon2
         zaves = zaves/fplon2
         ig    = 0
         do i = istart,istop
            ig           = ig + 1
            vb(i,k,js-1) = zavec*coslam(ig) + zaves*sinlam(ig)
         end do
      end do
C
C Fill southern lines beyond pole line.
C
      if( js-2 .ge. 1 )then
         do j = 1,js-2
            do k = 1,pkdim
cdir$ ivdep
               do i = istart,istart+plon2-1
                  vb(      i,k,j) = -vb(plon2+i,k,2*js-2-j)
                  vb(plon2+i,k,j) = -vb(      i,k,2*js-2-j)
               end do
            end do
         end do
      end if
C
      end subroutine extyv

c##############################################################################

      subroutine g2spos(dttmp   ,lam     ,phib    ,phi     ,cosphi  ,
     $                  sinphi  ,upr     ,vpr     ,lamgc   ,phigc   ,
     $                  lamsc   ,phisc   )
C-----------------------------------------------------------------------
C
C Transform position coordinates for a set of points, each of which is
C associated with a grid point in a global latitude slice, from local
C geodesic to spherical coordinates.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C-----------------------------------------------------------------------
      integer pgls
      parameter(pgls = plon*plev) ! length of latitude slice
      real(r8) fac                    ! 1. - eps, eps from machine precision
#if ( REALSIZE == 4 )
      parameter ( fac = 1. - 1.e-06 )
#else
      parameter ( fac = 1. - 1.e-12 )
#endif
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) dttmp,               ! time step
     $     lam(plon),           ! model longitude coordinates
     $     phib(platd),         ! extended grid latitude coordinates
     $     phi,                 ! current latitude coordinate (radians)
     $     cosphi,              ! cos of current latitude
     $     sinphi,              ! sin of current latitude
     $     upr  (pgls),         ! u-wind in geodesic coord
     $     vpr  (pgls),         ! v-wind in geodesic coord
     $     lamgc(pgls),         ! geodesic long coord. of dep. point
     $     phigc(pgls)          ! geodesic lat  coord. of dep. point
C
C Output arguments
C
      real(r8) lamsc(plon,plev),    ! spherical long coord. of dep. point
     $     phisc(pgls)          ! spherical lat  coord. of dep. point
C-----------------------------------------------------------------------
C
C  dttmp  Time step over which midpoint/endpoint trajectory is
C         calculated (seconds).
C  lam    Longitude coordinates of the global grid points in spherical
C         system.  The grid points in the global array are the reference
C         points for the local geodesic systems.
C  phib   Latitude values for the extended grid.
C  phi    Latitude coordinate (in the global grid) of the current
C         latitude slice.
C  cosphi cos( phi )
C  sinphi sin( phi )
C  upr    zonal      velocity at departure point in local geodesic coord
C  vpr    Meridional velocity at departure point in local geodesic coord
C  lamgc  Longitude coordinate of points in geodesic coordinates.
C  phigc  Latitude coordinate of points in geodesic coordinates.
C  lamsc  Longitude coordinate of points in spherical coordinates.
C  phisc  Latitude coordinate of points in spherical coordinates.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                   ! index
     $     ii,                  ! index
     $     k,                   ! index
     $     nval,                ! number of values returned from whenfgt
     $     indx(pgls)           ! index holder
      real(r8) pi,                  ! 4.*atan(1.)
     $     twopi,               ! 2.*pi
     $     pi2,                 ! pi/2
     $     sgnphi,              ! holds sign of phi
     $     sphigc,              ! sin(phigc)
     $     cphigc,              ! cos(phigc)
     $     clamgc,              ! cos(lamgc)
     $     slam2,               ! sin(lamgc)**2
     $     phipi2,              ! tmp variable
     $     slamgc(pgls),        ! sin(lamgc)
     $     dlam(pgls),          ! zonal extent of trajectory
     $     coeff,               ! tmp variable
     $     distmx,              ! max distance
     $     dist(pgls)           ! approx. distance traveled along traj.
C
C-----------------------------------------------------------------------
C
      pi     = 4.*atan(1._r8)
      twopi  = pi*2.
      pi2    = pi/2.
      coeff  = (1.1_r8*dttmp)**2
      distmx = (sign(pi2,phi) - phi)**2/coeff
      sgnphi = sign( 1._r8, phi )
C
      do i=1,pgls
         sphigc    = sin( phigc(i) )
         cphigc    = cos( phigc(i) )
         slamgc(i) = sin( lamgc(i) )
         clamgc    = cos( lamgc(i) )
         phisc(i) = asin((sphigc*cosphi + cphigc*sinphi*clamgc)*fac)
         if ( abs(phisc(i)) .ge. phib(j1+plat)*fac )
     $        phisc(i) = sign( phib(j1+plat),phisc(i) )*fac
         dlam(i) = asin((slamgc(i)*cphigc/cos(phisc(i)))*fac)
C
C Compute estimated trajectory distance based upon winds alone
C
         dist(i) = upr(i)**2 + vpr(i)**2
      end do
C
C Determine which trajectories may have crossed over pole
C
      call whenfgt(pgls,dist,1,distmx,indx,nval)
C
C Check that proper branch of arcsine is used for calculation of
C dlam for those trajectories which may have crossed over pole.
C
CDIR$ IVDEP
      do ii=1,nval
         i = indx(ii)
         slam2 = slamgc(i)**2
         phipi2 = asin((sqrt((slam2 - 1._r8)/(slam2 - 1._r8/cosphi**2)))*fac)
         if (sgnphi*phigc(i) .gt. phipi2) then
            dlam(i) = sign(pi,lamgc(i)) - dlam(i)
         end if
      end do
C
      do k=1,plev
         do i=1,plon
            lamsc(i,k) = lam(i) + dlam((k-1)*plon + i)
C
C Restrict longitude to be in the range [0, twopi).
C
            if( lamsc(i,k) .ge. twopi ) lamsc(i,k) = lamsc(i,k) - twopi
            if( lamsc(i,k) .lt. 0.0 )   lamsc(i,k) = lamsc(i,k) + twopi
         end do
      end do

      end subroutine g2spos

c##############################################################################

      subroutine herxin(pf      ,pkcnst  ,fb      ,fxl     ,fxr     ,
     $                  x       ,xdp     ,idp     ,jdp     ,fint
     $,                 nadv, x1, xint )
C-----------------------------------------------------------------------
C
C For each departure point in the latitude slice being forecasted,
C interpolate (using equally spaced Hermite cubic formulas) to its
C x value at each latitude required for later interpolation in the y
C direction.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      integer
     $     pf,                           ! dimension (number of fields)
     $     pkcnst                        ! dimension (number of fields)
     $,    nadv                          ! number of misc advected quantities
      real(r8) fb (plond,plev,pkcnst,platd), ! input constituents
     $     fxl(plond,plev,pf    ,platd), ! left  x-derivatives 
     $     fxr(plond,plev,pf    ,platd), ! right x-derivatives
     $     x(plond),                     ! longitudinal grid coordinates
     $     xdp(plon,plev)                ! departure point coordinates
      real(r8), dimension(:,:,:,:) ::
     $     x1                            ! misc. advected quantities (plond,plev,padv,platd)
      integer
     $     idp(plon,plev),               ! longitude index of dep pt.
     $     jdp(plon,plev)                ! latitude  index of dep pt.
C
C Output arguments
C
      real(r8) fint(plon,plev,ppdy,pf)       ! x-interpolants
      real(r8), dimension(:,:,:,:) ::
     $  xint                             ! x interpolants (plon,plev,ppdy,padv)
C
C-----------------------------------------------------------------------
C
C  pf      Number of fields being interpolated.
C  pkcnst  Dimensioning construct for 3-D arrays. (see common/com3d/)
C  fb      extended array of data to be interpolated.
C  fxl     x derivatives at the left edge of each interval containing 
C          the departure point
C  fxr     x derivatives at the right edge of each interval containing 
C          the departure point
C  x       Equally spaced x grid values in extended arrays.
C  xdp     xdp(i,k) is the x-coordinate (extended grid) of the
C          departure point that corresponds to global grid point (i,k)
C          in the latitude slice being forecasted.
C  idp     idp(i,k) is the index of the x-interval (extended grid) that
C          contains the departure point corresponding to global grid
C          point (i,k) in the latitude slice being forecasted.
C          Note that
C                x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
C  jdp     jdp(i,k) is the index of the y-interval (extended grid) that
C          contains the departure point corresponding to global grid
C          point (i,k) in the latitude slice being forecasted.
C          Suppose yb contains the y-coordinates of the extended array
C          and ydp(i,k) is the y-coordinate of the departure point
C          corresponding to grid point (i,k).  Then,
C                yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
C  fint    (fint(i,k,j,n),j=1,ppdy) contains the x interpolants at each
C          latitude needed for the y derivative estimates at the
C          endpoints of the interval that contains the departure point
C          for grid point (i,k).  The last index of fint allows for
C          interpolation of multiple fields.
C
C---------------------------Local workspace-----------------------------
C
      integer
     $     i,k,m                         ! indices
      real(r8) dx,                           ! x-increment
     $     rdx,                          ! 1./dx
     $     xl (plon,plev),               ! |
     $     xr (plon,plev),               ! |
     $     hl (plon,plev),               ! | --interpolation coefficients
     $     hr (plon,plev),               ! |
     $     dhl(plon,plev),               ! |
     $     dhr(plon,plev)                ! |
C
C-----------------------------------------------------------------------
C
      if(ppdy .ne. 4) then
         write(6,*) 'HERXIN:Fatal error: ppdy must be set to 4'
         call endrun
      end if
C
      dx = x(nxpt+2) - x(nxpt+1)
      rdx = 1._r8/dx
      do k=1,plev
         do i=1,plon
            xl (i,k) = ( x(idp(i,k)+1) - xdp(i,k) )*rdx
            xr (i,k) = 1._r8 - xl(i,k)
            hl (i,k) = ( 3._r8 - 2.*xl(i,k) )*xl(i,k)**2
            hr (i,k) = ( 3._r8 - 2.*xr(i,k) )*xr(i,k)**2
            dhl(i,k) = -dx*( xl(i,k) - 1._r8 )*xl(i,k)**2
            dhr(i,k) =  dx*( xr(i,k) - 1._r8 )*xr(i,k)**2
         end do
      end do
C
C x interpolation at each latitude needed for y interpolation.
C Once for each field.
C 
      do m = 1,pf
         do k = 1,plev
            do i = 1,plon
               fint(i,k,1,m) =
     $             fb (idp(i,k)  ,k,m,jdp(i,k)-1)*hl (i,k) +
     $             fb (idp(i,k)+1,k,m,jdp(i,k)-1)*hr (i,k) +
     $             fxl(idp(i,k)  ,k,m,jdp(i,k)-1)*dhl(i,k) +
     $             fxr(idp(i,k)  ,k,m,jdp(i,k)-1)*dhr(i,k)
               fint(i,k,2,m) =
     $             fb (idp(i,k)  ,k,m,jdp(i,k)  )*hl (i,k) +
     $             fb (idp(i,k)+1,k,m,jdp(i,k)  )*hr (i,k) +
     $             fxl(idp(i,k)  ,k,m,jdp(i,k)  )*dhl(i,k) +
     $             fxr(idp(i,k)  ,k,m,jdp(i,k)  )*dhr(i,k)
               fint(i,k,3,m) =
     $             fb (idp(i,k)  ,k,m,jdp(i,k)+1)*hl (i,k) +
     $             fb (idp(i,k)+1,k,m,jdp(i,k)+1)*hr (i,k) +
     $             fxl(idp(i,k)  ,k,m,jdp(i,k)+1)*dhl(i,k) +
     $             fxr(idp(i,k)  ,k,m,jdp(i,k)+1)*dhr(i,k)
               fint(i,k,4,m) =
     $             fb (idp(i,k)  ,k,m,jdp(i,k)+2)*hl (i,k) +
     $             fb (idp(i,k)+1,k,m,jdp(i,k)+2)*hr (i,k) +
     $             fxl(idp(i,k)  ,k,m,jdp(i,k)+2)*dhl(i,k) +
     $             fxr(idp(i,k)  ,k,m,jdp(i,k)+2)*dhr(i,k)
            end do
         end do
      end do

      if ( nadv .gt. 0 ) then
         do m = 1, nadv
            do k = 1,plev
               do i = 1,plon
                  xint(i,k,1,m) =
     $                 x1 (idp(i,k)  ,k,m,jdp(i,k)-1)*hl (i,k) +
     $                 x1 (idp(i,k)+1,k,m,jdp(i,k)-1)*hr (i,k) +
     $                 xxl(idp(i,k)  ,k,m,jdp(i,k)-1)*dhl(i,k) +
     $                 xxr(idp(i,k)  ,k,m,jdp(i,k)-1)*dhr(i,k)
                  xint(i,k,2,m) =
     $                 x1 (idp(i,k)  ,k,m,jdp(i,k)  )*hl (i,k) +
     $                 x1 (idp(i,k)+1,k,m,jdp(i,k)  )*hr (i,k) +
     $                 xxl(idp(i,k)  ,k,m,jdp(i,k)  )*dhl(i,k) +
     $                 xxr(idp(i,k)  ,k,m,jdp(i,k)  )*dhr(i,k)
                  xint(i,k,3,m) =
     $                 x1 (idp(i,k)  ,k,m,jdp(i,k)+1)*hl (i,k) +
     $                 x1 (idp(i,k)+1,k,m,jdp(i,k)+1)*hr (i,k) +
     $                 xxl(idp(i,k)  ,k,m,jdp(i,k)+1)*dhl(i,k) +
     $                 xxr(idp(i,k)  ,k,m,jdp(i,k)+1)*dhr(i,k)
                  xint(i,k,4,m) =
     $                 x1 (idp(i,k)  ,k,m,jdp(i,k)+2)*hl (i,k) +
     $                 x1 (idp(i,k)+1,k,m,jdp(i,k)+2)*hr (i,k) +
     $                 xxl(idp(i,k)  ,k,m,jdp(i,k)+2)*dhl(i,k) +
     $                 xxr(idp(i,k)  ,k,m,jdp(i,k)+2)*dhr(i,k)
               end do
            end do
         end do
      end if

      end subroutine herxin

c##############################################################################

      subroutine heryin(pf      ,fint    ,fyb     ,fyt     ,y       ,
     $                  dy      ,ydp     ,jdp     ,fdp
     $,                 nadv, xint, xyb, xyt, xhr )
C-----------------------------------------------------------------------
C
C For each departure point in the latitude slice to be forecasted,
C interpolate (using unequally spaced Hermite cubic formulas) the
C x interpolants to the y value of the departure point.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none

      integer pgls
      parameter( pgls = plon*plev ) ! length of latitude slice
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     pf                   ! dimension (number of fields)
     $,    nadv                 ! number of misc advected quantities
      real(r8) fint(pgls,ppdy,pf),  ! x-interpolants
     $     fyb (pgls,pf),       ! y-derivatives at bottom of interval
     $     fyt (pgls,pf),       ! y-derivatives at top    of interval
     $     xint(pgls,ppdy,padv),! x interpolants
     $     xyb (pgls,padv),     ! y-derivatives at bot of int.
     $     xyt (pgls,padv),     ! y-derivatives at top of int.
     $     y   (platd),         ! latitude grid coordinates
     $     dy  (platd),         ! intervals between latitude grid pts.
     $     ydp (pgls)           ! lat. coord of departure point.
      integer
     $     jdp (pgls)           ! lat. index of departure point.
C
C Output arguments
C
      real(r8) fdp (pgls,pf)        ! y-interpolants
     $,    xhr (pgls,padv)      ! horizontal interpolants
C
C-----------------------------------------------------------------------
C
C  pf      Number of fields being interpolated.
C  fint    (fint(i,k,j,m),j=ppdy/2,ppdy/2 + 1) contains the x
C          interpolants at the endpoints of the y-interval that
C          contains the departure point for grid point (i,k).  The last
C          index of fint allows for interpolation of multiple fields.
C          fint is generated by a call to herxin.
C  fyb     fyb(i,k,.) is the derivative at the "bottom" of the
C          y-interval that contains the departure point of grid
C          point (i,k).  fyb is generated by a call to cubydr.
C  fyt     fyt(i,k,.) is the derivative at the "top" of the y-interval
C          that contains the departure point of grid point (i,k).
C          fyt is generated by a call to cubydr.
C  y       y-coordinate (latitude) values in the extended array.
C  dy      Increment in the y-coordinate value for each interval in the
C          extended array.
C  ydp     ydp(i,k) is the y-coordinate of the departure point that
C          corresponds to global grid point (i,k) in the latitude slice
C          being forecasted.
C  jdp     jdp(i,k) is the index of the y-interval that contains the
C          departure point corresponding to global grid point (i,k) in
C          the latitude slice being forecasted.
C          Note that
C                y(jdp(i,k)) .le. ydp(i,k) .lt. y(jdp(i,k)+1) .
C  fdp     Horizontally interpolated field values at the departure point
C          for the latitude slice being forecasted.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                   ! index
     $     jb,                  ! index corresponding to bot of interval
     $     jt,                  ! index corresponding to top of interval
     $     m                    ! index
      real(r8) dyj(pgls),           ! latitude interval containing dep. pt.
     $     yb (pgls),           ! |
     $     yt (pgls),           ! |
     $     hb (pgls),           ! | -- interpolation coefficients
     $     ht (pgls),           ! |
     $     dhb(pgls),           ! |
     $     dht(pgls)            ! |
C
      jb = ppdy/2
      jt = jb + 1
C
      do i = 1,pgls
         dyj(i) = dy(jdp(i))
         yb (i) = ( y(jdp(i)+1) - ydp(i) )/dyj(i)
         yt (i) = 1._r8 - yb(i)
         hb (i) = ( 3._r8 - 2.*yb(i) )*yb(i)**2
         ht (i) = ( 3._r8 - 2.*yt(i) )*yt(i)**2
         dhb(i) = -dyj(i)*( yb(i) - 1._r8 )*yb(i)**2
         dht(i) =  dyj(i)*( yt(i) - 1._r8 )*yt(i)**2
      end do
C
C Loop over fields.
C
      do m = 1,pf
         do i = 1,pgls
            fdp(i,m) = fint(i,jb,m)*hb(i) + fyb(i,m)*dhb(i) +
     $                 fint(i,jt,m)*ht(i) + fyt(i,m)*dht(i)
         end do
      end do

      if ( nadv .gt. 0 ) then
         do m = 1, nadv
            do i = 1,pgls
               xhr(i,m) = xint(i,jb,m)*hb(i) + xyb(i,m)*dhb(i) +
     $                 xint(i,jt,m)*ht(i) + xyt(i,m)*dht(i)
            end do
         end do
      end if
C
      end subroutine heryin

c##############################################################################

      subroutine herzin(pkdim   ,pf      ,f       ,fst     ,fsb     ,
     $                  sig     ,dsig    ,sigdp   ,kdp     ,fdp
     $,                 nadv, xhr, xhst, xhsb, x2 )
C-----------------------------------------------------------------------
C
C Interpolate field on vertical slice to vertical departure point using
C Hermite cubic interpolation.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      integer
     $     pkdim,                ! vertical dimension
     $     pf                    ! dimension (number of fields)
     $,    nadv                  ! number of misc advected quantities
      real(r8) f    (plon,pkdim,pf), ! fields
     $     fst  (plon,pkdim,pf), ! z-derivatives at top edge of interval
     $     fsb  (plon,pkdim,pf), ! z-derivatives at bot edge of interval
     $     sig  (pkdim),         ! vertical grid coordinates
     $     dsig (pkdim),         ! intervals between vertical grid pts.
     $     sigdp(plon,plev)      ! vertical coord. of departure point
      real(r8), dimension(:,:,:) ::
     $  xhr   ! fields (plon,pkdim,padv)
     $, xhst  ! z-derivatives at top edge of interval (plon,pkdim,padv)
     $, xhsb  ! z-derivatives at bot edge of interval (plon,pkdim,padv)
      integer
     $     kdp  (plon,plev)      ! vertical index  of departure point
C
C Output arguments
C
      real(r8) fdp(plond,plev,pf)    ! z-interpolants
      real(r8), dimension(:,:,:) ::
     $  x2                       ! z-interpolants (plond,plev,padv)
C
C-----------------------------------------------------------------------
C
C  pkdim   Vertical dimension of vertical slice arrays.
C  pf      Number of fields being interpolated.
C  f       Vertical slice of data to be interpolated.
C  fst     z-derivatives at the top edge of each interval contained in f
C  fsb     z-derivatives at the bot edge of each interval contained in f
C  sig     Sigma values corresponding to the vertical grid
C  dsig    Increment in sigma value for each interval in vertical grid.
C  sigdp   Sigma value at the trajectory midpoint or endpoint for each
C          gridpoint in a vertical slice from the global grid.
C  kdp     Vertical index for each gridpoint.  This index points into a
C          vertical slice array whose vertical grid is given by sig.
C          E.g.,   sig(kdp(i,j)) .le. sigdp(i,j) .lt. sig(kdp(i,j)+1) .
C  fdp     Value of field at the trajectory midpoints or endpoints.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,k,m                 ! | -- indices
      real(r8) dzk,                  ! vert interval containing the dep. pt.
     $     zt,                   ! |
     $     zb,                   ! |
     $     ht (plon),            ! | -- interpolation coefficients
     $     hb (plon),            ! |
     $     dht(plon),            ! |
     $     dhb(plon)             ! |
C

      do k=1,plev
         do i=1,plon
            dzk = dsig(kdp(i,k))
            zt = ( sig(kdp(i,k)+1) - sigdp(i,k) )/dzk
            zb = 1._r8 - zt
            ht (i) = ( 3._r8 - 2._r8*zt )*zt**2
            hb (i) = ( 3._r8 - 2._r8*zb )*zb**2
            dht(i) = -dzk*( zt - 1._r8 )*zt**2
            dhb(i) =  dzk*( zb - 1._r8 )*zb**2
         end do
C
C Loop over fields.
C
         do m=1,pf
            do i=1,plon
               fdp(i,k,m) = f(i,kdp(i,k)  ,m)* ht(i) +
     $                      fst(i,kdp(i,k),m)*dht(i) +
     $                      f(i,kdp(i,k)+1,m)* hb(i) +
     $                      fsb(i,kdp(i,k),m)*dhb(i)
            end do
         end do

         if ( nadv .gt. 0 ) then
            do m=1,nadv
               do i=1,plon
                  x2(i,k,m) = xhr(i,kdp(i,k)  ,m)* ht(i) +
     $                        xhst(i,kdp(i,k),m)*dht(i) +
     $                        xhr(i,kdp(i,k)+1,m)* hb(i) +
     $                        xhsb(i,kdp(i,k),m)*dhb(i)
               end do
            end do
         end if

      end do
C
      end subroutine herzin

c##############################################################################

      subroutine hrintp(pf      ,pkcnst  ,fb      ,fxl     ,fxr     ,
     $                  x       ,y       ,dy      ,wdy     ,xdp     ,
     $                  ydp     ,idp     ,jdp     ,limitd  ,fint    ,
     $                  fyb     ,fyt     ,fdp
     $,                 nadv, x1, xint, xyb, xyt, xhr )
C-----------------------------------------------------------------------
C
C Interpolate 2-d field to departure point using tensor product
C Hermite cubic interpolation.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      integer
     $     pf,                           ! dimension (number of fields)
     $     pkcnst                        ! dimension (see ext. document)
     $,    nadv                          ! number of misc advected quantities
      real(r8) fb (plond,plev,pkcnst,platd), ! input fields
     $     fxl(plond,plev,pf    ,platd), ! x-derivs at left  of interval
     $     fxr(plond,plev,pf    ,platd), ! x-derivs at right of interval
     $     x1(plond,plev,padv,platd),    ! misc. advected quantities
     $     x  (plond),                   ! long. grid coordinates
     $     y  (platd),                   ! lat.  grid coordinates
     $     dy (platd),                   ! intervals betwn lat grid pts.
     $     wdy(4,2,platd),               ! lat. derivative weights
     $     xdp(plon,plev),               ! x-coord of dep. pt.
     $     ydp(plon,plev)                ! y-coord of dep. pt.
      integer
     $     idp(plon,plev),               ! i index of dep. pt.
     $     jdp(plon,plev)                ! j index of dep. pt.
      logical limitd                     ! flag for shape-preservation
C
C Output arguments
C
      real(r8) fint(plon,plev,ppdy,pf),      ! x interpolants
     $     fyb (plon,plev,pf),           ! y-derivatives at bot of int.
     $     fyt (plon,plev,pf),           ! y-derivatives at top of int.
     $     fdp (plon,plev,pf)            ! horizontal interpolants
     $,    xint(plon,plev,ppdy,padv),    ! x interpolants
     $     xyb (plon,plev,padv),         ! y-derivatives at bot of int.
     $     xyt (plon,plev,padv),         ! y-derivatives at top of int.
     $     xhr (plon,plev,padv)          ! horizontal interpolants
C
C-----------------------------------------------------------------------
C
C  pf      Number of fields being interpolated.
C  pkcnst  dimensioning construct for 3-D arrays. (see ext. document)
C  fb      Extended array of data to be interpolated.
C  fxl     x-derivatives at the left  edge of each interval containing
C          the departure point.
C  fxr     x-derivatives at the right edge of each interval containing
C          the departure point.
C  x       Equally spaced x grid values in extended arrays.
C  y       y-coordinate (latitude) values in the extended array.
C  dy      Increment in the y-coordinate value for each interval in the
C          extended array.
C  wdy     Weights for Lagrange cubic derivative estimates on the
C          unequally spaced y-grid.  If grid interval j (in extended
C          array is surrounded by a 4 point stencil, then the
C          derivative at the "bottom" of the interval uses the weights
C          wdy(1,1,j),wdy(2,1,j), wdy(3,1,j), and wdy(4,1,j).  The
C          derivative at the "top" of the interval uses wdy(1,2,j),
C          wdy(2,2,j), wdy(3,2,j) and wdy(4,2,j).
C  xdp     xdp(i,k) is the x-coordinate of the departure point that
C          corresponds to global grid point (i,k) in the latitude slice
C          being forecasted.
C  ydp     ydp(i,k) is the y-coordinate of the departure point that
C          corresponds to global grid point (i,k) in the latitude slice
C          being forecasted.
C  idp     idp(i,k) is the index of the x-interval that contains the
C          departure point corresponding to global grid point (i,k) in
C          the latitude slice being forecasted.
C          Note that
C                x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
C  jdp     jdp(i,k) is the index of the y-interval that contains the
C          departure point corresponding to global grid point (i,k) in
C          the latitude slice being forecasted.
C          Suppose yb contains the y-coordinates of the extended array
C          and ydp(i,k) is the y-coordinate of the departure point
C          corresponding to grid point (i,k).  Then,
C                yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
C  limitd  Logical flag to specify whether or not the y-derivatives will
C          be limited.
C  fint    WORK ARRAY, results not used on return
C  fyb     WORK ARRAY, results not used on return
C  fyt     WORK ARRAY, results not used on return
C  fdp     Value of field at the horizontal departure points.
C
C-----------------------------------------------------------------------
C
C Hermite cubic interpolation to the x-coordinate of each
C departure point at each y-coordinate required to compute the
C y-derivatives.
C
      call herxin(pf      ,pkcnst  ,fb      ,fxl     ,fxr     ,
     $            x       ,xdp     ,idp     ,jdp     ,fint
     $,           nadv, x1, xint )
C
C Compute y-derivatives.
C
      call cubydr(pf      ,fint    ,wdy     ,jdp     ,fyb     ,
     $            fyt
     $,           nadv, xint, xyb, xyt )

      if( limitd )then
         call limdy(pf      ,fint    ,dy      ,jdp     ,fyb     ,
     $              fyt     )
         if ( nadv .gt. 0 ) then
            call limdy( nadv, xint, dy, jdp, xyb, xyt )
         end if
      end if
C
C Hermite cubic interpolation in the y-coordinate.
C
      call heryin(pf      ,fint    ,fyb     ,fyt     ,y       ,
     $            dy      ,ydp     ,jdp     ,fdp
     $,           nadv, xint, xyb, xyt, xhr )
C
      end subroutine hrintp

c##############################################################################

      subroutine hrvint(pf      ,pkcnst  ,fb      ,fxl     ,fxr     ,
     $                  x       ,wdy     ,xdp     ,ydp     ,idp     ,
     $                  jdp     ,fint    ,fdp     )
C-----------------------------------------------------------------------
C
C Interpolate 2-d field to departure point using equivalent of tensor
C product Lagrange cubic interpolation. For economy, code actually
C uses Hermite cubic with Lagrange cubic derivative estimate in x
C and Lagrange cubic in y.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      integer
     $     pf,                           ! dimension (number of fields)
     $     pkcnst                        ! dimension (see ext. document)
      real(r8) fb (plond,plev,pkcnst,platd), ! input fields
     $     fxl(plond,plev,pf    ,platd), ! x-derivs at left  of int.
     $     fxr(plond,plev,pf    ,platd), ! x-derivs at right of int.
     $     x(plond),                     ! x-coordinates of grid
     $     wdy(4,2,platd),               ! latitude interp weights
     $     xdp(plon,plev),               ! x-coord of dep. pt.
     $     ydp(plon,plev)                ! y-coord of dep. pt.
      integer
     $     idp(plon,plev),               ! x-index of dep. pt.
     $     jdp(plon,plev)                ! y-index of dep. pt.
C
C Output arguments
C
      real(r8) fint(plon,plev,ppdy,pf),      ! x-interpolants
     $     fdp(plon,plev,pf)             ! final horizontal interpolants
C
C-----------------------------------------------------------------------
C
C  pf      Number of fields being interpolated.
C  pkcnst  Dimensioning construct for 3-D arrays. (see ext. document)
C  fb      Extended array of data to be interpolated.
C  fxl     x-derivatives at the left  edge of each interval containing
C          the departure point.
C  fxr     x-derivatives at the right edge of each interval containing
C          the departure point.
C  x       Equally spaced x grid values in extended arrays.
C  wdy     Weights for Lagrange cubic interpolation on the
C          unequally spaced y-grid.
C  xdp     xdp(i,k) is the x-coordinate of the departure point that
C          corresponds to global grid point (i,k) in the latitude slice
C          being forecasted.
C  ydp     ydp(i,k) is the y-coordinate of the departure point that
C          corresponds to global grid point (i,k) in the latitude slice
C          being forecasted.
C  idp     idp(i,k) is the index of the x-interval that contains the
C          departure point corresponding to global grid point (i,k) in
C          the latitude slice being forecasted.
C          Note that
C                x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
C  jdp     jdp(i,k) is the index of the y-interval that contains the
C          departure point corresponding to global grid point (i,k) in
C          the latitude slice being forecasted.
C          Suppose yb contains the y-coordinates of the extended array
C          and ydp(i,k) is the y-coordinate of the departure point
C          corresponding to grid point (i,k).  Then,
C                yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
C  fint    WORK ARRAY, results not used on return
C  fdp     Value of field at the departure points.
C

c Local variables:

      integer nadv
      real(r8) x1dum(1,1,1,1), xintdum(1,1,1,1)  ! dummies

C-----------------------------------------------------------------------
C
C Hermite cubic interpolation to the x-coordinate of each
C departure point at each y-coordinate required to compute the
C y-interpolants.
C
      nadv = 0

      call herxin(pf      ,pkcnst  ,fb      ,fxl     ,fxr     ,
     $            x       ,xdp     ,idp     ,jdp     ,fint
     $,           nadv, x1dum, xintdum )
C
C Lagrange cubic interpolation in y
C
      call lagyin(pf      ,fint    ,wdy     ,ydp     ,jdp     ,
     $            fdp     )
C
      end subroutine hrvint

c##############################################################################

      subroutine kdpfnd(pkdim   ,sig     ,sigdp   ,kdpmap  ,
     $                  kdp     )
C-----------------------------------------------------------------------
C
C Determine vertical departure point indices that point into a grid
C containing the full or half sigma levels.  Use an artificial evenly 
C spaced vertical grid to map into the true model levels.
C
C Note:
C Indices are computed assuming the the sigdp values have
C been constrained so that sig(1) .le. sigdp(i,j) .lt. sig(pkdim).
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
      integer pgls
      parameter ( pgls = plon*plev ) ! length of latitude slice
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     pkdim                ! dimension of "sig"
      real(r8) sig  (pkdim),        ! vertical grid coordinates
     $     sigdp(pgls)          ! vertical coords. of departure points
      integer
     $     kdpmap(pmap)         ! array of model grid indices which
C                               ! are mapped into the artificial grid.
C
C Output arguments
C
      integer kdp(pgls)         ! vertical index for each dep. pt.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,ii                 ! indices
      real(r8) rdel,                ! reciprocal of interval in artificial grid
     $     sig1ln               ! ln (sig(1))
C
C-----------------------------------------------------------------------
C
      rdel   = real(pmap,r8)/( log(sig(pkdim)) - log(sig(1)) )
      sig1ln = log( sig(1) )
C
      do i = 1,pgls
C 
C First guess of the departure point's location in the model grid
C
         ii = max0(1,min0(pmap,int((log(sigdp(i))-sig1ln)*rdel + 1._r8)))
         kdp(i) = kdpmap(ii)
C
C Determine if location is in next interval
C
         if(sigdp(i) .ge. sig(kdp(i)+1)) then
            kdp(i) = kdp(i) + 1
         end if
      end do
C
      end subroutine kdpfnd

c##############################################################################

      subroutine lagyin(pf      ,fint    ,wdy     ,ydp     ,jdp     ,
     $                  fdp     )
C-----------------------------------------------------------------------
C
C For each departure point in the latitude slice to be forecasted,
C interpolate (using unequally spaced Lagrange cubic formulas) the
C x interpolants to the y value of the departure point.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C-----------------------------------------------------------------------
      integer pgls              ! length of latitude slice
      parameter( pgls = plon*plev )
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     pf                   ! dimension (number of fields)
      real(r8) fint(pgls,ppdy,pf),  ! x-interpolants
     $     wdy(4,2,platd),      ! y-interpolation weights
     $     ydp(pgls)            ! y-coordinates of departure pts.
      integer
     $     jdp(pgls)            ! j-index of departure point coord.
C
C Output arguments
C
      real(r8) fdp(pgls,pf)         ! interpolants at the horiz. depart. pt.
C
C-----------------------------------------------------------------------
C
C  pf      Number of fields being interpolated.
C  fint    (fint(i,k,j,m),j=ppdy/2,ppdy/2 + 1) contains the x
C          interpolants at the endpoints of the y-interval that contains
C          the departure point for grid point (i,k).  The last index of
C          fint allows for interpolation of multiple fields.  fint is
C          generated by a call to herxin.
C  wdy     Grid values and weights for Lagrange cubic interpolation on
C          the unequally spaced y-grid.
C  ydp     ydp(i,k) is the y-coordinate of the departure point that
C          corresponds to global grid point (i,k) in the latitude slice
C          being forecasted.
C  jdp     jdp(i,k) is the index of the y-interval that contains the
C          departure point corresponding to global grid point (i,k) in
C          the latitude slice being forecasted.
C          Note that
C                y(jdp(i,k)) .le. ydp(i,k) .lt. y(jdp(i,k)+1) .
C  fdp     Horizontally interpolated field values at the departure point
C          for the latitude slice being forecasted.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,m                  ! indices
      real(r8) ymy1,                ! |
     $     ymy2,                ! |
     $     ymy3,                ! |
     $     ymy4,                ! |
     $     coef12,              ! |
     $     coef34,              ! | -- interpolation weights/coeffs.
     $     term1(pgls),         ! |
     $     term2(pgls),         ! |
     $     term3(pgls),         ! |
     $     term4(pgls)          ! |
C
C-----------------------------------------------------------------------
C
      if( ppdy .ne. 4) then
         write(6,*) 'LAGYIN:Error:  ppdy .ne. 4'
         call endrun
      end if
C
      do i = 1,pgls
         ymy3     = ydp(i) - wdy(3,1,jdp(i))
         ymy4     = ydp(i) - wdy(4,1,jdp(i))
         coef12   = ymy3*ymy4   
         ymy2     = ydp(i) - wdy(2,1,jdp(i))
         term1(i) = coef12*ymy2*wdy(1,2,jdp(i))
         ymy1     = ydp(i) - wdy(1,1,jdp(i))
         term2(i) = coef12*ymy1*wdy(2,2,jdp(i))
         coef34   = ymy1*ymy2   
         term3(i) = coef34*ymy4*wdy(3,2,jdp(i))
         term4(i) = coef34*ymy3*wdy(4,2,jdp(i))
      end do
C
C Loop over fields.
C
      do m = 1,pf
         do i = 1,pgls
            fdp(i,m) = fint(i,1,m)*term1(i) + fint(i,2,m)*term2(i) +
     $                 fint(i,3,m)*term3(i) + fint(i,4,m)*term4(i)
         end do
      end do
C
      end subroutine lagyin

c##############################################################################

      subroutine lcbas( grd, bas1, bas2 )
C-----------------------------------------------------------------------
C
C Evaluate the partial Lagrangian cubic basis functions (denominator
C only ) for the grid points and gather grid values
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) grd(4)               ! grid stencil
C
C Output arguments
C
      real(r8) bas1(4),             ! grid values on stencil
     $     bas2(4)              ! lagrangian basis functions
C
C---------------------------Local variables-----------------------------
C
      real(r8) x0mx1,               ! |
     $     x0mx2,               ! |
     $     x0mx3,               ! |- grid value differences used in
     $     x1mx2,               ! |  weights
     $     x1mx3,               ! |
     $     x2mx3                ! |
C
C-----------------------------------------------------------------------
C
      x0mx1   = grd(1) - grd(2)
      x0mx2   = grd(1) - grd(3)
      x0mx3   = grd(1) - grd(4)
      x1mx2   = grd(2) - grd(3)
      x1mx3   = grd(2) - grd(4)
      x2mx3   = grd(3) - grd(4)
C
      bas1(1) = grd(1)
      bas1(2) = grd(2)
      bas1(3) = grd(3)
      bas1(4) = grd(4)
C
      bas2(1) =  1._r8/ ( x0mx1 * x0mx2 * x0mx3 )
      bas2(2) = -1._r8/ ( x0mx1 * x1mx2 * x1mx3 )
      bas2(3) =  1._r8/ ( x0mx2 * x1mx2 * x2mx3 )
      bas2(4) = -1._r8/ ( x0mx3 * x1mx3 * x2mx3 )
C
      end subroutine lcbas

c##############################################################################

      subroutine lcdbas(grd,dbas2,dbas3)
C-----------------------------------------------------------------------
C
C Calculate weights used to evaluate derivative estimates at the
C inner grid points of a four point stencil based on Lagrange
C cubic polynomial through four unequally spaced points.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C

C-----------------------------------------------------------------------
      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) grd(4)               ! grid stencil
C
C Output arguments
C
      real(r8) dbas2(4),            ! derivatives at grid point 2.
     $     dbas3(4)             ! derivatives at grid point 3.
C
C-----------------------------------------------------------------------
C
C  grd    Coordinate values of four points in stencil.
C  dbas2  Derivatives of the four basis functions at grid point 2.
C  dbas3  Derivatives of the four basis functions at grid point 3.
C
C---------------------------Local variables-----------------------------
C
      real(r8) x1,                  ! |
     $     x2,                  ! |- grid values
     $     x3,                  ! |
     $     x4,                  ! |
     $     x1mx2,               !  |
     $     x1mx3,               !  |
     $     x1mx4,               !  |- differences of grid values
     $     x2mx3,               !  |
     $     x2mx4,               !  |
     $     x3mx4                !  |
C
C-----------------------------------------------------------------------
C
      x1 = grd(1)
      x2 = grd(2)
      x3 = grd(3)
      x4 = grd(4)
      x1mx2 = x1 - x2
      x1mx3 = x1 - x3
      x1mx4 = x1 - x4
      x2mx3 = x2 - x3
      x2mx4 = x2 - x4
      x3mx4 = x3 - x4
C
      dbas2(1) =   x2mx3 * x2mx4 / ( x1mx2 * x1mx3 * x1mx4 )
      dbas2(2) =   -1._r8/x1mx2 + 1._r8/x2mx3 + 1._r8/x2mx4
      dbas2(3) = - x1mx2 * x2mx4 / ( x1mx3 * x2mx3 * x3mx4 )
      dbas2(4) =   x1mx2 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )
C
      dbas3(1) = - x2mx3 * x3mx4 / ( x1mx2 * x1mx3 * x1mx4 )
      dbas3(2) =   x1mx3 * x3mx4 / ( x1mx2 * x2mx3 * x2mx4 )
      dbas3(3) =   -1._r8/x1mx3 - 1._r8/x2mx3 + 1._r8/x3mx4
      dbas3(4) = - x1mx3 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )
C
      end subroutine lcdbas

c##############################################################################

      subroutine limdx(pidim   ,ibeg    ,len     ,dx      ,f       ,
     $                 fxl     ,fxr     )
C-----------------------------------------------------------------------
C
C Limit the derivative estimates for data on an equally spaced grid
C so they satisfy the SCM0 condition, that is, the spline will be 
C monotonic, but only C0 continuous on the domain
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
c++bee
c      integer pbpts             ! (length of latitude slice)*fields
c      parameter(pbpts = plond*plev*pcnst) 
c--bee
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     pidim,               ! vector dimension
     $     ibeg,                ! index of vector to begin computation
     $     len                  ! length of vector to compute
      real(r8) dx,                  ! length of grid inteval
     $     f(pidim)             ! field
C
C Input/output arguments
C
      real(r8) fxl(pidim),          ! x-derivs at left  edge of interval
     $     fxr(pidim)           ! x-derivs at right edge of interval
C
C-----------------------------------------------------------------------
C
C  pidim   Length of f, fxl, and fxr.
C  ibeg    First interval of grid for which derivatives are computed.
C  len     Number of grid intervals for which derivatives are computed.
C          (There are pidim - 1 intervals between the pidim gridpoints
C          represented in f, fxl, and fxr.)
C  dx      Value of grid spacing.
C  f       Values on equally spaced grid from which derivatives fxl and
C          fxr were computed.
C  fxl     fxl(i) is the limited derivative at the left  edge of
C          interval
C  fxr     fxr(i) is the limited derivative at the right edge of
C          interval
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                   ! index
     $     iend                 ! index to end work on vector
      real(r8) rdx                  ! 1./dx
      real(r8) deli(pidim)          ! simple linear derivative
C
C-----------------------------------------------------------------------
C
C
      iend = ibeg + len - 1
      rdx  = 1._r8/dx
C
      do i = ibeg,iend
         deli(i) = ( f(i+1) - f(i) )*rdx
      end do
C
C Limiter
C
      call scm0( len, deli(ibeg), fxl(ibeg), fxr(ibeg) )
C
      end subroutine limdx

c##############################################################################

      subroutine limdy(pf      ,fint    ,dy      ,jdp     ,fyb     ,
     $                 fyt     )
C-----------------------------------------------------------------------
C
C Limit the y-derivative estimates so they satisy the SCM0 for the
C x-interpolated data corresponding to the departure points of a single
C latitude slice in the global grid, that is, they are monotonic, but
C spline has only C0 continuity
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      integer
     $     pf                        ! dimension (number of fields)
      real(r8) fint(plon,plev,ppdy,pf),  ! x-interpolants
     $     dy(platd)                 ! interval lengths in lat grid
      integer
     $     jdp(plon,plev)            ! j-index of coord. of dep. pt.
C
C Input/output arguments
C
      real(r8) fyb(plon,plev,pf),        ! y-derivatives at bot of interval
     $     fyt(plon,plev,pf)         ! y-derivatives at top of interval
C
C-----------------------------------------------------------------------
C
C  pf      Number of fields being interpolated.
C  fint    (fint(i,k,j,m),j=1,ppdy) contains the x interpolants at each
C          latitude needed for the y derivative estimates at the
C          endpoints of the interval that contains the departure point
C          for grid point (i,k).  The last index of fint allows for
C          interpolation of multiple fields.  fint is generated by a
C          call to herxin.
C  dy      Increment in the y-coordinate value for each interval in the
C          extended array.
C  jdp     jdp(i,k) is the index of the y-interval that contains the
C          departure point corresponding to global grid point (i,k) in
C          the latitude slice being forecasted.
C          Suppose yb contains the y-coordinates of the extended array
C          and ydp(i,k) is the y-coordinate of the departure point
C          corresponding to grid point (i,k).  Then,
C                yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
C  fyb     fyb(i,k,.) is the limited derivative at the bot of the y
C          interval that contains the departure point of global grid
C          point (i,k).
C  fyt     fyt(i,k,.) is the limited derivative at the top of the y
C          interval that contains the departure point of global grid
C          point (i,k).
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,k,m,                    ! indices
     $     jb,                       ! index for bottom of interval
     $     jt                        ! index for top    of interval
      real(r8) rdy (plon,plev),          ! 1./dy
     $     deli(plon,plev)           ! simple linear derivative
C
C-----------------------------------------------------------------------
C
      jb = ppdy/2
      jt = jb + 1
C
      do k = 1,plev
         do i = 1,plon
            rdy(i,k) = 1._r8/dy(jdp(i,k))
         end do
      end do
C
C Loop over fields.
C
      do m = 1,pf
         do k = 1,plev
            do i = 1,plon
               deli(i,k) = ( fint(i,k,jt,m) - fint(i,k,jb,m) )*rdy(i,k)
            end do
         end do
C
C Limiter
C
         call scm0( plon*plev, deli, fyb(1,1,m), fyt(1,1,m) )
      end do
C
      end subroutine limdy

c##############################################################################

      subroutine limdz(pf,f,dsig,fst,fsb)
C-----------------------------------------------------------------------
C
C Apply SCMO limiter to vertical derivative estimates on a vertical
C slice.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
      integer plevm1
      parameter( plevm1 = plev - 1 )
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $  pf                        ! number of input fields
      real(r8) f(plon,plev,pf),       ! input field
     $     dsig(plev)             ! size of vertical interval
C
C Input/output arguments
C
      real(r8) fst(plon,plev,pf),  ! z-derivative at top of interval
     $     fsb(plon,plev,pf)   ! z-derivative at bot of interval
C
C-----------------------------------------------------------------------
C
C  f       Field values used to compute the discrete differences for
C          each interval in the vertical grid.
C  dsig    Increment in the sigma-coordinate value for each interval.
C  fst     Limited derivative at the top of each interval.
C  fsb     Limited derivative at the bottom of each interval.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                     ! longitude   index
     $     k,                     ! vertical    index
     $     m                      ! constituent index
      real(r8) rdsig,                 ! 1./dsig
     $     deli(plon,plevm1)      ! simple linear derivative
C
C-----------------------------------------------------------------------
C
C Loop over fields.
C
      do m = 1,pf
         do k = 1,plevm1
            rdsig = 1._r8/dsig(k)
            do i = 1,plon
               deli(i,k) = ( f(i,k+1,m) - f(i,k,m) )*rdsig
            end do
         end do
         call scm0( plon*plevm1, deli, fst(1,1,m), fsb(1,1,m) )
      end do
C
      end subroutine limdz

c##############################################################################

      subroutine s2gphi(lam     ,cosphi  ,sinphi  ,lamsc   ,phisc   ,
     $                  phigc   )
C-----------------------------------------------------------------------
C
C Calculate transformed local geodesic latitude coordinates for a set
C of points, each of which is associated with a grid point in a global
C latitude slice. Transformation is spherical to local geodesic.
C (Williamson and Rasch, 1991)
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      real(r8) lam(plon),           ! long coordinates of model grid
     $     cosphi,              ! cos(latitude)
     $     sinphi,              ! sin(latitude)
     $     lamsc(plon,plev),    ! spher. long coords of dep points
     $     phisc(plon,plev)     ! spher. lat  coords of dep points
C
C Output arguments
C
      real(r8) phigc(plon,plev)     ! loc geod. lat coords of dep points
C
C-----------------------------------------------------------------------
C
C  lam    longitude coordinates of the global grid points in spherical
C         system. The grid points in the global array are the reference
C         points for the local geodesic systems.
C  cosphi cosine of the latitude of the global latitude slice.
C  sinphi sine of the latitude of the global latitude slice.
C  lamsc  longitude coordinate of dep. points in spherical coordinates.
C  phisc  latitude  coordinate of dep. points in spherical coordinates.
C  phigc  latitude  coordinate of dep. points in local geodesic coords.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,k                  ! longitude, level indices
      real(r8) sphisc,              ! |
     $     cphisc,              ! | -- temporary variables
     $     clamsc               ! |
C-----------------------------------------------------------------------
C
      do k = 1,plev
         do i = 1,plon
            sphisc = sin( phisc(i,k) )
            cphisc = cos( phisc(i,k) )
            clamsc = cos( lam(i) - lamsc(i,k) )
            phigc(i,k) = asin( sphisc*cosphi - cphisc*sinphi*clamsc )
         end do
      end do
C
      end subroutine s2gphi

c##############################################################################

      subroutine s2gvel(udp     ,vdp     ,lam     ,cosphi  ,sinphi  ,
     $                  lamdp   ,phidp   ,upr     ,vpr     )
C-----------------------------------------------------------------------
C
C Transform velocity components at departure points associated with a
C single latitude slice from spherical coordinates to local geodesic
C coordinates. (Williamson and Rasch, 1991)
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      real(r8) udp(plon,plev),      ! u in spherical coords.
     $     vdp(plon,plev),      ! v in spherical coords.
     $     lam(plon),           ! x-coordinates of model grid
     $     cosphi,              ! cos(latitude)
     $     sinphi,              ! sin(latitude)
     $     lamdp(plon,plev),    ! spherical longitude coord of dep pt.
     $     phidp(plon,plev)     ! spherical latitude  coord of dep pt.
C
C Output arguments
C
      real(r8) upr(plon,plev),      ! u in local geodesic coords.
     $     vpr(plon,plev)       ! v in local geodesic coords.
C
C-----------------------------------------------------------------------
C
C  udp    u-component of departure point velocity in spherical coords.
C  vdp    v-component of departure point velocity in spherical coords.
C  lam    Longitude of arrival point position (model grid point) in
C         spherical coordinates.
C  cosphi Cos of latitude of arrival point positions (model grid pt).
C  sinphi Sin of latitude of arrival point positions (model grid pt).
C  lamdp  Longitude of departure point position in spherical
C         coordinates.
C  phidp  Latitude  of departure point position in spherical
C         coordinates.
C  upr    u-component of departure point velocity in geodesic coords.
C  vpr    v-component of departure point velocity in geodesic coords.
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,k                  ! longitude, level indices
      real(r8) cdlam,               ! |
     $     clamp,               ! |
     $     cphid,               ! |
     $     cphip,               ! |
     $     dlam,                ! | -- temporary variables
     $     sdlam,               ! |
     $     slamp,               ! |
     $     sphid,               ! |
     $     sphip                ! |
C-----------------------------------------------------------------------
C
      do k = 1,plev
         do i = 1,plon
            dlam  = lam(i) - lamdp(i,k)
            sdlam = sin( dlam )
            cdlam = cos( dlam )
            sphid = sin( phidp(i,k) )
            cphid = cos( phidp(i,k) )
            sphip = sphid*cosphi - cphid*sinphi*cdlam
            cphip = cos( asin( sphip ) )
            slamp = -sdlam*cphid/cphip
            clamp = cos( asin( slamp ) )
C
            vpr(i,k) = (vdp(i,k)*(cphid*cosphi + sphid*sinphi*cdlam) -
     $                  udp(i,k)*sinphi*sdlam)/cphip
            upr(i,k) = (udp(i,k)*cdlam + vdp(i,k)*sphid*sdlam +
     $                  vpr(i,k)*slamp*sphip)/clamp
         end do
      end do
C
      end subroutine s2gvel

c##############################################################################

      subroutine scm0(n,deli,df1,df2)
C-----------------------------------------------------------------------
C
C Apply SCM0 limiter to derivative estimates.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      integer
     $     n                    ! length of vectors
      real(r8) deli(n)              ! discrete derivative
C
C Input/Output arguments
C
      real(r8) df1(n),              ! limited left -edge derivative
     $     df2(n)               ! limited right-edge derivative
C
C-----------------------------------------------------------------------
C
C  n      Dimension of input arrays.
C  deli   deli(i) is the discrete derivative on interval i, i.e.,
C         deli(i) = ( f(i+1) - f(i) )/( x(i+1) - x(i) ).
C  df1    df1(i) is the limited derivative at the left  edge of interval
C  df2    df2(i) is the limited derivative at the right edge of interval
C
C---------------------------Local variables/Parameters------------------
C
      real(r8) eps     ! machine precision parameter to prevent rounding
C                  ! negatives in interpolant
      real(r8) fac     ! factor applied in limiter
#if ( REALSIZE == 4 )
      parameter( eps = 1.e-06, fac = 3.*( 1. - eps ) )
#else
      parameter( eps = 1.e-12, fac = 3.*( 1. - eps ) )
#endif
      integer i                 ! index
      real(r8) tmp1,                ! derivative factor
     $     tmp2                 ! abs(tmp1)
C
C-----------------------------------------------------------------------
C
      do i = 1,n
         tmp1 = fac*deli(i)
         tmp2 = abs( tmp1 )
C
         if( deli(i)*df1(i)   .le. 0.0  ) df1(i) = 0.
         if( deli(i)*df2(i)   .le. 0.0  ) df2(i) = 0.
         if( abs(    df1(i) ) .gt. tmp2 ) df1(i) = tmp1
         if( abs(    df2(i) ) .gt. tmp2 ) df2(i) = tmp1
      end do
C
      end subroutine scm0

c##############################################################################

      subroutine sltini(sinlam  ,coslam  ,ub      ,
     $                  fb      ,uxl     ,uxr     ,qxl     ,qxr
     $,                 nadv, x1 )
C-----------------------------------------------------------------------
C
C Prepare the extended arrays for use in the SLT routines
C
C   1)  Fill latitude extensions.
C   2)  Fill longitude extensions.
C   3)  Compute x-derivatives
C
C Computational note: The latitude loop in this routine is multitasked
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
      integer puvpts,           ! number of u/v pts in lat slice
     $        pqpts             ! number of constituent pts in lat slice
     $,       pxpts             ! number of variable pts in lat slice
      parameter(puvpts = plond*plev,
     $          pqpts  = plond*plev*pcnst) 
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8)
     $     sinlam(plon),                     ! sin(lamda)
     $     coslam(plon)                      ! cos(lamda)

      integer
     $  nadv               ! number of advected variables (non-species)
C
C Input/Output arguments
C
      real(r8) ub(plond,plev,2    ,platd),   ! u-component of wind
     $     fb(plond,plev,pcnst,platd)    ! constituent fields
     $,    x1(plond,plev,padv,platd)     ! advected variables (non-species)
C
C Output arguments
C
      real(r8) uxl (plond,plev, 2    ,platd),    ! x-deriv of u/v on left
     $     uxr (plond,plev, 2    ,platd),    ! x-deriv of u/v on right
     $     qxl (plond,plev, pcnst,platd),    ! x-deriv of constituent
C                                            !  on Left
     $     qxr (plond,plev, pcnst,platd)     ! x-deriv of constituent
C                                            !  on Right
C-----------------------------------------------------------------------
C
C  sinlam  Sin of longitudes in global grid (model grid pts only).
C  coslam  Cos of longitudes in global grid (model grid pts only).
C  ub      u/v-velocity component with extensions filled.
C  fb      array of "pcnst" scalar species with extensions filled
C  uxl     uxl ( , , 1, ) = x-derivatives of u at the left  (west) edge
C          of given interval
C          uxl ( , , 2, ) = x-derivatives of v at the left  (west) edge
C          of given interval
C  uxr     uxr ( , , 1, ) = x-derivatives of u at the right (east) edge
C          of given interval
C          uxr ( , , 2, ) = x-derivatives of v at the right (east) edge
C          of given interval
C  qxl     x-derivatives of scalar species at the left  (west) edge
C          of given interval
C  qxr     x-derivatives of scalar species at the right (east) edge
C          of given interval
C
C---------------------------Local variables-----------------------------
C
      integer
     $     j,                  ! index
     $     m                   ! index
C
C-----------------------------------------------------------------------
C
      pxpts = plond*plev*nadv

C Fill latitude extensions beyond the southern- and northern-most
C latitudes in the global grid
C
      call extyv( 2,plev,coslam,sinlam,ub )
      call extyv( 2,plev,coslam,sinlam,ub(1,1,2,1) )
      do m = 1,pcnst
         call extys( pcnst,plev,fb(1,1,m,1) )
      end do
      if ( nadv .gt. 0 ) then
         do m = 1, nadv
            call extys( nadv,plev,x1(1,1,m,1) )
         end do
      end if
C
C Fill longitude extensions
C
      call extx( 2,plev,ub )
      call extx( 2,plev,ub(1,1,2,1) )
      do m = 1,pcnst
         call extx( pcnst,plev,fb(1,1,m,1) )
      end do
      if ( nadv .gt. 0 ) then
         do m = 1, nadv
            call extx( nadv,plev,x1(1,1,m,1) )
         end do
      end if
C
C Compute x-derivatives.
C
! OpenMP directive
!$OMP PARALLEL DO
      do j = 1,platd
         call cubxdr(puvpts,2,puvpts-3,dlam,ub(1,1,1,j),
     $               uxl(1,1,1,j),uxr(1,1,1,j))
         call cubxdr(puvpts,2,puvpts-3,dlam,ub(1,1,2,j),
     $               uxl(1,1,2,j),uxr(1,1,2,j))
         call cubxdr(pqpts,2,pqpts-3,dlam,fb(1,1,1,j),
     $               qxl(1,1,1,j),qxr(1,1,1,j))
         if ( nadv .gt. 0 ) then
            call cubxdr(pxpts,2,pxpts-3,dlam,x1(1,1,1,j),
     $                  xxl(1,1,1,j),xxr(1,1,1,j))
         end if
         if( plimdr )then
            call limdx(pqpts,2,pqpts-3,dlam,fb(1,1,1,j),
     $                 qxl(1,1,1,j),qxr(1,1,1,j))
            if ( nadv .gt. 0 ) then
               call limdx(pxpts,2,pxpts-3,dlam,x1(1,1,1,j),
     $                    xxl(1,1,1,j),xxr(1,1,1,j))
            end if
         end if
      end do
C
      end subroutine sltini

c##############################################################################

      subroutine sphdep(jcen    ,jgc     ,dt      ,iterdp  ,
     $                  locgeo  ,ub      ,uxl     ,uxr     ,lam     ,
     $                  phib    ,lbasiy  ,lammp   ,phimp   ,lamdp   ,
     $                  phidp   ,idp     ,jdp     )
C-----------------------------------------------------------------------
C
C Compute departure points for semi-Lagrangian transport on surface of
C sphere using midpoint quadrature.  Computations are done in:
C
C   1) "local geodesic"   coordinates for "locgeo" = .true.
C   2) "global spherical" coordinates for "locgeo" = .false.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      integer
     $     jcen,                            ! index of lat slice (extnd)
     $     jgc                              ! index of lat slice (model)
      real(r8) dt                               ! time step (seconds)
      integer
     $     iterdp                           ! number of iterations
      logical
     $     locgeo                           ! computation type flag
c++bee
      real(r8) ub (plond,plev, 2       ,platd), ! u-comp of wind
c--bee
     $     uxl(plond,plev, 2       ,platd), ! left  edge x-deriv of u/v 
     $     uxr(plond,plev, 2       ,platd), ! right edge x-deriv of u/v
     $     lam(plond),                      ! long. coord. of model grid
     $     phib(platd),                     ! lat.  coord. of model grid
     $     lbasiy(4,2,platd)                ! lat interpolation weights
C
C Input/Output arguments
C
      real(r8) lammp(plon,plev),                ! long coord of midpoint
     $     phimp(plon,plev)                 ! lat  coord of midpoint
C
C Output arguments
C
      real(r8) lamdp(plon,plev),                ! long coord of dep. point
     $     phidp(plon,plev)                 ! lat  coord of dep. point
      integer
     $     idp  (plon,plev),                ! long index of dep. point
     $     jdp  (plon,plev)                 ! lat  index of dep. point
C
C-----------------------------------------------------------------------
C
C  jcen    Index in extended grid corresponding to latitude being
C          forecast.
C  jgc     Index in model    grid corresponding to latitude being
C          forecast.
C  dt      Time interval that parameterizes the parcel trajectory.
C  iterdp  Number of iterations used for departure point calculation.
C  locgeo  Logical flag to indicate computation in "local geodesic" or
C          "global spherical" space.
C  ub      Longitudinal velocity components in spherical coordinates.
C  uxl     uxl ( , , 1, ) = x-derivatives of u at the left  (west) edge
C          of given interval
C          uxl ( , , 2, ) = x-derivatives of v at the left  (west) edge
C          of given interval
C  uxr     uxr ( , , 1, ) = x-derivatives of u at the right (east) edge
C          of given interval
C          uxr ( , , 2, ) = x-derivatives of v at the right (east) edge
C          of given interval
C  lam     Longitude values for the extended grid.
C  phib    Latitude  values for the extended grid.
C  lbasiy  Weights for Lagrange cubic interpolation on the unequally
C          spaced latitude grid.
C  lammp   Longitude coordinates of the trajectory mid-points of the
C          parcels that correspond to the global grid points contained
C          in the latitude slice being forecast.  On entry lammp
C          is an initial guess.
C  phimp   Latitude coordinates of the trajectory mid-points of the
C          parcels that correspond to the global grid points contained
C          in the latitude slice being forecast.  On entry phimp
C          is an initial guess.
C  lamdp   Longitude coordinates of the departure points that correspond
C          to the global grid points contained in the latitude slice
C          being forecast.  lamdp is constrained so that
C                     0.0 .le. lamdp(i) .lt. 2*pi .
C  phidp   Latitude coordinates of the departure points that correspond
C          to the global grid points contained in the latitude slice
C          being forecast.  If phidp is computed outside the latitudinal
C          domain of the extended grid, then an abort will be called by
C          subroutine "trjgl".
C  idp     Longitude index of departure points.  This index points into
C          the extended arrays, e.g.,
C              lam (idp(i,k)) .le. lamdp(i,k) .lt. lam (idp(i,k)+1).
C  jdp     Latitude  index of departure points.  This index points into
C          the extended arrays, e.g.,
C              phib(jdp(i,k)) .le. phidp(i,k) .lt. phib(jdp(i,k)+1).
C
C------------------------ local variables ------------------------------
C
      integer
     $     iter,                   ! |
     $     i,                      ! |--indices
     $     k                       ! |
      real(r8) finc,                   ! time step factor
     $     dttmp,                  ! time step (seconds)
     $     dlam,                   ! increment of grid in x-direction
     $     phicen,                 ! latitude coord of current lat slice
     $     cphic,                  ! cos(phicen)
     $     sphic,                  ! sin(phicen)
     $     upr  (plon,plev),       ! u in local geodesic coords
     $     vpr  (plon,plev),       ! v in local geodesic coords
     $     lampr(plon,plev),       ! relative long coord of dep pt
     $     phipr(plon,plev),       ! relative lat  coord of dep pt
     $     uvmp (plon,plev,2),     ! u/v (spherical) interpltd to dep pt
     $     fint (plon,plev,ppdy,2) ! u/v x-interpolants
C
C-----------------------------------------------------------------------
C
      dlam   = lam(nxpt+2) - lam(nxpt+1)
      phicen = phib(jcen)
      cphic  = cos( phicen )
      sphic  = sin( phicen )
C
C Convert latitude coordinates of trajectory midpoints from spherical
C to local geodesic basis.
C
      if( locgeo ) call s2gphi(lam(i1),cphic,sphic,lammp,phimp,phipr)
C
C Loop over departure point iterates.
C
      do 30 iter = 1,iterdp
C
C Compute midpoint indicies.
C
         call bandij(dlam    ,phib    ,lammp   ,phimp   ,idp     ,
     $               jdp     )
C
C Interpolate velocity fields to midpoint locations (tensor product
C Lagrange cubic interpolation).
C
c++bee
         call hrvint(2       ,2       ,ub      ,uxl     ,uxr     ,
c--bee
     $               lam     ,lbasiy  ,lammp   ,phimp   ,idp     ,
     $               jdp     ,fint    ,uvmp    )
C
C Put u/v on unit sphere
C
         do k = 1,plev
            do i = 1,plon
               uvmp(i,k,1) = uvmp(i,k,1)*ra
               uvmp(i,k,2) = uvmp(i,k,2)*ra
            end do
         end do
C
C For local geodesic:
C
C   a) Convert velocity coordinates at trajectory midpoints from
C      spherical coordinates to local geodesic coordinates,
C   b) Estimate midpoint parcel trajectory,
C   c) Convert back to spherical coordinates
C
C Else, for global spherical
C
C      Estimate midpoint trajectory with no conversions
C
         if ( locgeo ) then
            call s2gvel(uvmp(1,1,1),uvmp(1,1,2),lam(i1),cphic,sphic,
     $                  lammp      ,phimp      ,upr    ,vpr  )
            call trajmp(dt,upr,vpr,phipr,lampr)
            dttmp = 0.5*dt
            call g2spos(dttmp   ,lam(i1) ,phib    ,phicen  ,cphic   ,
     $                  sphic   ,upr     ,vpr     ,lampr   ,phipr   ,
     $                  lammp   ,phimp   )
         else
            call trjmps(dt,uvmp(1,1,1),uvmp(1,1,2),phimp,lampr,phipr)
            finc = 1._r8
            call trjgl(jgc     ,finc    ,phicen  ,lam(i1) ,phib    ,
     $                 lampr   ,phipr   ,lammp   ,phimp   )
         end if
   30 continue          ! End of iter=1,iterdp loop
C
C Compute departure points in geodesic coordinates, and convert back
C to spherical coordinates.
C
C Else, compute departure points directly in spherical coordinates
C
      if (locgeo) then
         do k = 1,plev
            do i = 1,plon
               lampr(i,k) = 2.*lampr(i,k)
               phipr(i,k) = 2.*phipr(i,k)
            end do
         end do
         dttmp = dt
         call g2spos(dttmp   ,lam(i1) ,phib    ,phicen  ,cphic   ,
     $               sphic   ,upr     ,vpr     ,lampr   ,phipr   ,
     $               lamdp   ,phidp   )
      else
         finc = 2.
         call trjgl(jgc     ,finc    ,phicen  ,lam(i1) ,phib    ,
     $              lampr   ,phipr   ,lamdp   ,phidp   )
      end if
C
C Compute departure point indicies.
C
      call bandij(dlam    ,phib    ,lamdp   ,phidp   ,idp     ,
     $            jdp     )
C
      end subroutine sphdep

c##############################################################################

      subroutine trajmp(dt      ,upr     ,vpr     ,phipr   ,lampr   )
C-----------------------------------------------------------------------
C
C Estimate mid-point of parcel trajectory (geodesic coordinates) based
C upon horizontal wind field.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
      integer pgls           ! length of latitude slice
      parameter( pgls = plon*plev )
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) dt,               ! time step (seconds)
     $     upr(pgls),        ! u-component of wind in local geodesic
     $     vpr(pgls)         ! v-component of wind in local geodesic
C
C Input/Output arguments
C
      real(r8) phipr(pgls)       ! latitude coord of trajectory mid-point
C
C Output arguments
C
      real(r8) lampr(pgls)       ! longitude coord of trajectory mid-point
C
C-----------------------------------------------------------------------
C
C  dt      Time interval that corresponds to the parcel trajectory.
C  upr     u-coordinate of velocity corresponding to the most recent
C          estimate of the trajectory mid-point (in geodesic system).
C  vpr     v-coordinate of velocity corresponding to the most recent
C          estimate of the trajectory mid-point (in geodesic system).
C  phipr   Phi value at trajectory mid-point (geodesic coordinates).
C          On entry this is the most recent estimate.
C  lampr   Lambda value at trajectory mid-point (geodesic coordinates).
C
C---------------------------Local variables-----------------------------
C
      integer i              ! index
C
C-----------------------------------------------------------------------
C
      do i = 1,pgls
         lampr(i) = -.5*dt* upr(i) / cos( phipr(i) )
         phipr(i) = -.5*dt* vpr(i)
      end do
C
      end subroutine trajmp

c##############################################################################

      subroutine trjgl(jgc     ,finc    ,phicen  ,lam     ,phib    ,
     $                 lampr   ,phipr   ,lamp    ,phip    )
C-----------------------------------------------------------------------
C
C Map relative trajectory mid/departure point coordinates to global
C latitude/longitude coordinates and test limits
C           
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
      integer pgls              ! number of points in latitude slice
      parameter( pgls = plon*plev )
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     jgc                  ! latitude index in model grid
      real(r8) finc,                ! number of time increments
     $     phicen,              ! current latitude value in extnded grid
     $     lam(plon),           ! longitude values for the extended grid
     $     phib(platd),         ! latitude  values for the extended grid
     $     lampr(plon,plev),    ! relative x-coordinate of departure pt.
     $     phipr(plon,plev)     ! relative y-coordinate of departure pt.
C
C Output arguments
C
      real(r8) lamp (plon,plev),    ! long coords of traj midpoints
     $     phip (plon,plev)     ! lat  coords of traj midpoints
C
C-----------------------------------------------------------------------
C
C  jgc     Index in model grid corresponding to latitude being forecast.
C  finc    Time step factor (1. for midpoint, 2. for dep. point)
C  phicen  Latitude value for current latitude being forecast.
C  lam     Longitude values for the extended grid.
C  phib    Latitude  values for the extended grid.
C  lampr   Longitude coordinates (relative to the arrival point) of the
C          trajectory mid-points of the parcels that correspond to the
C          global grid points contained in the latitude slice being
C          forecast.
C  phipr   Latitude coordinates (relative to the arrival point) of the
C          trajectory mid-points of the parcels that correspond to the
C          global grid points contained in the latitude slice being
C          forecast.
C  lamp    Longitude coordinates of the trajectory mid-points of the
C          parcels that correspond to the global grid points contained
C          in the latitude slice being forecast.
C  phip    Latitude  coordinates of the trajectory mid-points of the
C          parcels that correspond to the global grid points contained
C          in the latitude slice being forecast.
C
C--------------------------Local variables------------------------------
C
      integer
     $     i,                   ! longitude index
     $     k,                   ! level index
     $     imax,                ! |
     $     imin,                ! | 
     $     iimin,               ! | -- indices
     $     kkmin,               ! |
     $     iimax,               ! |
     $     kkmax                ! |
      real(r8) pi,                  ! 3.14.......
     $     twopi                ! 2*pi
C
C------------------------------Externals--------------------------------
C
      integer
     $     ismin,               ! returns index of array minimum
     $     ismax                ! returns index of array maximum
C
C-----------------------------------------------------------------------
C
      pi = 4.*atan(1._r8)
      twopi = pi*2.
C
      do k = 1,plev
         do i = 1,plon
            lamp(i,k) = lam(i) + finc*lampr(i,k)
            phip(i,k) = phicen + finc*phipr(i,k)
            if(lamp(i,k) .ge. twopi) lamp(i,k) = lamp(i,k) - twopi
            if(lamp(i,k) .lt.   0.0) lamp(i,k) = lamp(i,k) + twopi
         end do
      end do
C
C Test that the latitudinal extent of trajectory is NOT over the poles
C
      imax  = ismax( pgls, phip, 1 )
      imin  = ismin( pgls, phip, 1 )
      kkmax = (imax-1)/plon + 1
      iimax = imax - (kkmax-1)*plon
      kkmin = (imin-1)/plon + 1
      iimin = imin - (kkmin-1)*plon
      if ( phip(iimax,kkmax) .ge. phib(j1+plat) ) then
         write(6,*)'TRJGL: ****** MODEL IS BLOWING UP *********'
         write(6,9000) iimax,kkmax,jgc
         call endrun
      else if( phip(iimin,kkmin) .lt. phib(j1-1) ) then
         write(6,*)'TRJGL: ****** MODEL IS BLOWING UP *********'
         write(6,9000) iimin,kkmin,jgc
         call endrun
      end if
C
C Formats
C
 9000 format(//'Parcel associated with longitude ',i5,', level ',i5,
     $   ' and latitude ',i5,' is outside the model domain.')
      end subroutine trjgl

c##############################################################################

      subroutine trjmps(dt      ,upr     ,vpr     ,phimp   ,lampr   ,
     $                  phipr   )
C-----------------------------------------------------------------------
C
C Estimate mid-point interval of parcel trajectory (global spherical
C coordinates).
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
      integer pgls              ! length of latitude slice
      parameter( pgls = plon*plev )
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) dt,                  ! time step (seconds)
     $     upr  (pgls),         ! u-comp of wind at midpoint
     $     vpr  (pgls),         ! v-comp of wind at midpoint
     $     phimp(pgls)          ! lat coord at midpoint
C
C Output arguments
C
      real(r8) lampr(pgls),         ! relative long coord of midpoint
     $     phipr(pgls)          ! relative lat  coord of midpoint
C
C-----------------------------------------------------------------------
C
C  dt      Time interval that corresponds to the parcel trajectory.
C  upr     u-coordinate of velocity corresponding to the most recent
C          estimate of the trajectory mid-point.
C  vpr     v-coordinate of velocity corresponding to the most recent
C          estimate of the trajectory mid-point.
C  phimp   Phi value of trajectory midpoint (most recent estimate).
C  lampr   Longitude coordinate of trajectory mid-point relative to the
C          arrival point.
C  phipr   Latitude  coordinate of trajectory mid-point relative to the
C          arrival point.
C
C---------------------------Local variables-----------------------------
C
      integer i                 ! index
C
C-----------------------------------------------------------------------
C
      do i = 1,pgls
         lampr(i) = -.5*dt* upr(i) / cos( phimp(i) )
         phipr(i) = -.5*dt* vpr(i)
      end do
C
      end subroutine trjmps

c##############################################################################

      subroutine vdplim(pkdim,sig,sigdp)
C-----------------------------------------------------------------------
C
C Restrict vertical departure points to be between the top and bottom
C sigma levels of the "full-" or "half-" level grid
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
      integer pgls              ! length of latitude slice
      parameter( pgls = plon*plev)
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     pkdim                ! vertical dimension
      real(r8) sig(pkdim)           ! vertical coordinate of model grid
C
C Input/Output arguments
C
      real(r8) sigdp(pgls)          ! vertical coords. of departure points.
C
C-----------------------------------------------------------------------
C
C  pkdim   Vertical dimension of "sig"
C  sig     Sigma values at the "full" or "half" model levels
C  sigdp   Sigma value at the trajectory endpoint or midpoint for each
C          gridpoint in a vertical slice from the global grid.  This
C          routine restricts those departure points to within the
C          model's vertical grid.
C
C---------------------------Local variables-----------------------------
C
      real(r8) eps                  ! machine precision to restrict dep
C                               ! point to be inside last grid point
#if ( REALSIZE == 4 )
      parameter( eps = 1.e-06 )
#else
      parameter( eps = 1.e-12 )
#endif
      integer i                 ! index
C
C-----------------------------------------------------------------------
C
      do i = 1,pgls
         if (sigdp(i) .lt. sig(1)) then
            sigdp(i) = sig(1)
         end if
         if (sigdp(i) .ge. sig(pkdim)) then
            sigdp(i) = sig(pkdim) - eps
         end if
      end do
C
      end subroutine vdplim

c##############################################################################

      subroutine vrtdep(dt      ,iterdp  ,wb      ,wst     ,
     $                  wsb     ,sig     ,sigh    ,dsigh   ,kdpmpf  ,
     $                  kdpmph  ,sigmp   ,sigdp   ,kdp     )
C-----------------------------------------------------------------------
C
C Compute vertical departure point and departure point index.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C
C Input arguments
C
      real(r8) dt                   ! time step (seconds)
      integer
     $     iterdp               ! number of iterations
      real(r8) wb (plon,plevp),     ! vertical velocity
     $     wst(plon,plevp),     ! z-derivative of wb at top of interval
     $     wsb(plon,plevp),     ! z-derivative of wb at bot of interval
     $     sig  (plev ),        ! sigma values of model full levels
     $     sigh (plevp),        ! sigma values of model half levels
     $     dsigh(plevp)         ! increment between half levels
      integer
     $     kdpmpf(pmap),        ! artificial grid indices
     $     kdpmph(pmap)         ! artificial grid indices
C
C Input/Output arguments
C
      real(r8) sigmp(plon,plev)     ! vert coords of traj mid-points
C
C Output arguments
C
      real(r8) sigdp(plon,plev)     ! vert coords of traj departure points
      integer
     $     kdp(plon,plev)       ! vertical departure point indices
C
C-----------------------------------------------------------------------
C
C  dt      Time interval that parameterizes the parcel trajectory.
C  iterdp  Number of iterations used for departure point calculation.
C  wb      Vertical velocity component (sigma dot).
C  wst     z-derivs at the top edge of each interval contained in wb
C  wsb     z-derivs at the bot edge of each interval contained in wb
C  sig     Sigma values at the full-index levels.
C  sigh    Half-index sigma levels including sigh(1) = sigma(1/2) = 0.0
C          sigh(plev+1) = sigma(plev+1/2) = 1.0 .  Note that in general
C          sigh(k) .lt. sig(k)  where sig(k) is the sigma value at the
C          k_th full-index level.
C  dsigh   Increment in half-index sigma levels.
C  kdpmpf  Array of indices of the model full levels which are mapped
C          into an artificial evenly spaced vertical grid.  Used to aid
C          in search for vertical position of departure point 
C  kdpmph  Array of indices of the model half levels which are mapped
C          into an artificial evenly spaced vertical grid.  Used to aid
C          in search for vertical position of departure point 
C  sigmp   Sigma value at the trajectory midpoint for each gridpoint
C          in a vertical slice from the global grid.  On entry sigmp is
C          an initial guess.
C  sigdp   Sigma value at the trajectory endpoint for each gridpoint
C          in a vertical slice from the global grid.
C  kdp     Vertical index for each gridpoint.  This index points into a
C          vertical slice array whose vertical grid is given by sig.
C          E.g.,   sig(kdp(i,k)) .le. sigdp(i,k) .lt. sig(kdp(i,k)+1).
C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                   ! |
     $     iter,                ! |-- indices
     $     k                    ! |
      real(r8) wmp(plond,plev)      ! vert vel. at midpoint

      integer nadv
      real(r8) xhrdum(1,1,1), xhstdum(1,1,1), xhsbdum(1,1,1), x2dum(1,1,1)  ! dummies
C
C-----------------------------------------------------------------------
C
C Loop over departure point iterates.
C
      do iter = 1,iterdp
C
C Compute midpoint indices in half-index sigma-level arrays (use kdp
C as temporary storage).
C
         call kdpfnd(plevp   ,sigh    ,sigmp   ,kdpmph  ,
     $               kdp     )
C
C Interpolate sigma dot field to trajectory midpoints using Hermite
C cubic interpolant.
C
         nadv = 0
         call herzin(plevp   ,1       ,wb      ,wst     ,wsb     ,
     $               sigh    ,dsigh   ,sigmp   ,kdp     ,wmp
     $,              nadv, xhrdum, xhstdum, xhsbdum, x2dum )
C
C Update estimate of trajectory midpoint.
C
         do k = 1,plev
            do i = 1,plon
               sigmp(i,k) = sig(k) - .5*dt*wmp(i,k)
            end do
         end do
C
C Restrict vertical midpoints to be between the top and bottom half-
C index sigma levels.
C
         call vdplim(plevp,sigh,sigmp)
      end do
C
C Compute trajectory endpoints.
C
      do k = 1,plev
         do i = 1,plon
            sigdp(i,k) = sig(k) - dt*wmp(i,k)
         end do
      end do
C
C Restrict vertical departure points to be between the top and bottom
C full-index sigma levels.
C
      call vdplim(plev,sig,sigdp)
C
C Vertical indices for trajectory endpoints that point into full-index
C sigma level arrays.
C
      call kdpfnd(plev    ,sig     ,sigdp   ,kdpmpf  ,
     $            kdp     )
C
      end subroutine vrtdep

c##############################################################################

      subroutine vrtmap(pkdim   ,sigln   ,dsigln  ,kdpmap  )
C-----------------------------------------------------------------------
C
C Map indices of an artificial evenly spaced (in log) vertical grid to
C the indices of the log of the model vertical grid.  The resultant
C array of mapped indices will be used by "kdpfnd" to find the vertical
C location of any departure point relative to the model grid.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C
C-----------------------------------------------------------------------
      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer
     $     pkdim                ! dimension of "sigln" and "dsigln"
      real(DBLKIND) sigln (pkdim),       ! model levels (log(eta))
     $     dsigln(pkdim)        ! intervals between model levels (log)
C
C Output arguments
C
      integer kdpmap(pmap)      ! array of mapped indices
C
C---------------------------Local variables-----------------------------
C
      integer
     $     imin,                ! |
     $     k,                   ! |-- indices
     $     kk,                  ! |
     $     newmap               ! estimated value of "pmap"
      real(DBLKIND) del                 ! artificial grid interval
      real(r8)
     $     dp,                  ! artificial departure point
     $     eps                  ! epsilon factor

#if ( REALSIZE == 4 )
      parameter( eps = 1.e-06 )
#else
      parameter( eps = 1.e-12 )
#endif
C
C------------------------------Externals--------------------------------
      integer  ismin
C-----------------------------------------------------------------------
C
      del = ( sigln(pkdim) - sigln(1) )/real(pmap,r8)
      imin = ismin( pkdim-1,real(dsigln,r8), 1 )
c-drb      imin = ismin( pkdim-1,dsigln, 1 )

      if(del + eps .ge. dsigln(imin)) then
         newmap = ( sigln(pkdim) - sigln(1) )/dsigln(imin) + 1
         write(6,9000) pmap,newmap
         call endrun
      end if      
C
      kdpmap(1) = 1
      do kk = 2,pmap
         dp = sigln(1) + real(kk-1,r8)*del
         do k = 1,pkdim-1
            if(dp .gt. sigln(k) + eps) then
               kdpmap(kk) = k
            end if
         end do
      end do
C
 9000 format(' VRTMAP:  Not enough artificial grid',
     $       ' intervals.'/
     $       ' Currently, "pmap" is set to ',i20/
     $       ' Reset parameter "pmap" to at least ',i20)
      end subroutine vrtmap

c##############################################################################

      end module advect
