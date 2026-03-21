! $Name$
! $Date$

#define DEBUG

      module advect

! Interface routines for SPITFIRE advection.

      use precision

      implicit none
      save
      private
      integer, parameter:: ilook=64, jlook=1, klook=18
      public ::
     $  adv_type    ! returns 'spitfire'
     $, ini_adv     ! Initialize advection.
     $, adv_dr      ! Advection driver.
     $, radvrst     ! read restart data
     $, wadvrst     ! write restart data
     $, advfluxout  ! write advective fluxes to history output

c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

      character(len=*), parameter ::
     $  type = 'spitfire'

      real(r8)
     $  etaint(plevp)       ! eta values at layer interfaces
     $, deta(plev)          ! etaint(k+1) - etaint(k)
     $, etamid(plev)        ! eta values at layer midpoints
     $, coslam(plon)        ! cos lamda
     $, sinlam(plon)        ! sin lamda

      real(r8) dgw(plat), dxlon

#ifdef CALC_ADV_FLUX
c pcnst constituents + air + Q
      real(r8)
     $  xfluxo(plond,plev,pcnst+2,plat) ! fluxes on x cell walls (except wrap point)
     $, yfluxo(plond,plev,pcnst+2,plat) ! fluxes on y cell walls (except N. pole)
     $, zfluxo(plond,plev,pcnst+2,plat) ! vertical mass flux (except top interface)
#else
      real(r8)
     $  xfluxo(1,1,1,1), yfluxo(1,1,1,1), zfluxo(1,1,1,1)
#endif

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

! Initialize SPITFIRE advection module.

      use spitfire, only: initrgrd, inifilter, inifixflux

      implicit none

c Input arguments:
      real(r8)
     $  xdlam           ! longitude spacing (not used)
     $, xlam(plond)     ! longitudes in extended grid (radians)
     $, xphi(platd)     ! Gauss latitudes in extended SLT grid (radians)
     $, gw(plat)        ! gauss weights for latitude grid
     $, xetaint(plevp)  ! hybrid coords at layer interfaces
     $, xetamid(plev)   ! hybrid coords at layer midpoints
     $, hybi(plevp)     ! sigma component of hybrid coordinate at interfaces
      integer
     $  nadv            ! number of misc transported species (not used)

      real(r8) work

c Local variables:
      integer
     $  i, ix, k
c-----------------------------------------------------------------------

      do k = 1, plevp
         etaint(k) = xetaint(k)
      end do
      work = 0.
      do k = 1, plev
         etamid(k) = xetamid(k)
         deta(k) = etaint(k+1) - etaint(k)
         work = work + deta(k)
      end do
      write (6,*) ' deta sum in advect: ', work
      write (6,*) ' eta top, bot and diff ', 
     $     etaint(plevp), etaint(1), etaint(plevp)-etaint(1) 

      do i = 1, plon
         ix = nxpt + i
         sinlam(i) = sin( xlam(ix) )
         coslam(i) = cos( xlam(ix) )
      end do

c     Initialize grid
      call initrgrd( xetaint, hybi, gw, xphi(j1) )
      dgw = gw
      dxlon = xlam(2)-xlam(1)
      write (6,*) ' dgw and dxlon ', dgw, dxlon

c     Initialize filter
      call inifilter

c     Initialize routines to fix the horizontal winds
c     to give mass convergence similar to the sfc pres change
      call inifixflux

      end subroutine ini_adv

c##############################################################################

      subroutine adv_dr( nstep, dtime, ufld, wfld
     $,                  ps1, as1, ps2, as2
     $,                  hw1, nadv, x1, x2, mcar, mca, wfldr )

c Advection Driver.

      use spitfire, only: hwctw, hdep, hadvn, setv, vadvn
     $,                   fixhflux, netdcal3, fixvflux


      implicit none

c Input arguments:
      integer 
     $  nstep        ! current model timestep
     $, nadv         ! number of misc advected variables
      real(r8)
     $  dtime                       ! time step in seconds
     $, ufld(plond,plev,2,platd)    ! u,v wind components
     $, ps1(plond,plat)             ! dry surface pressure at begining of timestep
     $, as1(plond,plev,pcnst,platd) ! species before advection (dry mass m.r.)
     $, x1(plond,plev,padv,platd)   ! miscellaneous variables before advection

c Input/Output arguments:
      real(r8)
     $  wfld(plon,plevp,plat)       ! vertical wind ("eta-dot").  Is either
                                    ! input from dynamics files or calculated.
     $, ps2(plond,plat)             ! dry surface pressure at end of timestep.  When the
                                    ! mass/wind consistency mods are used the input values
                                    ! are overwritten with the ps from SPITFIRE.
c Output arguments:
      real(r8)
     $  as2(plond,plev,pcnst,platd) ! species after advection (dry mass m.r.)
     $, x2(plond,plev,padv,platd)   ! miscellaneous variables after advection
     $, hw1(pcnst)                  ! global mass before advection
     $, wfldr(plon,plevp,plat)      ! raw vertical wind ("eta-dot")
     $, mca(plon,plevp,plat)        ! convergence of mass above (Pa/s)
     $, mcar(plon,plevp,plat)       ! raw convergence of mass above (Pa/s)

c Local variables:
      character*2
     $  trnum                       ! tracer number for messages from negative fixer
      integer
     $  i, j, k, m
     $, jcen                        ! lat index (extended grid) of forecast
      real(r8)
     $  pmid1(plond,plev)           ! pressure at layer midpoints at beginning of step
     $, pmid2(plond,plev)           ! pressure at layer midpoints at end of step
     $, pint1(plond,plevp)          ! pressure at interfaces at beginning of timestep
     $, pint2(plond,plevp)          ! pressure at interfaces at end of timestep
     $, pdel1(plond,plev)           ! pressure diffs between interfaces at beginning of timestep
     $, pdel2(plond,plev)           ! pressure diffs between interfaces at end of timestep
     $, d(plon,plev)                ! divergence
     $, pstend(plond,plat)          ! surface pressure tendency due to dry air redistribution
      real(r8)
     $  hw1j(pcnst,plat)            ! latitude contribution to global mass
      real(r8)
     $  pde(plond,plat,plev)        ! d(p)/d(eta) at beginning of timestep
     $, pdn(plond,plat,plev)        ! d(p)/d(eta) at end of timestep
     $, dn1(plond,plat,plev)        ! tmp array for densities
     $, dn(plond,plat,plev)         ! tmp array for densities

      real(r8) xflux(plond+1,plat,plev) ! fluxes on x cell walls
      real(r8) yflux(plond,plat+1,plev) ! fluxes on y cell walls
      real(r8) zflux(plond,plevp,plat)  ! vertical mass flux
      real(r8) zflx2(plond,plevp,plat)  ! vertical mass flux
      real(r8) xfluxint(plond+1,plat) ! integrated fluxes on x cell walls
      real(r8) yfluxint(plond,plat+1) ! integreated fluxes on y cell walls

      real(r8) pstmp(plond,plat), psdev, psobt(plond,plat)
      integer idev, jdev

      real(r8) mpsobs, mpsold, mpsnew   ! mean (global surface pressures)
      real(r8) pi

c Use static allocation for all local 3D fields.
      save pde, pdn, dn1, dn, xflux, yflux, zflux, zflx2

c Externals:
      logical
     $  calc_etadot     ! estimate etadot internally
     $, adj_winds       ! adjust winds for consistency between advection of the dry
                        ! atmospheric mass and the input surface pressure tendency.
      real(r8)
     $  p_top

      integer count0, count1, count2, count3,count_rate

c------------------------------------------------------------------------------

       call system_clock(count0,count_rate)

       hw1j = 0.0

c     Extend velocity array.
      call uvfill( ufld )

         mpsold = 0.
         mpsobs = 0.
         do j = 1,plat
            do i = 1,plon
               mpsold = mpsold + ps1(i,j)*dgw(j)
               mpsobs = mpsobs + ps2(i,j)*dgw(j)
            end do
         end do
         pi = 4*atan(1.0_r8)
         mpsold = mpsold*dxlon/(4*pi)
         mpsobs = mpsobs*dxlon/(4*pi)
!        write (6,*) ' mean old and obs pressure ', mpsold, mpsobs

! OpenMP directive
!$OMP PARALLEL DO PRIVATE(JCEN, PINT1, PINT2, PMID1, PMID2, PDEL1, PDEL2, D)
      do j = 1, plat
         jcen = j + nxpt + jintmx

c        Compute global mass integrals "hw1"
         call plevs0( ps1(1,j), pint1, pmid1, pdel1 )
         call plevs0( ps2(1,j), pint2, pmid2, pdel2 )
         call gamdry2( j, pcnst, as1(i1,1,1,jcen), pdel1, hw1j(1,j) )

c        calculate the appropriate density in the eta coordinate system
c        (rho*dz/deta)
         do k = 1, plev
            do i = 1, plon
               pde(i,j,k) =  pdel1(i,k)/deta(k)
            end do
            pde(plon+1,j,k) = pde(1,j,k)
         end do
         if ( adj_winds() ) then
            pstend(1:plon,j) = (ps2(1:plon,j) - ps1(1:plon,j))/dtime
         else
            if ( calc_etadot() ) then
c              Approximate the divergence field.
               call dcal( j, ufld(1,1,1,jcen-1), ufld(1,1,2,jcen-1), d )
c              Calculate etadot.
               call etdcal( j, ufld(i1,1,1,jcen), ufld(i1,1,2,jcen), d
     $,                     ps1, pdel1, pmid1, wfld(1,1,j) )
            end if
         end if
      end do

#ifdef DEBUG
      write (6,*) ' desired values (n, n+1obs, pstend) ', ilook, jlook, ps1(ilook,jlook), ps2(ilook,jlook), pstend(ilook,jlook)
#endif

c     Move horizontal winds from cell centers to cell walls.
      call hwctw( ufld(i1,1,1,j1), ufld(i1,1,2,j1), dtime, nstep )

c     Calculate horizontal departure points.
      call hdep( dtime, nstep )

c     Horizontal transport on the density.
#ifdef DEBUG
      call pdnavg('advect: input sfc pressure ', pde, deta)
#endif
      call hadvn( pde, dtime, nstep, dn1, xflux, yflux )
#ifdef DEBUG
      call pdnavg('advect: first guess dn1 sfc pressure ', dn1, deta)
#endif
      do j = 1,plat
         k = 1
         do i = 1,plon
            mcar(i,k,j) = 0.
         end do
         do k = 1,plev
            do i = 1,plon
               mcar(i,k+1,j) = mcar(i,k,j) + (dn1(i,j,k)-pde(i,j,k))*deta(k)/dtime
            end do
         end do
      end do
#ifdef DEBUG      
      write (6,*) ' sfc mcar, ptend ', ilook, jlook, mcar(ilook,plevp,jlook), pstend(ilook,jlook)
#endif

c     unadjusted pressure tendencies
         pstmp = 0.0
         do k = 1,plev
            do j = 1, plat
               do i = 1, plon
                  pstmp(i,j)= pstmp(i,j) + dn1(i,j,k)*deta(k)
               end do
            end do
         end do
#ifdef DEBUG
         call pgavg('advect: first guess surface pressure ', pstmp)
#endif
         ! pstmp now has the new pressure
         pstmp(1:plon,:) = pstmp(1:plon,:) + p_top()

         mpsnew = 0.
         do j = 1,plat
            do i = 1,plon
               mpsnew = mpsnew + pstmp(i,j)*dgw(j)
            end do
         end do
         mpsnew = mpsnew*dxlon/(4*pi)
!        write (6,*) ' mean new pressure ', mpsnew

         ! pstmp now has the new tendencies
         pstmp(1:plon,:) = (pstmp(1:plon,:) - ps1(1:plon,:))/dtime
#ifdef DEBUG         
         write (6,*) ' current pstend ', pstmp(ilook,jlook)
#endif

c        pure diagnostic
c        Now given the horizontal mass transport we can calculate a vert
c        vel.  wfld contains a estimated vertical velocity which balances
c        the lateral fluxes and surface pressure tendency.  zflx2 contains
c        the vertical flux.
         call netdcal3( xflux, yflux, dn1, pstmp, dtime, wfldr, zflx2 )
c         write (6,*) ' wfldr(ilook,plev,jlook) ', wfldr(ilook,plev,jlook) 

#ifdef DEBUG
         idev = 0
         jdev = 0
         psdev = 0.
         do j = 1, plat
            do i = 1, plon
               if (abs(pstmp(i,j)-pstend(i,j)).gt.psdev) then
                  idev = i
                  jdev = j
                  psdev = abs(pstmp(i,j)-pstend(i,j))
               endif
            end do
         end do
         
         write (6,*) ' advect.F: max deviation from desired tend is ', idev, jdev, psdev
         write (6,*) ' model tend before corr, observedtend there is ', pstmp(idev,jdev), pstend(idev,jdev)
         write (6,*) ' old pressure ps1 is ', ps1(idev,jdev)
         write (6,*) ' observed new pressure ps2 there is ', ps2(idev,jdev)
         write (6,*) ' model new pressure before corr. is ', ps1(idev,jdev)+pstmp(idev,jdev)*dtime
         write (6,*) ' advect.F: ilook jlook ', ilook, jlook
         write (6,*) ' model tend before corr, observedtend there is ', pstmp(ilook,jlook), pstend(ilook,jlook)
         write (6,*) ' old pressure ps1 is ', ps1(ilook,jlook)
         write (6,*) ' observed new pressure ps2 there is ', ps2(ilook,jlook)
         write (6,*) ' model new pressure before corr. is ', ps1(ilook,jlook)+pstmp(ilook,jlook)*dtime
#endif

      if ( adj_winds() ) then
         call system_clock(count1,count_rate)
c         psobt = pstend
         psobt(1:plon,1:plat) = pstend(1:plon,1:plat)

#ifdef DEBUG
         yfluxint = 0.
         xfluxint = 0.
         do k = 1,plev
            write (6,*) ' spitfire fluxes before correction ', k, idev, jdev,
     $           xflux(idev,jdev,k), xflux(idev+1,jdev,k), yflux(idev,jdev,k), yflux(idev,jdev+1,k)
            do j = 1,plat
               do i = 1,plon
                  xfluxint(i,j) = xfluxint(i,j) + xflux(i,j,k)*deta(k)
                  yfluxint(i,j) = yfluxint(i,j) + yflux(i,j,k)*deta(k)
               end do
            end do
         end do
         write (6,*) ' integrated spitfire fluxes before correction ', idev, jdev,
     $        xfluxint(idev,jdev), xfluxint(idev+1,jdev), yfluxint(idev,jdev), yfluxint(idev,jdev+1)
#endif

c        Check the predicted sfc pressure change compared to the interpolated
c        value and adjust winds where necessary to enforce consistency.
         call fixhflux( xflux, yflux, pstend, deta, dtime, pde, dn1 )
#ifdef DEBUG
!         call pgavg('advect: surface pressure tend after fixhflux ', pstend)
#endif

c        Recalculate departure points and horizontal transport.
         call hdep( dtime, nstep )
         call hadvn( pde, dtime, nstep, dn1, xflux, yflux )
!!!! DO IT AGAIN TO SEE WHAT HAPPENS
c        Check the predicted sfc pressure change compared to the interpolated
c        value and adjust winds where necessary to enforce consistency.
         call fixhflux( xflux, yflux, pstend, deta, dtime, pde, dn1 )
c        Recalculate departure points and horizontal transport.
         call hdep( dtime, nstep )
         call hadvn( pde, dtime, nstep, dn1, xflux, yflux )
         do j = 1,plat
            k = 1
            do i = 1,plon
               mca(i,k,j) = 0.
            end do
            do k = 1,plev
               do i = 1,plon
                  mca(i,k+1,j) = mca(i,k,j) + (dn1(i,j,k)-pde(i,j,k))*deta(k)/dtime
               end do
            end do
         end do

#ifdef DEBUG
         write (6,*) ' after fixhflux sfc mca, ptend at ', ilook, jlook, mca(ilook,plevp,jlook), pstend(ilook,jlook)

c         call pgavg('advect: surface pressure tend after second dep calc ', mca(:,plevp,:))
         yfluxint = 0.
         xfluxint = 0.
         do k = 1,plev
!            write (6,*) ' spitfire fluxes after correction ', k,
!     $           xflux(11,2,k), xflux(12,2,k), yflux(11,2,k), yflux(11,3,k)
            do j = 1,plat
               do i = 1,plon
                  xfluxint(i,j) = xfluxint(i,j) + xflux(i,j,k)*deta(k)
                  yfluxint(i,j) = yfluxint(i,j) + yflux(i,j,k)*deta(k)
               end do
            end do
         end do
         write (6,*) ' integrated spitfire fluxes after correction ',
     $        xfluxint(idev,jdev), xfluxint(idev+1,jdev), yfluxint(idev,jdev), yfluxint(idev,jdev+1)
#endif

c        Replace the end of timestep ps and pstend with the value calculated in SPITFIRE.
         ps2 = 0.0
         do k = 1,plev
            do j = 1, plat
               do i = 1, plon
                  ps2(i,j)= ps2(i,j) + dn1(i,j,k)*deta(k)
               end do
            end do
         end do
         ps2(1:plon,:) = ps2(1:plon,:) + p_top()
         pstend(1:plon,:) = (ps2(1:plon,:) - ps1(1:plon,:))/dtime


#ifdef DEBUG
         call pgavg('advect: surface pressure after hor correction ', ps2)

         write (6,*) ' pstend after correction ilook, jlook, ', ilook, jlook, pstend(ilook,jlook)
         write (6,*) ' ps2 after correction ilook, jlook, ', ilook, jlook, ps2(ilook,jlook)
         write (6,*) ' new pres and tend after corr from hfluxes at ',idev,jdev,
     $        ps2(idev,jdev), pstend(idev,jdev)
#endif

c        Now given the horizontal mass transport we can calculate a vert
c        vel.  wfld contains a estimated vertical velocity which balances
c        the lateral fluxes and surface pressure tendency.  zflx2 contains
c        the vertical flux.
         call netdcal3( xflux, yflux, dn1, pstend, dtime, wfld, zflx2 )
c         write (6,*) ' first guess wfld(ilook,plev,jlook) ', wfld(ilook,plev,jlook) 
c         write (6,*) ' zflux needed for mass balance ', zflx2(ilook,plev,jlook) 
         call system_clock(count2)
         write (6,*) ' elapsed time in hor wind adjustment is ', 1.*(count2-count1)/count_rate,' secs '
      end if

c     Calculate vertical departure points.
      call setv( wfld, dtime )

c     Vertical advection on the density.
      call vadvn( dn1, dtime, pdn, zflux )
c      write (6,*) ' zflux resulting from first guess ', zflux(ilook,plev,jlook)

      if ( adj_winds() ) then
         call system_clock(count1,count_rate)

c        Now that we know the tendencies and fluxes from spitfire adjust
c        the vertical velocities so that the spitfire fluxes match those
c        implied by the input fields.
         call fixvflux (wfld, zflx2, zflux)

c         write (6,*) ' revised wfld ', wfld(ilook,plev,jlook)
c        Calculate new vertical departure points.
         call setv( wfld, dtime )

c        Vertical advection on the density.
         call vadvn( dn1, dtime, pdn, zflux )
c         write (6,*) ' zflux final ', zflux(ilook,plev,jlook)
c        Replace the end of timestep ps with the value calculated in SPITFIRE.
         ps2 = 0.0
         do k = 1,plev
            do j = 1, plat
               do i = 1, plon
                  ps2(i,j) = ps2(i,j) + pdn(i,j,k)*deta(k)
               end do
            end do
         end do
         ps2(1:plon,:) = ps2(1:plon,:) + p_top()
         pstend(1:plon,:) = (ps2(1:plon,:) - ps1(1:plon,:))/dtime
         call system_clock(count2)
         write (6,*) ' elapsed time in vert wind adjustment is ', 1.*(count2-count1)/count_rate,' secs'
      end if
#ifdef DEBUG
      call pgavg('advect: surface pressure after vertical correction ', ps2)
      write (6,*) ' advect.F: new pres, tend after all fluxes now', idev, jdev, ps2(idev,jdev), pstend(idev,jdev)
      write (6,*) ' advect.F: sfc mca, ptend ', ilook, jlook, mca(ilook,plevp,jlook), pstend(ilook,jlook)
#endif



      pstmp(1:plon,:) = ps1(1:plon,:) + psobt(1:plon,:)*dtime

      idev = 0
      jdev = 0
      psdev = 0.
      do j = 1, plat
         do i = 1, plon
            if (abs(pstmp(i,j)-ps2(i,j)).gt.psdev) then
               idev = i
               jdev = j
               psdev = abs(pstmp(i,j)-ps2(i,j))
            endif
         end do
      end do
      write (6,*) ' worst agreement in mass/wind sfc Ps is at ', idev, jdev, psdev, ' Pa'
!     if we cant keep the surface pressure close to the obs issue a message
      if (psdev.gt.100.) then
        write (6,*) ' warning mass/wind consistency is not so good ', psdev
        call endrun()
      endif

#ifdef CALC_ADV_FLUX
      do j = 1, plat
         do k = 1, plev
            do i = 1, plon
               xfluxo(i,k,pcnst+1,j) = xflux(i,j,k)
               yfluxo(i,k,pcnst+1,j) = yflux(i,j,k)
               zfluxo(i,k,pcnst+1,j) = zflux(i,k+1,j)
            end do
         end do
      end do
#endif

!      write (6,*) ' possibly fixed vertical velocities ', wfld(ilook,plev,jlook), wfld(ilook,plevp,jlook)

#ifdef DEBUG

! remember, psobt now has the observed tendencies
         idev = 0
         jdev = 0
         psdev = 0.;
         do j = 1, plat
            do i = 1, plon
               if (abs(psobt(i,j)-pstend(i,j)).gt.psdev) then
                  idev = i
                  jdev = j
                  psdev = abs(psobt(i,j)-pstend(i,j))
               endif
            end do
         end do
         
      write (6,*) 'rescan of max deviation from observed tend is now at ', idev, jdev, psdev
      write (6,*) ' obtend, newtend there is ', psobt(idev,jdev), pstend(idev,jdev)
      write (6,*) ' ps1, ps2, psob there is now', ps1(idev,jdev), ps2(idev,jdev), 
     $     ps1(idev,jdev) + psobt(idev,jdev)*dtime
!      if (nstep.ge.4) then
!         write (6,*) ' early quit '
!         call abort()
!      endif
#endif

c     ***********************
c     * Advect the tracers. *
c     ***********************


!      write (6,*) ' tracers '
      do m = 1,pcnst
         write(trnum,'(I2)') m


c        Convert tracer mixing ratio to a density.
         do j = 1, plat
            jcen = j + nxpt + jintmx
            do k = 1, plev
               do i = 1, plon
                  dn(i,j,k) = as1(i1+i-1,k,m,jcen)*pde(i,j,k)
               end do
               dn(plon+1,j,k) = dn(1,j,k)
            end do
         end do

!         write (6,*) ' adv_dr bh ', dn(ilook,2,18), m
c        Horizontal advection.
         call hadvn( dn, dtime, nstep, dn1, xflux, yflux)

!         write (6,*) ' adv_dr ah ', dn1(ilook,2,18), m

c        Vertical advection.
         call vadvn( dn1, dtime, dn, zflux )

!        write (6,*) ' adv_dr av ', dn(ilook,2,18), m

#ifdef CALC_ADV_FLUX
         do j = 1, plat
            do k = 1, plev
               do i = 1, plon
                  xfluxo(i,k,m,j) = xflux(i,j,k)
                  yfluxo(i,k,m,j) = yflux(i,j,k)
                  zfluxo(i,k,m,j) = zflux(i,k+1,j)
               end do
            end do
         end do
#endif

c        Convert density back to mixing ratio.
         do j = 1, plat
            jcen = j + nxpt + jintmx
            do k = 1, plev
               do i = 1, plon
                  if (pdn(i,j,k).eq.0.) then
                     write (6,*) 'advect: WARNING!!!! PSEUDO AIR DENSITY ZERO ', i, j, k, 
     $                    ' fluxes are ',
     $                    xflux(i,j,k), xflux(i+1,j,k),
     $                    yflux(i,j,k), yflux(i,j+1,k),
     $                    zflux(i,k,j), zflux(i,k+1,j),
     $                    ' previous pseudo-density ',
     $                    pde(i,j,k)
                     as2(i1+i-1,k,m,jcen) = 0.
                  else
                     as2(i1+i-1,k,m,jcen) = dn(i,j,k)/pdn(i,j,k)
                  endif
               end do
            end do
         end do

      end do

c     Compute global integrals
      do m = 1, pcnst
         hw1(m) = 0.
      end do
      do m = 1, pcnst
         do j = 1, plat
            hw1(m) = hw1(m) + hw1j(m,j)
         end do
      end do

c     ********************************************
c     * Update miscellaneous advected variables. *
c     ********************************************

      if ( nadv .gt. 0 ) then
!        write (6,*) ' adv vars '
         do m = 1, nadv
c           Weight scalar field with density.
            do j = 1, plat
               jcen = j + nxpt + jintmx
               do k = 1, plev
                  do i = 1, plon
                     dn(i,j,k) = x1(i1+i-1,k,m,jcen)*pde(i,j,k)
                  end do
                  dn(plon+1,j,k) = dn(1,j,k)
               end do
            end do

!         write (6,*) ' adv_dr bh ', dn(ilook,2,18), m

            call hadvn( dn, dtime, nstep, dn1, xflux, yflux)
         
!         write (6,*) ' adv_dr ah ', dn1(ilook,2,18), m
            call vadvn( dn1, dtime, dn, zflux )

!         write (6,*) ' adv_dr av ', dn(ilook,2,18), m
#ifdef CALC_ADV_FLUX
            if ( m .eq. 1 ) then   ! specific humidity flux
               do j = 1, plat
                  do k = 1, plev
                     do i = 1, plon
                        xfluxo(i,k,pcnst+2,j) = xflux(i,j,k)
                        yfluxo(i,k,pcnst+2,j) = yflux(i,j,k)
                        zfluxo(i,k,pcnst+2,j) = zflux(i,k+1,j)
                     end do
                  end do
               end do
            end if
#endif

c           Remove density weighting.
            do j = 1, plat
               jcen = j + nxpt + jintmx
               do k = 1, plev
                  do i = 1, plon
                  if (pdn(i,j,k).eq.0.) then
                     write (6,*) 'advect: WARNING!!!! x2 PSEUDO AIR DENSITY ZERO ', i, j, k, 
     $                    ' fluxes are ',
     $                    xflux(i,j,k), xflux(i+1,j,k),
     $                    yflux(i,j,k), yflux(i,j+1,k),
     $                    zflux(i,j,k), zflux(i,j,k+1),
     $                    ' previous pseudo-density ',
     $                    pde(i,j,k)
                     x2(i1+i-1,k,m,jcen) = 0.
                  else
                     x2(i1+i-1,k,m,jcen) = dn(i,j,k)/pdn(i,j,k)
                  endif
                  end do
               end do
            end do
         end do
      end if

      call system_clock(count3,count_rate)
      write (6,*) ' elapsed time in adv_dr is ', 1.*(count3-count0)/count_rate,' secs'

      end subroutine adv_dr

c##############################################################################

      subroutine uvfill( uv )

C   1)  Fill latitude extensions.
C   2)  Fill longitude extensions.

      implicit none

c Input/Output argument:
      real(r8) uv(plond,plev,2,platd)   ! u and v components of wind
c-----------------------------------------------------------------------

c     Fill latitude extensions beyond the southern- and northern-most
c     latitudes in the global grid

      call extyv( 2, plev, 1, uv )
      call extyv( 2, plev, 2, uv )

c     Fill longitude extensions

      call extx( 2, plev, 1, uv )
      call extx( 2, plev, 2, uv )

      end subroutine uvfill

c##############################################################################

      subroutine extx( pkcnst, pkdim, pnter, fb )

c Copy data to the longitude extensions of the extended array

      implicit none

      integer i2pi              ! start of eastern long. extension
      parameter( i2pi = nxpt+plon+1 ) 

c Input arguments:
      integer
     $     pkcnst,              ! dimension construct for 3-D arrays
     $     pkdim                ! vertical dimension
      integer pnter

c Input/Output argument:
      real(r8) fb(plond,pkdim,pkcnst,platd) ! array to be extended

c Local variables:
      integer i, j, k
c-----------------------------------------------------------------------

c     Fill west edge points.

      if(nxpt .ge. 1) then
         do j=1,platd
            do k=1,pkdim
               do i=1,nxpt
                  fb(i,k,pnter,j) = fb(i+plon,k,pnter,j)
               end do
            end do
         end do
      end if

c     Fill east edge points

      do j=1,platd
         do k=1,pkdim
            do i=i2pi,plond
               fb(i,k,pnter,j) = fb(i-plon,k,pnter,j)
            end do
         end do
      end do

      end subroutine extx

c##############################################################################

      subroutine extyv( pkcnst, pkdim, pnter, vb )

c Fill latitude extensions of a vector component extended array.
c This is done in 2 steps:

c   1) interpolate to the pole points;
c      use coefficients for zonal wave number 1 on the Gaussian
c      latitude closest to the pole.
c   2) add latitude lines beyond the poles.

      implicit none

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

c Input arguments:
      integer
     $     pkcnst,              ! dimensioning construct for 3-D arrays
     $     pkdim                ! vertical dimension
     $,    pnter                ! index into u or v

c Input/Output argument:
      real(r8) vb(plond,pkdim,pkcnst,platd) ! one component of vector field

c Local variables
      integer
     $     i,                   ! index
     $     ig,                  ! index
     $     j,                   ! index
     $     k                    ! index
      real(r8) zavec,               ! accumulator for wavenumber 1
     $     zaves,               ! accumulator for wavenumber 1
     $     fplon2               ! float(plon2)
c-----------------------------------------------------------------------

      fplon2 = float(plon2)

c Fill north pole line.

      do k = 1,pkdim
         zavec = 0.0
         zaves = 0.0
         ig = 0
         do i = istart,istop
            ig    = ig + 1
            zavec = zavec + vb(i,k,pnter,jn  )*coslam(ig)
            zaves = zaves + vb(i,k,pnter,jn  )*sinlam(ig)
         end do
         zavec = zavec/fplon2
         zaves = zaves/fplon2
         ig    = 0
         do i = istart,istop
            ig           = ig + 1
            vb(i,k,pnter,jn+1) = zavec*coslam(ig) + zaves*sinlam(ig)
         end do
      end do

c Fill northern lines beyond pole line.

      if( jn+2 .le. platd )then
         do j = jn+2,platd
            do k = 1,pkdim
cdir$ ivdep
               do i = istart,istart+plon2-1
                  vb(      i,k,pnter,j) = -vb(plon2+i,k,pnter,2*jn+2-j)
                  vb(plon2+i,k,pnter,j) = -vb(      i,k,pnter,2*jn+2-j)
               end do
            end do
         end do
      end if
 
c Fill south pole line.

      do k = 1,pkdim
         zavec = 0.0
         zaves = 0.0
         ig = 0
         do i = istart,istop
            ig    = ig + 1
            zavec = zavec + vb(i,k,pnter,js  )*coslam(ig)
            zaves = zaves + vb(i,k,pnter,js  )*sinlam(ig)
         end do
         zavec = zavec/fplon2
         zaves = zaves/fplon2
         ig    = 0
         do i = istart,istop
            ig           = ig + 1
            vb(i,k,pnter,js-1) = zavec*coslam(ig) + zaves*sinlam(ig)
         end do
      end do

c Fill southern lines beyond pole line.

      if( js-2 .ge. 1 )then
         do j = 1,js-2
            do k = 1,pkdim
cdir$ ivdep
               do i = istart,istart+plon2-1
                  vb(      i,k,pnter,j) = -vb(plon2+i,k,pnter,2*js-2-j)
                  vb(plon2+i,k,pnter,j) = -vb(      i,k,pnter,2*js-2-j)
               end do
            end do
         end do
      end if
C
      end subroutine extyv

c##############################################################################

      subroutine advfluxout( lat, obuf )

c Write advective history to history/output buffers.

      use histout, only: outfld

      implicit none
c-----------------------------------------------------------------------
#include <tracnam.h>
c-----------------------------------------------------------------------

c Input argument:
      integer
     $  lat         ! latitude index for this slice

c Output arguments:
      real(r8)
     $  obuf(*)     ! output buffer

c Local variables:
      integer m
c-----------------------------------------------------------------------

      do m = 1, pcnst+2
         call outfld( xflxnam(m), xfluxo(1,1,m,lat), plond, lat, obuf )
         call outfld( yflxnam(m), yfluxo(1,1,m,lat), plond, lat, obuf )
         call outfld( zflxnam(m), zfluxo(1,1,m,lat), plond, lat, obuf )
      end do

      end subroutine advfluxout

c##############################################################################

      subroutine radvrst( lun )

c Read restart data for advection module.

      implicit none

c Input argument:

      integer
     $  lun    ! logical unit for fortran i/o.
c-----------------------------------------------------------------------

      end subroutine radvrst

c##############################################################################

      subroutine wadvrst( lun )

c Write restart data from advection module.

      implicit none

c Input argument:

      integer
     $  lun    ! logical unit for fortran write.
c-----------------------------------------------------------------------

      end subroutine wadvrst

c##############################################################################

      end module advect
