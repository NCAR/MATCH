! $Name$
! $Date$

! define for more info on fixflux
#define DEBUG
#undef DEBUG

! define for info on spitfire
#undef DEBUGSPITFIRE

! define for level info on fixflux
#define DEBUGLEV
#undef DEBUGLEV

! define REMEMBER to remember previous streamfunction 
! as first guess for this time
#undef REMEMBER

#define LIMITON
#define MONO
#undef FLUXFORM2
#define NWORK 1000

      module spitfire

! SPITFIRE advective transport algorithm.

      use precision

      implicit none
      save
      private
      public ::
     $  initrgrd
     $, inifilter
     $, inifixflux
     $, hwctw
     $, hdep
     $, hadvn
     $, setv
     $, vadvn
     $, fixhflux
     $, netdcal3
     $, fixvflux
     $, chi2fl
     $, laplc

       integer iter
       integer nitsm
       integer icount
       integer nitin
       real(r8) eer,eem

c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

      integer platext
      parameter (platext = 2*plat + 1)
      real(DBLKIND) ywext(platext+5)
      real(DBLKIND) dymin
      integer pjemax
      parameter (pjemax = 40000)
      integer jeven(pjemax)
      integer njeven
      real(DBLKIND) dzmin
      integer pkemax
      parameter (pkemax = 330000)
      integer keven(pkemax)
      integer nkeven

      integer nptslat(plat), bini(plond,plat)
      integer idep(plond+1,plat,plev)
      integer idom(plond+1,plat,plev)
      integer kdep(plond,plevp,plat+1)
      integer jdep(platext,plon/2,plev)
      integer jdom(platext,plon/2,plev)

      real(DBLKIND) hybi(plevp)         ! sigma component of hybrid coordinate at interfaces
      real(DBLKIND) deta(plev)          ! etaint(k+1) - etaint(k)
      real(DBLKIND) zw(plond,plevp)
      real(DBLKIND) xc(plond)
      real(DBLKIND) xw(plond+1)
      real(DBLKIND) yc(plat)
      real(DBLKIND) dmu(plat)
      real(DBLKIND) yw(plat+1)
      real(DBLKIND) cyw(plat+1)
      real(DBLKIND) ylat(plat)
      real(DBLKIND) cylat(plat)
      real(DBLKIND) dxlon
      real(DBLKIND) arad

      real(DBLKIND) uwdp(plond+1,plat,plev)
      real(DBLKIND) vwdp(platext,plon/2,plev)
      real(DBLKIND) wwdp(plond,plevp,plat+1)
      real(DBLKIND) uw(plond+1,plat,plev)
      real(DBLKIND) vw(plond,plat+1,plev)
      real(r8) debugdt

      real(DBLKIND) :: tinynum = 2*epsilon(1._DBLKIND)
      real(r8) :: tinyr8 = 2*epsilon(1._r8)


      real(DBLKIND) :: one = 1.
      real(DBLKIND) :: two = 2.
      real(DBLKIND) :: zero = 0.

      integer :: debugit = 0
      integer :: ilook = 64, jlook=1, klook = 18

c##############################################################################
      contains
c##############################################################################

      subroutine initrgrd( sigh, xhybi, gw, lat )


c     author: Phil Rasch

      implicit none

c     input arguments

      real(r8) sigh(plevp)                  ! sigma or eta coordinates
     $, xhybi(plevp)                    ! sigma component of hybrid coordinate at interfaces
     $, gw(plat)                        ! gaussian weights 
     $, lat(plat)                       ! gauss lats

c     local variables
      integer i, j, k
      real(DBLKIND) pi                           ! work var
      real(DBLKIND) yeven                        ! work var
!      real(DBLKIND) cosc(plat), sinc(plat), colat(plat)
!     $        , wos2(plat)


c     make sure the work arrays are big enough
      i = max(plond,2*plat+5)
      if (i .gt. NWORK) then
         write (6,*) 'initgrd: increase nwork to at least ', i
         stop
      endif

      pi = 4.*atan(1.0)


      arad = 1.
      arad = 6.37e6
      write (6,*) ' radius of earth is assumed to be ', arad, ' meters'

      do k = 1,plevp
         hybi(k) = xhybi(k)
      end do
      call initvt (sigh)
      
c     calculate our own gaus lats, since NCEP returns only 5 digits
c     call lggaus( plat, cosc, gw, sinc, colat, wos2 )

c     grid for horizontal transforms
      yw(1) = -1.
      do j = 1,plat
         yw(j+1) = yw(j) + gw(j)
         dmu(j) = gw(j)
c        lat(j) = pi/2.-colat(j)
      end do
      yw(1) = -pi/2.
      yw(plat+1) = pi/2.
      cyw(1) = 0.
      cyw(plat+1) = 0.
      do j = 2,plat
         yw(j) = asin(yw(j))
         cyw(j) = cos(yw(j))
      end do
      do j = 1,plat
c        ylat(j) = (yw(j+1)+yw(j))*0.5
         ylat(j) = lat(j)
         cylat(j) = cos(ylat(j))
         yc(j) = ylat(j)
      end do

      dxlon = 2.*pi/plon
      do i = 1,plon+1
c        xw(i) = (i-1)*dxlon
         xw(i) = (i-1.5)*dxlon
      end do

c     cell centered quantities
      do i = 1,plon
         xc(i) = (xw(i) + xw(i+1))*0.5
      end do
      xc(plon+1) = xc(plon) + dxlon

      do j = 1,plat+1
#ifdef FLUXFORM2
         ywext(j) = yw(j)
#else    
         ywext(j) = sin(yw(j))
#endif   
      end do
      do j = plat+2,platext
         ywext(j) = ywext(j-plat) + (ywext(plat+1)-ywext(1))
      end do

c      write (6,*) 'ywext initialized '
c      do j = 1,platext
c      write (6,*) j, ywext(j)
c      end do

      dymin = 1.e36
      do j = 1,platext-1
         dymin = min(dymin,abs(ywext(j+1)-ywext(j)))
      end do

      njeven = abs(ywext(platext)-ywext(1))/dymin + 2
      write (6,*) 'even y grid requires nypts = ', njeven, dymin
      if (njeven.gt.pjemax) then
         write (6,*) ' the even grid jeven isnt big enough '
         write (6,*) ' reset to greater than ', njeven
         stop 10
      endif
      dymin = (ywext(platext)-ywext(1))/(njeven-1)
      write (6,*) ' dymin reset to ', dymin
      do j = 1,njeven
         yeven = min( (j-1)*dymin + ywext(1), ywext(platext) )
         jeven(j) = 0
         do i = 1,platext-1
            if (yeven.ge.ywext(i).and.yeven.le.ywext(i+1)) then
               jeven(j) = i
            endif
         end do
         if (jeven(j).eq.0.) then
            write (6,*) 'couldnt map yeven ', j, yeven
            stop
         else
            i = jeven(j)
c           write (6,*) 'even grid ', j, ' at y ', yeven,
c     $           ' is on uneven grid interval ', 
c     $           jeven(j), ywext(i), ywext(i+1)
         endif
      end do

      end subroutine initrgrd

c##############################################################################

      subroutine initvt (sigh)

c     author: Phil Rasch

c     initialize the vertical grid

      implicit none

c     input variables
      real(r8) sigh(plevp)                  ! vertical coords of cell walls
      real(DBLKIND) sigh2(plevp)                  ! vertical coords of cell walls
      

c     local variables
      integer i
      integer k
      real(DBLKIND) zeven

      sigh2 = (1._DBLKIND)*sigh

      do k = 1,plev
         deta(k) = sigh2(k+1) - sigh2(k)
      end do
 
c     grid for vertical transforms
      write (6,*) ' number of interfaces ', plevp
      do k = 1,plevp
         write (6,*) ' initvt: sigh, sigh2 ', k, sigh(k), sigh2(k)
         do i = 1,plon
            zw(i,k) = sigh2(k)
         end do
      end do

      dzmin = 1.e36
      do k = 1,plev
c        write (6,*) ' dzmin dz ', dzmin, abs(sigh(k+1)-sigh(k))
         dzmin = min(dzmin,abs(sigh2(k+1)-sigh2(k)))
      end do
      write (6,*) ' smallest z grid increment is ', dzmin

      nkeven = abs(sigh2(plevp)-sigh2(1))/dzmin + 2
      write (6,*) 'regular z grid requires nypts = ', nkeven, dzmin
      if (nkeven.gt.pkemax) then
         write (6,*) ' the even grid keven isnt big enough '
         write (6,*) ' reset to greater than ', nkeven
         stop 20
      endif

      dzmin = (sigh2(plevp)-sigh2(1))/(nkeven-1)
      write (6,*) ' dzmin reset to ', dzmin
      do k = 1,nkeven

         zeven = (k-1)*dzmin + sigh2(1)
         keven(k) = 0
         do i = 1,plev
#ifdef PGF90
! put this write statement in with Portland Group Compiler to prevent optimizer
! breaking the following if check            
            write (6,*) 'zeven-sigh2(i) ', zeven-sigh2(i),zeven-sigh2(i+1)
#endif
            if (zeven.ge.sigh2(i).and.zeven.le.sigh2(i+1)) then
               keven(k) = i
            endif
         end do

         if (keven(k).eq.0.) then
            write (6,*) 'couldnt map zeven ', k, zeven
            stop
         else
            i = keven(k)
c           write (6,*) 'even grid ', k, ' at z ', zeven,
c           $           ' is on uneven grid interval ', 
c           $           i, sigh(i), sigh(i+1)
         endif
      end do

      write (6,*) ' initvt complete '

      end subroutine initvt

c##############################################################################

      subroutine getdyngrd( xwall, ywall )

c Return dynamics grid variables.

      implicit none

c Output arguments:
      real(DBLKIND)
     $  xwall(plon+1)       ! longitudes in cell walls (radians)
     $, ywall(plat+1)       ! latitudes in cell walls (radians)

c Local variables:
      integer
     $  i, j
c-----------------------------------------------------------------------
      do i = 1, plon+1
         xwall(i) = xw(i)
      end do
      do j = 1, plat+1
         ywall(j) = yw(j)
      end do

      end subroutine getdyngrd

c##############################################################################

      subroutine hadvn (den, deltat, n, denout, xflux, yflux)

c     horizontal advection algorithm
c     author: Phil Rasch
c     modified 28 December: modifications for mass vs sfc pres consistency

      implicit none

c     input
      real(r8) deltat                       ! unused, left for argument compatibility
      integer n                         ! time step index used for toggling operators
      
      real(r8) den(plond,plat,plev)         ! density on input

c     output
      real(r8) denout(plond,plat,plev)      ! density on output
      real(r8) yflux(plond,plat+1,plev)     ! fluxes on y cell walls
      real(r8) xflux(plond+1,plat,plev)     ! fluxes on x cell walls

c     local variables
      integer i
      integer j
      integer k
      real(r8) denx(plond,plat)             ! work variable
      real(r8) deny(plond,plat)             ! work variable
      real(r8) dent(plond,plat)             ! work variable
      real(r8) xtend(plond,plat)            ! work variable
      real(r8) ytend(plond,plat)            ! work variable
      real(r8) rdiff
      integer iprob, jprob
      integer ii

c     .....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
c     .... psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
c     .... velw1.....velw2.....velw3.....velw4.....velw5.....velw6
c     ..........xc1.......xc2........xc3......xc4........xc5.......
c     ..........dx1.......dx2........dx3......dx4........dx5.......
c     .........phi1......phi2.......phi3.....phi4.......phi5.......

      debugdt = deltat

! OpenMP directive
!$OMP PARALLEL DO PRIVATE( denx, deny, dent, xtend, ytend, iprob, jprob, rdiff )
      do k = 1,plev

#ifdef DEBUGSPITFIRE
         write (6,*) 'hadvn: range0 ', minval(den(:plon,:plat,k)), maxval(den(:plon,:plat,k))
#endif

c        toggle order of x and y operations each time step
         if (mod(n,2).eq.1) then

c           y advection
            call yupdate(den(1,1,k), deny, ytend, yflux(1,1,k), k)
#ifdef DEBUGSPITFIRE
            write (6,*) 'hadvn: range0.5 ', minval(deny(:plon,:plat)), maxval(deny(:plon,:plat))
#endif
            

            do j = 1,plat
               do i = 1,plon+1
c                 very bad advective approximation
                  dent(i,j) = deny(i,j) + den(i,j,k)
     $                 *deltat*(vw(i,j+1,k)-vw(i,j,k))/dmu(j)
                  if (dent(i,j).lt.0.) then
#ifdef DEBUGSPITFIRE
                     write (6,*) 'SPITFIRE bad adv app going negative ', i,j,k,den(i,j,k), deny(i,j), dent(i,j)
                     write (6,*) 'SPITFIRE vw (n, s) ', vw(i,j+1,k), vw(i,j,k)
                     write (6,*) 'SPITFIRE app cour num (n, s) ', vw(i,j+1,k)*deltat/dmu(j),
     $                                                    vw(i,j,k)*deltat/dmu(j)
                     write (6,*) 'SPITFIRE diff cour num ', deltat*(vw(i,j+1,k)-vw(i,j,k))/dmu(j)
#endif
                     dent(i,j) = max(0._r8,dent(i,j))
                  endif

               end do
            end do

#ifdef DEBUGSPITFIRE
            write (6,*) 'hadvn: range1 ', minval(dent(:plon,:plat)), maxval(dent(:plon,:plat))
#endif
c           x advection
            call xupdate(dent, denx, xtend, xflux(1,1,k), k)
#ifdef DEBUGSPITFIRE
            write (6,*) 'hadvn: range1.5 ', minval(denx(:plon,:plat)), maxval(denx(:plon,:plat))
            if  (minval(denx(:plon,:plat)).lt.0.) then
               do j = 1,plat
                  do i = 1,plon
                     if (denx(i,j).lt.0.) then
                        write (6,*) ' fluxes at an offending latitude ', i, j, denx(i,j),
     $                       minval(xflux(:plon+1,j,k)), maxval(xflux(:plon+1,j,k))
                        write (6,*) ' tinynum ', 2.*epsilon(1._r8)
                        do ii = 1,plon
                           write (6,*) ' flux ', ii, xflux(ii,j,k), xflux(ii+1,j,k)
                        end do
!                       call endrun()
                     end if
                  end do
               end do
            endif
#endif

            do j = 1,plat
               do i = 1,plon+1
                  denout(i,j,k) = deny(i,j) + xtend(i,j)
               end do
            end do
c REMOVE
c            write (6,*) 'hadvn odd: ', k, den(96,56,k), dent(96,56), denx(96,56), deny(96,56), denout(96,56,k)

#ifdef DEBUGSPITFIRE
            write (6,*) 'hadvn: range2 ', minval(denout(:plon,:plat,k)), maxval(denout(:plon,:plat,k))
#endif
         else

c           x advection
            call xupdate(den(1,1,k), denx, xtend, xflux(1,1,k), k)
#ifdef DEBUGSPITFIRE
            write (6,*) 'hadvn: range2.5 ', minval(denx(:plon,:plat)), maxval(denx(:plon,:plat))
#endif

            do j = 1,plat
               do i = 1,plon
c                 very bad advective approximation
                  dent(i,j) = denx(i,j) + den(i,j,k)
     $                 *deltat*(uw(i+1,j,k)-uw(i,j,k))/dxlon
                  if (dent(i,j).lt.0.) then
#ifdef DEBUGSPITFIRE
                     write (6,*) 'SPITFIRE bad adv app going negative ', i,j,k,den(i,j,k), denx(i,j), dent(i,j)
                     write (6,*) 'SPITFIRE vw (n, s) ', uw(i+1,j,k), uw(i,j,k)
                     write (6,*) 'SPITFIRE app cour num (e, w) ', uw(i+1,j,k)*deltat/dxlon,
     $                                                    uw(i,j,k)*deltat/dxlon
                     write (6,*) 'SPITFIRE diff cour num ', deltat*(uw(i+1,j,k)-uw(i,j,k))/dxlon
#endif
                  endif
                  dent(i,j) = max(0._r8,dent(i,j))
               end do
               dent(plon+1,j) = dent(1,j)
            end do

#ifdef DEBUGSPITFIRE
            write (6,*) 'hadvn: range3 ', minval(dent(:plon,:plat)), maxval(dent(:plon,:plat))
#endif

c           y advection
            call yupdate(dent, deny, ytend, yflux(1,1,k), k)
#ifdef DEBUGSPITFIRE
            write (6,*) 'hadvn: range3.5 ', minval(deny(:plon,:plat)), maxval(deny(:plon,:plat))
#endif

            do j = 1,plat
               do i = 1,plon+1
                  denout(i,j,k) = denx(i,j) + ytend(i,j)
               end do
            end do

#ifdef DEBUGSPITFIRE
            write (6,*) 'hadvn: range4 ', minval(denout(:plon,:plat,k)), maxval(denout(:plon,:plat,k))
#endif

         endif

      end do

      end subroutine hadvn

c##############################################################################

      subroutine hdep (deltat, n)

c     calculate velocity information for all tracers
      
c     calculate the departure points (uwdp, vwdp)
c     their intervals (idep, jdep)
c     and
c     whether they wrap outside the domain

c.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
c.... velw1.....velw2.....velw3.....velw4.....velw5.....velw6
c..........xc1.......xc2........xc3......xc4........xc5.......
c..........dx1.......dx2........dx3......dx4........dx5.......

      implicit none

c     input
      integer n
      real(r8) deltat

c     local

      integer i                         ! work variable
      integer j                         ! work variable
      integer k                         ! work variable

      real(DBLKIND) vwext(platext+5)             ! work variable

! OpenMP directive
!$OMP PARALLEL DO PRIVATE( vwext )
      do k = 1,plev


c        take the wall centered winds, and calculate their departure interval
c        and increment

c        first do the y cell walls
         do i = 1,plon/2

c           put v winds into the reshaped arrays
CDIR$ INLINE
            call getys2 (vw(1,1,k), vwext, platext, i)

c           calculate the departure point for particles located at each cell and
c           the interval it occurs in
#ifdef DEBUGPOINT
            if (i.eq.ilook.and.k.eq.klook) then
               debugit = 1
            endif
#endif
CDIR$ INLINE
            call getf2 (platext, ywext, vwext, deltat,
     $           2, jeven, njeven, 
     $           jdep(1,i,k), jdom(1,i,k), vwdp(1,i,k))

#ifdef DEBUGPOINT
            debugit = 0

            if (i.eq.ilook.and.k.eq.klook) then
               write (6,*) ' hdep: ilook, jlook,klook vw, vwext, jdep, jdom, vwdp', 
     $              ilook, jlook, klook,
     $              vw(ilook,jlook,klook),  jdep(jlook,ilook,klook),  jdom(jlook,ilook,klook),
     $              vwdp(jlook,ilook,klook)
!              write (6,*) ' vwext ', vwext(1:20)
            endif
#endif

         end do

c        now do the x cell walls
         do j = 1,plat

c           calculate the departure point for particles located at each cell and
c           the interval it occurs in
CDIR$ INLINE
            call getf2 (plon+1, xw, uw(1,j,k), deltat,
     $           1, jeven, njeven, 
     $           idep(1,j,k), idom(1,j,k), uwdp(1,j,k))

         end do

      end do

#ifdef DEBUGPOINT
      i = 1
      j = 1
      do k = 1,plev
         write (6,*) ' hdep: advect: uw, w,e ', k, uw(i+1,j,k), uw(i,j,k)
         write (6,*) ' hdep: advect: uwdp, w,e ', k, uwdp(i+1,j,k), uwdp(i,j,k)
      end do
#endif

      end subroutine hdep

c##############################################################################

      subroutine hwctw (u3, v3, deltat, n)

c     move horizontal winds from cell centers to cell walls
      

c.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
c.... velw1.....velw2.....velw3.....velw4.....velw5.....velw6
c..........xc1.......xc2........xc3......xc4........xc5.......
c..........dx1.......dx2........dx3......dx4........dx5.......

      implicit none

c     input
      real(r8) u3(plond,plevd,plat)
      real(r8) v3(plond,plevd,plat)
      integer n
      real(r8) deltat

c     local

      integer i                         ! work variable
      integer j                         ! work variable
      integer k                         ! work variable

      real(DBLKIND) c1, c2, c3                   ! work variable
      real(r8) cour, cour2, cmax, dcour
      integer imax, jmax, kmax

      c1 = 0.5/arad 

! OpenMP directive
!$OMP PARALLEL DO PRIVATE( c2, c3 )
      do k = 1,plev

c        map from cell centered velocities on ccm grid to cell wall velocities
c        convert from traditional velocities to lamdadot and mudot (mu=sin(phi))
c        or      from traditional velocities to lamdadot and phidot 
         do j = 1,plat
            c2 = cos(yw(j+1))*c1
            c3 = c1/cos((yw(j+1)+yw(j))*0.5)
            do i = 1,plon
               uw(i+1,j,k) = (u3(i,k,j)+u3(i+1,k,j))*c3
               if (j.lt.plat) then
#ifdef FLUXFORM2  
                  vw(i,j+1,k) = (v3(i,k,j)+v3(i,k,j+1))*c1
#else             
                  vw(i,j+1,k) = (v3(i,k,j)+v3(i,k,j+1))*c2
#endif            
               endif
            end do
            uw(1,j,k) =  (u3(1,k,j)+u3(plon,k,j))*c3
            uw(plon+1,j,k) = uw(1,j,k)
            vw(plon+1,j,k) = vw(1,j,k)
            do i = 1,plon+1
               vw(i,1,k) = 0.
               vw(i,plat+1,k) = 0.
            end do
         end do

      end do

!     now look for courant number problems
      imax = 1
      jmax = 1
      kmax = 1
      cmax = 0.
      do k = 1,plev
         do j = 1,plat
            do i = 1,plon
               cour = uw(i,j,k)*deltat/dxlon
               cour2 = uw(i+1,j,k)*deltat/dxlon	 
               dcour = abs(cour-cour2)
               if (abs(dcour).gt.cmax) then
                  imax = i
                  jmax = j
                  kmax = k
                  cmax = dcour
               endif
            end do
         end do
      end do
      if (cmax > 1.0) then
         write (6,*) ' EXPECT NEGATIVES zonal differential courant number exceed unity --- max val and ijk indices are ', 
     $        cmax, imax, jmax, kmax
      endif

      end subroutine hwctw

c##############################################################################

      subroutine cfint1x (x, f, fdot, ndata, npts, psistar, 
     $     jdep2, vwdp2)

c     interpolate a function specified by (x, f, fdot) at points vwdp2
c     return values of psistar

c     assume periodic in the first dimension

      implicit none

c     input
      integer ndata
      real(DBLKIND) x(ndata)                     ! coords of function
      real(DBLKIND) f(ndata)                     ! value of function at x
      real(DBLKIND) fdot(ndata)                  ! deriv of function at x

      integer npts                      ! number of points to interpolate to
      integer jdep2(npts)               ! interval where each interpolation is done
      real(DBLKIND) vwdp2(npts)                  ! coord to construct the interpolant for

c     output
      real(DBLKIND) psistar(npts)                ! interpolated value

c     local workspace
      integer i                         ! work variable
      integer ii                        ! work variable
      real(DBLKIND) s                            ! work variable
      real(DBLKIND) c2                           ! work variable
      real(DBLKIND) c3                           ! work variable
      real(DBLKIND) xx                           ! work variable
#ifdef LIMITON
      real(DBLKIND) psi1, psi2, psi3, psim
      real(DBLKIND) cfnew
#endif
      real(DBLKIND) cfint
      integer nwork
      parameter (nwork = NWORK)
      real(DBLKIND) fext(-4:nwork)
      real(DBLKIND) dx(nwork)
      integer ip1
      integer im1
      integer next

c     the minmod function 
      real(DBLKIND) a, b, c
      real(DBLKIND) minmod
      real(DBLKIND) medan
      minmod(a,b) = 0.5*(sign(one,a) + sign(one,b))*min(abs(a),abs(b))
      medan(a,b,c) = a + minmod(b-a,c-a)


c     construct an extended version of the grid
      do ii = 1,ndata-1
         fext(ii) = f(ii)
         dx(ii) = (x(ii+1)-x(ii))
      end do
      fext(ndata) = f(ndata)
      dx(ndata) = dx(1)

      next = 3
cdir$ ivdep
      do ii = 0,next
         fext(ii+ndata+1) = f(ii+2)  + (f(ndata)-f(1))
         fext(-ii)     = f(ndata-ii-1) - (f(ndata)-f(1))
      end do

c     now do the interpolation at each point
      do ii = 1,npts

c        i is the departure point interval
         i = jdep2(ii)
         xx = (vwdp2(ii)-x(i))
         ip1 = mod(ndata+i-2,ndata-1) + 2
c        ip2 = mod(ndata+i-1,ndata-1) + 2
         im1 = mod(ndata+i-4,ndata-1) + 2

c        now interpolate
         s = (f(ip1)-f(i))/dx(i)
         c2 = (3*s-2*fdot(i)-fdot(ip1))/dx(i)
         c3 = (fdot(i)+fdot(ip1)-2*s)/dx(i)**2
c        fxdot =  (3*c3*xx + 2*c2)*xx + fdot(i)
c        fxdd = 6*c3*xx + 2*c2
         cfint = ((c3*xx + c2)*xx + fdot(i))*xx + f(i)

#ifdef LIMITON
c        leonards limiter
         psi1 = fext(i)   + (fext(i+1)-fext(i))*xx/dx(i)
         psi2 = fext(i)   + (fext(i)-fext(i-1))*xx/dx(im1)
         psi3 = fext(i+1) - (fext(i+2)-fext(i+1))*(dx(i)-xx)/dx(ip1)
         psim = medan(psi1, psi2, psi3)
         cfnew = medan(cfint, psi1, psim)
         cfint = cfnew
#endif   
!-drb         cfint = medan(cfint,f(i),f(ip1))
         cfint = max(min(f(i),f(i+1)),cfint)
         cfint = min(max(f(i),f(i+1)),cfint)

         psistar(ii) = cfint

      end do

      end subroutine cfint1x

c##############################################################################

      subroutine cfdot1dp (x, f, npts, fdot)

c     calculate the derivative for the interpolating polynomial
c     periodic BCs
c     for extending the values,
c     we assume x(1)-x(0) = x(npts)-x(npts-1),
c     f(1)-f(0) = f(npts)-f(npts-1)

      implicit none

c     input
      integer npts
      real(DBLKIND) x(npts)                      ! nodes where data are
      real(DBLKIND) f(npts)                      ! values of data

c     output
      real(DBLKIND) fdot(npts)                   ! derivative of f at x

c     assumed variable distribution
c     x1.......x2.......x3.......x4.......x5.......x6     1,npts points
c     f1.......f2.......f3.......f4.......f5.......f6     1,npts points
c     ...sh1.......sh2......sh3......sh4......sh5....     1,npm1 points
c     .........d2.......d3.......d4.......d5.........     2,npm1 points
c     .........s2.......s3.......s4.......s5.........     2,npm1 points
c     .............dh2......dh3......dh4.............     2,npm2 points
c     .............eh2......eh3......eh4.............     2,npm2 points
c     ..................e3.......e4..................     3,npm2 points
c     .................ppl3......ppl4................     3,npm2 points
c     .................ppr3......ppr4................     3,npm2 points
c     .................t3........t4..................     3,npm2 points
c     ................fdot3.....fdot4................     3,npm2 points


c     work variables

      integer i
      integer npm1
      integer npm2
      integer nwork

      parameter (nwork = NWORK)
      real(DBLKIND) s(nwork)                     ! first divided differences at nodes
      real(DBLKIND) sh(nwork)                    ! first divided differences between nodes
      real(DBLKIND) d(nwork)                     ! second divided differences at nodes
      real(DBLKIND) dh(nwork)                    ! second divided differences between nodes
      real(DBLKIND) e(nwork)                     ! third divided differences at nodes
      real(DBLKIND) eh(nwork)                    ! third divided differences between nodes
      real(DBLKIND) pp                           ! p prime
      real(DBLKIND) ppl                          ! p prime on left
      real(DBLKIND) ppr                          ! p prime on right
      real(DBLKIND) qpl
      real(DBLKIND) qpr
      real(DBLKIND) ttt
      real(DBLKIND) t
      real(DBLKIND) tmin
      real(DBLKIND) tmax

      integer ii, im1, ip1
      real(DBLKIND) dxh(nwork)
      real(DBLKIND) dx2(nwork)
      real(DBLKIND) dx3h(nwork)
      real(DBLKIND) fd1, fd2

c     the minmod function 
      real(DBLKIND) minmod
      real(DBLKIND) medan
      real(DBLKIND) a, b, c
      minmod(a,b) = 0.5*(sign(one,a) + sign(one,b))*min(abs(a),abs(b))
      medan(a,b,c) = a + minmod(b-a,c-a)

      npm1 = npts-1
      npm2 = npts-2

c     a bunch of divided differences are related functions
      do i = 1,npm1
         dxh(i) = x(i+1)-x(i)
         sh(i) = (f(i+1)-f(i))/dxh(i)
      end do

      do i = 2,npm1
         dx2(i) = (dxh(i)+dxh(i-1))     ! equals x(i+1)-x(i-1)
         d(i) = (sh(i)-sh(i-1))/dx2(i)
         s(i) = minmod(sh(i),sh(i-1))
      end do

      do i = 2,npm2
         dx3h(i) = (x(i+2)-x(i-1))
         eh(i) = (d(i+1)-d(i))/dx3h(i)
         dh(i) = minmod(d(i),d(i+1))
      end do

      do i = 3,npm2
         e(i) = minmod(eh(i),eh(i-1))
      end do

      sh(npts) = sh(1)
      dxh(npts) = dxh(1)
      dx2(1) = dxh(1) + dxh(npm1)
      dx2(npts) = dx2(1)
      dx3h(npts) = dxh(npm1) + dxh(npts) + dxh(2)
      dx3h(1) = dx3h(npts)
      dx3h(npm1) = dxh(npts-2) + dxh(npm1) + dxh(npts)
      d(npts) = (sh(npts)-sh(npm1))/dx2(npts)
      d(1) = d(npts)
      s(npts) = minmod(sh(npts),sh(npm1))
      s(1) = s(npts)
      eh(npm1) = (d(npts)-d(npm1))/dx3h(npm1)
      dh(npm1) = minmod(d(npm1),d(npts))
      eh(npts) = (d(2)-d(npts))/dx3h(npts)
      eh(1) = eh(npts)
      dh(npts) = minmod(d(npts),d(2))
      dh(1) = dh(npts)
      e(2) = minmod(eh(2), eh(1))
      e(1) = minmod(eh(1), eh(npm1))
      e(npts) = e(1)
      e(npm1) = minmod(eh(npm1),eh(npts-2))

c     now calculate the derivatives in interior
      do i = 3,npm2

         im1 = i-1
         ip1 = i+1

c        p prime on left side (what huynh calls i-0.5)
         ppl=sh(im1) + dh(im1)*dxh(im1)
c        p prime on right side (what huynh calls i+0.5)
         ppr=sh(i)   - dh(i)  *dxh(i)
         t = minmod(ppl,ppr)
c        tmax = sign(one,t)*max(3.*abs(s(i)),1.5*abs(t))
c        fdot(i) = 0.5*(ppl+ppr) !  derivative M3-A
c        fdot(i) = minmod(fdot(i),tmax)


c        derivative from parabola thru f(im1), f(i), and f(ip1)
         pp = sh(im1) + d(i)*dxh(im1)
         fd1 = pp - dxh(im1)*dxh(i)
     $        *(eh(im1)*dx2(ip1) + eh(i)*dx2(im1))
     $        /(dx2(ip1)+dx2(im1))      !(x(i+2)-x(i-2))
c        $                 /(x(i+2)-x(i-2))

#ifdef MONO
         qpl = sh(im1) 
     $        + (dxh(im1))*minmod(d(im1) + e(im1)*dx2(im1),
     $        d(i)   - e(i)  *dxh(i))
         qpr = sh(i)
     $        + (dxh(i))  *minmod(d(i)   + e(i)  *dxh(i),
     $        d(ip1) - e(ip1)*dx2(ip1))
         ttt = minmod(qpl, qpr)
         tmin = min(zero,3*s(i),1.5*t,ttt)
         tmax = max(zero,3*s(i),1.5*t,ttt)
         fd2 = medan(fd1, qpl, qpr)
         fdot(i) = fd2 + minmod(tmin-fd2, tmax-fd2)
#else    
         fdot(i) = fd1
#endif   

      end do

      do ii = -1, 2
         i = mod(npts+ii-1,npts) + 1
         im1 = mod(npts+i-3,npts-1) + 1
         ip1 = mod(npts+i-1,npts-1) + 1

c        p prime on left side (what huynh calls i-0.5)
         ppl=sh(im1) + dh(im1)*dxh(im1)
c        p prime on right side (what huynh calls i+0.5)
         ppr=sh(i)   - dh(i)  *dxh(i)
         t = minmod(ppl,ppr)
c        tmax = sign(one,t)*max(3.*abs(s(i)),1.5*abs(t))
c        fdot(i) = 0.5*(ppl+ppr) !  derivative M3-A
c        fdot(i) = minmod(fdot(i),tmax)


c        derivative from parabola thru f(im1), f(i), and f(ip1)
         pp = sh(im1) + d(i)*dxh(im1)
         fd1 = pp - dxh(im1)*dxh(i)
     $        *(eh(im1)*dx2(ip1) + eh(i)*dx2(im1))
     $        /(dx2(ip1)+dx2(im1))      !(x(i+2)-x(i-2))

#ifdef MONO
         qpl = sh(im1) 
     $        + (dxh(im1))*minmod(d(im1) + e(im1)*dx2(im1),
     $        d(i)   - e(i)  *dxh(i))
         qpr = sh(i)
     $        + (dxh(i))  *minmod(d(i)   + e(i)  *dxh(i),
     $        d(ip1) - e(ip1)*dx2(ip1))
         ttt = minmod(qpl, qpr)
         tmin = min(zero,3*s(i),1.5*t,ttt)
         tmax = max(zero,3*s(i),1.5*t,ttt)
         fd2 = medan(fd1, qpl, qpr)
         fdot(i) = fd2 + minmod(tmin-fd2, tmax-fd2)
#else    
         fdot(i) = fd1
#endif   

      end do

      end subroutine cfdot1dp


c##############################################################################

      subroutine getflxa (npts, xw, psi, flux,
     $     jdep2, jdom2,vwdp2)


c     calculate the horizontal mass fluxes at cell walls

c     xw are coords of cell walls
c     psi is integral of mean tracer densities at cell walls

c     assume number of cell centers is npts-1 
c     number of cell walls is npts

c     tracers are updated from the output from this code as
c     do i = 1,npts-1
c     q(i) = q(i) - (flux(i+1)-flux(i))/(xw(i+1)-xw(i))
c     end do

c.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
c.... psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
c.... velw1.....velw2.....velw3.....velw4.....velw5.....velw6
c.........phi1......phi2.......phi3.....phi4.......phi5.......

      implicit none

c     input
      integer npts                      ! number of points
      real(DBLKIND) xw(npts)                     ! cell walls
      real(DBLKIND) psi(npts)                    ! integral of densities
      integer jdep2(npts)               ! departure point interval
      integer jdom2(npts)               ! whether the departure point wrapped the domain
      real(DBLKIND) vwdp2(npts)                  ! departure point

c     output
      real(r8) flux(npts)                   ! flux over the time step


c     local vars

      integer i

      integer nwork
      parameter (nwork=NWORK)
      real(DBLKIND) fdot(nwork)                  ! derivative of psi
      real(DBLKIND) psistar(nwork)               ! value of interpolated psi at dep points
      real(DBLKIND) omsm
      real(r8) w1, w2

      omsm  = 1.-2*epsilon(one)


c     calculate derivatives
CDIR$ INLINE
      call cfdot1dp (xw, psi, npts, fdot)

c     interpolate to departure points
CDIR$ INLINE
      call cfint1x (xw, psi, fdot, npts, 
     $     npts, psistar,
     $     jdep2, vwdp2)

#ifdef DEBUGSPITFIRE
! sanity checks
      do i = 1,npts-1
         w1 = xw(i+1)-xw(i)
         if (w1.le.0.) then
            write (6,*) 'SPITFIRE  xw is not monotonic ', i, xw(i), xw(i+1)
            call endrun()
         endif
         w1 = psi(i+1)-psi(i)
         if (w1.lt.0.) then
            write(6,*) 'SPITFIRE  psi is not monotonic ', i, psi(i), psi(i+1)
            call endrun()
         endif
         w1 = vwdp2(i) - xw(jdep2(i))
         if (w1.lt.0.) then
            write (6,*) 'SPITFIRE  vwdp2 check 1 ', i, jdep2(i), vwdp2(i), xw(jdep2(i))
            call endrun()
         endif
         w1 = xw(jdep2(i)+1) - vwdp2(i)  
         if (w1.lt.0.) then
            write (6,*) 'SPITFIRE  vwdp2 check 2 ', i, jdep2(i), vwdp2(i), xw(jdep2(i)+1)
            call endrun()
         endif
         w1 = psistar(i) - psi(jdep2(i))
         if (w1.lt.0.) then
            write (6,*) 'SPITFIRE  psistar check 1 ', i, jdep2(i), psistar(i), psi(jdep2(i))
!-drb            call endrun()
         endif
         w1 = psi(jdep2(i)+1) - psistar(i)  
         if (w1.lt.0.) then
            write (6,*) 'SPITFIRE  psistar check 2 at wall point', i, xw(i)
            write (6,*) ' departure point is ', vwdp2(i)
            write (6,*) ' should be on interval ', jdep2(i)
            write (6,*) ' between walls ', xw(jdep2(i)),' and ', xw(jdep2(i)+1)
            write (6,*) ' the interpolated value is ', psistar(i)
            write (6,*) ' which should be between ', psi(jdep2(i)),' and ', psi(jdep2(i)+1)
!-drb            call endrun()
         endif
      end do
#endif

c     fluxes
      do i = 1,npts
         flux(i) = (psi(i)-psistar(i)) 
     $        - jdom2(i)*psi(npts)      ! periodicity correction
         flux(i) = flux(i)*omsm
      end do
      
      end subroutine getflxa

c##############################################################################

      subroutine getf2 (npts, xw, vel, deltat, 
     $     depflag, jeven, njeven, idep, idom, xwd)

c     calculate the wall departure point and intervals

c     xw are coords of cell walls
c     vel are velocites at cell walls
c     deltat is time step

c     assume number of cell centers is npts-1 
c     number of cell walls is npts


c     .....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
c     ....velw1.....velw2.....velw3.....velw4.....velw5.....velw6

      implicit none

c     input
      integer npts                      ! number of cell walls
      real(DBLKIND) vel(npts)                    ! velocity at cell walls
      real(DBLKIND) xw(npts)                     ! coords of cell walls
      integer depflag                   ! even or uneven grid flag
      integer njeven                    ! number of points in even grid
      integer jeven(*)                  ! even to uneven grid mapping index
      real(r8) deltat                       ! time step

c     output
      integer idep(npts)                ! departure point interval
      integer idom(npts)                ! flag the wrapping of dep points
      real(DBLKIND) xwd(npts)                    ! departure point

c     local
      integer i
      integer ii
      integer j

      real(DBLKIND) xint
      real(DBLKIND) tmp
      real(DBLKIND) dxm

c     xwd is the departure point for a 
c     an arrival point originating at the wall

c     idom flags whether the departure point is outside and 
c     greater than the max value of the interval (1), or
c     less than the min value of the interval (-1)
      xint = xw(npts)-xw(1)
      do i = 1,npts
         xwd(i) = xw(i)-vel(i)*deltat
         idom(i) = 0
         if (xwd(i).lt.xw(1)) then
            idom(i) = -1
            tmp = xwd(i)+xint
c           xwd(i) = tmp
            xwd(i) = min(tmp,xw(npts))
         endif
         if (xwd(i).gt.xw(npts)) then
            idom(i) = 1
            tmp = xwd(i)-xint
c           xwd(i) = tmp
            xwd(i) = max(tmp,xw(1))
         endif
c        write (6,*) ' xw, xwd, vel, idom ', i, xw(i), xwd(i), vel(i), idom(i)
      end do

      if (depflag .eq.1) then
c        guess the interval assuming x is distributed uniformly
         dxm = xint/(npts-1)
         do ii = 1,npts-1
c           index i is the interval of the departure point
c            i = int((xwd(ii)-xw(1) - tinynum)/dxm)+1
            i = int(((xwd(ii)-xw(1))*(1. - tinynum))/dxm)+1
            if (i.ge.npts .or. i.lt.1) then
               write (6,*) ' something is fishy ', npts, i, ii, xw(1), dxm
               write (6,*) ' xw(ii), xwd(i) ', xw(ii), xwd(ii)
               write (6,*) ' before int calc ', (xwd(ii)-xw(1))/dxm
               write (6,*) ' tinynum ', tinynum
               tmp = ((xwd(ii)-xw(1))*(1. - tinynum))/dxm
               write (6,*) ' float of i including tinynum ', tmp
               write (6,*) ' int of float of i including tinynum+1 ', 
     $              int(tmp)+1;
               stop
            endif
            if ((xwd(ii)-xw(i))*(xw(i+1)-xwd(ii)).ge.0) then
               idep(ii) = i
            else if (xwd(ii).lt.xw(i)) then
               idep(ii) = i-1
            else if (xwd(ii).gt.xw(i+1)) then
               idep(ii) = i+1
            else
               write (6,*) ' spitfire: should never get here ', npts, i, ii, xw(1), dxm
               write (6,*) ' xw(ii), xwd(i) ', xw(ii), xwd(ii)
               write (6,*) ' before int calc ', (xwd(ii)-xw(1))/dxm
               stop
            endif
         end do
         idep(npts) = idep(1)
      else if (depflag.eq.2) then
c        guess the interval using a precalculated mapping from an even distribution to 
c        the actual grid
         dxm = xint/(njeven-1)
         do ii = 1,npts-1
            j = int((xwd(ii)-xw(1))/dxm)+1
            i = jeven(j)
            if ((xwd(ii)-xw(i))*(xw(i+1)-xwd(ii)).ge.0) then
               idep(ii) = i
            else if (xwd(ii).lt.xw(i)) then
               idep(ii) = i-1
            else if (xwd(ii).gt.xw(i+1)) then
               idep(ii) = i+1
            endif
         end do
         idep(npts) = idep(1)
      else
         write (6,*) ' cfint1x: illegal value for depflag ', depflag
         stop
      endif

c     now check to see whether we actually found it
      j = 0
      do ii = 1,npts-1
         i = idep(ii)
         if ((xwd(ii)-xw(i))*(xw(i+1)-xwd(ii)).lt.0) then
            j = ii
         endif
#ifdef DEBUGPOINT
         if (debugit .eq. 1) then
            write (6,*) ' wall point ', ii, xw(ii)
            write (6,*) ' has departure point ', ii, xwd(ii)
            write (6,*) ' which should be in interval ', i,
     $           ' defined by endpoints ', xw(i), xw(i+1)
            write (6,*) ' implied and true velocity is ', 
     $           (xw(ii) - xwd(ii))/deltat, vel(ii)
         endif
#endif
      end do
      if (j.ne.0) then
         ii = j
         i = idep(ii)
         write (6,*) ' cfint1x: somethings wrong at dep calc ', ii, i
         write (6,*) ' xwd(ii) ', xwd(ii)
         write (6,*) ' we thought it should be between ',
     $        xw(i),' and ', xw(i+1)
         stop 11
      endif

      end subroutine getf2

c##############################################################################
      
      subroutine getyslice (den, denext, platext, i)

c     move the densities from a normal lon/lat grid to the rotated grid

      implicit none

c     input
      integer platext                   ! number of points on rotated grid
      real(r8) den(plond, plat)             ! input densities
      integer i                         ! longitude to extrace

c     output
      real(r8) denext(platext+5)            ! rotated densities

c     local
      integer ii 
      integer j

      ii = mod(i+plon/2-1,plon) + 1     ! x index on other side of globe
cdir$ ivdep
      do j = 1,plat
         denext(j) = den(i,j)
         denext(2*plat-j+1) = den(ii,j)
      end do

      denext(2*plat+1) = denext(1)

      end subroutine getyslice

c##############################################################################

      subroutine getys2 (vw, vwext, platext, i)

c     move the velocities from normal lon/lat grid to the rotated grid

      implicit none

c     input
      integer platext
      real(DBLKIND) vw(plond, plat+1)
      integer i

c     output
      real(DBLKIND) vwext(platext+5)

c     local
      integer ii 
      integer j

      ii = mod(i+plon/2-1,plon) + 1

cdir$ ivdep
      do j = 1,plat
         vwext(j) = vw(i,j)
         vwext(2*plat-j+2) = -vw(ii,j)
      end do

      vwext(plat+1) = vw(i,plat+1)
      vwext(2*plat+1) = vwext(1)

      end subroutine getys2

c##############################################################################

      subroutine putyslice(yfluxext, yflux, platext, i)

c     map the fluxes from the rotated grid back to the lon/lat grid

      implicit none 

c     input
      integer i
      integer platext
      real(r8) yfluxext(platext+5)

c     output
      real(r8) yflux(plond,plat+1)

c     local
      integer ii
      integer j

      ii = mod(i+plon/2-1,plon) + 1

      do j = 2,plat
         yflux(i,j) = yfluxext(j)  
         yflux(ii,j) = -yfluxext(2*plat-j+2) 
      end do
      yflux(i,plat+1) =  yfluxext(plat+1) 
      yflux(ii,plat+1) = -yfluxext(plat+1)
      yflux(i,1) =  yfluxext(1) 
      yflux(ii,1) =  -yfluxext(1) 

      end subroutine putyslice

c##############################################################################

      subroutine yupdate(den, deny, ytend, yflux, k)

c     calculate the new densities after advection in y direction

      implicit none

c     input
      real(r8) den(plond,plat)              ! input density
      integer k                         ! lev to advect

c     output
      real(r8) deny(plond,plat)             ! output density
      real(r8) ytend(plond,plat)
      real(r8) yflux(plond,plat+1)

c     local
      integer i, ii
      integer j, jj
      real(r8) yfluxext(platext+5)
      real(r8) den2(platext+5,plon/2)
      real(DBLKIND) deni(platext+5,plon/2)
      real(r8) dena(platext+5,plon/2)
      real(r8) denar(plond,plat)            
      real(r8) tmp

      do i = 1,plon/2

c     move densities to rotated grid
CDIR$ INLINE
         call getyslice (den, den2(1,i), platext, i)

      end do

c     compute the y integral of the densities
      do i = 1,plon/2
         deni(1,i) = 0.
      end do
      do j = 2,platext
         do i = 1,plon/2
            deni(j,i) = deni(j-1,i) + den2(j-1,i)*(ywext(j)-ywext(j-1))
            dena(j-1,i) = (deni(j,i)-deni(j-1,i))/(ywext(j)-ywext(j-1)) ! this is the density calculated from the 
                                                                        ! integral function so it will not generate
                                                                        ! negatives, but will have precision loss
         end do
      end do

c     calculate fluxes
      do i = 1,plon/2
CDIR$ INLINE
         call getflxa(platext, ywext, deni(1,i),
     $        yfluxext,
     $        jdep(1,i,k), jdom(1,i,k),vwdp(1,i,k))


c        move fluxes back to regular lon/lat grid
CDIR$ INLINE
         call putyslice(yfluxext, yflux, platext, i)

      end do

!     unscramble the approximate densities
      do j = 1,plat
         do i = 1,plon/2
            ii = mod(i+plon/2-1,plon) + 1
            jj = 2*plat-j+1
            denar(i,j) = dena(j,i)
            denar(ii,j) = dena(jj,i)
         end do
      end do

c     add wrap point
cdir$ ivdep
      do j = 1,plat+1
         yflux(plon+1,j) = yflux(1,j)
#ifdef FLUXFORM2
         do i = 1,plon+1
            yflux(i,j) = yflux(i,j)*cos(yw(j))
         end do
#endif   
      end do


c     calculate tendencies and update densities
      do j = 1,plat

         do i = 1,plon
            tmp = - (yflux(i,j+1)- yflux(i,j))/dmu(j)
            ytend(i,j) = tmp
!           deny(i,j) = den(i,j) + tmp
            deny(i,j) = denar(i,j) + tmp
         end do

         deny(plon+1,j) = deny(1,j)
         ytend(plon+1,j) = ytend(1,j)

      end do

      end subroutine yupdate

c##############################################################################

      subroutine xupdate(den, denx, xtend, xflux, k)

c     calculate the new densities after advection in x direction

      implicit none

c     input
      real(r8) den(plond,plat)
      integer k

c     output
      real(r8) denx(plond,plat)
      real(r8) xtend(plond,plat)
      real(r8) xflux(plond+1,plat)

c     local
      integer i
      integer j

      real(DBLKIND) deni(plond,plat)
      real(DBLKIND) tmp
      real(r8) fluxmax
      integer ii

c     construct x update and calculate tendencies

c     calculate the x integral of the density
      do j = 1,plat
         deni(1,j) = 0.
      end do

      do i = 2,plon+1
         do j = 1,plat
            deni(i,j) = deni(i-1,j) + den(i-1,j)*(xw(i)-xw(i-1))
         end do
      end do

c     calculate the x fluxes

      do j = 1, plat
CDIR$ INLINE
         call getflxa(plon+1, xw, deni(1,j),
     $        xflux(1,j),
     $        idep(1,j,k), idom(1,j,k),uwdp(1,j,k))


         xflux(plon+1,j) = xflux(1,j)
#ifdef OFFIT
!        make a pass through the fluxes zeroing out the ones that are zero to roundoff
         fluxmax = 0.
         do i = 1,plon+1
            fluxmax = abs(xflux(i,j))
         end do
         fluxmax = fluxmax*tinynum

         do i = 1,plon+1
            if (abs(xflux(i,j)).lt.fluxmax) then
               xflux(i,j) = 0.
            endif
         end do
#endif

      end do

c     update the fields
      do j = 1,plat
         do i = 1,plon
            tmp = - (xflux(i+1,j) - xflux(i,j))/dxlon
            xtend(i,j) = tmp
!normal            denx(i,j) = den(i,j) + tmp
            denx(i,j) = (deni(i+1,j)-deni(i,j))/(xw(i+1)-xw(i)) + tmp
            if (abs(denx(i,j)).lt.abs(den(i,j))*tinyr8) denx(i,j) = 0._r8
#ifdef DEBUGSPITFIRE
            if (denx(i,j).lt.0._r8) then
               write (6,*) 'SPITFIRE  xupdate: problem ', i, j, k, denx(i,j), xflux(i+1,j), xflux(i,j), dx
               write (6,*) ' old density ', den(i,j)
               write (6,*) ' uw at right and left walls ', uw(i+1,j,k), uw(i,j,k)
               write (6,*) ' xw at right and left walls ', xw(i+1), xw(i)
               write (6,*) ' uwdp at right and left walls ', uwdp(i+1,j,k), uwdp(i,j,k)
               write (6,*) ' stupid calculation of departure point ', xw(i+1)-uw(i+1,j,k)*debugdt,
     $                     xw(i)-uw(i,j,k)*debugdt
               write (6,*) ' some controlling fields = i, den, deni, uwdp, idep, idom '
               do ii = 1,plon
                  write (6,*) ii, den(ii,j), deni(ii,j), uwdp(ii,j,k), idep(ii,j,k), idom(ii,j,k)
               end do
!               write (6,*) 'SPITFIRE debug calling abort'
!               call endrun()
            endif
#endif
         end do
         denx(plon+1,j) = denx(1,j)
         xtend(plon+1,j) = xtend(1,j)
      end do
      
      end subroutine xupdate

c##############################################################################

      subroutine inifilter()

      implicit none

      integer i, j, n

      real(DBLKIND) reduce
      integer npole
      integer nfilt
      integer nmin
      integer pow2(16)
c      data pow2 /1, 2, 4, 8, 16, 32, 64, 128 /
 
      reduce = .3  ! factor by which to reduce the stencil (0 (max reduction) to 1. (standar))

      do n = 1,16
         pow2(n) = 2**(n-1)
      end do

      do j = 1,plat
         nptslat(j) = max(1,min(plon,
     $                    int(reduce/(cos(ylat(j))+1.e-10))))
c        write (6,*) ' deb1 nptslat ', j, nptslat(j)
         if (nptslat(j).gt.1) then
            npole = abs(nptslat(j)-pow2(1))
            nmin = 1
            do n = 1,16
               if (abs(nptslat(j)-pow2(n)).le.npole) then
                  npole = abs(nptslat(j)-pow2(n))
                  nmin = n
               endif
            end do
            nptslat(j) = pow2(nmin)
         endif
      end do

      nptslat(1) = plon
      nptslat(plat) = plon
      do j = 1,plat
         nfilt = nptslat(j)
         do i = 1,plon
            bini(i,j) = (i+nfilt-1)/nfilt
         end do
         if (nptslat(j).ne.1) then
         write (6,*) nptslat(j), ' points in filter at lat ', j
         endif
      end do
 
      end subroutine inifilter

c##############################################################################

      subroutine filter (field)

      implicit none

      real(r8) field(plond,plat)
      integer ifirst
      real(r8) t1

      integer i
      integer j

      real(r8) temp(plond)

c     update solution
      
      do j = 1,plat
         if (nptslat(j).gt.1) then
            do i = 1,plon
               temp(i) = 0.
            end do
            do i = 1,plon
               temp(bini(i,j)) = temp(bini(i,j)) + field(i,j)
            end do
            do i = 1,plon
               field(i,j) = temp(bini(i,j))/nptslat(j)
            end do
            field(plon+1,j) = field(1,j)
         endif
      end do

      end subroutine filter

c##############################################################################

      subroutine vadvn (dn, dt, dnout, flux)

      implicit none

c     input
      real(r8) dn(plond, plat, plev)
      real(r8) dt

c     output
      real(r8) dnout(plond, plat, plev)
      real(r8) flux(plond,plevp,plat)

c     local
      real(r8) dnxz(plond,plev)
      real(DBLKIND) psi(plond,plevp)
      integer i
      integer k
      integer lat
      real(DBLKIND) tvcour, bvcour
      real(DBLKIND) rdiff
      integer iprob, jprob

! OpenMP directive
!$OMP PARALLEL DO PRIVATE( dnxz, iprob, jprob, rdiff , psi, tvcour, bvcour)
      do lat=1,plat

c        move 3d fields to standard lon/height slice
         do k = 1,plev
            do i = 1,plon
               dnxz(i,k) = dn(i,lat,k)
            end do
            dnxz(plon+1,k) = dnxz(1,k)
         end do

c        calculate the mass moving past each face
         call getflxv(zw, dnxz, flux(1,1,lat),
     $        wwdp(1,1,lat), kdep(1,1,lat), psi)

c        update the species array for vertical transport
         do k = 1,plev
            do i = 1,plon      
               dnout(i,lat,k) = dn(i,lat,k) 
     $              - (flux(i,k+1,lat)-flux(i,k,lat))
     $              /(zw(i,k+1)-zw(i,k))
               tvcour = (zw(i,k)-wwdp(i,k,lat))/(zw(i,k+1)-zw(i,k))
               bvcour = (zw(i,k+1)-wwdp(i,k+1,lat))/(zw(i,k+1)-zw(i,k))
               if (abs(tvcour-bvcour).gt.1.0) then
c-drb (doesn't break)               if ((abs(tvcour-bvcour).gt.1.0) .and. (k.ne.42))then
                  write (6,*) 'vadvn: vertical courant numbers top, bot, diff,', i,k,lat,
     $                 tvcour, bvcour, tvcour-bvcour
                  write (6,*) ' zw, wwdp, i,k ', zw(i,k), wwdp(i,k,lat)
                  write (6,*) ' zw, wwdp, i,k+1 ', zw(i,k+1), wwdp(i,k+1,lat)
               endif
            end do
            dnout(plon+1,lat,k) = dnout(1,lat,k)
         end do                         ! do k = 1,plev
         iprob = 0
         jprob = 0
         do k = 1,plev
            do i = 1,plon
               rdiff = (dn(i,lat,k)-dnout(i,lat,k))
     $              /(abs(dn(i,lat,k))+abs(dnout(i,lat,k))+1.e-36)
               if (rdiff.gt.2.) then
                  iprob = i
                  jprob = k
               endif
            end do
         end do
         if (iprob.ne.0) then
            write (6,*) 'vadv: large change ',iprob, lat, jprob, rdiff,
     $           dn(iprob, k,jprob), dnout(iprob,k,jprob)
            call endrun
         endif
      end do                            ! do lat = 1,plat
#ifdef DEBUGSPITFIRE
      write (6,*) 'SPITFIRE  vadvn: range1 ', minval(dnout(:plon,:plat,:plev)), maxval(dnout(:plon,:plat,:plev))
#endif

      end subroutine vadvn

c##############################################################################

      subroutine getflxv (zw2, phi, flux, dp, kdep2, psi)

c     author: Phil Rasch

c     calculate the vertical mass fluxes 

c     zw2 are coords of cell walls
c     phi are mean tracer densities in cell (nominally at cell centers)

c     tracers are updated from the output from this code as
c     do k = 1,plev
c     do i = 1,plon      
c     q(i,k) = q(i,k) - (flux(i,k+1)-flux(i,k))/(zw2(i,k+1)-zw2(i,k))
c     end do
c     end do

c     .....zw21.......zw22.......zw23.......zw24.......zw25.......zw26
c     .... psiw1......psiw2......psiw3......psiw4......psiw5......psiw6
c     .... velw1......velw2......velw3......velw4......velw5......velw6
c     ...........phi1.......phi2........phi3......phi4........phi5.......

      implicit none

c     input
      real(DBLKIND) zw2(plond,plevp)             ! vertical wall coords
      real(r8) phi(plond,plevp-1)           ! densities
      real(DBLKIND) dp(plond,plevp)              ! departure point of cell walls
      integer kdep2(plond,plevp)        ! interval of departure point

c     output
      real(r8) flux(plond,plevp) 

c     local
      real(DBLKIND) psi(plond,plevp)
      real(DBLKIND) fdot(plond,plevp)

      real(DBLKIND) psistar(plond)
      real(r8) phinew
      integer i
      integer k


      do i = 1,plon
c        integral of phi
         psi(i,1) = 0.
c        fluxes at boundaries
         flux(i,1) = 0
         flux(i,plevp) = 0.
      end do

c     integral function
      do k = 2,plevp
         do i = 1,plon
            psi(i,k) = phi(i,k-1)*(zw2(i,k)-zw2(i,k-1)) + psi(i,k-1)
         end do
      end do

c     calculate the derivatives for the interpolating polynomial
      call cfdotmc (zw2, psi, fdot)


      do k = 2,plev
CDIR$ INLINE
         call cfintv(plond, plevp, plon, zw2, psi, fdot, dp(1,k), 
     $        psistar, kdep2(1,k))
         do i = 1,plon
            flux(i,k) = (psi(i,k)-psistar(i))
         end do
      end do
      
!      useful for double checking calculation
!      do k = 1,plev
!         do i = 1,plon
!            phinew = phi(i,k) - (flux(i,k+1)-flux(i,k))/(zw2(i,k+1)-zw2(i,k))
!            if (phinew.lt.0.) then
!               write (6,*) ' getflxv: phinew,phiold ', i,k, phinew, phi(i,k),
!     $              flux(i,k+1), flux(i,k), zw2(i,k+1), zw(i,k)
!               write (6,*) ' alternate phi ', (psi(i,k+1)-psi(i,k))/(zw2(i,k+1)-zw2(i,k))
!              call endrun()
!            endif
!         end do
!      end do

      end subroutine getflxv

c##############################################################################

      subroutine cfdotmc (x, f, fdot)


c     author: Phil Rasch

c     calculate the derivative for the interpolating polynomial
c     multi column version

      implicit none

c     input
      real(DBLKIND) x(plond, plevp)
      real(DBLKIND) f(plond, plevp)

c     output
      real(DBLKIND) fdot(plond, plevp)           ! derivative at nodes

c     assumed variable distribution
c     x1.......x2.......x3.......x4.......x5.......x6     1,plevp points
c     f1.......f2.......f3.......f4.......f5.......f6     1,plevp points
c     ...sh1.......sh2......sh3......sh4......sh5....     1,plev points
c     .........d2.......d3.......d4.......d5.........     2,plev points
c     .........s2.......s3.......s4.......s5.........     2,plev points
c     .............dh2......dh3......dh4.............     2,plev-1 points
c     .............eh2......eh3......eh4.............     2,plev-1 points
c     ..................e3.......e4..................     3,plev-1 points
c     .................ppl3......ppl4................     3,plev-1 points
c     .................ppr3......ppr4................     3,plev-1 points
c     .................t3........t4..................     3,plev-1 points
c     ................fdot3.....fdot4................     3,plev-1 points


c     work variables


      integer i
      integer k

      real(DBLKIND) a                            ! work var
      real(DBLKIND) b                            ! work var
      real(DBLKIND) c                            ! work var
      real(DBLKIND) s(plond,plevp)               ! first divided differences at nodes
      real(DBLKIND) sh(plond,plevp)              ! first divided differences between nodes
      real(DBLKIND) d(plond,plevp)               ! second divided differences at nodes
      real(DBLKIND) dh(plond,plevp)              ! second divided differences between nodes
      real(DBLKIND) e(plond,plevp)               ! third divided differences at nodes
      real(DBLKIND) eh(plond,plevp)              ! third divided differences between nodes
      real(DBLKIND) pp                           ! p prime
      real(DBLKIND) ppl(plond,plevp)             ! p prime on left
      real(DBLKIND) ppr(plond,plevp)             ! p prime on right
      real(DBLKIND) qpl
      real(DBLKIND) qpr
      real(DBLKIND) ttt
      real(DBLKIND) t
      real(DBLKIND) tmin
      real(DBLKIND) tmax
      real(DBLKIND) delxh(plond,plevp)


c     the minmod function 
      real(DBLKIND) minmod
      real(DBLKIND) medan
      minmod(a,b) = 0.5*(sign(one,a) + sign(one,b))*min(abs(a),abs(b))
      medan(a,b,c) = a + minmod(b-a,c-a)

      do k = 1,plev


c        first divided differences between nodes
         do i = 1, plon
            delxh(i,k) = (x(i,k+1)-x(i,k))
            sh(i,k) = (f(i,k+1)-f(i,k))/delxh(i,k)
         end do

c        first and second divided differences at nodes
         if (k.ge.2) then
            do i = 1,plon
               d(i,k) = (sh(i,k)-sh(i,k-1))/(x(i,k+1)-x(i,k-1))
               s(i,k) = minmod(sh(i,k),sh(i,k-1))
            end do
         endif
      end do

      if (plev.eq.1) then
         do i = 1,plon
            fdot(i,1) = 0.
            fdot(i,2) = 0.
         end do
         return
      endif

c     second and third divided diffs between nodes
      do k = 2,plev-1
         do i = 1, plon
            eh(i,k) = (d(i,k+1)-d(i,k))/(x(i,k+2)-x(i,k-1))
            dh(i,k) = minmod(d(i,k),d(i,k+1))
         end do
      end do

c     treat the boundaries
      do i = 1,plon
         e(i,2) = eh(i,2)
         e(i,plev) = eh(i,plev-1)
c        outside level
         fdot(i,1) = sh(i,1) - d(i,2)*delxh(i,1)
     $        - eh(i,2)*delxh(i,1)*(x(i,1)-x(i,3))
         fdot(i,1) = minmod(fdot(i,1),3*sh(i,1))
         fdot(i,plevp) = sh(i,plev) + d(i,plev)*delxh(i,plev) 
     $        + eh(i,plev-1)*delxh(i,plev)*(x(i,plevp)-x(i,plev-1))
         fdot(i,plevp) = minmod(fdot(i,plevp),3*sh(i,plev))
c        one in from boundary
         fdot(i,2) = sh(i,1) + d(i,2)*delxh(i,1)
     $        - eh(i,2)*delxh(i,1)*delxh(i,2)
         fdot(i,2) = minmod(fdot(i,2),3*s(i,2))
         fdot(i,plev) = sh(i,plev) - d(i,plev)*delxh(i,plev) 
     $        - eh(i,plev-1)*delxh(i,plev)*delxh(i,plev-1)
         fdot(i,plev) = minmod(fdot(i,plev),3*s(i,plev))
      end do


      do k = 3,plev-1
         do i = 1,plon
            e(i,k) = minmod(eh(i,k),eh(i,k-1))
         end do
      end do



      do k = 3,plev-1

         do i = 1,plon

c           p prime at k-0.5
            ppl(i,k)=sh(i,k-1) + dh(i,k-1)*delxh(i,k-1)  
c           p prime at k+0.5
            ppr(i,k)=sh(i,k)   - dh(i,k)  *delxh(i,k)

            t = minmod(ppl(i,k),ppr(i,k))

c           derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
            pp = sh(i,k-1) + d(i,k)*delxh(i,k-1) 

c           quartic estimate of fdot
            fdot(i,k) = pp 
     $           - delxh(i,k-1)*delxh(i,k)
     $           *(  eh(i,k-1)*(x(i,k+2)-x(i,k  )) 
     $           + eh(i,k  )*(x(i,k  )-x(i,k-2))
     $           )/(x(i,k+2)-x(i,k-2))

c           now limit it
            qpl = sh(i,k-1) 
     $        + delxh(i,k-1)*minmod(d(i,k-1)+e(i,k-1)*(x(i,k)-x(i,k-2)),
     $        d(i,k)  -e(i,k)*delxh(i,k))
            qpr = sh(i,k)
     $        + delxh(i,k  )*minmod(d(i,k)  +e(i,k)*delxh(i,k-1),
     $        d(i,k+1)+e(i,k+1)*(x(i,k)-x(i,k+2)))

            fdot(i,k) = medan(fdot(i,k), qpl, qpr)

            ttt = minmod(qpl, qpr)
            tmin = min(zero,3*s(i,k),1.5*t,ttt)
            tmax = max(zero,3*s(i,k),1.5*t,ttt)

            fdot(i,k) = fdot(i,k) 
     $           + minmod(tmin-fdot(i,k), tmax-fdot(i,k))

         end do

      end do

      end subroutine cfdotmc

c##############################################################################

      subroutine cfintv (plond, plevp, plon, x, f, fdot, xin, 
     $     psistar, kdep2)


c     author: Phil Rasch

c     do plon simultaneous interpolations over vectors of length plevp
c     in arrays dimensioned plond,plev

      implicit none

c     input
      integer plevp
      integer plond
      integer plon
      real(DBLKIND) x(plond, plevp)
      real(DBLKIND) f(plond, plevp)
      real(DBLKIND) fdot(plond, plevp)
      real(DBLKIND) xin(plond)
      integer kdep2(plond)

c     output
c     real(DBLKIND) fxdot(plond)
c     real(DBLKIND) fxdd(plond)
      real(DBLKIND) psistar(plond)

c     local
      integer i
      integer k

      real(DBLKIND) dx
      real(DBLKIND) s
      real(DBLKIND) c2
      real(DBLKIND) c3
      real(DBLKIND) xx
#ifdef LIMITON
      real(DBLKIND) psi1, psi2, psi3, psim
      real(DBLKIND) cfnew
#endif
      real(DBLKIND) cfint

      integer nwork
      parameter (nwork=NWORK)
      real(DBLKIND) xins


c     the minmod function 
      real(DBLKIND) a, b, c
      real(DBLKIND) minmod
      real(DBLKIND) medan
      minmod(a,b) = 0.5*(sign(one,a) + sign(one,b))*min(abs(a),abs(b))
      medan(a,b,c) = a + minmod(b-a,c-a)

      do i = 1,plon
         xins = max(x(i,1), min(xin(i), x(i,plevp)))
         k = kdep2(i)
         dx = (x(i,k+1)-x(i,k))
         s = (f(i,k+1)-f(i,k))/dx
         c2 = (3*s-2*fdot(i,k)-fdot(i,k+1))/dx
         c3 = (fdot(i,k)+fdot(i,k+1)-2*s)/dx**2
         xx = (xins-x(i,k))
c        fxdot(i) =  (3*c3*xx + 2*c2)*xx + fdot(i,k)
c        fxdd(i) = 6*c3*xx + 2*c2
         cfint = ((c3*xx + c2)*xx + fdot(i,k))*xx + f(i,k)
         psistar(i) = cfint
#ifdef LIMITON
c        leonards limiter
c        limit the interpolant
         psi1 = f(i,k)+(f(i,k+1)-f(i,k))*xx/dx
         if (k.eq.1) then
            psi2 = f(i,1)
         else
            psi2 = f(i,k) + (f(i,k)-f(i,k-1))*xx/(x(i,k)-x(i,k-1))
         endif
         if (k+1.eq.plevp) then
            psi3 = f(i,plevp)
         else
            psi3 = f(i,k+1) 
     $           - (f(i,k+2)-f(i,k+1))*(dx-xx)/(x(i,k+2)-x(i,k+1))
         endif
         psim = medan(psi1, psi2, psi3)
         cfnew = medan(cfint, psi1, psim)
c        if (abs(cfnew-cfint)/(abs(cfnew)+abs(cfint)+1.e-36)
c        $        .gt..03) then
c        write (6,*) ' cfintv limiting important ', cfint, cfnew
c        endif
         psistar(i) = cfnew
#endif   
      end do

      end subroutine cfintv

c##############################################################################

      subroutine setv (etadot, deltat)


c     author: Phil Rasch

c     calculate the departure points of vertical cell walls
c     and their intervals given vertical velocities and time step

      implicit none

c     input
c**********NOTE THAT THIS IS NOT IN STANDARD ORDER, OR DIMENSIONED PLOND!!!
      real(r8) etadot(plon,plevp,plat)      ! vertical velocities
      real(r8) deltat                       ! time step

c     local
      integer i
      integer j
      integer jj
      integer k
      integer kk
      real(DBLKIND) xins

! OpenMP directive
!$OMP PARALLEL DO PRIVATE( jj, kk, xins )
      do j = 1,plat

         do i = 1,plon
            wwdp(i,1,j) = 0.
            wwdp(i,plevp,j) = 1.
         end do

c        calculate departure points of cell walls
         do k = 2,plev
            do i = 1,plon
c              
c              the departure point
               wwdp(i,k,j) = zw(i,k)-etadot(i,k,j)*deltat
c-drb               write(6,*)'DRBDBG spitfire setv i,k,j,wwdp',i,k,j,wwdp(i,k,j)

c              guess the interval 
c              using a precalculated mapping from an even distribution to 
c              the actual grid, then refine it
c              xins = medan(zw(i,1), wwdp(i,k,j), zw(i,plevp))
               xins = max(zw(i,1), min(wwdp(i,k,j),zw(i,plevp)))
               kdep(i,k,j) = 0
               jj = int((xins-zw(i,1))/dzmin)+1
               kk = keven(jj)
               if (xins.lt.zw(i,kk)) then
                  kdep(i,k,j) = kk-1
               else if (xins.gt.zw(i,kk+1)) then
                  kdep(i,k,j) = kk+1
               else
                  kdep(i,k,j) = kk
               endif

            end do
         end do

      end do

      end subroutine setv

c##############################################################################

      subroutine fixhflux (fx, fy, pstend, deta, deltat, pde, dn1)

c     author: Phil Rasch (Dec 1998)

c---------------------------------------------------------------------
c DESCRIPTION
c     assume a continuity equation such that

c     dpsobs = -del . M

c     M = vert integral of mass fluxes at each level
c     dpsobs = surface pressure change in column

c     now assume 
c     M = Mmodel - Mdiff
c     where
c        Mmodel are the model predicted mass fluxes
c        Mdiff are the correction terms to get consistency
c     and assume
c     -Mdiff = del(chi)                                (G)

c     then
c     dpsobs = -del . (Mmodel + del(chi))
c            = -del . Mmodel - del^2(chi)

c     so 
c     del^2(chi) = -dpsobs - del.Mmodel = dpsdif       (H)

c     or
c     -del . Mdiff = dpsdif
c
c     so correction is done by.....
c 1) first calculate dpsobs
c 2) then calculate -del . Mmodel which we call dpsold
c 3) then calculate dpsdif
c 4) then solve helmholtz eqn (H) for chi
c 5) then calculate -Mdiff  from (G)
c 6) Mxnew and Mynew are merely M - Mdiff
c
c ----------------------------------------------------------------------

      implicit none

c include <fixflux.h>

c     input
      real(r8) deta(plev)                   ! thickness of layer in eta
      real(r8) pstend(plond, plat)          ! observed sfc pressure tendency
      real(r8) deltat                       ! time step
      real(r8) pde(plond,plat,plev)
      real(r8) dn1(plond,plat,plev)

c     input/output
      real(r8) fx(plond+1, plat, plev)      ! model x fluxes
      real(r8) fy(plond, plat+1, plev)       ! model y fluxes

c     work arrays
      
      integer i, ii
      integer j
      integer k

      real(r8) dpsobs(plon, plat)           ! observed sfc pressure change
      real(r8) dpsdeb(plon, plat)           ! observed sfc pressure change
      real(r8) dpsold(plon, plat)           ! original modeled sfc pressure change
      real(r8) dpsnew(plon, plat)           ! revised model sfc pressure change
      real(r8) dpsnobs(plon, plat)           ! revised model sfc pressure change
      real(r8) dpsnmod(plon, plat)           ! revised model sfc pressure change
      real(r8) dpsrmod(plon, plat)           ! revised model sfc pressure change
      real(r8) dpschk(plon, plat)            ! sfc pressure change for check
      real(r8) dpsdif(plon, plat)           ! difference in sfc pressure change
      real(r8) fxint(plon+1, plat)          ! vertically integrated x fluxes
      real(r8) fyint(plon, plat+1)          ! vertically integrated y fluxes
      real(r8) fxaint(plon+1, plat)          ! vertically integrated abs of x fluxes
      real(r8) fyaint(plon, plat+1)          ! vertically integrated abs of y fluxes
      real(r8) mxobs(plon+1, plat)             ! revised model x fluxes
      real(r8) myobs(plon, plat+1)             ! revised model y fluxes
      real(r8) mxrold(plon+1, plat)             ! revised model x fluxes
      real(r8) myrold(plon, plat+1)             ! revised model y fluxes
      real(r8) mxdold(plon+1, plat)             ! revised model x fluxes
      real(r8) mydold(plon, plat+1)             ! revised model y fluxes
      real(r8) mxdif(plon+1, plat)             ! revised model x fluxes
      real(r8) mydif(plon, plat+1)             ! revised model y fluxes
      real(r8) mxnew(plon+1, plat)             ! revised model x fluxes
      real(r8) mynew(plon, plat+1)             ! revised model y fluxes

      real(r8) dmufd(plat)
      real(r8) t1
      real(r8) pi
      real(r8) fxnew(plon+1,plat,plev), fynew(plon,plat+1,plev)
      real(r8) fxinew, fyinew
      real(r8) vfrac, sumdeta2
      real(r8) uvm(plon+1,plat,plev), vvm(plon,plat+1,plev)
      real(r8) dpstmp, dpsdev
      integer idev, jdev

c     Move this stack storage to static to avoid stack size limits.
      save fxnew, fynew, uvm, vvm

c     Do steps 1 and 2 simultaneously
c     step 1
c     calculate dpsobs = observed ps tendency
c     step 2
c     calculate dpsold = -del . Mmodel
      do j = 1,plat
         do i = 1,plon
            dpsdeb(i,j) = 0.
            dpsobs(i,j) = pstend(i,j)*deltat    ! step 1
            fxint(i,j)  = 0.
            fxaint(i,j)  = 0.
            fyint(i,j)  = 0.
            fyaint(i,j)  = 0.
         end do
      end do
      do j = 1,plat
         fxint(plon+1,j) = 0.
         fxaint(plon+1,j) = 0.
      end do
      do i = 1,plon
         fyint(i,plat+1) = 0.
         fyaint(i,plat+1) = 0.
      end do
               
      do k = 1,plev
         do j = 1,plat+1
            do i = 1,plon
               if (j.le.plat) then
                  fxint(i,j) = fxint(i,j) + fx(i,j,k)*deta(k)
                  fxaint(i,j) = fxaint(i,j) + abs(fx(i,j,k)*deta(k))
               endif
               fyint(i,j) = fyint(i,j) + fy(i,j,k)*deta(k)
               fyaint(i,j) = fyaint(i,j) + abs(fy(i,j,k)*deta(k))
            end do
         end do
      end do
      do j = 1,plat
         fxint(plon+1,j) = fxint(1,j)
         fxaint(plon+1,j) = fxaint(1,j)
      end do

c     dpsdeb is the pressure change calculated from spitfire
c     which should match the pressure change dpsold calculated below
      do k = 1,plev
         do j = 1,plat
            do i = 1,plon
               dpsdeb(i,j) = dpsdeb(i,j)+(dn1(i,j,k)-pde(i,j,k))*deta(k)
            end do
         end do
      end do

c     step 2
      dpsdev = 0.
      idev = 0
      jdev = 0
      do j = 1,plat
         do i = 1,plon
            dpsold(i,j) = - ( (fxint(i+1,j)-fxint(i,j))/dxlon
     $                       +(fyint(i,j+1)-fyint(i,j))/dmu(j)
     $                      )
            dpstmp = abs(dpsobs(i,j)-dpsold(i,j))
            if (dpstmp.gt.dpsdev) then
               idev = i
               jdev = j
               dpsdev = dpstmp 
            endif
         end do
      end do
#ifdef DEBUG
      i = idev
      j = jdev
      write (6,*) 'FIXHFLUX max dps = obstend - oldmodtend at ',
     $     idev, jdev, dpsdev/deltat, dpsobs(i,j)/deltat,dpsold(i,j)/deltat,
     $     dpsobs(i,j)/deltat - dpsold(i,j)/deltat
      write (6,*) ' fxint and dx, ', fxint(i+1,j), fxint(i,j), dx
      write (6,*) ' fyint and dy, ', fyint(i,j+1), fyint(i,j), dmu(j)
#endif

c     step 3
c     now calculate the difference between the assumed pressure change
c     and the modelled pressure change
      do j = 1,plat
         do i = 1,plon
            dpsdif(i,j) =  (-dpsobs(i,j) + dpsold(i,j))
         end do
      end do

c
c     step 4 and 5, compute the correction fluxes
      call dsc2fl (dpsdif, dxlon, dmu, mxdif, mydif)
c      call dsc2fl (dpsobs, dx, dmu, mxobs, myobs)
c      call dsc2fl (dpsold, dx, dmu, mxdold, mydold)

c
c     step 6 update the fluxes
      do j = 1,plat+1
         if (j.lt.plat+1) then
            do i = 1,plon+1
c               mxrold(i,j) = fxint(i,j) - mxdold(i,j)
               mxnew(i,j) = fxint(i,j) - mxdif(i,j)
            end do
         endif
         do i = 1,plon
c            myrold(i,j) = fyint(i,j) - mydold(i,j)
            mynew(i,j) = fyint(i,j) - mydif(i,j)
         end do
      end do
         
c     dpsnew contains the predicted pressure change from the new mass fluxes
c     dpsnew should be very close to dpsobs
      dpsdev = 0.
      do j = 1,plat
         do i = 1,plon
!            dpsnobs(i,j) =              ! pressure change from divergent obs fluxes
!     $              - (mxobs(i+1,j)-mxobs(i,j))/(dx)
!     $              - (myobs(i,j+1)-myobs(i,j))/dmu(j)
            dpsnew(i,j) =              ! revised pressure change 
     $              - (mxnew(i+1,j)-mxnew(i,j))/(dxlon)
     $              - (mynew(i,j+1)-mynew(i,j))/dmu(j)
            dpschk(i,j) =              ! revised pressure change 
     $              - (mxdif(i+1,j)-mxdif(i,j))/(dxlon)
     $              - (mydif(i,j+1)-mydif(i,j))/dmu(j)
c            dpsnmod(i,j) =              ! divergent model pressure change
c     $              - (mxdold(i+1,j)-mxdold(i,j))/(dx)
c     $              - (mydold(i,j+1)-mydold(i,j))/dmufd(j)
c            dpsrmod(i,j) =              ! rotational pressure change (should be zero)
c     $              - (mxrold(i+1,j)-mxrold(i,j))/(dx)
c     $              - (myrold(i,j+1)-myrold(i,j))/dmufd(j)
         end do
      end do
#ifdef DEBUG
      i = idev
      j = jdev
      write (6,*) 'FIXHFLUX comparison of dpsdif and dpschk at worst pt ', dpsdif(i,j), dpschk(i,j)
      write (6,*) 'FIXHFLUX obstend, proposed newtend at ',
     $     idev, jdev, dpsobs(i,j)/deltat,dpsnew(i,j)/deltat
      write (6,*) 'FIXHFLUX correction fluxes are ', mxnew(i,j), mxnew(i+1,j),  mynew(i,j), mynew(i,j+1)
#endif

c     now correct the horizontal velocities 
c     we assume that the mass flux is proportionate to the velocity
c     and partition it vertically according to the relative magnitude of the 
c     flux at each level to the total flux

      sumdeta2 = 0.
      do k = 1, plev
         sumdeta2 = sumdeta2 + deta(k)**2
      end do

      do k = 1,plev
         vfrac = deta(k)/sumdeta2
#ifdef DEBUGLEV
         i = ilook
         j = jlook
         write (6,*) ' k, vfrac is ', ilook, jlook, k, vfrac
         write (6,*) ' velocity before correction ',
     $        uw(i,j,k), uw(i+1,j,k), vw(i,j,k), vw(i,j+1,k)
         write (6,*) 'cour nums  ', uw(i+1,j,k)*deltat/dx, uw(i,j,k)*deltat/dx,
     $        vw(i,j+1,k)*deltat/dmu(j), vw(i,j,k)*deltat/dmu(j)
         write (6,*) ' fluxes before correction ',
     $        fx(i,j,k), fx(i+1,j,k), fy(i,j,k), fy(i,j+1,k)
#endif

        do j = 1, plat
            do i = 1, plon
               fxnew(i,j,k) = fx(i,j,k) - mxdif(i,j)*vfrac        ! new mass fluxes

               if ( fx(i,j,k) .ne. 0.0 ) then
                  uvm(i,j,k) = fxnew(i,j,k)/fx(i,j,k)
                  uw(i,j,k) = uw(i,j,k)*uvm(i,j,k)
               else
                  uvm(i,j,k) = 1.0
                  if ( fxnew(i,j,k) > 0.0 ) then
!pjr                     ii = modulo( i-1, plon ) + 1
                     ii = modulo( plon+i-2, plon ) + 1
                     uw(i,j,k) = fxnew(i,j,k)/pde(ii,j,k)/deltat
                  else
                     uw(i,j,k) = fxnew(i,j,k)/pde(i,j,k)/deltat
                  end if
               end if
c REMOVE
c               if (i.eq.1.and.j.eq.76) then
c                  ii = modulo( i-1, plon ) + 1
c                  write (6,*) 'fixflux: k, uwnew, ii, pde(ii), pde(i) ', k, uw(i,j,k), ii, pde(ii,j,k), pde(i,j,k)
c               endif
            end do
            fxnew(plon+1,j,k) = fxnew(1,j,k)
            uvm(plon+1,j,k) = uvm(1,j,k)
            uw(plon+1,j,k) = uw(1,j,k)
         end do

         do j = 2, plat
            do i = 1, plon
               fynew(i,j,k) = fy(i,j,k) - mydif(i,j)*vfrac        ! new mass fluxes
               if ( fy(i,j,k) .ne. 0.0 ) then
                  vvm(i,j,k) = fynew(i,j,k)/fy(i,j,k)
                  vw(i,j,k) = vw(i,j,k)*vvm(i,j,k)
               else
                  vvm(i,j,k) = 1.0
                  if ( fynew(i,j,k) > 0.0 ) then
                     vw(i,j,k) = fynew(i,j,k)/pde(i,j-1,k)/deltat
                  else
                     vw(i,j,k) = fynew(i,j,k)/pde(i,j,k)/deltat
                  end if
               end if
            end do
         end do

#ifdef DEBUGLEV
         i = ilook
         j = jlook
         write (6,*) ' velocity after correction ', i, j, k,
     $        uw(i,j,k), uw(i+1,j,k), vw(i,j,k), vw(i,j+1,k)
         write (6,*) 'cour nums  ', uw(i+1,j,k)*deltat/dx, uw(i,j,k)*deltat/dx,
     $        vw(i,j+1,k)*deltat/dmu(j), vw(i,j,k)*deltat/dmu(j)
         write (6,*) ' fluxes after correction ', i, j, k,
     $        fxnew(i,j,k), fxnew(i+1,j,k), fynew(i,j,k), fynew(i,j+1,k)
#endif

      end do

#ifdef DEBUG
      write (6,*) ' final scan for remaining worst deviation '
      dpsdev = 0.
      idev = 0
      jdev = 0
      do j = 1,plat
         do i = 1,plon
            dpstmp = abs(dpsobs(i,j)-dpsnew(i,j))
            if (dpstmp.gt.dpsdev) then
               idev = i
               jdev = j
               dpsdev = dpstmp 
            endif
         end do
      end do
      i = idev
      j = jdev
      write (6,*) ' max dps, obstend, newmodtend at ',
     $     idev, jdev, dpsdev/deltat, dpsobs(i,j)/deltat,dpsnew(i,j)/deltat
      write (6,*) ' new fluxes are ', mxnew(i,j), mxnew(i+1,j),  mynew(i,j), mynew(i,j+1)
#endif

      end subroutine fixhflux

c##############################################################################

      subroutine netdcal3( fx, fy, dpde, pstend, dtime, etadot, zflux )

c     Calculate a etadot using an integrated flux form of the mass continuity eqn.
c     The continuity eqn is 
c
c     d dpde                       d (edot dpde)
c     ------ + del . (V dpde) = -  --- 
c     d t                          d eta
c      
c     so
c                                                                           
c     (edot dpde) (bot-top) = - vint(del . V dpde) - pdot(bot-top)

      implicit none

c Input argments
      real(r8) fx(plond+1,plat,plev)        ! fluxes on x cell walls
      real(r8) fy(plond,plat+1,plev)        ! fluxes on y cell walls
      real(r8) dpde(plond,plat,plev)        ! d(p)/d(eta) after horizontal advection
      real(r8) pstend(plond,plat)           ! surface pressure tendency (dry)
      real(r8) dtime                        ! timestep in seconds
c Output arguments
      real(r8) etadot(plon,plevp,plat)      ! vertical velocity in eta coordinates
      real(r8) zflux(plond,plevp,plat)

c Local variables
      integer i, j, k, kk
      real(DBLKIND) fxint(plon+1,plat,plev), fyint(plon,plat+1,plev)
      real(DBLKIND) etdotdpde, divf, dpdot
c---------------------------------------------------------------------------

c Integrate fluxes from the top to each interface

      fxint(:,:,1) = fx(1:plon+1,:,1) * deta(1)
      do k = 2, plev
         do j = 1, plat
            do i = 1, plon+1
               fxint(i,j,k) = fxint(i,j,k-1) + fx(i,j,k)*deta(k)
            end do
         end do
      end do

      fyint(:,:,1) = fy(1:plon,:,1) * deta(1)
      do k = 2, plev
         do j = 1, plat+1
            do i = 1, plon
               fyint(i,j,k) = fyint(i,j,k-1) + fy(i,j,k)*deta(k)
            end do
         end do
      end do

c     solve for edot dpde (bottom to top)
c     convert to etadot by dividing by dp/deta
      do j = 1,plat

         do i = 1,plon
            etadot(i,1,j) = 0.
            etadot(i,plevp,j) = 0.
            zflux(i,1,j) = 0.
            zflux(i,plevp,j) = 0.
         end do

         do k = 2, plev
            do i = 1, plon
               divf = (fxint(i+1,j,k-1)-fxint(i,j,k-1))/dxlon +
     $                (fyint(i,j+1,k-1)-fyint(i,j,k-1))/dmu(j)

               dpdot = pstend(i,j)*hybi(k)

               etdotdpde = -divf/dtime - dpdot

               if ( etdotdpde .ge. 0. ) then
                  etadot(i,k,j) = etdotdpde / dpde(i,j,k-1)
               else
                  etadot(i,k,j) = etdotdpde / dpde(i,j,k)
               end if

               zflux(i,k,j) = etdotdpde*dtime
            end do
         end do
      end do

      end subroutine netdcal3

c##############################################################################

      subroutine fixvflux (wfld, fluxo, flux)


c     author: Phil Rasch (Dec 1998)

      implicit none

      real(r8) wfld(plon,plevp,plat)        ! vertical wind ("eta-dot").  Is either
      real(r8) fluxo(plond,plevp,plat)       ! optimal vertical mass flux
      real(r8) flux(plond,plevp,plat)        ! vertical mass flux from spitfire

      integer i, j, k
      real(DBLKIND) corm

      do j = 1,plat
        do k = 1,plev
           do i = 1,plon
             corm = fluxo(i,k,j)/(flux(i,k,j)+sign(1.e-36_r8,flux(i,k,j)))
             corm = min(max(corm,-two),two)
#ifdef DEBUG
             if (i.eq.22.and.j.eq.24.and.k.eq.plev) then
                write (6,*) 'FIXVFLUX fluxo(i,k,j) ', fluxo(i,k,j)
                write (6,*) 'FIXVFLUX flux(i,k,j) ', flux(i,k,j)
                write (6,*) 'FIXVFLUX corm ', corm
             endif            
#endif 
             wfld(i,k,j) = wfld(i,k,j)*corm
           end do
        end do
      end do

      end subroutine fixvflux

c##############################################################################

      subroutine chi2fl (chig, mx, my)

c     author: Phil Rasch (Dec 1998)

c---------------------------------------------------------------------
c DESCRIPTION

c     chi, find the fluxes fx and fy such that
c               d fx      d fy
c      dsc = -  ----   -  ----
c               d lam     d mu

c ----------------------------------------------------------------------

      implicit none

c     input
      real(r8) chig(plon, plat)

c     output
      real(r8) mx(plon+1, plat)      ! model x fluxes
      real(r8) my(plon, plat+1)       ! model y fluxes

c     work arrays
      
      integer i
      integer j
      integer k
      integer m
      integer n

      real(DBLKIND) den
      real(DBLKIND) vwcalc
      real(DBLKIND) mean

      real(DBLKIND) dmufd(plat)

c     x fluxes
      do j = 1,plat
         den = dxlon*cos(ylat(j))**2
         do i = 1,plon-1
            mx(i+1,j) = -(chig(i+1,j)-chig(i,j))/den
c     $           + fxint(i+1,j)
         end do
         mx(plon+1,j) = -(chig(1,j)-chig(plon,j))/den
c     $           + fxint(plon+1,j)
         mx(1,j) =  mx(plon+1,j)
      end do


c     y fluxes
c     interior
      do j = 1,plat-1
         do i = 1,plon
            vwcalc = -(chig(i,j+1)-chig(i,j))/(ylat(j+1)-ylat(j))
            my(i,j+1) = vwcalc*cos(yw(j+1))
         end do
      end do

c     poles
      do i = 1,plon
         my(i,plat+1) = 0.
         my(i,1) = 0.
      end do

c     now double check to see whether we have achieved our goal

c                1.           d^2 (a)              d                 d (a)

      end subroutine chi2fl

c##############################################################################

      subroutine laplc(chi, dpsnew)

c     author: Phil Rasch (Dec 1998)

      implicit none

c     input/output
      real(r8) chi(plon, plat)

c     work arrays
      real(DBLKIND) mx(plon+1, plat)      ! model x fluxes
      real(DBLKIND) my(plon, plat+1)       ! model y fluxes

      integer i
      integer j
      integer k
      integer m
      integer n
      integer nmin
      integer nmax

      real(r8) dpsnew(plon, plat)           ! revised model sfc pressure change
      real(DBLKIND) dpstmp(plon, plat)           ! revised model sfc pressure change
      real(DBLKIND) denx
      real(DBLKIND) conn
      real(DBLKIND) cons
      real(DBLKIND) myn
      real(DBLKIND) mys

c                1.           d^2 (a)              d                 d (a)
c     del^2(a) = --          ----       + ------------- [ cos(phi)   ------ ]
c                cos(phi)^2   dlam^2      cos(phi) d phi             d phi

c     call laplc_v1(chi, dpstmp)

c     interior points
      do j = 2,plat-1
         denx = (dxlon*cylat(j))**2
         conn = cyw(j+1)/(dmu(j)*(ylat(j+1)-ylat(j)))
         cons = cyw(j)/(dmu(j)*(ylat(j)-ylat(j-1)))
         do i = 2,plon-1
            myn = (chi(i,j+1)-chi(i,j))*conn
            mys = (chi(i,j)-chi(i,j-1))*cons
            dpsnew(i,j) = 
     $               + (chi(i+1,j)-2*chi(i,j)+chi(i-1,j))/(denx)
     $               + (myn-mys)
         end do
         i = 1
         myn = (chi(i,j+1)-chi(i,j))*conn
         mys = (chi(i,j)-chi(i,j-1))*cons
         dpsnew(i,j) = 
     $        + (chi(i+1,j)-2*chi(i,j)+chi(plon,j))/(denx)
     $        + (myn-mys)
         i = plon
         myn = (chi(i,j+1)-chi(i,j))*conn
         mys = (chi(i,j)-chi(i,j-1))*cons
         dpsnew(i,j) = 
     $        + (chi(1,j)-2*chi(i,j)+chi(i-1,j))/(denx)
     $        + (myn-mys)
      end do

      j = plat
      denx = (dxlon*cylat(j))**2
      cons = cyw(j)/(dmu(j)*(ylat(j)-ylat(j-1)))
      myn = 0.
      do i = 2,plon-1
         mys =  (chi(i,j)-chi(i,j-1))*cons
         dpsnew(i,j) = 
     $        + (chi(i+1,j)-2*chi(i,j)+chi(i-1,j))/(denx)
     $        + (myn-mys)
      end do
      i = 1
      mys =  (chi(i,j)-chi(i,j-1))*cons
      dpsnew(i,j) = 
     $     + (chi(i+1,j)-2*chi(i,j)+chi(plon,j))/(denx)
     $     + (myn-mys)
      i = plon
      mys =  (chi(i,j)-chi(i,j-1))*cons
      dpsnew(i,j) = 
     $     + (chi(1,j)-2*chi(i,j)+chi(i-1,j))/(denx)
     $     + (myn-mys)

      j = 1
      denx = (dxlon*cylat(j))**2
      conn = cyw(j+1)/(dmu(j)*(ylat(j+1)-ylat(j)))
      mys =  0.
      do i = 2,plon-1
         myn = (chi(i,j+1)-chi(i,j))*conn
         dpsnew(i,j) = 
     $        + (chi(i+1,j)-2*chi(i,j)+chi(i-1,j))/(denx)
     $        + (myn-mys)
      end do
      i = 1
      myn = (chi(i,j+1)-chi(i,j))*conn
      dpsnew(i,j) = 
     $     + (chi(i+1,j)-2*chi(i,j)+chi(plon,j))/(denx)
     $     + (myn-mys)
      myn = (chi(i,j+1)-chi(i,j))*conn
      i = plon
      dpsnew(i,j) = 
     $     + (chi(1,j)-2*chi(i,j)+chi(i-1,j))/(denx)
     $     + (myn-mys)

      end subroutine laplc

c##############################################################################

      subroutine inifixflux
 

c     author: Phil Rasch (Dec 1998)

c     initialize the fixflux routine

      use lgtran
      implicit none
#include <pmgrid.h>
#include <fixflux.h>


c     -initialize the transforms
      write (6,*) 'initializing transforms '
      write (6,*) ' nlon, nlat ', plon, plat 
      write (6,*) ' with spectral truncation params, mf, ms, ma ',
     $     mf, ms, ma
      write (6,*) ' space allocated for integer array ', idim
      write (6,*) ' space allocated for real arrays ', rdim + wdim
      if (plon > mf*3 .or. plat > int(ms*1.5)) then
         write (6,*) ' warning, gross mismatch of grid and spectral parameters' 
         call abort()
      endif
      call lginit( plon, plat, mf, ms, ma, idim, rdim, wdim
     $     , iwff, rwff, wkalff, rlff)
      write (6,*) ' complete '

      piff=acos(-1.)

      return
      end subroutine inifixflux

c##############################################################################

      subroutine dsc2fl (dscin, dx, dmu, mx, my)

c     author: Phil Rasch (Dec 1998)

c---------------------------------------------------------------------
c DESCRIPTION

c     given del^(chi), find the fluxes fx and fy such that
c               d fx      d fy
c      dsc = -  ----   -  ----
c               d lam     d mu

c ----------------------------------------------------------------------

      use lgtran
      implicit none
#include <pmgrid.h>
#include <fixflux.h>

c     input
      real(r8) dscin(plon, plat)
      real(DBLKIND) dx
      real(DBLKIND) dmu(plat)

c     output
      real(r8) mx(plon+1, plat)      ! model x fluxes
      real(r8) my(plon, plat+1)       ! model y fluxes

c     work arrays
      
      integer i
      integer j
      integer k
      integer m
      integer n
      integer nmin
      integer nmax
      integer nmm

      real(r8) dpsnew(plon, plat)           ! revised model sfc pressure change
      real(r8) dpsnew2(plon, plat)           ! revised model sfc pressure change
      real(r8) dsc(plon, plat)           ! difference in sfc pressure change
      real(r8) chig(plon,plat)               ! potential function
      real(r8) dschig(plon,plat)               ! del^2 of potential function
      real(DBLKIND) den
      real(r8) vwcalc
      real(r8) mean
      complex(r8) chis(0:mf,0:ms)           ! spectral coefs of chi
      complex(r8) chisn(0:mf,0:ms)           ! spectral coefs of chi new
      real(r8) rchis(2,0:mf,0:ms)  
      equivalence (chis, rchis)
      real(r8) maxcoef
      real(r8) err, maxerr
      real(r8) dmufd(plat)
      integer itr
      real(r8) eps0
      real(r8) t1
      common/itero/ niter,nitsm,icount,nitin,eer,eem
      real(r8) eer, eem
      integer niter, nitsm, icount, nitin
      integer firsttime
      data firsttime /0/
      save chig

      do j = 1,plat
         dmufd(j) = dmu(j)
      end do
c      dmufd(1) = (yw(2)-yw(1))*cos(ylat(1))
c      dmufd(plat) = dmufd(1)
      
c     since we are returning a flux divergence we cant resolve the 
c     mean field, so remove it for convenience

      mean = 0.
      do j = 1,plat
         do i = 1,plon
            mean = mean + dscin(i,j)*dmu(j)
         end do
      end do

      mean = mean/(2.*plon)
c     write (6,*) ' mean of dsc ', mean

      do j = 1,plat
         do i = 1,plon
            dsc(i,j) = dscin(i,j) - mean
         end do
c         i = 1
c         write (6,*) ' dsc no mean ', j, dsc(i,j)
      end do
      if (firsttime == 0) then
         do j = 1,plat
            do i = 1,plon
               chig(i,j) = 0.
            end do
         end do
! enable the next line if you want to remember previous solution
! rather than use spectral first guess
#ifdef REMEMBER
        firsttime = 1
#endif
      else
         go to 100
      endif

c     now solve for chi

c     -transform RHS to spectral space
      call lggts ( dsc, rchis, 1, iwff, rwff, workff )

#ifdef DEBSPC
c     reverse to identify the truncation
      call lgstg (chis, dschig, 1, iwff, rwff, workff )

      call lggts ( dschig, chisn, 1, iwff, rwff, workff )
#endif

c     invert laplacian in spectral space
#ifdef DEBSPC
      maxerr = 0.
      do m = 0,mf
         nmax = min(ms+ma,m+ms)
         nmin = max(1,m)
         do n = nmin,nmax
            nmm = n-m
            maxcoef = max(maxcoef,abs(chis(m,nmm)),abs(chisn(m,nmm)))
         end do
      end do
#endif
      do m = 0,mf
         nmax = min(ms+ma,m+ms)
         nmin = max(1,m)
c         write (6,*) ' rescale for m = ', m
c         write (6,*) ' n goes from ',nmin,' to ', nmax
         do n = nmin,nmax
            nmm = n-m
#ifdef DEBSPC
            err = abs(chis(m,nmm)-chisn(m,nmm))
     $           /(maxcoef + 1.e-36)
c     $           /(abs(chis(m,nmm))+abs(chisn(m,nmm)) + 1.e-36)
            maxerr = max(err,maxerr)
            if (err.gt.1.e-5) then
               write (6,*) 'bigerr ', 
     $              err, m, n, chis(m,nmm), chisn(m,nmm), maxerr
               stop
            endif
#endif
            chis(m,nmm) = -chis(m,nmm)/float(n*(n+1))
c            chis(m,nmm) = -chis(m,nmm)*float(n*(n+1))
         end do
      end do
#ifdef DEBSPC
      write (6,*) ' max relative error in transforms ', maxerr
#endif
      chis(0,0) = 0.
c     write (6,*) ' chis ', chis(1,1), chis(mf/2,mf/2)

C     
C     -transform back to grid space
      call lgstg (rchis, chig, 1, iwff, rwff, workff )

c     we now know chi at cell centers on interior of grid
c     calculate del(chi) to get fluxes

      call chi2fl(chig, mx, my)
c     choose sign of fluxes so that

c     now double check to see whether we have achieved our goal

c                1.           d^2 (a)              d                 d (a)
c     del^2(a) = --          ----       + ------------- [ cos(phi)   ------ ]
c                cos(phi)^2   dlam^2      cos(phi) d phi             d phi

c      write (6,*) ' dsc2fl: rhsorig, dscspectral '
      do j = 1,plat
         do i = 1,plon
            dpsnew(i,j) = 
     $              - (mx(i+1,j)-mx(i,j))/(dx)
     $              - (my(i,j+1)-my(i,j))/dmufd(j)
         end do
c         i = 1
c         write (6,66) j, dsc(i,j), dpsnew(i,j)
c   66    format (i5,5f15.5)
      end do

#ifdef DEBUG
      write (6,*) ' max dev between spectral and rhs ', 
     $     maxval(abs(dpsnew-dsc))
#endif

  100 continue
      itr = 100  ! max number of times that we do the iteration
      itr = 30  ! max number of times that we do the iteration
!      itr = 300
      eps0 = 1.e-30
c      write (6,*) ' at call gcrk '
      call gcrk(chig,dsc,itr,eps0)
c     niter is actual number of interations
c     nitsm is total number of times it was ever called
c     icount number of times gcrk was called
c     nitin is not used
c     eer is 2 norm error (initially 1)
c     eem is inf norm error (initially 1)
      if (max(eer,eem).gt.0.1) then
         write (6,*) ' slow convergence of gcrk, RMS and MAX error', eer, eem
      endif
      call laplc(chig, dpsnew2)

c      write (6,*) ' improved: dscorig, dpsnew, dpsnew2 '
c      do j = 1,plat
c         i = 1
c         write (6,66) j, dsc(i,j), dpsnew(i,j), dpsnew2(i,j)
c      end do

c     and translate back to fluxes
      call chi2fl(chig, mx, my)

      do j = 1,plat
         do i = 1,plon
            dpsnew(i,j) = 
     $              - (mx(i+1,j)-mx(i,j))/(dx)
     $              - (my(i,j+1)-my(i,j))/dmufd(j)
         end do
      end do

#ifdef DEBUG
      write (6,*) ' max dev between gcrk and rhs ', 
     $     maxval(abs(dpsnew-dsc)),
     $     maxval(abs(dpsnew2-dsc))
#endif

      return
      end subroutine dsc2fl

c##############################################################################

      subroutine gcrk(p,rhs,itr,eps0)

c     solves a helmholtz equation using conjugate gradients techniques

c     author: Piotr Smolarkiewicz (Dec 1998)

#include <pmgrid.h>

      real(r8) p(plon*plat),rhs(plon*plat)
      integer itr

      integer nn
      parameter (nn=plon*plat)
      common// r(nn),q(nn),ar(nn)
      real(r8) r, q, ar
c     niter is actual number of interations
c     nitsm is total number of times it was ever called
c     icount number of times gcrk was called
c     nitin is not used
c     eer is 2 norm error (initially 1)
c     eem is inf norm error (initially 1)
      integer niter, nitsm, icount, nitin
      real(r8) eer, eem
      common/itero/ niter,nitsm,icount,nitin,eer,eem
      integer lord
      parameter (lord=3)
      real(r8) x(nn,lord),ax(nn,lord),ax2(lord),axar(lord),del(lord)
      real(r8) eps0
convergence test modes **************************************************
      logical ctest
#ifdef DEBUG
      data ctest/.true./
#else
      data ctest/.false./
#endif
      integer nplt
      parameter (nplt=100)
      real(r8) err(0:nplt),xitr(0:nplt)
      integer ner
      real(r8) snorm
      real(r8) eps, eer0, eem0, rl20, rax
      integer nml, nlc, itmn, iprc
      integer k, i, ier, it, ii
      real(r8) epa, beta, dvmx, rl2, plmx
      if(ctest) then
!     itr=900/lord
      ner=9
      snorm=1./float(nn)
      eps0=1.e-10
      endif
convergence test modes **************************************************

      eps=eps0
      nml=nn
c     epa=1.e-60
      epa=1.e-36
      nlc=0
      itmn=3
      iprc=1                            ! 1 for preconditioning
      iprc=0                            ! 0 for no preconditioning

      do k=1,nml
        r(k)=0.
       ar(k)=0.
        q(k)=0.
      enddo
      do i=1,lord
       do k=1,nml
         x(k,i)=0.
        ax(k,i)=0.
       enddo
      enddo
      call laplc(p,r)
       do k=1,nml
         r(k)=r(k)-rhs(k)
       enddo
        call precon(r,q,ar,iprc)
      eer0=0.
      eem0=-1.e15
      rl20=0.
      do k=1,nml
      eer0=eer0+q(k)**2
      eem0=max(eem0,abs(q(k)))
      rl20=rl20+r(k)**2
      enddo
      eer0=max(eer0,epa)
      eem0=max(eem0,epa)
      rl20=max(rl20,epa)
convergence test modes **************************************************
      if(ctest) then
      do ier=0,nplt
      err(ier)=eps
      enddo
      eer=-1.e15
      do k=1,nml
         eer=max(eer,abs(r(k)))
      end do
      err(0)=eer
      print 300,  err(0)
  300 format(4x,e11.4,' residual error at it=1')
      endif
convergence test modes **************************************************
       do k=1,nml
        x(k,1)=q(k)
       enddo
      call laplc(x(1,1),ax(1,1))

      do 100 it=1,itr
       do i=1,lord
        ax2(i)=0.
        rax=0.
         do k=1,nml
          rax=rax+r(k)*ax(k,i)
          ax2(i)=ax2(i)+ax(k,i)*ax(k,i)
         enddo
        ax2(i)=max(epa,ax2(i))
        beta=-rax/ax2(i)
        dvmx=-1.e15
        rl2=0.
c        relerr = 0.
         do k=1,nml
          p(k)=p(k)+beta* x(k,i)
          r(k)=r(k)+beta*ax(k,i)
          dvmx=max(dvmx,abs(r(k)))
          rl2=rl2+r(k)*r(k)
c          relerr = relerr + abs(r(k)/(abs(rhs(k))xxxxxxxxxxxxxxxxxxxx
         enddo
c        write (6,*) ' dvmx, rl2 ', it, dvmx, rl2
       if(dvmx.le.eps.and.it.ge.itmn) go to 200
       if(rl2.ge.rl20.and..not.ctest) go to 200
          rl20=max(rl2,epa)
        call precon(r,q,ar,iprc)
        call laplc(q,ar)
        nlc=nlc+1
         do ii=1,i
          axar(ii)=0.
           do k=1,nml
            axar(ii)=axar(ii)+ax(k,ii)*ar(k)
           enddo
          del(ii)=-axar(ii)/ax2(ii)
cccc      del(ii)=max(del(ii),0.2)
         enddo
        if(i.lt.lord) then
          do k=1,nml
            x(k,i+1)= q(k)
           ax(k,i+1)=ar(k)
          enddo
           do ii=1,i
            do k=1,nml
              x(k,i+1)= x(k,i+1)+del(ii)* x(k,ii)
             ax(k,i+1)=ax(k,i+1)+del(ii)*ax(k,ii)
            enddo
           enddo
        else
          do k=1,nml
            x(k,1)= q(k)+del(1)* x(k,1)
           ax(k,1)=ar(k)+del(1)*ax(k,1)
          enddo
           do ii=2,i
            do k=1,nml
              x(k,1 )= x(k,1)+del(ii)* x(k,ii)
              x(k,ii)=0.
             ax(k,1 )=ax(k,1)+del(ii)*ax(k,ii)
             ax(k,ii)=0.
            enddo
           enddo
        endif
convergence test modes **************************************************
      if(ctest) then
      if(nlc/ner*ner.eq.nlc) then
      ier=nlc/ner
      eer=-1.e15
      do k=1,nml
         eer=max(eer,abs(r(k)))
      end do
      err(ier)=eer
      endif
      endif
convergence test modes **************************************************
       enddo
  100 continue
  200 continue
      eer=0.
      eem=-1.e15
      do k=1,nml
      eer=eer+q(k)**2 ! 2 norm
      eem=max(eem,abs(q(k))) ! inf norm
      enddo
      eer=sqrt(eer/eer0)
      eem=eem/eem0
      niter=nlc
      nitsm=nitsm+niter
      icount=icount+1

convergence test modes **************************************************
      if(ctest) then
         write (6,*) ' conj grad error '
      print 301, (err(ier),ier=1,nplt,1)
  301 format(4x,5e11.4)
!      do ier=0,nplt
!      xitr(ier)=ier*ner
!      err(ier)=log10(err(ier)/err(0))
!      end do
!      plmx=float(itr*lord)
c      call set(.1,.9,.1,.9,0.,plmx,-10.,2.,1)
c      call labmod('(i4)','(f5.0)',4,4,2,2,20,20,0)
c      call periml(3,10,7,2)
c      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)
c      call curved(xitr,err,nplt+1)
c      i1=int(102.4+409.6)
c      call wtstr(cpux(i1),cpuy(50),'niter',3,0,0)
c      call wtstr(cpux(17),cpuy(i1),'log e',3,90,0)
c      call frame
      endif
convergence test modes **************************************************
      return
      end subroutine gcrk

c##############################################################################

      subroutine precon(rhs,p,r,iflg)

c     author: Phil Rasch (Dec 1998)

#include <pmgrid.h>
#include <fixflux.h>

      integer nn
      parameter (nn=plon*plat)
!     real(r8) rhs(*),p(*),r(*)
      real(r8) rhs(nn),p(nn),r(*)
      complex(DBLKIND) chis(0:mf,0:ms)           ! spectral coefs of chi
c     work arrays
      real(DBLKIND) mx(plon+1, plat)      ! model x fluxes
      real(DBLKIND) my(plon, plat+1)       ! model y fluxes
      real(DBLKIND) t1
      integer iflg
      integer i

      if(iflg.eq.0) then
       do i=1,nn
        p(i)=rhs(i)
       enddo
      else
         do i = 1,nn
            p(i) = 0.
         end do
         call iterchi (p, mx, my, rhs)
      endif

      return
      end subroutine precon

c##############################################################################

      subroutine iterchi (chi, mx, my, rhs)

c     author: Phil Rasch (Dec 1998)

c---------------------------------------------------------------------
c DESCRIPTION

c     given chi, and an RHS use an iteration to improve the solution


c      del . del(chi) = RHS

c ----------------------------------------------------------------------

      implicit none
#include <pmgrid.h>
#include <fixflux.h>

c     input
      real(r8) rhs(plon, plat)

c     input/output
      real(r8) chi(plon, plat)

c     work arrays
      real(DBLKIND) mx(plon+1, plat)      ! model x fluxes
      real(DBLKIND) my(plon, plat+1)       ! model y fluxes

      integer i
      integer j
      integer k
      integer m
      integer n
      integer nmin
      integer nmax

      real(r8) dpsnew(plon, plat)           ! revised model sfc pressure change
      real(DBLKIND) dsc(plon, plat)           ! difference in sfc pressure change
      real(DBLKIND) dschig(plon,plat)               ! del^2 of potential function
      real(DBLKIND) den
      real(DBLKIND) vwcalc
      real(DBLKIND) mean
      real(DBLKIND) err, maxerr, maxlast
      integer iter
      real(DBLKIND) dtau

c     choose sign of fluxes so that

c     - del . fluxes = del^2(chi)


c      write (6,*) ' iterchi: rhs '
c      do j = 1,plat
c         i = 1
c         write (6,66) j, rhs(i,j)
c      end do

      dtau = 2.e-7 

      do iter = 1,4

         call laplc(chi, dpsnew)

         do j = 1,plat
            do i = 1,plon
               err = dpsnew(i,j)-rhs(i,j)
               chi(i,j) = chi(i,j) + dtau*err
            end do
         end do
      
      end do

      return
      end subroutine iterchi

c##############################################################################

      end module spitfire
