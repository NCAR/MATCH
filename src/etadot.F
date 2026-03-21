c $Id$
      subroutine inietd( rearth, xdlam, xnprlev, xhybi, xhybd, xetamid
     $,                  xfilter, pi, xphi )

c Initialize /etadot0/.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <etadot.com>
c-----------------------------------------------------------------------

c Input arguments:

      integer
     $  xnprlev        ! index of lowest pure pressure interface

      logical
     $  xfilter        ! true => apply filter to calculated etadot values

      real(r8)
     $  rearth         ! radius of earth (m)
     $, xdlam          ! longitude grid spacing (radians)
     $, xhybi(plevp)   ! ps component of hybrid coordinate - interfaces
     $, xhybd(plev)    ! difference of the hybi's
     $, xetamid(plev)  ! hybrid coord at layer midpoints
     $, pi             ! 3.14159...
     $, xphi(platd)    ! Latitude in radians.

c Local variables:

      integer
     $  j, j0, jcen, k                 ! indices
     $, nint

      real(r8)
     $  dloneq, dlonmx, dlon
c-----------------------------------------------------------------------

      arad = rearth
      dlam = xdlam
      nprlev = xnprlev
      filter = xfilter
      do k = 1, plevp
         hybi(k) = xhybi(k)
      end do
      do k = 1, plev
         hybd(k) = xhybd(k)
         etamid(k) = xetamid(k)
      end do
      do j = 1, platd
         phi(j) = xphi(j)
      end do

      if ( .not. filter ) return

c     longitude interval in latitude circle closest to equator
c     (on unit sphere).

      j0 = plat/2 + 1
      jcen = j0 + j1m
      dloneq = 2.*pi/plat
      dlonmx =  dloneq * cos( phi(jcen) )

c     At each lat calc longitude interval.  If the interval is less than
c     2/3 of dlonmx then double the number of intervals and make sure that
c     the sum of those intervals is at least 2/3 of dlonmx.  Repeat until
c     this condition is met.  The resulting number of intervals divided by 
c     2 (integer arithmetic) is the number of points on each side of the
c     midpoint to be used in the filter.

      nint = 1
      nfltpts(j0) = 2*(nint/2) + 1
      do j = j0+1, plat
         jcen = j + j1m
         dlon = dloneq * cos( phi(jcen) )
   10    continue
         if ( nint*dlon  .lt.  2.*dlonmx/3. ) then
            nint = nint*2
            go to 10
         end if
         nfltpts(j) = min( 2*(nint/2) + 1, plon )
      end do
      do j = 1, plat/2
         nfltpts(j) = nfltpts(plat-j+1)
      end do

      return
      end

c##############################################################################

      subroutine etdcal( lat, u, v, d, ps, pdel, pmid, etadot )

c Calculate etadot as is done in the ccm2 calculation.

      use precision

      implicit none
c---------------------------------------------------------------------
#include <pmgrid.h>
c---------------------------------------------------------------------
#include <etadot.com>
c---------------------------------------------------------------------

c Input arguments:

      integer
     $  lat                  ! lat index

      real(r8)
     $  u(plond,plev)      ! zonal wind 
     $, v(plond,plev)      ! meridional wind 
     $, d(plon,plev)       ! divergence
     $, ps(plond,plat)     ! surface pressure
     $, pdel(plond,plev)   ! layer thicknesses (pressure)
     $, pmid(plond,plev)   ! pressure at full levels

c Output arguments:

      real(r8)
     $  etadot(plon,plevp) ! vertical velocity in eta coordinates

c Local variables:

      integer
     $  i, k                     ! indices
     $, jcen                     ! latitude index for extended phi array
     $, nsum, istart, iend, ii   ! for filter

      real(r8)
     $  edotdpde(plon,plevp) ! etadot*dp/deta
     $, ddpk(plon)           ! partial sum of d*delta p
     $, ddpn(plon)           ! complete sum of d*delta p
     $, dpsl(plon)           ! longitudinal component of grad ln(ps)
     $, dpsm(plon)           ! latitudinal component of grad ln(ps)
     $, vpdsn(plon)          ! complete sum V dot grad(ln(ps)) delta b
     $, vpdsk(plon)          ! partial sum  V dot grad(ln(ps)) delta b
     $, vkdp                 ! (v(k)*grad(lnps)*db(k))
     $, tmp                  ! temporary variable
     $, psl(plon)            ! log of sfc pressure (center lat)
     $, psln(plon)           ! log of sfc pressure (north lat)
     $, psls(plon)           ! log of sfc pressure (south lat)
     $, pavg                 ! sfc pressure at poles
     $, sum(plon)            ! tmp array for filter

c Functions:

      integer imap, in, len
      imap(in, len) = mod(in+len*2-1, len)+1
c---------------------------------------------------------------------

      jcen = lat + j1m

c     calculate log of sfc pressure in 3 lat band
c     for later use in 2nd order differences to the gradient terms
c     also initialize some accumulators

c     central band
      do i = 1,plon
         psl(i) = log(ps(i,lat))
         vpdsn(i) = 0.
         vpdsk(i) = 0.
         ddpn(i) = 0.
         ddpk(i) = 0.
      end do

c     south band
      if (lat.eq.1) then
         pavg = 0.
         do i = 1,plon
            pavg = pavg + ps(i,1)
         end do
         pavg = pavg/float(plon)
         do i = 1,plon
            psls(i) = log(pavg)
         end do
      else
         do i = 1,plon
            psls(i) = log(ps(i,lat-1))
         end do
      endif

c     north band
      if (lat.eq.plat) then
         pavg = 0.
         do i = 1,plon
            pavg = pavg + ps(i,plat)
         end do
         pavg = pavg/float(plon)
         do i = 1,plon
            psln(i) = log(pavg)
         end do
      else
         do i = 1,plon
            psln(i) = log(ps(i,lat+1))
         end do
      endif
c
c     calculate second order diffs to gradiant terms
c
      do i = 1,plon
         dpsl(i) = (psl(imap(i+1,plon)) - psl(imap(i-1,plon))
     $             ) / (2.*arad*dlam*cos(phi(jcen)))
         dpsm(i) = (psln(i) - psls(i)
     $             ) / (arad*(phi(jcen+1)-phi(jcen-1)))
      end do
C
C     sum(plev)(d(k)*dp(k))
C
      do k=1,plev
         do i=1,plon
            ddpn(i) = ddpn(i) + d(i,k)*pdel(i,k)
         end do 
      end do
C
C     sum(plev)(v(k)*grad(lnps)*db(k))
C
      do k=nprlev,plev
         do i=1,plon
            vkdp = (u(i,k)*dpsl(i) + v(i,k)*dpsm(i))*ps(i,lat)
            vpdsn(i) = vpdsn(i) + vkdp*hybd(k)
         end do
      end do
C
C     Compute etadot (dp/deta) (k+1/2).  Note: sum(k)(d(lat)*dp(lat)) required in
C     pressure region. sum(k)(d(lat)*dp(lat)) and sum(k)(v(lat)*grad(pstar)*db(lat)) 
C     required in hybrid region
C
      do i=1,plon
         edotdpde(i,1) = 0.
         edotdpde(i,plevp) = 0.
      end do
      do k=1,nprlev-1
         do i=1,plon
            ddpk(i) = ddpk(i) + d(i,k)*pdel(i,k)
            edotdpde(i,k+1) = -ddpk(i)
         end do
      end do
      do k=nprlev,plev-1
         do i=1,plon
            ddpk(i) = ddpk(i) + d(i,k)*pdel(i,k)
            vkdp = (u(i,k)*dpsl(i) + v(i,k)*dpsm(i))*ps(i,lat)
            vpdsk(i) = vpdsk(i) + vkdp*hybd(k)
            edotdpde(i,k+1) = -ddpk(i) - vpdsk(i) +
     $                        hybi(k+1)*(ddpn(i)+vpdsn(i))
         end do
      end do
C
C     Convert eta-dot(dp/deta) to eta-dot (top and bottom = 0.)
C
      do i=1,plon
         etadot(i,1) = 0.
         etadot(i,plevp) = 0.
      end do
      do k=2,plev
         tmp = etamid(k) - etamid(k-1)
         do i=1,plon
            etadot(i,k) = edotdpde(i,k)*tmp/(pmid(i,k) - pmid(i,k-1))
         end do
      end do

c     Filter eta-dot if requested.
      if ( .not. filter ) return

c     If no filtering at this latitude exit
      nsum = nfltpts(lat)
      if ( nsum .eq. 1 ) return

      if ( nsum .eq. plon ) then
c        All points in filter.  Only need calc one average.
         do k = 2, plev
            tmp = 0.0
            do i = 1, plon
               tmp = tmp + etadot(i,k)
            end do
            tmp = tmp / plon
            do i = 1, plon
               etadot(i,k) = tmp
            end do
         end do
      else
c        Filter eta-dot using a running average.
         do k = 2, plev
            do i = 1, plon
               istart = i - (nsum-1)/2
               iend = i + (nsum-1)/2
               sum(i) = 0.0
               do ii = istart, iend
                  sum(i) = sum(i) + etadot(imap(ii,plon),k)
               end do
               sum(i) = sum(i) / float(nsum)
            end do
            do i = 1, plon
               etadot(i,k) = sum(i)
            end do
         end do
      end if

      return
      end


