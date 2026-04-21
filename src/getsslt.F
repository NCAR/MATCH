#include <params.h>
      subroutine getsslt (lat, rlat, ioro, u10, zm, zi, rpdel, sslt)

c     compute the sea salt optical depth, give the 10m wind speed
c     author: Phil RaschC     modified to remove all instances of tausslt   7/19/00

      use precision
      implicit none

#include <pmgrid.h>
#include <ssltcom.h>
#include <control.com>

      integer lat
      real(r8) rlat             ! [rdn] latitude of this row
      integer ioro(plond)
      real(r8) u10(plond)
      real(r8) zm(plond,plev)
      real(r8) zi(plond,plevp)
      real(r8) rpdel(plond,plev)
      integer ii
      integer npts
      integer indx(plond)

      real(r8) sslt(plond,plev)

      real(r8) zmin
      real(r8) u
      real(r8) z
      real(r8) dz
      real(r8) extc
      integer i
      integer k

      real(r8) alpha, beta, gamma
      real(r8) zib
      real(r8) h
      integer kstar
      real(r8) a(plond)
      real(r8) :: gravit = 9.80616

      real(r8) :: one = 1.0

      real(r8) lat_deg          ! latitude in degrees
      real(r8) rgn_scale        ! per-band multiplier (1.0 outside all bands)
      real(r8) tot_scale        ! global * regional composed multiplier
      integer ib

      ii = 1
      do i = 1,plon
         if (ioro(i) .eq. 0) then
            indx(ii) = i
            ii = ii + 1
         endif
      end do
      npts = ii-1
         


      do k = 1,plev
         do i = 1,plon
            sslt(i,k) = 0.
         end do
      end do

c     note that our burden is only from the sfc to about 500m 
C          (for the NCEP grid)
c     many explicit models find significant salt burdens above this level
c     to maybe 5km, although the studies of blanchard & woodcock show
c     concentrations at least an order of magnitude smaller than the
c     500m level values, and their fit is identified as being 
C          appropriate for 1 to 300m.
C

      do k = plev,1,-1
         zmin = 1.e36
c
c        write (6,*) ' points with salt, ug/m3, zt, zb, ug/m2 at lev ', k
c
         alpha = 5.e-6          ! gm/m3
         beta = 6.3e-6          ! /m
         do ii = 1,npts
            i = indx(ii)
            z = max(min(300._r8,zm(i,k)),one)
            zmin = min(zmin,zm(i,k))
            u = max(min(14._r8,u10(i)),one)
            dz = zi(i,k)-zi(i,k+1)
            gamma = 0.21-0.39*log10(u)
c
c           the next line follows the profile of blanchard and woodcock
c           sslt(i,k) = alpha*(beta*z)**gamma ! grams/m3
c           write (6,*) ' sslt midpoint (micrograms/m3)', 
c    $                  sslt(i,k)*1.e6, z
c            
            zib = max(zi(i,k+1),one)
c
c           the next line is the mean value over the model layers 
c           of the profile of blanchard and woodcock
c
            sslt(i,k) = alpha/(dz*(gamma+1.)*beta)
     $           *( (beta*zi(i,k))**(gamma+1.)
     $           -(beta*zib)**(gamma+1.)
     $           )
         end do
c
c        write (6,*) ' zmin for level k ', zmin, k
c        if all layers are 300 meters at this eta level, break
c
         if (zmin.gt.300) go to 100
      end do
 100  continue
      
      if (.true.) then
c
c        above this level we just make the sea salt get small fast up to 5 km
c        then forget it
c
         h = 500.      ! 500 meter efolding depth (sort of from fig 7 from b&w)
         kstar = k
         do ii = 1,npts
            i = indx(ii)
            a(i) = sslt(i,kstar)*exp(zm(i,kstar)/h)
         end do
         do k = kstar-1,1,-1
            zmin = 1.e36
c
c           write (6,*) ' sea salt spec ', k
c
            do ii = 1,npts
               i = indx(ii)
               zmin = min(zmin,zm(i,k))
               dz = zi(i,k)-zi(i,k+1)
               sslt(i,k) = a(i)*exp(-zm(i,k)/h)
c
c              write (6,67) i, sslt(i,k)*1.e6,zm(i,k),
c     $                     sslt(i,kstar)*1.e6,zm(i,kstar)
c   67         format (i5,4f10.3)
c
C
            end do
            if (zmin.gt.5.e3) go to 200
         end do
 200     continue
      endif
      
      
c++bee, convert sslt concentration from g/m3 --> kg/kg
      do k = 1,plev
         do i = 1,plon
            dz = zi(i,k)-zi(i,k+1)
            sslt(i,k) = sslt(i,k)*.001*dz*rpdel(i,k)*gravit
         end do
      end do
c--bee
c     Apply namelist sea salt scaling: global sslt_scale * per-latitude-band
c     factor from sslt_bands (if any band covers this row's latitude).
      rgn_scale = 1.0
      if (n_sslt_bands > 0) then
         lat_deg = rlat * 57.29577951308232_r8
         do ib = 1, n_sslt_bands
            if (lat_deg .ge. sslt_bands(1,ib) .and.
     $          lat_deg .le. sslt_bands(2,ib)) then
               rgn_scale = sslt_bands(3,ib)
               exit
            end if
         end do
      end if
      tot_scale = sslt_scale * rgn_scale
      if (tot_scale .ne. 1.0) then
         do k = 1,plev
            do i = 1,plon
               sslt(i,k) = tot_scale * sslt(i,k)
            end do
         end do
      end if
C++   CTB
      
      do k = 1,plev
         do i = 1,plon
            sslt_profiles(i, k, lat ) = sslt(i,k)

         end do
      end do   
C-- CTB
      return
      end


