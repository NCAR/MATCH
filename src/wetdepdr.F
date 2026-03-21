c $Id$

      subroutine wetdepdr( lat, dtime, t, sh, pmid, pdel, as
     $,                    cmfdqr, conicw, cwat, cldt, cldc, cldv
     $,                    cme, evapr, prain, rain, fracis
     $,                    obuf )

c Apply wet deposition to advected species.

      use histout, only: outfld
      use wetdep, only: wetdepa, wetdepg
c++caer
      use caer, only: docaer, caer_ncomp, caer_idx1
c--caer

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
c++caer
c----------------------------------------------------------------------
#include <tracnam.h>
c----------------------------------------------------------------------
c--caer

c Input arguments:
      integer
     $  lat         ! latitude index for this slice
      real(r8)
     $  dtime                  ! timestep
     $, t(plond,plev)          ! temperature (K)
     $, sh(plond,plev)         ! specific humidity (kg/kg)
     $, pmid(plond,plev)       ! pressure at layer midpoints
     $, pdel(plond,plev)       ! pressure difference across layers
     $, cmfdqr(plond,plev)     ! rate of production of convective precip
     $, conicw(plond,plev)     ! convective cloud water
!    $, icwmr(plond,plev)      ! in cloud water mixing ratio (convective)
     $, cwat(plond,plev)       ! cloud water (kg/kg)
     $, cldt(plond,plevp)      ! total cloud fraction
     $, cldc(plond,plev)       ! convective cloud fraction
     $, cldv(plond,plev)       ! fraction occupied by rain or cloud water 
     $, cme(plond,plev)        ! rate of cond-evap within the cloud
     $, evapr(plond,plev)      ! rate of evaporation of precip
     $, prain(plond,plev)      ! rate of production of stratiform precip
     $, rain(plond,plev)       ! mixing ratio of rain (kg/kg)

c Input/Output arguments:
      real(r8)
     $  as(plond,plev,pcnst)   ! advected species

c Output arguments:
      real(r8) 
     $  fracis(plond,plev,pcnst) ! fraction of transported species in gas phase 
     $, obuf(*)

c Local variables:
      integer
     $  i, k, m
      real(r8)
     $  molwt                    ! molecular weight
     $, solconst(plond,plev)     ! Henry's law coefficient
     $, solconst1
     $, scavt(plond,plev)
     $, iscavt(plond,plev)
!    $, icwmr2(plond,plev)
c++caer
      integer mm, nphob
c--caer
c-----------------------------------------------------------------------

c     The in-cloud water mixing ratio for zhang and hack schemes are added
c     together inside wetdepg (for now).  Use the input arg for icwmr which
c     represents the sum of the zhang and hack contributions
c     and zero the icwmr2 arg.
c      do k = 1,plev
c         do i = 1,plon
c            icwmr2(i,k) = 0.0
c         end do
c      end do
c     Scavenging of gas phase constituents by henry's law.
c      m = ?
c      molwt = ?
c      do k = 1,plev
c         do i = 1,plon
c            solconst(i,k) = ?
c         end do
c      end do
c      call wetdepg( lat, t, pmid, sh, pdel, cldt, cldc
c     $,             cmfdqr, prain, evapr, rain
c     $,             cwat, as(1,1,m), dtime, molwt, solconst
c     $,             scavt, iscavt, cldv, icwmr, icwmr2, fracis(1,1,m) )
c      do k = 1,plev
c         do i = 1,plon
c            as(i,k,m) = as(i,k,m) + scavt(i,k)*dtime
c         end do
c      end do
c      call outfld( '?WET', scavt, plond, lat, obuf )

c     Wet deposition of a highly soluble aerosol.
c      m = ?
c      molwt = ?
c      solconst1 = ?
c      call wetdep( lat, t, pmid, sh, pdel, cldt, cldc
c     $,             cmfdqr, prain, cme, evapr
c     $,             cwat, as(1,1,m), dtime, molwt,
c     $,             scavt, iscavt, cldv, fracis(1,1,m) )
c      do k = 1,plev
c         do i = 1,plon
c            as(i,k,m) = as(i,k,m) + scavt(i,k)*dtime
c         end do
c      end do
c      call outfld( '?WET', scavt, plond, lat, obuf )

c++caer
      molwt = 0.                ! not used
      solconst1 = 0.            ! not used

      if ( docaer() ) then
c        Wet deposition of hydrophilic carbon aerosol species.
         nphob = caer_ncomp()/2
         do m = 1, nphob
            mm = caer_idx1() + nphob + m - 1
            call wetdepa( lat, t, pmid, sh, pdel
     $,                cldt, cldc, cmfdqr, conicw, prain, cme
     $,                evapr, cwat, as(:,:,mm), dtime
     $,                scavt, iscavt, cldv, fracis(:,:,mm), 0.3_r8 )
            do k = 1,plev
               do i = 1,plon
                  as(i,k,mm) = as(i,k,mm) + scavt(i,k)*dtime
               end do
            end do
            call outfld( trim(tracnam(mm))//'WET', scavt, plond
     $,                  lat, obuf )
         end do
      end if
c--caer

      return
      end

c##############################################################################

