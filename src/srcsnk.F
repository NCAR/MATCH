c $Id$
      subroutine srcsnkdr( lat, dtime, calday, ioro, q, as, obuf )

c Atmospheric source/sink terms.  (Simple chemistry schemes go here.)

      use histout, only: outfld

#ifdef CAER
c++caer
      use caer, only: docaer, caer_ncomp, caer_idx1, caercv
c--caer
#endif

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <tracnam.h>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  lat         ! latitude index
     $, ioro(plond) ! nint(oro(i))
      real(r8)
     $  dtime             ! timestep size (s)
     $, calday            ! Julian day + fraction (Greenwich)
     $, q(plond,plev)     ! specific humidity

c Input/Output arguments:
      real(r8)
     $  as(plond,plev,pcnst)   ! advected species

c Output argument:
      real(r8)
     $  obuf(*)                ! output buffer

c Local variables:
      integer i, k, m
      real(r8)
     $  srcsnk(plond,plev,pcnst)
#ifdef CAER
c++caer
      integer mm, nphob
c--caer
#endif
c-----------------------------------------------------------------------

      call setr( plond*plev*pcnst, 0._r8, srcsnk )

#ifdef CAER
c Add code to calculate source/sink contributions to tracers here.
c++caer
      if ( docaer() ) then
c        Conversion of hydrophobic to hydrophilic carbon aerosol.
         nphob = caer_ncomp()/2
         do m = 1, nphob
            mm = caer_idx1() + m - 1
            call caercv( as(1,1,mm), dtime, srcsnk(1,1,mm)
     $,                  srcsnk(1,1,mm+nphob) )
            call outfld( trim(tracnam(mm))//'SRC', srcsnk(1,1,mm), plond
     $,                  lat, obuf )
         end do
      end if
c--caer
#endif

c     Add contribution of src/snk tendency to constituent fields

      do m = 1, pcnst
         do k = 1, plev
            do i = 1, plon
               as(i,k,m) = as(i,k,m) + dtime*srcsnk(i,k,m)
            end do
         end do
      end do

      return
      end
