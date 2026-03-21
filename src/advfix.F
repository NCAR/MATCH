c $Id$
      subroutine iniadvfx( xlimfix, xetamid )

c Initialize /advfix0/.


      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <advfix.com>
c-----------------------------------------------------------------------

c Input arguments:

      logical
     $  xlimfix         ! true => limit the fixer

      real(r8)
     $  xetamid(plev)   ! hybrid coord at layer midpoints

c Local variables:

      integer k
c-----------------------------------------------------------------------

      limfix = xlimfix

      do k = 1, plev
         etamid(k) = xetamid(k)
      end do

      return
      end

c##############################################################################

      subroutine advfix( ncnst, dtime, lat, qin, qout )


      use precision

      implicit none
c-----------------------------------------------------------------------
c
c "Fix" the constituent fields (mass conservation)
c
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <advfix.com>
c-----------------------------------------------------------------------
c
c ... On input
c
      integer ncnst               ! number of fields
      real(r8) dtime                  ! Time step (seconds)
      integer lat                 ! Current latitude
      real(r8) qin (plond,plev,ncnst) ! input  constituent field
c
c ... Input/output
c
      real(r8) qout(plond,plev,ncnst) ! output constituent field 
c                                 ! on input : "unfixed"
c                                 ! on output: "fixed"

c----------------------------- local variables -------------------------
      integer i,k,m
      integer ifcnt
      real(r8) worst, wm, corm, absf

      real(r8) qf3m(plond,plev,ncnst) ! Fixer tendency
c----------------------------------------------------------------------

      do m = 1,ncnst 
         do k = 1,plev

            if ( .not. limfix ) then

               do i = 1,plon
                  qf3m (i,k,m) = alpha(m) * qout(i,k,m) *
     $                        etamid(k)*abs( qout(i,k,m) - qin(i,k,m) )
                  qout (i,k,m) = qout(i,k,m) + qf3m(i,k,m)
                  qf3m (i,k,m) = qf3m(i,k,m)/dtime
               enddo

            else    ! limit fixer

               do i = 1, plon
                  qf3m(i,k,m) = alpha(m) * 
     $                        etamid(k)*abs( qout(i,k,m) - qin(i,k,m) )
               enddo

               ifcnt = 0
               worst = 0
               wm = 0
               corm = 0.1

               do i = 1,plon
                  absf = abs( qf3m(i,k,m) )
                  if ( absf .gt. corm ) then
                     ifcnt = ifcnt + 1
                     worst = max( absf, worst )
                     wm = wm + absf
                     qf3m(i,k,m) = sign( corm, qf3m(i,k,m) )
                  endif
                  qf3m(i,k,m) = qout(i,k,m) * qf3m(i,k,m) / dtime
                  qout(i,k,m) = qout(i,k,m) + dtime*qf3m(i,k,m)
               end do

               wm = wm/max(ifcnt,1)
               if (ifcnt.gt.0.) then
                  write (6,*) ' fixer exceeded ', corm, ifcnt, wm,
     $                 k, m, lat, worst
               endif

            endif
         enddo
      enddo

      return
      end

c##############################################################################

      subroutine advfixcf( ncnst, ps, q1, q2, hw1 )

c Compute global mass integrals required to evaluate alpha for the advection
c mass fixer.


      use precision

      implicit none

c------------------------ parameters -----------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <advfix.com>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  ncnst               ! number of fields
      real(r8)
     $  ps(plond,plat)      ! surface pressure (time interpolant)
     $, q1(plond,plev,ncnst,platd)  ! advected species before adv
     $, q2(plond,plev,ncnst,platd)  ! advected species after adv
     $, hw1(ncnst)          ! constituent mass before adv

c Local variables:
      integer j, m
      real(r8)
     $  pint(plond,plevp)   ! pressure at layer interfaces
     $, pmid(plond,plev)    ! pressure at layer midpoints
     $, pdel(plond,plev)    ! pressure difference across layers
      real(r8)
     $  hw2j(ncnst,plat)    ! integration accumulator
     $, hw3j(ncnst,plat)    ! integration accumulator
     $, hw2(ncnst)          ! constituent mass after adv
     $, hw3(ncnst)          ! global integral of q2*|q2 - q1|**(3/2)
c-----------------------------------------------------------------------

! OpenMP directive
!$OMP PARALLEL DO PRIVATE(PINT, PMID, PDEL)
      do j = 1, plat

c        Compute global mass integrals "hw2" and "hw3"

         call plevs0( ps(1,j), pint, pmid, pdel )
         call qmassa( j, ncnst, q2(1,1,1,j), pdel, hw2j(1,j) )
         call qmassd( j, ncnst, q1(1,1,1,j), q2(1,1,1,j), pdel
     $,               hw3j(1,j) )
      end do

c     Sum global mass integrals (hw1 is summed in scanadv)

      do m = 1, ncnst
         hw2(m) = 0.
         hw3(m) = 0.
         do j = 1, plat
            hw2(m) = hw2(m) + hw2j(m,j)
            hw3(m) = hw3(m) + hw3j(m,j)
         end do
         if ( hw3(m) .ne. 0.0 ) then 
            alpha(m) = (hw1(m) - hw2(m))/hw3(m)
         else
            alpha(m) = 1.0
         end if
      end do

      return
      end

c##############################################################################

      subroutine advfixcf2( ncnst, ps, sh, as1, as2, hw1 )

c Compute global mass integrals required to evaluate alpha for the advection
c mass fixer.
c Tracer concentration units are dry mixing ratio (kg tracer / kg dry air).


      use precision

      implicit none

c------------------------ parameters -----------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <advfix.com>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  ncnst               ! number of fields
      real(r8)
     $  ps(plond,plat)       ! surface pressure
     $, sh(plond,plev,plat)  ! specific humidity
     $, as1(plond,plev,ncnst,platd)  ! advected species before adv
     $, as2(plond,plev,ncnst,platd)  ! advected species after adv
     $, hw1(ncnst)          ! constituent mass before adv

c Local variables:
      integer j, m
      real(r8)
     $  pint(plond,plevp)   ! pressure at layer interfaces
     $, pmid(plond,plev)    ! pressure at layer midpoints
     $, pdel(plond,plev)    ! pressure difference across layers
      real(r8)
     $  hw2j(ncnst,plat)    ! integration accumulator
     $, hw3j(ncnst,plat)    ! integration accumulator
     $, hw2(ncnst)          ! constituent mass after adv
     $, hw3(ncnst)          ! global integral of as2*|as2 - as1|**(3/2)
c-----------------------------------------------------------------------

! OpenMP directive
!$OMP PARALLEL DO PRIVATE(PINT, PMID, PDEL)
      do j = 1, plat

c        Compute global mass integrals "hw2" and "hw3"

         call plevs0( ps(1,j), pint, pmid, pdel )
         call gamdry( j, ncnst, as2(1,1,1,j), sh(1,1,j), pdel
     $,               hw2j(1,j) )
         call fixint( j, ncnst, as1(1,1,1,j), as2(1,1,1,j)
     $,               sh(1,1,j), pdel, hw3j(1,j) )
      end do

c     Sum global mass integrals (hw1 is summed in scanadv)

      do m = 1, ncnst
         hw2(m) = 0.
         hw3(m) = 0.
         do j = 1, plat
            hw2(m) = hw2(m) + hw2j(m,j)
            hw3(m) = hw3(m) + hw3j(m,j)
         end do
         if ( hw3(m) .ne. 0.0 ) then 
            alpha(m) = (hw1(m) - hw2(m))/hw3(m)
         else
            alpha(m) = 1.0
         end if
      end do

      return
      end

c##############################################################################

      subroutine advfixcf3( ncnst, ps, as1, as2, hw1 )

c Compute global mass integrals required to evaluate alpha for the advection
c mass fixer.
c Tracer concentration units are dry mixing ratio (kg tracer / kg dry air).
c Use dry surface pressure.


      use precision

      implicit none

c------------------------ parameters -----------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------
#include <advfix.com>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  ncnst               ! number of fields
      real(r8)
     $  ps(plond,plat)      ! dry surface pressure
     $, as1(plond,plev,ncnst,platd)  ! advected species before adv
     $, as2(plond,plev,ncnst,platd)  ! advected species after adv
     $, hw1(ncnst)          ! constituent mass before adv

c Local variables:
      integer j, m
      real(r8)
     $  pint(plond,plevp)   ! pressure at layer interfaces
     $, pmid(plond,plev)    ! pressure at layer midpoints
     $, pdel(plond,plev)    ! pressure difference across layers
      real(r8)
     $  hw2j(ncnst,plat)    ! integration accumulator
     $, hw3j(ncnst,plat)    ! integration accumulator
     $, hw2(ncnst)          ! constituent mass after adv
     $, hw3(ncnst)          ! global integral of as2*|as2 - as1|**(3/2)
c-----------------------------------------------------------------------

! OpenMP directive
!$OMP PARALLEL DO PRIVATE(PINT, PMID, PDEL)
      do j = 1, plat

c        Compute global mass integrals "hw2" and "hw3"

         call plevs0( ps(1,j), pint, pmid, pdel )
         call gamdry2( j, ncnst, as2(1,1,1,j), pdel, hw2j(1,j) )
         call fixint2( j, ncnst, as1(1,1,1,j), as2(1,1,1,j), pdel
     $,                hw3j(1,j) )
      end do

c     Sum global mass integrals (hw1 is summed in scanadv)

      do m = 1, ncnst
         hw2(m) = 0.
         hw3(m) = 0.
         do j = 1, plat
            hw2(m) = hw2(m) + hw2j(m,j)
            hw3(m) = hw3(m) + hw3j(m,j)
         end do
         if ( hw3(m) .ne. 0.0 ) then 
            alpha(m) = (hw1(m) - hw2(m))/hw3(m)
         else
            alpha(m) = 1.0
         end if
      end do

      return
      end
