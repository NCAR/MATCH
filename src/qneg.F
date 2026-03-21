c $Id$
      subroutine iniqneg( xqmin )


      use precision

      implicit none
c----------------------------------------------------------------------------
#include <pmgrid.h>
c----------------------------------------------------------------------------
#include <qneg.com>
c----------------------------------------------------------------------------

c Input arguments:

      real(r8) 
     $      xqmin(pcnst)  ! Minimum mixing ratio for constituents.

c Local variables:

      integer m
c----------------------------------------------------------------------------

      do m = 1, pcnst
         qmin(m) = xqmin(m)
      end do

      return
      end

c##############################################################################

      subroutine qneg3(subnam  ,lat     ,q       )
C-----------------------------------------------------------------------
C
C Check moisture and tracers for minimum value, reset any below
C minimum value to minimum value and return information to allow
C warning message to be printed. The global average is NOT preserved.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, August 1992
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c

      use precision

#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <qneg.com>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      character*(*) subnam       ! name of calling routine
      integer lat              ! latitude index
C
C Input/Output arguments
C
      real(r8) q(plond,plev,pcnst) ! moisture/tracer field
C
C---------------------------Local workspace-----------------------------
C
      integer indx(plond)      ! array of indices of points < qmin
      integer nval             ! number of points < qmin for 1 level
      integer nvals            ! number of values found < qmin
      integer i,ii,k           ! longitude, level indices
      integer m                ! constituent index
      logical found            ! true => at least 1 minimum violator found
      integer iw,kw            ! i,k indices of worst violator
      real(r8) worst               ! biggest violator

#undef DEBUG
#ifdef DEBUG
      integer ilook,jlook,klook
#endif
C
C-----------------------------------------------------------------------
C
#ifdef DEBUG
      ilook = 99
      jlook = 72
      klook = 23
#endif

#ifdef DEBUG
         if ( lat.eq.jlook ) write(6,*)'qneg3 0.0 ',q(ilook,klook,pcnst),qmin(pcnst)
#endif

      do m=1,pcnst

!        write (6,*) 'pre-qneg3 range:',subnam, m, lat,minval(q(1:plon,:,m)), maxval(q(1:plon,:,m))

         nvals = 0
         found = .false.
         worst = 1.e35
C
C Test all field values for being less than minimum value. Set q = qmin
C for all such points. Trace offenders and identify worst one.
C
         do k=1,plev
            call whenflt(plon,q(1,k,m),1,qmin(m),indx,nval)
            if (nval.gt.0) then
               found = .true.
               nvals = nvals + nval
               do ii=1,nval
                  i = indx(ii)
                  if (q(i,k,m).lt.worst) then
                     worst = q(i,k,m)
                     kw = k
                     iw = i
                  end if
                  q(i,k,m) = qmin(m)
               end do
            end if
         end do
         if ( found .and. abs(worst) .gt. qmin(m) ) then
            write(6,9000)subnam,m,lat,nvals,qmin(m),worst,iw,kw
         end if
      end do
C

      return
 9000 format(' QNEG3 from ',a8,':m=',i3,' lat=',i3,
     $       ' Min. mixing ratio violated at ',i4,' points.  Reset to ',
     $       1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)
      end

c##############################################################################

      subroutine shneg( mes, lat, q )

C Check fields and reset any negative values to zero.


      use precision

      implicit none
#if ( defined RS6000 )

      use precision

      implicit automatic (a-z)
#endif
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      character mes*(*)        ! message
      integer lat              ! latitude index

C Input/Output arguments

      real(r8) q(plond,plev)       ! specific humidity
C
C---------------------------Local workspace-----------------------------
C
      integer indx(plond)      ! array of indices of points < qmin
      integer nval             ! number of points < qmin for 1 level
      integer nvals            ! number of values found < qmin
      integer i,ii,k           ! longitude, level indices
      logical found            ! true => at least 1 minimum violator found
      integer iw,kw            ! i,k indices of worst violator
      real(r8) worst               ! biggest violator
      real(r8) qmin
      parameter( qmin = 1.e-12 ) ! from CCM2
C-----------------------------------------------------------------------

      nvals = 0
      found = .false.
      worst = 1.e35
      do k=1,plev
         call whenflt(plon,q(1,k),1,qmin,indx,nval)
         if (nval.gt.0) then
            found = .true.
            nvals = nvals + nval
            do ii=1,nval
               i = indx(ii)
               if (q(i,k).lt.worst) then
                  worst = q(i,k)
                  kw = k
                  iw = i
               end if
               q(i,k) = qmin
            end do
         end if
      end do
      if ( found  .and. mes.ne.'VDIFF:sh') then
         write(6,9000)mes,lat,nvals,qmin,worst,iw,kw
      end if
C
      return
 9000 format('shneg: ',a8,' lat=',i3,
     $       ' neg values at ',i4,' points.  Reset to ',
     $       1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)
      end

