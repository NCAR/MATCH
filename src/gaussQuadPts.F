      subroutine gauaw(oa      ,ow      ,k)
C-----------------------------------------------------------------------
C Jan 2001
C modified to return a different precision from that actually used to calculate
c the gaussian wts by Phil Rasch
C
C Calculate sine of latitudes a(k) and weights w(k) for the gaussian
C quadrature. The algorithm is described in Davis and Rabinowitz,
C Journal of Research of the NBS, V 56, Jan 1956.
C The zeros of the  bessel function j0, which are obtained 
C from bsslzr, are used as a first guess for the abscissa.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      L. Bath, June 1992
C Reviewed:          D. Williamson, J. Hack, August 1992
C
C-----------------------------------------------------------------------

      use precision

      implicit none
#include <pmgrid.h>
c      integer, parameter :: DBLKIND = selected_real_kind( p = 12 )
C------------------------------Arguments--------------------------------
C
C Input argument
C
      integer k         !  number of latitudes pole to pole
C
C Output arguments
C
      real(r8) ::
     $  oa(k)            !  sine of latitudes
     $, ow(k)            !  weights
C
C---------------------------Local variables-----------------------------
C
      real(DBLKIND) ::
     $  a(k)            !  sine of latitudes
     $, w(k)            !  weights
C
      real(DBLKIND) ::
     $     pi,          ! value of pi
     $     eps,         ! convergence criterion
     $     c,           ! constant combination
     $     fk,          ! real k
     $     xz,          ! abscissa estimate
     $     pkm1,        ! |
     $     pkm2,        ! |-polynomials
     $     pkmrk,       ! |
     $     pk,          ! |
     $     sp,          ! current iteration latitude increment
     $     avsp,        ! |sp|
     $     fn           ! real n
      integer kk,       ! k/2 (number of latitudes in hemisphere)
     $        is,       ! latitude index
     $        iter,     ! iteration counter
     $        n,l       ! indices
C
C-----------------------------------------------------------------------
C
      eps = 1.e-13
      eps = 2*epsilon(pi)
      write (6,*) ' gausswts: eps is ', eps
      pi = 4.*atan(1.)
C
C The value eps, used for convergence tests in the iterations, 
C can be changed.  Newton iteration is used to find the abscissas.
C
      c = (1.-(2./pi)**2)*0.25
      fk = k
      kk = k/2
      call bsslzr(a,kk)
      do 30 is=1,kk
         xz = cos(a(is)/sqrt((fk+0.5)**2+c))
C
C This is the first approximation to xz
C
         iter = 0
   10    pkm2 = 1.
         pkm1 = xz
         iter = iter + 1
         if (iter.gt.10) then
C
C Error exit
C
            write(6,*)'GAUAW:Error exit,no convergence in 10 iterations'
            stop
         end if
C
C Computation of the legendre polynomial
C
         do 20 n=2,k
            fn = n
            pk = ((2.*fn-1.)*xz*pkm1-(fn-1.)*pkm2)/fn
            pkm2 = pkm1
            pkm1 = pk
   20    continue
         pkm1 = pkm2
         pkmrk = (fk*(pkm1-xz*pk))/(1.-xz**2)
         sp = pk/pkmrk
         xz = xz - sp
         avsp = abs(sp)
         if (avsp.gt.eps) go to 10
         a(is) = xz
         w(is) = (2.*(1.-xz**2))/(fk*pkm1)**2
   30 continue
      if (k.ne.kk*2) then
C
C For odd k computation of weight at the equator
C
         a(kk+1) = 0.
         pk = 2./fk**2
         do 40 n=2,k,2
            fn = n
            pk = pk*fn**2/(fn-1.)**2
   40    continue
         w(kk+1) = pk
      end if
C
C Complete the sets of abscissas and weights, using the symmetry.
C
      do 60 n=1,kk
         l = k + 1 - n
         a(l) = -a(n)
         w(l) = w(n)
   60 continue

      do n = 1,k
         ow(n) = w(n)
         oa(n) = a(n)
      end do

      return
      end

c##############################################################################

      subroutine bsslzr(bes,n)
C-----------------------------------------------------------------------
C
C Return n zeros (or if n>50, approximate zeros), of the Bessel function
C j0,in the array bes. The first 50 zeros will be given exactly, and the
C remaining zeros are computed by extrapolation,and therefore not exact.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, D. Williamson, August 1992
C
C-----------------------------------------------------------------------
      implicit none
      integer, parameter :: DBLKIND = selected_real_kind( p = 12 )
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer n              ! Number of zeros to return
C
C Output arguments
C
      real(DBLKIND) ::
     $  bes(n)               ! Array containing zeros of j0
C
C---------------------------Local workspace-----------------------------
C
      integer j,nn           ! loop indices
      real(DBLKIND) ::
     $  bz(50)               ! table of first 50 zeros
     $, pi                   ! 3.14.......
      save bz                ! ensure re-entrancy
C
      data bz           / 2.4048255577,   5.5200781103,
     $    8.6537279129,  11.7915344391,  14.9309177086,  18.0710639679,
     $   21.2116366299,  24.3524715308,  27.4934791320,  30.6346064684,
     $   33.7758202136,  36.9170983537,  40.0584257646,  43.1997917132,
     $   46.3411883717,  49.4826098974,  52.6240518411,  55.7655107550,
     $   58.9069839261,  62.0484691902,  65.1899648002,  68.3314693299,
     $   71.4729816036,  74.6145006437,  77.7560256304,  80.8975558711,
     $   84.0390907769,  87.1806298436,  90.3221726372,  93.4637187819,
     $   96.6052679510,  99.7468198587, 102.8883742542, 106.0299309165,
     $  109.1714896498, 112.3130502805, 115.4546126537, 118.5961766309,
     $  121.7377420880, 124.8793089132, 128.0208770059, 131.1624462752,
     $  134.3040166383, 137.4455880203, 140.5871603528, 143.7287335737,
     $  146.8703076258, 150.0118824570, 153.1534580192, 156.2950342685/
C
      pi = 4.*atan(1.)
      nn = n
      if (n.gt.50) then
         bes(50) = bz(50)
         do j=51,n
            bes(j) = bes(j-1) + pi
         end do
         nn = 49
      end if
      do j=1,nn
         bes(j) = bz(j)
      end do
      return
      end

