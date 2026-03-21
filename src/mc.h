c $Id$

      common /mc/ rhonot, t0, cldmin, small, 
     $            c, d, nos, pi, prhonos, thrpd, 
     $            rhos, gam3pd, gam4pd, rhow, rhoi,
     $            esi, esw,
     $            mcon01, mcon02, mcon03, mcon04,
     $            mcon05, mcon06, mcon07, mcon08,
     $            landm

      real(r8) rhonot
      real(r8) t0
      real(r8) cldmin
      real(r8) small
      real(r8) c
      real(r8) d
      real(r8) esi
      real(r8) esw
      real(r8) nos
      real(r8) pi
      real(r8) prhonos
      real(r8) thrpd
      real(r8) gam3pd
      real(r8) gam4pd
      real(r8) rhoi
      real(r8) rhos
      real(r8) rhow
      real(r8) mcon01
      real(r8) mcon02
      real(r8) mcon03
      real(r8) mcon04
      real(r8) mcon05
      real(r8) mcon06
      real(r8) mcon07
      real(r8) mcon08
      real(r8) landm(plond, plat)

      integer ktop  ! top level for pcw calculation

      common / mc_I / ktop

