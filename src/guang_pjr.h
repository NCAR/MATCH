      real(r8) a
      real(r8) b
      real(r8) eps1
      real(r8) c1
      real(r8) c2
      real(r8) c3
      real(r8) tfreez
      real(r8) grav
      real(r8) rgrav
      real(r8) cpg
      real(r8) rgas
      real(r8) tau
      real(r8) dd
      real(r8) rl                   ! latent heat of vap
      real(r8) qmin
      common /guang_pjr_r/a,b,eps1,c1,c2,c3,tfreez,
     $     grav, rgrav, cpg, rgas, tau, dd, rl, qmin

      logical trigon
      common /guang_pjr_l/ trigon

      integer msg
#ifdef DEBCONV
      integer latlook
      integer ilook
#endif
      common /guang_pjr_i/ msg
#ifdef DEBCONV
     $,  latlook, ilook
#endif

