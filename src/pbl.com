c $Id$
c
c pbl constants
c
      real(r8)
     $     betam,  ! cnstnt in wind gradient expression
     $     betas,  ! cnstnt in surface layer gradient expression
     $     betah,  ! cnstnt in temp gradient expression 
     $     fak,    ! constant in srfc tmp excess         
     $     g,      ! gravitational acceleration
     $     onet,   ! 1/3 power in wind gradient expression
     $     fakn,   ! constant in turbulent prandtl number
     $     ricr,   ! critical richardson number
     $     sffrac, ! surface layer fraction of boundary layer
     $     vk,     ! von karman's constant
     $     ccon,   ! fak * srfface * vk
     $     binm,   ! betam * srffac
     $     binh    ! betah * srffac
c
      common /pbl/
     $            betam , betas, betah,   fak,    g,
     $            onet  ,  fakn,  ricr,sffrac,   vk,
     $            ccon  ,  binm,  binh
c
