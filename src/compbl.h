c $Id$
C
C Pbl constants
C
      common /compbl/ betam   ,betas   ,betah   ,fak     ,g       ,
     $                onet    ,fakn    ,ricr    ,sffrac  ,vk      ,
     $                ccon    ,binm    ,binh
C
      real(r8) betam,  ! constant in wind gradient expression
     $     betas,  ! constant in surface layer gradient expression
     $     betah,  ! constant in temperature gradient expression 
     $     fak,    ! constant in surface temperature excess         
     $     g,      ! gravitational acceleration
     $     onet,   ! 1/3 power in wind gradient expression
     $     fakn,   ! constant in turbulent prandtl number
     $     ricr,   ! critical richardson number
     $     sffrac, ! surface layer fraction of boundary layer
     $     vk,     ! von karman's constant
     $     ccon,   ! fak * sffrac * vk
     $     binm,   ! betam * sffrac
     $     binh    ! betah * sffrac
C
