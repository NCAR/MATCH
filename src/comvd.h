c $Id$
C
C Constants used in vertical diffusion and pbl
C
      common/comvd/cpair,cpvir,gravit,rair,zkmin,ml2(plevp),
c++bee
     $ qmincg(pcnst),
c--bee
     $             ntopfl,npbl
C
      real(r8) cpair,    ! specific heat of dry air
     $     cpvir,    ! derived constant for cp moist air
     $     gravit,   ! acceleration due to gravity
     $     rair,     ! gas const for dry air
     $     zkmin,    ! minimum kneutral*f(ri)
     $     ml2       ! mixing lengths squared
c++bee
     $,    qmincg    ! Min. constituent concentration counter-gradient term
c--bee
      integer
     $     ntopfl,   ! top level to which vertical diffusion is applied.
     $     npbl      ! maximum number of levels in pbl from surface
C
