c $Id$
      common / surf0 / gravit, rair

      real(r8)
     $  gravit       ! Gravitational acceleration
     $, rair         ! Gas constant for dry air
