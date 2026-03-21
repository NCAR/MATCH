c $Id$
C Data for virtual temperature calculation

      common /virtem0/ zvir

      real(r8)
     $  zvir         ! rh2o/rair - 1.

