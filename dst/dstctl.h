c     $Header$ -*-fortran-*- 
      
c     Purpose: dstctl.h contains common blocks which store model control variables
      
c     These variables are initialized in aer()
      
c     Usage:
c     #include <dstctl.h> /* Model control variables */
      
c     CCM/MATCH control variables
      integer mcdate            ! [day] Current date in YYMMDD format
      integer mcsec             ! [s] Seconds past mcdate at 0Z
      integer nbdate            ! [day] Simulation start date in YYMMDD format
      integer nbsec             ! [s] Simulation start second relative to nbdate
      integer ndcur             ! [day] Current day number of simulation
      integer nscur             ! [s] Seconds relative to ndcur
      integer nstep             ! [idx] Current timestep number
      common /aer_ctl_cmn/      ! aer_ctl_cmn is initialized in aer()
     $     mcdate,              ! [day] Current date in YYMMDD format
     $     mcsec,               ! [s] Seconds past mcdate at 0Z
     $     nbdate,              ! [day] Simulation start date in YYMMDD format
     $     nbsec,               ! [s] Simulation start second relative to nbdate
     $     ndcur,               ! [day] Current day number of simulation
     $     nscur,               ! [s] Seconds relative to ndcur
     $     nstep                ! [idx] Current timestep number

      
