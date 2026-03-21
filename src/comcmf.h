c
c $Id$
c $Author$
c
C
C Common block for moist convective mass flux procedure
C
      common/comcmf/cp      ,hlat    ,grav    ,c0      ,betamn  ,
     $              rhlat   ,rcp     ,rgrav   ,cmftau  ,rhoh2o  ,
     $              rgas    ,dzmin   ,limcnv  ,iloc    ,jloc    ,
     $              nsloc   ,tinyalt    ,eps     ,tpmax   ,shpmax  ,
     $              rlxclm
C
      real(r8) cp          ! specific heat of dry air
      real(r8) hlat        ! latent heat of vaporization
      real(r8) grav        ! gravitational constant       
      real(r8) c0          ! rain water autoconversion coefficient
      real(r8) betamn      ! minimum overshoot parameter
      real(r8) rhlat       ! reciprocal of hlat
      real(r8) rcp         ! reciprocal of cp
      real(r8) rgrav       ! reciprocal of grav
      real(r8) cmftau      ! characteristic adjustment time scale
      real(r8) rhoh2o      ! density of liquid water (STP)
      real(r8) rgas        ! gas constant for dry air
      real(r8) dzmin       ! minimum convective depth for precipitation
      real(r8) tinyalt        ! arbitrary small num used in transport estimates
      real(r8) eps         ! convergence criteria (machine dependent)
      real(r8) tpmax       ! maximum acceptable t perturbation (degrees C)
      real(r8) shpmax      ! maximum acceptable q perturbation (g/g)           
C
      integer iloc        ! longitude location for diagnostics
      integer jloc        ! latitude  location for diagnostics
      integer nsloc       ! nstep for which to produce diagnostics
      integer limcnv      ! top interface level limit for convection
C
      logical rlxclm      ! logical to relax column versus cloud triplet
C
 
