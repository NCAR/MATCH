c     $Header$ -*-fortran-*-

c     Purpose: Common block dstlsm stores LSM info needed by dust parameterization
c     dstlsm is initialized in CCM:control/ccm3() by dstbdlsm()
c     dstlsm most emphatically does NOT need params.h nor lsmpar.h to work

c     Usage:
c     #include <dstlsm.h> /* LSM block data */

c#ifndef DSTLSM_H
c#define DSTLSM_H

c     Contents of lsmpar.h
      integer mband   !number of solar radiation bands: vis, nir
      integer msl     !number of soil layers
      integer mst     !number of "soil" types (soil, ice, 2 lakes, wetland)
      integer mvt     !number of plant types
      integer msc     !number of soil color types
 
      parameter (mband=2, msl=6, mst=5, mvt=14, msc=9)

c     Contents of vegtyp.h with cover <-> pln_frc and plant <-> pln_typ
c     and vegcon_i <-> dstvegcon_i and vegcon_r <-> dstvegcon_r
      common / dstvegtyp_i / pln_typ ! dstvegtyp_i is initialized in dstbdlsm()
      common / dstvegtyp_r / pln_frc ! dstvegtyp_r is initialized in dstbdlsm()
 
      real(r8) pln_frc(0:28,3)      ! [frc] Weight of corresponding plant type (sums to 1.0)
      integer pln_typ(0:28,3)   ! [idx] LSM plant type (1..14 = nbr_LSM_pln_typ)

* ------------------------ code history ---------------------------
* source file:       vegcon.h
* purpose:           vegetation type constants 
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
* -----------------------------------------------------------------

      common /dstvegcon_i/ nic, noveg

      integer nic       !value for irrigated crop 
      integer noveg     !value for not vegetated 

      common /dstvegcon_r/ vw, rdp, ch2op, dleaf
 
      real(r8) vw(mvt)      !btran exponent: [(h2osoi-watdry)/(watopt-watdry)]**vw 
      real(r8) rdp(mvt)     !defines root fraction decrease with depth
      real(r8) ch2op(mvt)   !maximum intercepted h2o per unit lai+sai (mm) 
      real(r8) dleaf(mvt)   !characteristic leaf dimension (m) 

      common /dstvegcon_r/ c3psn  , kc25   , akc   , ko25  , ako  ,
     &                  vcmx25 , avcmx  , bp    , mp    , qe25 ,
     &                  aqe    , rmf25  , rms25 , rmr25 , arm  ,
     &                  dmcf   , folnmx , tmin

      real(r8) c3psn(mvt)   !photosynthetic pathway: 0. = c4, 1. = c3
      real(r8) kc25(mvt)    !co2 michaelis-menten constant at 25c (pa) 
      real(r8) akc(mvt)     !q10 for kc25 
      real(r8) ko25(mvt)    !o2 michaelis-menten constant at 25c (pa) 
      real(r8) ako(mvt)     !q10 for ko25 
      real(r8) vcmx25(mvt)  !maximum rate of carboxylation at 25c (umol co2/m**2/s)
      real(r8) avcmx(mvt)   !q10 for vcmx25 
      real(r8) bp(mvt)      !minimum leaf conductance (umol/m**2/s) 
      real(r8) mp(mvt)      !slope of conductance-to-photosynthesis relationship 
      real(r8) qe25(mvt)    !quantum efficiency at 25c (umol co2 / umol photon)
      real(r8) aqe(mvt)     !q10 for qe25 
      real(r8) rmf25(mvt)   !leaf maintenance respiration at 25c (umol co2/m**2/s)
      real(r8) rms25(mvt)   !stem maintenance respiration at 25c (umol co2/kg bio/s)
      real(r8) rmr25(mvt)   !root maintenance respiration at 25c (umol co2/kg bio/s)
      real(r8) arm(mvt)     !q10 for maintenance respiration
      real(r8) dmcf(mvt)    !co2-to-biomass conversion factor (ug biomass/umol co2) 
      real(r8) folnmx(mvt)  !foliage nitrogen concentration when f(n)=1 (%) 
      real(r8) tmin(mvt)    !minimum temperature for photosynthesis (kelvin) 

      common /dstvegcon_r/ xl, rhol, rhos, taul, taus

      real(r8) xl(mvt)          !leaf/stem orientation index
      real(r8) rhol(mvt,mband)  !leaf reflectance: 1=vis, 2=nir 
      real(r8) rhos(mvt,mband)  !stem reflectance: 1=vis, 2=nir 
      real(r8) taul(mvt,mband)  !leaf transmittance: 1=vis, 2=nir 
      real(r8) taus(mvt,mband)  !stem transmittance: 1=vis, 2=nir 

      common /dstvegcon_r/ binvvt, z0mvt, zpdvt  , stembvt, rootbvt,
     &                  folnvt, mrp  , soilcvt, cwpvt

      real(r8) binvvt(mvt)  !1/vkc*ln(z0m/z0h)
      real(r8) z0mvt(mvt)   !momentum roughness length (m) 
      real(r8) zpdvt(mvt)   !displacement height (m) 
      real(r8) stembvt(mvt) !stem biomass (kg /m**2)
      real(r8) rootbvt(mvt) !root biomass (kg /m**2)
      real(r8) folnvt(mvt)  !foliage nitrogen concentration (%)
      real(r8) mrp(mvt)     !microbial respiration parameter (umol co2 /kg c/ s)
      real(r8) soilcvt(mvt) !soil carbon (kg c /m**2)
      real(r8) cwpvt(mvt)   !empirical canopy wind parameter

      common /dstvegcon_r/ tai, gai, hvt, hvb

      real(r8) tai(mvt,12)  !monthly leaf area index + stem area index, one-sided
      real(r8) gai(mvt,12)  !monthly leaf area index, one-sided
      real(r8) hvt(mvt)     !top of canopy (m)
      real(r8) hvb(mvt)     !bottom of canopy (m)

* ------------------------ end vegcon.h ---------------------------

c#endif /* ! DSTLSM_H */ 
