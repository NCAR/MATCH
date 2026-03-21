module wetdep

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Wet deposition routines for both aerosols and gas phase constituents.
   ! 
   ! Author: module coded by B. Eaton.
   !-----------------------------------------------------------------------

   use precision
   use pmgrid

   implicit none
   save
   private
   public :: &
      wetdepa,     &! scavenging codes for very soluble aerosols
      wetdepg       ! scavenging of gas phase constituents by henry's law

#ifdef DEBUG
	 integer:: ilook = 142
	 integer:: jlook = 71
	 integer:: klook = 19
#endif

!##############################################################################
contains
!##############################################################################

   subroutine wetdepa( lat, t, p, q, pdel, &
                       cldt, cldc, cmfdqr, conicw, precs, conds, &
                       evaps, cwat, tracer, deltat, &
                       scavt, iscavt, cldv, fracis, sol_fact )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! scavenging code for very soluble aerosols
      ! 
      ! Author: P. Rasch
      !-----------------------------------------------------------------------

      implicit none

      integer, intent(in) ::&
         lat                    ! lat index
      real(r8), intent(in) ::&
         t(plond,plev),        &! temperature
         p(plond,plev),        &! pressure
         q(plond,plev),        &! moisture
         pdel(plond,plev),     &! pressure thikness
         cldt(plond,plevp),    &! total cloud fraction
         cldc(plond,plev),     &! convective cloud fraction
         cmfdqr(plond,plev),   &! rate of production of convective precip
         conicw(plond,plev),   &! convective cloud water
         cwat(plond,plev),     &! cloud water amount 
         precs(plond,plev),    &! rate of production of stratiform precip
         conds(plond,plev),    &! rate of production of condensate
         evaps(plond,plev),    &! rate of evaporation of precip
         cldv(plond,plev),     &! total cloud fraction
         deltat,               &! time step
         tracer(plond,plev),   &! trace species
         sol_fact               ! solubility factor (fraction of aer scavenged)

         
      real(r8), intent(out) ::&
         scavt(plond,plev),    &! scavenging tend 
         iscavt(plond,plev),   &! incloud scavenging tends
         fracis(plond,plev)     ! fraction of species not scavenged

      ! local variables

      integer i                 ! x index
      integer k                 ! z index

      real(r8) adjfac               ! factor stolen from cmfmca
      real(r8) aqfrac               ! fraction of tracer in aqueous phase
      real(r8) cmftau               ! time scale stolen from cmfmca
      real(r8) cwatc                ! local convective total water amount 
      real(r8) cwats                ! local stratiform total water amount 
      real(r8) cwatp                ! local water amount falling from above precip
      real(r8) fracev               ! fraction of precip from above that is evaporating
      real(r8) fracp                ! fraction of cloud water converted to precip
      real(r8) gafrac               ! fraction of tracer in gas phasea
      real(r8) grav                 ! gravity
      real(r8) hconst               ! henry's law solubility constant when equation is expressed
                                ! in terms of mixing ratios
      real(r8) mpla                 ! moles / liter H2O entering the layer from above
      real(r8) mplb                 ! moles / liter H2O leaving the layer below
      real(r8) omsm                 ! 1 - (a small number)
      real(r8) part                 !  partial pressure of tracer in atmospheres
      real(r8) patm                 ! total pressure in atmospheres
      real(r8) pdog                 ! work variable (pdel/grav)
      real(r8) precabc(plond)       ! conv precip from above (work array)
      real(r8) precabs(plond)       ! strat precip from above (work array)
      real(r8) precbl               ! precip falling out of level (work array)
      real(r8) precmin              ! minimum convective precip causing scavenging
      real(r8) rat                  ! ratio of amount available to amount removed
      real(r8) rhoh2o               ! density of water
      real(r8) scavab(plond)        ! scavenged tracer flux from above (work array)
      real(r8) scavabc(plond)       ! scavenged tracer flux from above (work array)
      real(r8) srcc                 ! tend for convective rain
      real(r8) srcs                 ! tend for stratiform rain
      real(r8) srct                 ! work variable
      real(r8) tracab(plond)        ! column integrated tracer amount
!      real(r8) vfall                ! fall speed of precip
      real(r8) fins                 ! fraction of rem. rate by strat rain
      real(r8) finc                 ! fraction of rem. rate by conv. rain
      real(r8) srcs1                ! work variable
      real(r8) srcs2                ! work variable
      real(r8) tc                   ! temp in celcius
      real(r8) weight               ! fraction of condensate which is ice
      real(r8) cldmabs(plond)       ! maximum cloud at or above this level
      real(r8) cldmabc(plond)       ! maximum cloud at or above this level
      real(r8) odds                 ! limit on removal rate (proportional to prec)
      real(r8) dblchek
      ! ------------------------------------------------------------------------

!      omsm = 1.-1.e-10          ! used to prevent roundoff errors below zero
      omsm = 1.-2*epsilon(1._r8) ! used to prevent roundoff errors below zero
      precmin =  0.1/8.64e4      ! set critical value to 0.1 mm/day in kg/m2/s

      cmftau = 3600.
      grav = 9.81               ! gravity
      rhoh2o = 1000.            ! density of water
      adjfac = deltat/(max(deltat,cmftau)) ! adjustment factor from hack scheme

      ! assume 4 m/s fall speed currently (should be improved)
!      vfall = 4.

      ! this section of code is for highly soluble aerosols,
      ! the assumption is that within the cloud that
      ! all the tracer is in the cloud water
      !
      ! for both convective and stratiform clouds, 
      ! the fraction of cloud water converted to precip defines
      ! the amount of tracer which is pulled out.
      !

      do i = 1,plond
         precabs(i) = 0
         precabc(i) = 0
         scavab(i) = 0
         scavabc(i) = 0
         tracab(i) = 0
         cldmabs(i) = 0
         cldmabc(i) = 0
      end do

      do k = 1,plev
         do i = 1,plon
            tc     = t(i,k) - 273.16
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice
            weight = 0.                                 ! assume no ice

            pdog = pdel(i,k)/grav

            ! calculate the fraction of strat precip from above 
            !                 which evaporates within this layer
            fracev = evaps(i,k)*pdel(i,k)/grav &
                     /max(1.e-12_r8,precabs(i))

            ! now do the convective scavenging

            ! set odds proportional to fraction of the grid box that is swept by the 
            ! precipitation =precabc/rhoh20*(area of sphere projected on plane
            !                                /volume of sphere)*deltat
            ! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
            ! unless the fraction of the area that is cloud is less than odds, in which
            ! case use the cloud fraction (assumes precabs is in kg/m2/s)
            ! is really: precabs*3/4/1000./1e-3*deltat
            ! here I use .1 from Balkanski
            !
            ! use a local rate of convective rain production for incloud scav
            !odds=max(min(1._r8, &
            !     cmfdqr(i,k)*pdel(i,k)/grav*0.1_r8*deltat),0._r8)
            !++mcb -- change cldc to cldt; change cldt to cldv (9/17/96)
            !            srcs1 =  cldt(i,k)*odds*tracer(i,k)*(1.-weight) &
            ! srcs1 =  cldv(i,k)*odds*tracer(i,k)*(1.-weight) &
            !srcs1 =  cldc(i,k)*odds*tracer(i,k)*(1.-weight) &
            !         /deltat 

            ! fraction of convective cloud water converted to rain
            fracp = cmfdqr(i,k)*deltat/max(1.e-8_r8,conicw(i,k))
            ! note cmfdrq can be negative from evap of rain, so constrain it
            fracp = max(min(1._r8,fracp),0._r8)
            ! remove that amount from within the convective area
!           srcs1 = cldc(i,k)*fracp*tracer(i,k)*(1.-weight)/deltat ! liquid only
!           srcs1 = cldc(i,k)*fracp*tracer(i,k)/deltat             ! any condensation
!           srcs1 = 0.
            srcs1 = sol_fact*cldt(i,k)*fracp*tracer(i,k)*(1.-weight)/deltat ! liquid only

            !--mcb

            ! scavenge below cloud

            !            cldmabc(i) = max(cldc(i,k),cldmabc(i))
            !            cldmabc(i) = max(cldt(i,k),cldmabc(i))
            cldmabc(i) = max(cldv(i,k),cldmabc(i))
            cldmabc(i) = cldv(i,k)

            odds=max( &
                 min(1._r8,precabc(i)/max(cldmabc(i),1.e-5_r8) &
                 *0.1_r8*deltat),0._r8)
            srcs2 = sol_fact*cldmabc(i)*odds*tracer(i,k)/deltat

            srcc = srcs1 + srcs2  ! convective tend by both processes
            finc = srcs1/(srcc + 1.e-36_r8)

            ! now do the stratiform scavenging

            ! incloud scavenging

            ! fracp is the fraction of cloud water converted to precip
            fracp =  precs(i,k)*deltat/max(cwat(i,k),1.e-12_r8)
            fracp = max(0._r8,min(1._r8,fracp))
!           fracp = 0.     ! for debug

            ! assume the corresponding amnt of tracer is removed
            !++mcb -- remove cldc; change cldt to cldv 
            !            srcs1 = (cldt(i,k)-cldc(i,k))*fracp*tracer(i,k)/deltat
            !            srcs1 = cldv(i,k)*fracp*tracer(i,k)/deltat &
!           srcs1 = cldt(i,k)*fracp*tracer(i,k)/deltat            ! all condensate
            srcs1 = sol_fact*cldt(i,k)*fracp*tracer(i,k)/deltat*(1.-weight)! only the liquid phase



            ! below cloud scavenging

!           volume undergoing below cloud scavenging
            cldmabs(i) = cldv(i,k)   ! precipitating volume
!           cldmabs(i) = cldt(i,k)   ! local cloud volume

            odds = precabs(i)/max(cldmabs(i),1.e-5_r8)*0.1_r8*deltat
            odds = max(min(1._r8,odds),0._r8)
            srcs2 =sol_fact*(cldmabs(i)*odds) *tracer(i,k)/deltat

            srcs = srcs1 + srcs2             ! total stratiform scavenging
            fins=srcs1/(srcs + 1.e-36_r8)    ! fraction taken by incloud processes

            ! make sure we dont take out more than is there
            ! ratio of amount available to amount removed
            rat = tracer(i,k)/max(deltat*(srcc+srcs),1.e-36_r8)
            if (rat.lt.1.) then
               srcs = srcs*rat
               srcc = srcc*rat
            endif
            srct = (srcc+srcs)*omsm

            
            ! fraction that is not removed within the cloud
            ! (assumed to be interstitial, and subject to convective transport)
            fracp = deltat*srct/max(cldmabs(i)*tracer(i,k),1.e-36_r8)  ! amount removed
            fracp = max(0._r8,min(1._r8,fracp))
            fracis(i,k) = 1. - fracp

            ! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
            scavt(i,k) = -srct + fracev*scavab(i)*grav/pdel(i,k)
            iscavt(i,k) = -(srcc*finc + srcs*fins)*omsm

            dblchek = tracer(i,k) + deltat*scavt(i,k)
            if (dblchek .lt. 0.) then
               write (6,*) ' wetdapa: negative value ', i, k, tracer(i,k), &
                    dblchek, scavt(i,k), srct, rat, fracev
            endif

            ! now keep track of scavenged mass and precip
            scavab(i) = scavab(i)*(1-fracev) + srcs*pdel(i,k)/grav
            precabs(i) = precabs(i) + (precs(i,k) - evaps(i,k))*pdel(i,k)/grav
            scavabc(i) = scavabc(i) + srcc*pdel(i,k)/grav
            precabc(i) = precabc(i) + (cmfdqr(i,k))*pdel(i,k)/grav
            tracab(i) = tracab(i) + tracer(i,k)*pdel(i,k)/grav

         end do
      end do

   end subroutine wetdepa

!##############################################################################

   subroutine wetdepg( lat, t, p, q, pdel, &
                       cldt, cldc, cmfdqr, precs, evaps, &
                       rain, cwat, tracer, deltat, molwt, &
                       solconst, scavt, iscavt, cldv, icwmr1, &
                       icwmr2, fracis )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! scavenging of gas phase constituents by henry's law
      ! 
      ! Author: P. Rasch
      !-----------------------------------------------------------------------

      implicit none

      integer, intent(in) ::&
         lat                    ! lat index
      real(r8), intent(in) ::&
         t(plond,plev),        &! temperature
         p(plond,plev),        &! pressure
         q(plond,plev),        &! moisture
         pdel(plond,plev),     &! pressure thikness
         cldt(plond,plev),     &! total cloud fraction
         cldc(plond,plev),     &! convective cloud fraction
         cmfdqr(plond,plev),   &! rate of production of convective precip
         rain (plond,plev),    &! total rainwater mixing ratio
         cwat(plond,plev),     &! cloud water amount 
         precs(plond,plev),    &! rate of production of stratiform precip
         evaps(plond,plev),    &! rate of evaporation of precip
         cldv(plond,plev),     &! estimate of local volume occupied by clouds
         icwmr1 (plond,plev),  &! in cloud water mixing ration for zhang scheme
         icwmr2 (plond,plev),  &! in cloud water mixing ration for hack  scheme
         deltat,               &! time step
         tracer(plond,plev),   &! trace species
         molwt                  ! molecular weights

      real(DBLKIND) &
         solconst(plond,plev)   ! Henry's law coefficient

      real(r8), intent(out) ::&
         scavt(plond,plev),    &! scavenging tend 
         iscavt(plond,plev),   &! incloud scavenging tends
         fracis(plond, plev)    ! fraction of constituent that is insoluble

      ! local variables

      integer i                 ! x index
      integer k                 ! z index

      real(r8) adjfac               ! factor stolen from cmfmca
      real(r8) aqfrac               ! fraction of tracer in aqueous phase
      real(r8) cmftau               ! time scale stolen from cmfmca
      real(r8) cwatc                ! local convective total water amount 
      real(r8) cwats                ! local stratiform total water amount 
      real(r8) cwatl                ! local cloud liq water amount 
      real(r8) cwatp                ! local water amount falling from above precip
      real(r8) cwatpl               ! local water amount falling from above precip (liq)
      real(r8) cwatt                ! local sum of strat + conv total water amount 
      real(r8) cwatti               ! cwatt/cldv = cloudy grid volume mixing ratio
      real(r8) fracev               ! fraction of precip from above that is evaporating
      real(r8) fracp                ! fraction of cloud water converted to precip
      real(r8) gafrac               ! fraction of tracer in gas phasea
      real(r8) grav                 ! gravity
      real(r8) hconst               ! henry's law solubility constant when equation is expressed
                                ! in terms of mixing ratios
      real(r8) molwta               ! molecular weight of air
      real(r8) mpla                 ! moles / liter H2O entering the layer from above
      real(r8) mplb                 ! moles / liter H2O leaving the layer below
      real(r8) omsm                 ! 1 - (a small number)
      real(r8) part                 !  partial pressure of tracer in atmospheres
      real(r8) patm                 ! total pressure in atmospheres
      real(r8) pdog                 ! work variable (pdel/grav)
      real(r8) precab(plond)        ! precip from above (work array)
      real(r8) precbl               ! precip work variable
      real(r8) precxx               ! precip work variable
      real(r8) precic               ! precip work variable
      real(r8) rat                  ! ratio of amount available to amount removed
      real(r8) rhoh2o               ! density of water
      real(r8) scavab(plond)        ! scavenged tracer flux from above (work array)
      real(r8) scavabc(plond)       ! scavenged tracer flux from above (work array)
      !      real(r8) vfall                ! fall speed of precip
      real(r8) scavmax              ! an estimate of the max tracer avail for removal
      real(r8) scavbl               ! flux removed at bottom of layer
      real(r8) fins                 ! in cloud fraction removed by strat rain
      real(r8) finc                 ! in cloud fraction removed by conv rain
      real(r8) rate                 ! max removal rate estimate
      real(r8) scavlimt             ! limiting value 1
      real(r8) scavt1               ! limiting value 2
      real(r8) scavin               ! scavenging by incloud processes
      real(r8) scavbc               ! scavenging by below cloud processes
      real(r8) tc
      real(r8) weight               ! ice fraction
      real(r8) wtpl                 ! work variable
      real(r8) cldmabs(plond)       ! maximum cloud at or above this level
      real(r8) cldmabc(plond)       ! maximum cloud at or above this level
      !-----------------------------------------------------------

      omsm = 1.-2*epsilon(1._r8)   ! used to prevent roundoff errors below zero

      cmftau = 3600.
      grav = 9.81               ! gravity
      rhoh2o = 1000.            ! density of water
      molwta = 28.97            ! molecular weight dry air gm/mole
      adjfac = deltat/(max(deltat,cmftau)) ! adjustment factor from hack scheme

      ! assume 4 m/s fall speed currently (should be improved)
      !      vfall = 4.

      ! zero accumulators
      do i = 1,plond
         precab(i) = 1.e-36_r8
         scavab(i) = 0.
         cldmabs(i) = 0.
      end do


      do k = 1,plev
         do i = 1,plon

            tc     = t(i,k) - 273.16
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice

            cldmabs(i) = max(cldmabs(i),cldt(i,k))

            ! partitioning coefs for gas and aqueous phase
            !              take as a cloud water amount, the sum of the stratiform amount
            !              plus the convective rain water amount 

            ! convective amnt is just the local precip rate from the hack scheme
            !              since there is no storage of water, this ignores that falling from above
            !            cwatc = cmfdqr(i,k)*deltat/adjfac
            !++mcb -- test cwatc
            cwatc = (icwmr1(i,k) + icwmr2(i,k)) * (1.-weight)
            !--mcb 

            ! strat cloud water amount and also ignore the part falling from above
            cwats = cwat(i,k) 

            ! cloud water as liq
            !++mcb -- add cwatc later (in cwatti)
            !            cwatl = (1.-weight)*(cwatc+cwats)
            cwatl = (1.-weight)*cwats
            ! cloud water as ice
            !*not used        cwati = weight*(cwatc+cwats)

            ! total suspended condensate as liquid
            cwatt = cwatl + rain(i,k)

            ! incloud version 
            !++mcb -- add cwatc here
            cwatti = cwatt/max(cldv(i,k), 0.00001_r8) + cwatc

            ! partitioning terms
            patm = p(i,k)/1.013e5 ! pressure in atmospheres
            hconst = molwta*patm*solconst(i,k)*cwatti/rhoh2o
            aqfrac = hconst/(1.+hconst)
            gafrac = 1/(1.+hconst)
            fracis(i,k) = gafrac


            ! partial pressure of the tracer in the gridbox in atmospheres
            part = patm*gafrac*tracer(i,k)*molwta/molwt

            ! use henrys law to give moles tracer /liter of water
            ! in this volume 
            ! then convert to kg tracer /liter of water (kg tracer / kg water)
            mplb = solconst(i,k)*part*molwt/1000.


            pdog = pdel(i,k)/grav

            ! this part of precip will be carried downward but at a new molarity of mpl 
            precic = pdog*(precs(i,k) + cmfdqr(i,k))

            ! we cant take out more than entered, plus that available in the cloud
            !                  scavmax = scavab(i)+tracer(i,k)*cldt(i,k)/deltat*pdog
            scavmax = scavab(i)+tracer(i,k)*cldv(i,k)/deltat*pdog

            ! flux of tracer by incloud processes
            scavin = precic*(1.-weight)*mplb

            ! fraction of precip which entered above that leaves below
            precxx = precab(i)-pdog*evaps(i,k)

            ! flux of tracer by below cloud processes
            !++mcb -- removed wtpl because it is now not assigned and previously
            !          when it was assigned it was unnecessary:  if(tc.gt.0)wtpl=1
            if (tc.gt.0) then
               !               scavbc = precxx*wtpl*mplb ! if liquid
               scavbc = precxx*mplb ! if liquid
            else
               scavbc = scavab(i)*precxx/(precab(i)) ! if ice
            endif

            scavbl = min(scavbc + scavin, scavmax)

            ! first guess assuming that henries law works
            scavt1 = (scavab(i)-scavbl)/pdog*omsm

#ifdef DEBUG
            if (lat.eq.jlook.and.i.eq.ilook.and.k.eq.klook) then
               write (6,*) ' wetdepg: blim is ', scavt1, tracer(i,k),cldv(i,k), tracer(i,k), &
                    tracer(i,k)*cldv(i,k)/deltat
            endif
#endif

            ! pjr this should not be required, but we put it in to make sure we cant remove too much
            ! remember, scavt1 is generally negative (indicating removal)
            scavt1 = max(scavt1,-tracer(i,k)*cldv(i,k)/deltat)

#ifdef DEBUG
            if (lat.eq.jlook.and.i.eq.ilook.and.k.eq.klook) then
               write (6,*) ' wetdepg: alim is ', scavt1
            endif
#endif

            !++mcb -- remove this limitation for gas species
            !c use the dana and hales or balkanski limit on scavenging
            !c            rate = precab(i)*0.1
            !            rate = (precic + precxx)*0.1
            !            scavlimt = -tracer(i,k)*cldv(i,k)
            !     $           *rate/(1.+rate*deltat)

            !            scavt(i,k) = max(scavt1, scavlimt)

            ! instead just set scavt to scavt1
            scavt(i,k) = scavt1
            !--mcb

            ! now update the amount leaving the layer
            scavbl = scavab(i) - scavt(i,k)*pdog 

            ! in cloud amount is that formed locally over the total flux out bottom
            fins = scavin/(scavin + scavbc + 1.e-36_r8)
            iscavt(i,k) = scavt(i,k)*fins

            scavab(i) = scavbl
            precab(i) = max(precxx + precic,1.e-36_r8)
            
         end do
      end do
      
   end subroutine wetdepg

!##############################################################################

end module wetdep
