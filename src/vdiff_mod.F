      module vdiff_mod

      use precision

      real(r8) :: zero = 0.

      contains

c $Id$
#include <params.h>
c++bee, 11 August 1998
c Vertical diffusion/PBL codes from ccm3.6.8 modified for MATCH as follows:
c * put qmincg into /comvd/
c * remove assumption that water vapor is the first constituent in the
c   constituent and surface flux arrays.  pass in arrays for water vapor
c   in separate arguments.  Change loops over the non-water species to
c   "do m=1,pcnst" instead of "do m=2,pcnst".
c * use MATCH's version of virtem.
c * could speed up and reduce memory requirements by removing calcs for
c   diffusion of momentum coefficients which aren't being used (in vdiff).
c--bee

c##############################################################################

      subroutine vdinti(cpwvx   ,cpairx  ,gravx   ,rairx   
c++bee
     $, qmin, hyam, hybm )
c--bee
C-----------------------------------------------------------------------
C
C Initialization of time independent fields for vertical diffusion.
C Call initialization routine for boundary layer scheme.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      J. Rosinski, June 1992
C Reviewed:          P. Rasch, B. Boville, August 1992
C Reviewed:          P. Rasch, April 1996
C Reviewed:          B. Boville, April 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comvd.h>
C-----------------------------------------------------------------------
c++bee calc hypm locally
cinclude <comhyb.h>
c--bee
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) cpwvx        ! spec. heat of water vapor at const. pressure
      real(r8) cpairx       ! specific heat of dry air
      real(r8) gravx        ! acceleration due to gravity
      real(r8) rairx        ! gas constant for dry air
c++bee
      real(r8)
     $  qmin(pcnst)    ! Global minimum constituent concentration
     $, hyam(plev)     ! p0 component of hybrid coordinate - midpoints
     $, hybm(plev)     ! ps component of hybrid coordinate - midpoints
c--bee
C
C---------------------------Local workspace-----------------------------
C
      integer k         ! vertical loop index
c++bee
      integer m
      real(r8)
     $  ps0             ! Base state surface pressure (pascals)
     $, hypm(plev)
c--bee
C
C-----------------------------------------------------------------------
C
C Hard-wired numbers.
C zkmin = minimum k = kneutral*f(ri)
C
      zkmin = 0.01
C
C Set physical constants for vertical diffusion and pbl:
C
      cpair  = cpairx
      gravit = gravx
      rair   = rairx
C
C Derived constants
C ntopfl = top level to which v-diff is applied
C npbl = max number of levels (from bottom) in pbl
C
c
c Limit pbl height to regions below 400 mb
c
c++bee
      ps0 = 1.0e5
      do k=1,plev
        hypm(k) = (hyam(k) + hybm(k))*ps0
      end do
c--bee
      npbl = 0
      do k=plev,1,-1
        if (hypm(k).ge.4.e4) then
          npbl = npbl + 1
        end if
      end do
      npbl = max(npbl,1)
      if (masterproc) then
        write(6,*)'VDINTI: PBL height will be limited to bottom ',npbl,
     $' model levels. Top is ',hypm(plevp-npbl),' pascals'
      end if

      ntopfl = 1
      if (plev.eq.1) ntopfl = 0
      cpvir = cpwvx/cpairx - 1.
C
C Set the square of the mixing lengths.
C
      ml2(1) = 0.
      do k=2,plev
        ml2(k) = 30.0**2
      end do
      ml2(plevp) = 0.
c++bee
c Set the minimum mixing ratio for the counter-gradient term.  Normally this
c should be the same as qmin.

      do m = 1, pcnst
         qmincg(m) = qmin(m)
      end do
c--bee
C
C Initialize pbl variables
C
      call pbinti(gravx)
C
      return
      end subroutine vdinti

c##############################################################################

      subroutine pbinti(gravx)
C-----------------------------------------------------------------------
C
C Initialize time independent variables of pbl package.
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Boville
C Standardized:      J. Rosinski, June 1992
C Reviewed:          B. Boville, P. Rasch, August 1992
C Reviewed:          B. Boville, P. Rasch, April 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c
#include <implicit.h>
C------------------------------Commons----------------------------------
#include <compbl.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) gravx     !  acceleration of gravity
C
C Basic constants
C
      g    = gravx
      onet = 1./3.
c++bee
c The following PBL constants are initialized in blkdat in CCM3.
      vk = 0.40
      fak = 8.5
      fakn = 7.2
      betam = 15.0
      betas = 5.0
      betah = 15.0
      sffrac = 0.1
      ricr = 0.3
c--bee
C
C Derived constants
C
      ccon = fak*sffrac*vk
      binm = betam*sffrac
      binh = betah*sffrac
C
      return
      end subroutine pbinti

c##############################################################################

      subroutine vdiff(lat     ,um1     ,vm1     ,tm1     ,qm1     ,
     $                 pmidm1  ,pintm1  ,rpdel   ,rpdeli  ,ztodt   ,
     $                 thm     ,zm      ,taux    ,tauy    ,shflx   ,
     $                 cflx    ,up1     ,vp1     ,thp     ,qp1     ,
     $                 pblh    ,ustar   ,kvh     ,kvm     ,tpert   ,
     $                 qpert   ,cgs
c++bee
     $, shm1, shp1, wvflx, wpert, obklen )
c--bee
C-----------------------------------------------------------------------
C
C Driver routine to compute vertical diffusion of momentum,
C moisture, trace constituents and potential temperature.
C
C Free atmosphere diffusivities are computed first; then modified
C by the boundary layer scheme; then passed to individual
C parameterizations mvdiff, qvdiff
C
C The free atmosphere diffusivities are based on standard mixing length 
C forms for the neutral diffusivity multiplied by functns of Richardson 
C number. K = l^2 * |dV/dz| * f(Ri). The same functions are used for 
C momentum, potential temperature, and constitutents.
C The stable Richardson num function (Ri>0) is taken from Holtslag and 
C Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*Ri*(1 + 8*Ri))
C The unstable Richardson number function (Ri<0) is taken from  CCM1.
C f = sqrt(1 - 18*Ri)
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          P. Rasch, B. Boville, August 1992
C Reviewed:          P. Rasch, March 1996
C Reviewed:          B. Boville, April 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comvd.h>
C-----------------------------------------------------------------------
c++bee,  qmincg moved to /comvd/
cinclude <comqmin.h>
c--bee
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat                 ! latitude index

      real(r8) um1(plond,plev)        ! u wind input
      real(r8) vm1(plond,plev)        ! v wind input
      real(r8) tm1(plond,plev)        ! temperature input
      real(r8) qm1(plond,plev,pcnst)  ! moisture and trace constituent input
      real(r8) pmidm1(plond,plev)     ! midpoint pressures
      real(r8) pintm1(plond,plevp)    ! interface pressures
      real(r8) rpdel(plond,plev)      ! 1./pdel  (thickness bet interfaces)
      real(r8) rpdeli(plond,plev)     ! 1./pdeli (thickness bet midpoints)
      real(r8) ztodt                  ! 2 delta-t
      real(r8) thm(plond,plev)        ! potential temperature
      real(r8) zm(plond,plev)         ! midpoint geoptl height above sfc
      real(r8) taux(plond)            ! x surface stress (n)
      real(r8) tauy(plond)            ! y surface stress (n)
      real(r8) shflx(plond)           ! surface sensible heat flux (w/m2)
      real(r8) cflx(plond,pcnst)      ! surface constituent flux (kg/m2/s)
c++bee
      real(r8) shm1(plond,plev)       ! specific humidity (kg/kg)
      real(r8) wvflx(plond)           ! water vapor flux (kg/m2/s)
c--bee
C
C Output arguments
C
      real(r8) up1(plond,plev)        ! u-wind after vertical diffusion
      real(r8) vp1(plond,plev)        ! v-wind after vertical diffusion
      real(r8) thp(plond,plev)        ! pot temp after vert. diffusion
      real(r8) qp1(plond,plev,pcnst)  ! moist, tracers after vert. diff
      real(r8) pblh(plond)            ! planetary boundary layer height
      real(r8) ustar(plond)           ! surface friction velocity
      real(r8) kvh(plond,plevp)       ! coefficient for heat and tracers
      real(r8) kvm(plond,plevp)       ! coefficient for momentum
      real(r8) tpert(plond)           ! convective temperature excess
      real(r8) qpert(plond)           ! convective humidity excess
      real(r8) cgs(plond,plevp)       ! counter-grad star (cg/flux)
c++bee
      real(r8) shp1(plond,plev)       ! specific humidity (kg/kg)
      real(r8) wpert(plond)
      real(r8) obklen(plond)          ! Obukhov length
c--bee
C
C---------------------------Local workspace-----------------------------
C
      real(r8) cah(plond,plev)        ! -upper diag for heat and constituts
      real(r8) cam(plond,plev)        ! -upper diagonal for momentum
      real(r8) cch(plond,plev)        ! -lower diag for heat and constits
      real(r8) ccm(plond,plev)        ! -lower diagonal for momentum
      real(r8) cgh(plond,plevp)       ! countergradient term for heat
      real(r8) cgq(plond,plevp,pcnst) ! countergrad term for constituent
c++bee
      real(r8) cgsh(plond,plevp)      ! countergrad term for sh
c--bee
      real(r8) dvdz2                  ! (du/dz)**2 + (dv/dz)**2
      real(r8) dz                     ! delta z between midpoints
      real(r8) fstab                  ! stable f(ri)
      real(r8) funst                  ! unstable f(ri)
      real(r8) kvf(plond,plevp)       ! free atmosphere kv at interfaces
      real(r8) rinub                  ! richardson no=(g/theta)(dtheta/dz)/
C                                                 (du/dz**2+dv/dz**2)
      real(r8) sstab                  ! static stability = g/th  * dth/dz
      real(r8) potbar(plond,plevp)    ! pintm1(k)/(.5*(tm1(k)+tm1(k-1))
      real(r8) tmp1(plond)            ! temporary storage
      real(r8) tmp2                   ! temporary storage
      real(r8) rcpair                 ! 1./cpair
      real(r8) ztodtgor               ! ztodt*gravit/rair
      real(r8) gorsq                  ! (gravit/rair)**2
      real(r8) dubot(plond)           ! lowest layer u change from stress
      real(r8) dvbot(plond)           ! lowest layer v change from stress
      real(r8) dtbot(plond)           ! lowest layer t change from heat flx
      real(r8) dqbot(plond,pcnst)     ! lowest layer q change from const flx
c++bee
      real(r8) dshbot(plond)          ! lowest layer sh change from wvflx
c--bee
      real(r8) thx(plond,plev)        ! temperature input + counter gradient
      real(r8) thv(plond,plev)        ! virtual potential temperature
      real(r8) qmx(plond,plev,pcnst)  ! constituents input + counter grad
c++bee
      real(r8) shmx(plond,plev)       ! sh input + counter grad
c--bee
      real(r8) zeh(plond,plev)        ! term in tri-diag. matrix system (T & Q)
      real(r8) zem(plond,plev)        ! term in tri-diag. matrix system (momentum)
      real(r8) termh(plond,plev)      ! 1./(1. + cah(k) + cch(k) - cch(k)*zeh(k-1))
      real(r8) termm(plond,plev)      ! 1./(1. + cam(k) + ccm(k) - ccm(k)*zem(k-1))
      real(r8) kvn                    ! neutral Kv
      real(r8) ksx(plond)             ! effective surface drag factor (x)
      real(r8) ksy(plond)             ! effective surface drag factor (y)
      real(r8) sufac(plond)           ! lowest layer u implicit stress factor
      real(r8) svfac(plond)           ! lowest layer v implicit stress factor

      integer indx(plon)          ! array of indices of potential q<0
      integer ilogic(plon)        ! 1 => adjust vertical profile
      integer nval                ! num of values which meet criteria
      integer ii                  ! longitude index of found points
      integer i                   ! longitude index
      integer k                   ! vertical index
      integer m                   ! constituent index
      integer ktopbl(plond)       ! index of first midpoint inside pbl
      integer ktopblmn            ! min value of ktopbl
C
C-----------------------------------------------------------------------
C
C Convert the surface fluxes to lowest level tendencies.
C Stresses are converted to effective drag coefficients if these are >0
C
      rcpair = 1./cpair
      do i=1,plon
        tmp1(i) = ztodt*gravit*rpdel(i,plev)

c++bee, commented out stuff is only needed for diffusion of momentum
c        ksx(i) = -taux(i) / um1(i,plev)
c        if (ksx(i) .gt. 0.) then
c           sufac(i) = tmp1(i) * ksx(i)
c           dubot(i) = 0.
c        else
c           ksx(i)   = 0.
c           sufac(i) = 0.
c           dubot(i) = tmp1(i) * taux(i)
c        end if

c        ksy(i) = -tauy(i) / vm1(i,plev)
c        if (ksy(i) .gt. 0.) then
c           svfac(i) = tmp1(i) * ksy(i)
c           dvbot(i) = 0.
c        else
c           ksy(i)   = 0.
c           svfac(i) = 0.
c           dvbot(i) = tmp1(i) * tauy(i)
c        end if

        dshbot(i)    = wvflx(i)*tmp1(i)
c--bee
        dtbot(i)     = shflx(i)*tmp1(i)*rcpair
        kvf(i,plevp) = 0.0
      end do
c++bee
      do m=1,pcnst
c--bee
        do i=1,plon
          dqbot(i,m) = cflx(i,m)*tmp1(i)
        end do
      end do
C
C Set the vertical diffusion coefficient above the top diffusion level
C
      do k=1,ntopfl
        do i=1,plon
          kvf(i,k) = 0.0
        end do
      end do
C
C Compute virtual potential temperature for use in static stability 
C calculation.  0.61 is 1. - R(water vapor)/R(dry air).  Use 0.61 instead
C of a computed variable in order to obtain an identical simulation to
C Case 414.
C
c++bee, use MATCH version of virtem
c      call virtem(thm     ,qm1     ,0.61    ,thv     )
      call virtem(thm     ,shm1     ,thv     )
c--bee
C
C Compute the free atmosphere vertical diffusion coefficients
C kvh = kvq = kvm. 
C
      do k=ntopfl,plev-1
        do i=1,plon
C
C Vertical shear squared, min value of (delta v)**2 prevents zero shear.
C
          dvdz2 = (um1(i,k)-um1(i,k+1))**2 + (vm1(i,k)-vm1(i,k+1))**2
          dvdz2 = max(dvdz2,1.e-10_r8)
          dz    = zm(i,k) - zm(i,k+1)
          dvdz2 = dvdz2/(dz**2)
C
C Static stability (use virtual potential temperature)
C
          sstab = gravit*2.0*( thv(i,k) - thv(i,k+1))/
     $                       ((thv(i,k) + thv(i,k+1))*dz)
C
C Richardson number, stable and unstable modifying functions
C
          rinub = sstab/dvdz2
          fstab = 1.0/(1.0 + 10.0*rinub*(1.0 + 8.0*rinub))
          funst = max(1. - 18.*rinub,zero)
C
C Select the appropriate function of the richardson number
C
          if (rinub.lt.0.0) fstab = sqrt(funst)
C
C Neutral diffusion coefficient
C compute mixing length (z), where z is the interface height estimated
C with an 8 km scale height.
C
          kvn = ml2(k)*sqrt(dvdz2)
C
C Full diffusion coefficient (modified by f(ri)),
C
          kvf(i,k+1) = max(zkmin,kvn*fstab)
        end do
      end do
C
C Determine the boundary layer kvh (=kvq), kvm, 
C counter gradient terms (cgh, cgq, cgs)
C boundary layer height (pblh) and 
C the perturbation temperature and moisture (tpert and qpert)
C The free atmosphere kv is returned above the boundary layer top.
C
c++bee
      call pbldif(thm     ,shm1    ,zm      ,um1     ,vm1     ,
c--bee
     $            tm1     ,pmidm1  ,kvf     ,cflx    ,shflx   ,
     $            taux    ,tauy    ,ustar   ,kvm     ,kvh     ,
     $            cgh     ,cgq     ,cgs     ,pblh    ,tpert   ,
     $            qpert   ,ktopbl  ,ktopblmn
c++bee
     $, wvflx, cgsh, wpert, obklen )
c--bee
C
C Add the counter grad terms to potential temp, specific humidity
C and other constituents in the bdry layer. Note, ktopblmn gives the 
c minimum vertical index of the first midpoint within the boundary layer.
C
C first set values above boundary layer
C
      do k=1,ktopblmn-2
        do i=1,plon
          thx(i,k)   = thm(i,k)
c++bee
          shmx(i,k) = shm1(i,k)
        end do
        do m=1,pcnst
c--bee
          do i=1,plon
            qmx(i,k,m) = qm1(i,k,m)
          end do
        end do
      end do
      do k=2,plev
        do i=1,plon
          potbar(i,k) = pintm1(i,k)/(0.5*(tm1(i,k) + tm1(i,k-1)))
        end do
      end do
      do i=1,plon
        potbar(i,plevp) = pintm1(i,plevp)/tm1(i,plev)
      end do
C
C now focus on the boundary layer
C
      ztodtgor = ztodt*gravit/rair
      do k=ktopblmn-1,plev
        do i=1,plon
          tmp1(i) = ztodtgor*rpdel(i,k)
          thx(i,k) = thm(i,k) + tmp1(i)*
     $           (potbar(i,k+1)*kvh(i,k+1)*cgh(i,k+1) -
     $            potbar(i,k  )*kvh(i,k  )*cgh(i,k  ))
c++bee
          shmx(i,k) = shm1(i,k) + tmp1(i) *
     $           (potbar(i,k+1)*kvh(i,k+1)*cgsh(i,k+1) -
     $            potbar(i,k  )*kvh(i,k  )*cgsh(i,k  ))
        end do
        do m=1,pcnst
c--bee
          do i=1,plon
            qmx(i,k,m) = qm1(i,k,m) + tmp1(i)*
     $           (potbar(i,k+1)*kvh(i,k+1)*cgq(i,k+1,m) -
     $            potbar(i,k  )*kvh(i,k  )*cgq(i,k  ,m))
          end do
        end do
      end do
C
C Check for neg q's in each constituent and put the original vertical
C profile back if a neg value is found. A neg value implies that the
C quasi-equilibrium conditions assumed for the countergradient term are
C strongly violated.
C Original code rewritten by Rosinski 7/8/91 to vectorize in longitude.
C
      do m=1,pcnst
        do i=1,plon
          ilogic(i) = 0
        end do
        do k=ktopblmn-1,plev
          do i=1,plon
            if (qmx(i,k,m).lt.qmincg(m)) ilogic(i) = 1
          end do
        end do
C
C Find long indices of those columns for which negatives were found
C
        call wheneq(plon,ilogic,1,1,indx,nval)
C
C Replace those columns with original values
C
        if (nval.gt.0) then
          do k=ktopblmn-1,plev
            do ii=1,nval
              i=indx(ii)
              qmx(i,k,m) = qm1(i,k,m)
            end do
          end do
        end if
      end do
c++bee
c Repeat above for shmx.
      do i=1,plon
         ilogic(i) = 0
      end do
      do k=ktopblmn-1,plev
         do i=1,plon
c           1.e-12 is the value of qmin (=qmincg) used in CCM2.
            if (shmx(i,k).lt.1.e-12) ilogic(i) = 1
         end do
      end do

C     Find long indices of those columns for which negatives were found
      call wheneq(plon,ilogic,1,1,indx,nval)

C     Replace those columns with original values
      if(nval.gt.0) then
         do k=ktopblmn-1,plev
            do ii=1,nval
               i=indx(ii)
               shmx(i,k) = shm1(i,k)
            end do
         end do
      end if
c--bee
C
C Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
C tridiagonal diffusion matrix. the diagonal elements are a combination of 
C ca and cc; they are not explicitly provided to the solver
C
      gorsq = (gravit/rair)**2
      do k=ntopfl,plev-1
        do i=1,plon
          tmp2 = ztodt*gorsq*rpdeli(i,k)*(potbar(i,k+1)**2)
          cah(i,k  ) = kvh(i,k+1)*tmp2*rpdel(i,k  )
          cam(i,k  ) = kvm(i,k+1)*tmp2*rpdel(i,k  )
          cch(i,k+1) = kvh(i,k+1)*tmp2*rpdel(i,k+1)
          ccm(i,k+1) = kvm(i,k+1)*tmp2*rpdel(i,k+1)
        end do
      end do
C
C The last element of the upper diagonal is zero.
C
      do i=1,plon
        cah(i,plev) = 0.
        cam(i,plev) = 0.
      end do
C
C Calculate e(k) for heat & momentum vertical diffusion.  This term is 
C required in solution of tridiagonal matrix defined by implicit diffusion eqn.
C
      do i=1,plon
        termh(i,ntopfl) = 1./(1. + cah(i,ntopfl))
        termm(i,ntopfl) = 1./(1. + cam(i,ntopfl))
        zeh(i,ntopfl) = cah(i,ntopfl)*termh(i,ntopfl)
        zem(i,ntopfl) = cam(i,ntopfl)*termm(i,ntopfl)
      end do
      do k=ntopfl+1,plev-1
        do i=1,plon
          termh(i,k) = 1./
     $           (1. + cah(i,k) + cch(i,k) - cch(i,k)*zeh(i,k-1))
          termm(i,k) = 1./
     $           (1. + cam(i,k) + ccm(i,k) - ccm(i,k)*zem(i,k-1))
          zeh(i,k) = cah(i,k)*termh(i,k)
          zem(i,k) = cam(i,k)*termm(i,k)
        end do
      end do
C
C Diffuse momentum
C
c++bee
c      call mvdiff(um1     ,vm1     ,dubot   ,dvbot   ,
c     $            sufac   ,svfac   ,ccm     ,
c     $            zem     ,termm   ,up1     ,vp1     )
c+
c Determine the difference between the implicit stress and the
c externally specified stress. Apply over boundary layer depth.
c-
c$$$      print *, 'vdiff: lat=, ktopblmn', lat, ktopblmn
c      do i = 1, plon
c         if (ksx(i) .gt. 0.) sufac(i) = (taux(i) + ksx(i)*up1(i,plev))
c     $        * ztodt * gravit / (pintm1(i,plevp) - pintm1(i,ktopbl(i)))
c         if (ksy(i) .gt. 0.) svfac(i) = (tauy(i) + ksy(i)*vp1(i,plev))
c     $        * ztodt * gravit / (pintm1(i,plevp) - pintm1(i,ktopbl(i)))
c$$$         print 1, i, lat
c$$$     $        , taux(i), -ksx(i)*up1(i,plev)
c$$$     $        , sufac(i)*86400./ztodt, ksx(i), um1(i,plev),up1(i,plev)
c$$$    1 format (1x,i4,i4, 1P, 4E20.10, 2f10.4)
c      end do
c      do k = plev, ktopblmn, -1
c         do i = 1, plon
c            if (k .ge. ktopbl(i)) then
c               up1(i,k) = up1(i,k) + sufac(i)
c               vp1(i,k) = vp1(i,k) + svfac(i)
c            end if
c         end do
c      end do
c--bee
C
C Diffuse constituents
C
      call qvdiff(pcnst   ,qmx     ,dqbot   ,cch     ,zeh     ,
     $            termh   ,qp1     )
C
C Identify and correct constituents exceeding user defined bounds
C 
      call qneg3('VDIFF   ',lat,qp1(1,1,1))
c++bee
c Diffuse sh

      call qvdiff(1   ,shmx     ,dshbot   ,cch     ,zeh     ,
     $            termh   ,shp1     )

c correct sh

      call shneg( 'VDIFF:sh', lat, shp1 )
c--bee
C
C Diffuse potential temperature
C
      call qvdiff(1       ,thx     ,dtbot   ,cch     ,zeh     ,
     $            termh   ,thp     )
C
      return
      end subroutine vdiff

c##############################################################################

      subroutine pbldif(th      ,q       ,z       ,u       ,v       ,
     $                  t       ,pmid    ,kvf     ,cflx    ,shflx   ,
     $                  taux    ,tauy    ,ustar   ,kvm     ,kvh     ,
     $                  cgh     ,cgq     ,cgs     ,pblh    ,tpert   ,
     $                  qpert   ,ktopbl  ,ktopblmn
c++bee
     $, wvflx, cgsh, wpert, obklen )
c--bee
C------------------------------------------------------------------------
C 
C Atmospheric boundary layer computation.
C
C Nonlocal scheme that determines eddy diffusivities based on a
C diagnosed boundary layer height and a turbulent velocity scale;
C also, countergradient effects for heat and moisture, and constituents
C are included, along with temperature and humidity perturbations which 
C measure the strength of convective thermals in the lower part of the 
C atmospheric boundary layer.
C
C For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
C Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
C Model. J. Clim., vol. 6., p. 1825--1842.
c
c Updated by Holtslag and Hack to exclude the surface layer from the
c definition of the boundary layer Richardson number. Ri is now defined
c across the outer layer of the pbl (between the top of the surface
c layer and the pbl top) instead of the full pbl (between the surface and
c the pbl top). For simiplicity, the surface layer is assumed to be the
c region below the first model level (otherwise the boundary layer depth 
c determination would require iteration).
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Boville
C Standardized:      J. Rosinski, June 1992
C Reviewed:          B. Boville, P. Rasch, August 1992
C Reviewed:          B. Boville, P. Rasch, April 1996
C
C Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
C >>>>>>>>>  (Use ricr = 0.3 in this formulation)
C
C-----------------------------------------------------------------------
c
c $Id$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <compbl.h>
C-----------------------------------------------------------------------
#include <comvd.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8) th(plond,plev)          ! potential temperature [K]
c++bee, Note: non-water tracers not used in this routine.
      real(r8) q(plond,plev)           ! specific humidity [kg/kg]
c--bee
      real(r8) z(plond,plev)           ! height above surface [m]
      real(r8) u(plond,plev)           ! windspeed x-direction [m/s]
      real(r8) v(plond,plev)           ! windspeed y-direction [m/s]
      real(r8) t(plond,plev)           ! temperature (used for density)
      real(r8) pmid(plond,plev)        ! midpoint pressures
      real(r8) kvf(plond,plevp)        ! free atmospheric eddy diffsvty [m2/s]
      real(r8) cflx(plond,pcnst)       ! surface constituent flux (kg/m2/s)
c++bee
      real(r8) wvflx(plond)            ! water vapor flux (kg/m2/s)
c--bee
      real(r8) shflx(plond)            ! surface heat flux (W/m2)
      real(r8) taux(plond)             ! surface u stress (N)
      real(r8) tauy(plond)             ! surface v stress (N)
C
C Output arguments
C
      real(r8) ustar(plond)            ! surface friction velocity [m/s]
      real(r8) kvm(plond,plevp)        ! eddy diffusivity for momentum [m2/s]
      real(r8) kvh(plond,plevp)        ! eddy diffusivity for heat [m2/s]
      real(r8) cgh(plond,plevp)        ! counter-gradient term for heat [K/m]
      real(r8) cgq(plond,plevp,pcnst)  ! counter-gradient term for constituents
c++bee
      real(r8) cgsh(plond,plevp)       ! counter-gradient term for sh
c--bee
      real(r8) cgs(plond,plevp)        ! counter-gradient star (cg/flux)
      real(r8) pblh(plond)             ! boundary-layer height [m]
      real(r8) tpert(plond)            ! convective temperature excess
      real(r8) qpert(plond)            ! convective humidity excess

      integer ktopbl(plond)        ! index of first midpoint inside pbl
      integer ktopblmn             ! min value of ktopbl
c++pjr
      real(r8) wpert(plond)            ! convective velocity scale
      real(r8) obklen(plond)           ! Obukhov length
c--pjr
C
C---------------------------Local parameters----------------------------
C
      real(r8) tiny                    ! lower bound for wind magnitude
      parameter (tiny=1.e-36)
C
C---------------------------Local workspace-----------------------------
C
      integer i                    ! longitude index
      integer k                    ! level index
      integer m                    ! constituent index

      real(r8) heatv(plond)            ! surface virtual heat flux
      real(r8) thvsrf(plond)           ! sfc (bottom) level virtual temperature
      real(r8) thvref(plond)           ! reference level virtual temperature
      real(r8) tkv                     ! model level potential temperature
      real(r8) therm(plond)            ! thermal virtual temperature excess
      real(r8) phiminv(plond)          ! inverse phi function for momentum
      real(r8) phihinv(plond)          ! inverse phi function for heat 
      real(r8) wm(plond)               ! turbulent velocity scale for momentum
      real(r8) vvk                     ! velocity magnitude squared
      real(r8) zm(plond)               ! current level height
      real(r8) zp(plond)               ! current level height + one level up
      real(r8) khfs(plond)             ! surface kinematic heat flux [mK/s]
      real(r8) kqfs(plond,pcnst)       ! sfc kinematic constituent flux [m/s]
c++bee
      real(r8) kshfs(plond)            ! sfc kinematic moisture flux [m/s]
c--bee
      real(r8) zmzp                    ! level height halfway between zm and zp
c      real(r8) vv2(plond)              ! velocity magnitude squared 
      real(r8) rino(plond,plev)        ! bulk Richardson no. from level to ref lev
      real(r8) tlv(plond)              ! ref. level pot tmp + tmp excess
      real(r8) fak1(plond)             ! k*ustar*pblh
      real(r8) fak2(plond)             ! k*wm*pblh
      real(r8) fak3(plond)             ! fakn*wstr/wm 
      real(r8) pblk(plond)             ! level eddy diffusivity for momentum
      real(r8) pr(plond)               ! Prandtl number for eddy diffusivities
      real(r8) zl(plond)               ! zmzp / Obukhov length
      real(r8) zh(plond)               ! zmzp / pblh
      real(r8) zzh(plond)              ! (1-(zmzp/pblh))**2
      real(r8) wstr(plond)             ! w*, convective velocity scale
      real(r8) rrho(plond)             ! 1./bottom level density (temporary)
c++bee, move to output args
c      real(r8) obklen(plond)           ! Obukhov length
      real(r8) ustr                    ! unbounded ustar
      real(r8) term                    ! intermediate calculation
      real(r8) fac                     ! interpolation factor
      real(r8) pblmin                  ! min pbl height due to mechanical mixing

      logical unstbl(plond)        ! pts w/unstbl pbl (positive virtual ht flx)
      logical stblev(plond)        ! stable pbl with levels within pbl
      logical unslev(plond)        ! unstbl pbl with levels within pbl
      logical unssrf(plond)        ! unstb pbl w/lvls within srf pbl lyr
      logical unsout(plond)        ! unstb pbl w/lvls in outer pbl lyr
      logical check(plond)         ! True=>chk if Richardson no.>critcal
C
C Compute kinematic surface fluxes
C
      do i=1,plon
        rrho(i) = rair*t(i,plev)/pmid(i,plev)
        ustr = sqrt(sqrt(taux(i)**2 + tauy(i)**2)*rrho(i))
        ustar(i) = max(ustr,0.01_r8)
        khfs(i) = shflx(i)*rrho(i)/cpair
c++bee
        kshfs(i) = wvflx(i)*rrho(i)
c--bee
      end do
c++bee
      do m=1,pcnst
c--bee
        do i=1,plon
          kqfs(i,m)= cflx(i,m)*rrho(i)
        end do
      end do
C
C Initialize output arrays with free atmosphere values
C
      do k=1,plevp
        do i=1,plon
          kvm(i,k) = kvf(i,k)
          kvh(i,k) = kvf(i,k)
          cgh(i,k) = 0.0
c++bee
          cgsh(i,k) = 0.0
c--bee
          cgs(i,k)   = 0.0
        end do
      end do
c++bee
      do m=1,pcnst
c--bee
        do k=1,plevp
          do i=1,plon
            cgq(i,k,m) = 0.0
          end do
        end do
      end do
C
C Compute various arrays for use later:
C
      do i=1,plon
c++bee
        thvsrf(i) = th(i,plev)*(1.0 + 0.61*q(i,plev))
        heatv(i)  = khfs(i) + 0.61*th(i,plev)*kshfs(i)
c--bee
        wm(i)     = 0.
        therm(i)  = 0.
        qpert(i)  = 0.
        tpert(i)  = 0.
        fak3(i)   = 0.  
        zh(i)     = 0.  
        obklen(i) = -thvsrf(i)*ustar(i)**3/
     $             (g*vk*(heatv(i) + sign(1.e-10_r8,heatv(i))))
      end do
C
C
C >>>> Define first a new factor fac=100 for use in Richarson number
C      Calculate virtual potential temperature first level
C      and initialize pbl height to z1
C
      fac = 100.
C
      do i=1,plon
c++bee
        thvref(i) = th(i,plev)*(1.0 + 0.61*q(i,plev))
c--bee
        pblh(i) = z(i,plev)
        check(i) = .true.
C
C Initialization of lowest level Ri number 
C (neglected in initial Holtslag implementation)
C
C--      vv2(i) = fac*ustar(i)**2
C--      vv2(i) = max(vv2(i),tiny)
C--      rino(i,plev) = g*(thvsrf(i) - thvref(i))*z(i,plev)/
C--  $                  (thvref(i)*vv2(i))
        rino(i,plev) = 0.0
      end do
C
C PBL height calculation:
C Search for level of pbl. Scan upward until the Richardson number between
C the first level and the current level exceeds the "critical" value.
C
      do k=plev-1,plev-npbl+1,-1
        do i=1,plon
          if (check(i)) then
            vvk = (u(i,k) - u(i,plev))**2 + (v(i,k) - v(i,plev))**2
     $           + fac*ustar(i)**2
            vvk = max(vvk,tiny)
c++bee
            tkv = th(i,k)*(1. + 0.61*q(i,k))
c--bee
            rino(i,k) = g*(tkv - thvref(i))*(z(i,k)-z(i,plev))
     $                 /(thvref(i)*vvk)
            if(rino(i,k).ge.ricr) then
              pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/
     $                 (rino(i,k) - rino(i,k+1))*(z(i,k) - z(i,k+1))
              check(i) = .false.
            end if
          end if
        end do
      end do
C
C Set pbl height to maximum value where computation exceeds number of
C layers allowed
C
      do i=1,plon
        if (check(i)) pblh(i) = z(i,plevp-npbl)
      end do
C
C Improve estimate of pbl height for the unstable points.
C Find unstable points (virtual heat flux is positive):
C
      do i=1,plon
        if (heatv(i) .gt. 0.) then
          unstbl(i) = .true.
          check(i) = .true.
        else
          unstbl(i) = .false.
          check(i) = .false.
        end if   
      end do
C
C For the unstable case, compute velocity scale and the
C convective temperature excess:
C
      do i=1,plon
        if (check(i)) then
          phiminv(i) = (1. - binm*pblh(i)/obklen(i))**onet
          wm(i)= ustar(i)*phiminv(i)
          therm(i) = heatv(i)*fak/wm(i)       
C--         rino(i,plev) = -g*therm(i)*z(i,plev)/(thvref(i)*vv2(i))
          rino(i,plev) = 0.0
          tlv(i) = thvref(i) + therm(i)
        end if
      end do
C
C Improve pblh estimate for unstable conditions using the
C convective temperature excess:
C
      do k=plev-1,plev-npbl+1,-1
        do i=1,plon
          if (check(i)) then
            vvk = (u(i,k) - u(i,plev))**2 + (v(i,k) - v(i,plev))**2
     $           + fac*ustar(i)**2
            vvk = max(vvk,tiny)
c++bee
            tkv = th(i,k)*(1. + 0.61*q(i,k))
c--bee
            rino(i,k) = g*(tkv - tlv(i))*(z(i,k)-z(i,plev))
     $                 /(thvref(i)*vvk)
            if(rino(i,k).ge.ricr) then
              pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/
     $                 (rino(i,k) - rino(i,k+1))*(z(i,k) - z(i,k+1))
              check(i) = .false.
            end if
          end if
        end do
      end do
C
C Points for which pblh exceeds number of pbl layers allowed;
C set to maximum
C
      do i=1,plon
        if (check(i)) pblh(i) = z(i,plevp-npbl)
      end do
C
C PBL height must be greater than some minimum mechanical mixing depth
C Several investigators have proposed minimum mechanical mixing depth
C relationships as a function of the local friction velocity, u*.  We 
C make use of a linear relationship of the form h = c u* where c=700.
C The scaling arguments that give rise to this relationship most often 
C represent the coefficient c as some constant over the local coriolis
C parameter.  Here we make use of the experimental results of Koracin 
C and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
C where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
C latitude value for f so that c = 0.07/f = 700.
C
      do i=1,plon
        pblmin  = 700.0*ustar(i)
        pblh(i) = max(pblh(i),pblmin)
      end do
C
C pblh is now available; do preparation for diffusivity calculation:
C
      do i=1,plon
        pblk(i) = 0.0
        fak1(i) = ustar(i)*pblh(i)*vk
C
C Do additional preparation for unstable cases only, set temperature
C and moisture perturbations depending on stability.
C
        if (unstbl(i)) then
          phiminv(i) = (1. - binm*pblh(i)/obklen(i))**onet
          phihinv(i) = sqrt(1. - binh*pblh(i)/obklen(i))
          wm(i)      = ustar(i)*phiminv(i)
          fak2(i)    = wm(i)*pblh(i)*vk
          wstr(i)    = (heatv(i)*g*pblh(i)/thvref(i))**onet 
          fak3(i)    = fakn*wstr(i)/wm(i)
          tpert(i)   = max(khfs(i)*fak/wm(i),zero)   
c++bee
          qpert(i)   = max(kshfs(i)*fak/wm(i),zero)    
c--bee
c++pjr
c-v1          wpert(i)   = wstr(i)
c         from h&b evaluated at z = 0.1h
          wpert(i)   = ustar(i)*phihinv(i)
c          write (6,*) ' unstable: ', i, wpert(i), wstr(i)
c--pjr
        else
          tpert(i)   = max(khfs(i)*fak/ustar(i),zero) 
c++bee
          qpert(i)   = max(kshfs(i)*fak/ustar(i),zero) 
c--bee
c++pjr
c-v1          wpert(i)   = 0.
c-v2 from h&b wm = wt = ustar/phih evaluated at z = obklen
c          wpert(i)   = ustar(i)/6.
c
c-v3 from h&b wm = wt = ustar/phih evaluated at z = 0
          wpert(i)   = ustar(i)
          
c          write (6,*) ' stable: ', i, wpert(i)
c--pjr
        end if
      end do
C
C Main level loop to compute the diffusivities and 
C counter-gradient terms:
C
      do 1000 k=plev,plev-npbl+2,-1
C
C Find levels within boundary layer:
C
        do i=1,plon
          unslev(i) = .false.
          stblev(i) = .false.
          zm(i) = z(i,k)
          zp(i) = z(i,k-1)
          if (zkmin.eq.0.0 .and. zp(i).gt.pblh(i)) zp(i) = pblh(i)
          if (zm(i) .lt. pblh(i)) then
            ktopbl(i) = k
            zmzp = 0.5*(zm(i) + zp(i))
            zh(i) = zmzp/pblh(i)
            zl(i) = zmzp/obklen(i)
            zzh(i) = 0.
            if (zh(i).le.1.0) zzh(i) = (1. - zh(i))**2
C
C stblev for points zm < plbh and stable and neutral
C unslev for points zm < plbh and unstable
C
            if (unstbl(i)) then
              unslev(i) = .true.
            else
              stblev(i) = .true.
            end if
          end if
        end do
C
C Stable and neutral points; set diffusivities; counter-gradient
C terms zero for stable case:
C
        do i=1,plon
          if (stblev(i)) then
            if (zl(i).le.1.) then
              pblk(i) = fak1(i)*zh(i)*zzh(i)/(1. + betas*zl(i))
            else
              pblk(i) = fak1(i)*zh(i)*zzh(i)/(betas + zl(i))
            end if
            kvm(i,k) = max(pblk(i),kvf(i,k))
            kvh(i,k) = kvm(i,k)
          end if
        end do
C
C unssrf, unstable within surface layer of pbl
C unsout, unstable within outer   layer of pbl
C
        do i=1,plon
          unssrf(i) = .false.
          unsout(i) = .false.
          if (unslev(i)) then
            if (zh(i).lt.sffrac) then
              unssrf(i) = .true.
            else
              unsout(i) = .true.
            end if
          end if
        end do
C
C Unstable for surface layer; counter-gradient terms zero
C
        do i=1,plon
          if (unssrf(i)) then
            term = (1. - betam*zl(i))**onet
            pblk(i) = fak1(i)*zh(i)*zzh(i)*term
            pr(i) = term/sqrt(1. - betah*zl(i))
          end if
        end do
C
C Unstable for outer layer; counter-gradient terms non-zero:
C
        do i=1,plon
          if (unsout(i)) then
            pblk(i) = fak2(i)*zh(i)*zzh(i)
            cgs(i,k) = fak3(i)/(pblh(i)*wm(i))
            cgh(i,k) = khfs(i)*cgs(i,k)
            pr(i) = phiminv(i)/phihinv(i) + ccon*fak3(i)/fak
c++bee
            cgsh(i,k) = kshfs(i)*cgs(i,k)
c--bee
          end if
        end do
c++bee
        do m=1,pcnst
c--bee
          do i=1,plon
            if (unsout(i)) cgq(i,k,m) = kqfs(i,m)*cgs(i,k)
          end do
        end do
C
C For all unstable layers, set diffusivities
C
        do i=1,plon
          if (unslev(i)) then
            kvm(i,k) = max(pblk(i),kvf(i,k))
            kvh(i,k) = max(pblk(i)/pr(i),kvf(i,k))
          end if
        end do
 1000 continue           ! end of level loop
c+
c Check whether last allowed midpoint is within pbl, determine ktopblmn
c-
      ktopblmn = plev
      k = plev-npbl+1
      do i = 1, plon
         if (z(i,k) .lt. pblh(i)) ktopbl(i) = k
         ktopblmn = min(ktopblmn, ktopbl(i))
      end do
C
      return
      end subroutine pbldif

c##############################################################################

      subroutine qvdiff(ncnst   ,qm1     ,qflx    ,cc      ,ze      ,
     $                  term    ,qp1     )
C-----------------------------------------------------------------------
C
C Solve vertical diffusion eqtn for constituent with explicit srfc flux.
C Procedure for solution of the implicit equation follows Richtmyer and 
C Morton (1967,pp 198-199).
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          P. Rasch, B. Boville, August 1992
C Reviewed:          P. Rasch, April 1996
C Reviewed:          B. Boville, May 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comvd.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer ncnst               ! number of constituents being diffused

      real(r8) qm1(plond,plev,ncnst)  ! initial constituent
      real(r8) qflx(plond,ncnst)      ! sfc q flux into lowest model level
      real(r8) cc(plond,plev)         ! -lower diag coeff.of tri-diag matrix
      real(r8) ze(plond,plev)         ! term in tri-diag. matrix system
      real(r8) term(plond,plev)       ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
C
C Output arguments
C
      real(r8) qp1(plond,plev,ncnst)  ! final constituent
C
C---------------------------Local workspace-----------------------------
C
      real(r8) zfq(plond,plev,pcnst)  ! terms appear in soln of tri-diag sys
      real(r8) tmp1d(plond)           ! temporary workspace (1d array)

      integer i,k                 ! longitude, vertical indices
      integer m                   ! constituent index
C
C-----------------------------------------------------------------------
C
C Calculate fq(k).  Terms fq(k) and e(k) are required in solution of 
C tridiagonal matrix defined by implicit diffusion eqn.
C Note that only levels ntopfl through plev need be solved for.
C No vertical diffusion is applied above this level
C
      do m=1,ncnst
        do i=1,plon
          zfq(i,ntopfl,m) = qm1(i,ntopfl,m)*term(i,ntopfl)
        end do
        do k=ntopfl+1,plev-1
          do i=1,plon
            zfq(i,k,m) = (qm1(i,k,m) + cc(i,k)*zfq(i,k-1,m))*term(i,k)
          end do
        end do
      end do
C
C Bottom level: (includes  surface fluxes)
C
      do i=1,plon
        tmp1d(i) = 1./(1. + cc(i,plev) - cc(i,plev)*ze(i,plev-1))
        ze(i,plev) = 0.
      end do
      do m=1,ncnst
        do i=1,plon
          zfq(i,plev,m) = (qm1(i,plev,m) + qflx(i,m) +
     $                     cc(i,plev)*zfq(i,plev-1,m))*tmp1d(i)
        end do
      end do
C
C Perform back substitution
C
      do m=1,ncnst
        do i=1,plon
          qp1(i,plev,m) = zfq(i,plev,m)
        end do
        do k=plev-1,ntopfl,-1
          do i=1,plon
            qp1(i,k,m) = zfq(i,k,m) + ze(i,k)*qp1(i,k+1,m)
          end do
        end do
      end do
      return
      end subroutine qvdiff

c $Id$

c Interface routine for CCM3 vertical diffusion and boundary layer
c parameterizations.
c Also includes routines for "archived" calculations.

      subroutine vdiffdr( divdiff, arvdiff
     $,                   lat, um1, vm1, tadv, t1, as1
     $,                   pmid, pint, rpdel, rpdeli, dtime
     $,                   thm, zm, taux, tauy, hflx
     $,                   sflx, shm, shflx, cgs, kvh
     $,                   obuf, thp, as2, shp
     $,                   pblh, tpert, qpert, wpert3d
     $,                   ustar, obklen
     $)

c Vertical diffusion calculations on a latitude slice of data.

      use histout, only: outfld

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments:

      logical
     $  divdiff     ! t => diagnosed vertical diffusion
     $, arvdiff     ! t => archived vertical diffusion (from CCM2 run)
      integer
     $  lat         ! latitude index
      real(r8)
     $  um1(plond,plev)   ! u wind
     $, vm1(plond,plev)   ! v wind
     $, tadv(plond,plev)  ! temperature corresponding to thm (K)
     $, t1(plond,plev)    ! temperature from input data interpolated to
c                           end of current timestep, (K)
     $, as1(plond,plev,pcnst) ! advected species
     $, pmid(plond,plev)  ! pressure at layer midpoints
     $, pint(plond,plevp) ! pressure at layer interfaces
     $, rpdel(plond,plev) ! 1/pdel
     $, rpdeli(plond,plev)! 1./(pmid(k+1)-pmid(k))
     $, dtime             ! timestep size (s)
     $, thm(plond,plev)   ! potential temperature
     $, zm(plond,plev)    ! potential height above surface at midpoints
     $, taux(plond)       ! x surface stress (n)
     $, tauy(plond)       ! y surface stress (n)
     $, hflx(plond)       ! surface sensible heat flux (w/m2)
     $, sflx(plond,pcnst) ! surface fluxes for advected species
     $, shm(plond,plev)   ! specific humidity (kg/kg)
     $, shflx(plond)      ! surface water vapor flux (kg/m2/s)

c Input/Output arguments:
      real(r8)
c cgs and kvh are input for archived vdiff and output for diagnosed vdiff
     $  cgs(plond,plevp)       ! counter-gradient coefficient
     $, kvh(plond,plevp)       ! vertical diffusion coefficient

c Output arguments:
      real(r8)
     $  thp(plond,plev)       ! potential temperature
     $, as2(plond,plev,pcnst) ! advected species
     $, shp(plond,plev)       ! specific humidity (kg/kg)
     $, pblh(plond)           ! planetary boundary layer height
     $, tpert(plond)          ! convective temperature excess
     $, qpert(plond)          ! convective humidity excess
     $, wpert3d(plond,plev)
     $, ustar(plond)          ! surface friction velocity
     $, obklen(plond)         ! Obukhov length
     $, obuf(*)               ! output buffer

c Local variables:
c     Diagnosed Vdiff
      integer i, k
      real(r8)
     $  up1(plond,plev)         ! u-wind after vertical diffusion
     $, vp1(plond,plev)         ! v-wind after vertical diffusion
     $, kvm(plond,plevp)        ! momentum k
     $, wpert(plond)            ! vertical velocity perturbation
     $, absvel
c-----------------------------------------------------------------------

      if ( divdiff ) then

         call vdiff( lat, um1, vm1, tadv, as1
     $,              pmid, pint, rpdel, rpdeli, dtime
     $,              thm, zm, taux, tauy, hflx
     $,              sflx, up1, vp1, thp, as2
     $,              pblh, ustar, kvh, kvm, tpert
     $,              qpert, cgs
     $,              shm, shp, shflx, wpert, obklen )

c        construct a 3d vertical velocity perturbation for the convection
c        this part is associated with the turbulent boundary layer and 
c        flow over topography
         do k = 1,plev
            do i = 1,plon
               absvel = abs(um1(i,plev)) + abs(vm1(i,plev))

               wpert3d(i,k) = 
c
c              assume the pbl pert penetrates 2 pbl heights
     $            wpert(i)*max(zero,1.-zm(i,k)/(2*pblh(i)))
c
c              assume the flow over oro penetrates max(2 pbl, sgh)
c     $          + absvel*sgh(i)*2.e-5*
c     $            max(zero,1.-zm(i,k)/max(2.*pblh(i),sgh(i)))
            end do
         end do

         call outfld( 'KVH', kvh(1,2), plond, lat, obuf )
         call outfld( 'CGS', cgs(1,2), plond, lat, obuf )
         call outfld( 'PBLH', pblh,    plond, lat, obuf )
         call outfld( 'TPERT', tpert,  plond, lat, obuf )
         call outfld( 'QPERT', qpert,  plond, lat, obuf )

      else if ( arvdiff ) then

c        Vertical diffusion using archived values of cgs and kvh.

         call vdiffar( lat, t1, as1
     $,                pmid, pint, rpdel, rpdeli, dtime
     $,                sflx, as2
     $,                kvh, cgs )

      end if

      return
      end subroutine vdiffdr

c##############################################################################

      subroutine vdiffar(lat     ,tm1     ,qm1     ,
     $                 pmidm1  ,pintm1  ,rpdel   ,rpdeli  ,ztodt   ,
     $                 cflx    ,qp1     ,
     $                 kvh     ,cgs )
C-----------------------------------------------------------------------
C
C Driver routine to compute vertical diffusion of trace constituents
C using archived coefficients for cgs and kvh.  This is a gutted
C version of vdiff.
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comvd.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat                 ! latitude index
      real(r8)
     $     tm1(plond,plev),       ! temperature input
     $     qm1(plond,plev,pcnst), ! moisture and trace constituent input
     $     pmidm1(plond,plev),    ! midpoint pressures
     $     pintm1(plond,plevp),   ! interface pressures
     $     rpdel(plond,plev),     ! 1./pdel  (thickness bet interfaces)
     $     rpdeli(plond,plev),    ! 1./pdeli (thickness bet midpoints)
     $     ztodt,                 ! 2 delta-t
     $     cflx(plond,pcnst),     ! surface constituent flux (kg/m2/s)
     $     kvh(plond,plevp),      ! coefficient for heat and tracers
     $     cgs(plond,plevp)       ! counter-grad star (cg/flux)
C
C Output arguments
C
      real(r8)
     $     qp1(plond,plev,pcnst)  ! moist, tracers after vert. diff
C
C---------------------------Local workspace-----------------------------
C
      real(r8) cah(plond,plev),       ! -upper diag for heat and constituts
     $     cch(plond,plev),       ! -lower diag for heat and constits
     $     cgq(plond,plevp,pcnst),! countergrad term for constituent
     $     potbar(plond,plevp),   ! pintm1(k)/(.5*(tm1(k)+tm1(k-1))
     $     tmp1(plond),           ! temporary storage
     $     tmp2,                  ! temporary storage
     $     ztodtgor,              ! ztodt*gravit/rair
     $     gorsq,                 ! (gravit/rair)**2
     $     dqbot(plond,pcnst),    ! lowest layer q change from const flx
     $     qmx(plond,plev,pcnst), ! constituents input + counter grad
     $     zeh(plond,plev),       ! term in tri-diag. matrix system (T & Q)
     $     termh(plond,plev)      ! 1./(1. + cah(k) + cch(k) - cch(k)*zeh(k-1))
      integer indx(plon),         ! array of indices of potential q<0
     $        ilogic(plon),       ! 1 => adjust vertical profile
     $        nval,               ! num of values which meet criteria
     $        ii                  ! longitude index of found points
      integer i,                  ! longitude index
     $        k,                  ! vertical index
     $        m                   ! constituent index
C
C------------------------------Externals--------------------------------
C
      external 
c    $         qvdiff,            ! constituent v-diff
     $         wheneq,            ! find indices where criteria met
     $         qneg3
c     $         pbldifar
C
C-----------------------------------------------------------------------
C
C Convert the surface fluxes to lowest level tendencies
C
      do i=1,plon
         tmp1(i)      = ztodt*gravit*rpdel(i,plev)
      end do
      do m=1,pcnst
         do i=1,plon
            dqbot(i,m) = cflx(i,m)*tmp1(i)
         end do
      end do

C counter gradient terms:
C
      call pbldifar(
     $            tm1     ,pmidm1  ,cflx    ,cgs     ,cgq )
C
C Add the counter grad terms to potential temp, specific humidity
C and other constituents in the bdry layer. Note, npbl gives the max
C num of levels which are permitted to be within the boundary layer.
C
C first set values above boundary layer
C
      do k=1,plev-npbl
         do m=1,pcnst
            do i=1,plon
               qmx(i,k,m) = qm1(i,k,m)
            end do
         end do
      end do
C
      do k=2,plev
         do i=1,plon
            potbar(i,k) = pintm1(i,k)/(0.5*(tm1(i,k) + tm1(i,k-1)))
         end do
      end do
      do i=1,plon
         potbar(i,plevp) = pintm1(i,plevp)/tm1(i,plev)
      end do
C
C now focus on the boundary layer
C
      ztodtgor = ztodt*gravit/rair
      do k=plev-npbl+1,plev
         do i=1,plon
            tmp1(i) = ztodtgor*rpdel(i,k)
         end do
         do m=1,pcnst
            do i=1,plon
               qmx(i,k,m) = qm1(i,k,m) + tmp1(i)*
     $           (potbar(i,k+1)*kvh(i,k+1)*cgq(i,k+1,m) -
     $            potbar(i,k  )*kvh(i,k  )*cgq(i,k  ,m))
            end do
         end do
      end do
C
C Check for neg q's in each constituent and put the original vertical
C profile back if a neg value is found. A neg value implies that the
C quasi-equilibrium conditions assumed for the countergradient term are
C strongly violated.
C Original code rewritten by Rosinski 7/8/91 to vectorize in longitude.
C
      do m=1,pcnst
         do i=1,plon
            ilogic(i) = 0
         end do
         do k=plev-npbl+1,plev
            do i=1,plon
               if (qmx(i,k,m).lt.qmincg(m)) ilogic(i) = 1
            end do
         end do
C
C Find long indices of those columns for which negatives were found
C
         call wheneq(plon,ilogic,1,1,indx,nval)
C
C Replace those columns with original values
C
         if(nval.gt.0) then
            do k=plev-npbl+1,plev
               do ii=1,nval
                  i=indx(ii)
                  qmx(i,k,m) = qm1(i,k,m)
               end do
            end do
         end if
      end do
C
C Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs
C of the tridiagonal diffusion matrix. the diagonal elements are a
C combination of ca and cc; they are not explicitly provided to the 
C solver
C
      gorsq = (gravit/rair)**2
      do k=ntopfl,plev-1
         do i=1,plon
            tmp2 = ztodt*gorsq*rpdeli(i,k)*(potbar(i,k+1)**2)
            cah(i,k  ) = kvh(i,k+1)*tmp2*rpdel(i,k  )
            cch(i,k+1) = kvh(i,k+1)*tmp2*rpdel(i,k+1)
         end do
      end do
C
C The last element of the upper diagonal is zero.
C
      do i=1,plon
         cah(i,plev) = 0.
      end do
C
C Calculate e(k) for heat vertical diffusion.  This term is 
C required in solution of tridiagonal matrix defined by implicit diffusion eqn.
C
      do i=1,plon
         termh(i,ntopfl) = 1./(1. + cah(i,ntopfl))
         zeh(i,ntopfl) = cah(i,ntopfl)*termh(i,ntopfl)
      end do
      do k=ntopfl+1,plev-1
         do i=1,plon
            termh(i,k) = 1./
     $           (1. + cah(i,k) + cch(i,k) - cch(i,k)*zeh(i,k-1))
            zeh(i,k) = cah(i,k)*termh(i,k)
         end do
      end do
C
C Diffuse constituents
C
      call qvdiff(pcnst   ,qmx     ,dqbot   ,cch     ,zeh     ,
     $            termh   ,qp1     )
C
C Identify and correct constituents exceeding user defined bounds
C 
      call qneg3('VDIFF   ',lat,qp1(1,1,1))

      return
      end subroutine vdiffar

c##############################################################################

      subroutine pbldifar(
     $                  t       ,pmid    ,cflx    ,cgs     ,cgq )

c Gutted version of pbldif which only calculates cgq given cgs.
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comvd.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(r8)
     $     t(plond,plev),          ! temperature (used for density)
     $     pmid(plond,plev),       ! midpoint pressures
     $     cflx(plond,pcnst),       ! surface constituent flux (kg/m2/s)
     $     cgs(plond,plevp)        ! counter-gradient star (cg/flux)
C
C Output arguments
C
      real(r8)
     $     cgq(plond,plevp,pcnst)  ! counter-gradient term for constituents
C
C---------------------------Local workspace-----------------------------
C
      integer i,                   ! longitude index
     $        k,                   ! level index
     $        m                    ! constituent index
      real(r8)
     $     rrho(plond),            ! 1./bottom level density
     $     kqfs(plond,pcnst)       ! sfc kinematic constituent flux [m/s]
C------------------------------------------------------------------------
C
C Compute kinematic surface fluxes
C
      do i=1,plon
         rrho(i) = rair*t(i,plev)/pmid(i,plev)
      end do
      do m=1,pcnst
         do i=1,plon
            kqfs(i,m)= cflx(i,m)*rrho(i)
         end do
      end do
C
C Initialize output arrays with free atmosphere values
C
      do m=1,pcnst
         do k=1,plevp
            do i=1,plon
               cgq(i,k,m) = 0.0
            end do
         end do
      end do
C
C Compute the counter-gradient terms:
C
      do k=plev,plev-npbl+2,-1
         do m=1,pcnst
            do i=1,plon
               cgq(i,k,m) = kqfs(i,m)*cgs(i,k)
            end do
         end do
      end do
C
      return
      end subroutine pbldifar

c##############################################################################

      end module vdiff_mod
