      subroutine optdep( name, mr, pdel, q, rh, od )

c Optical depth calculation.
      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input

      character*(*) name   ! keyword for aerosol
      real(r8)
     $  mr(plond,plev)     ! aerosol mixing ratio (dry)
     $, pdel(plond,plev)   ! pressure difference across each layer (Pa)
     $, q(plond,plev)      ! specific humidity (kg/kg)
     $, rh(plond,plev)     ! relative humidity (fraction)
c Output
      real(r8)
     $  od(plond)          ! optical depth

c Local
      integer i, k
      real(r8)
     $  aer_ext
     $, spcext(plond,plev)
     $, x(plond,plev)

c External
      real(r8) ext_coeff
c-----------------------------------------------------------------------

c set rh < 0.99 to avoid infinite optical depths for hydroscopic growth species


      do k = 1, plev
         do i = 1, plon
            if ( rh(i,k).gt.0.99 ) rh(i,k) = 0.99
         end do
      end do

      if ( trim(name) .eq. 'SO4' ) then

         call so4_extinction( rh, spcext )
         do k = 1, plev
            do i = 1, plon
               x(i,k) = mr(i,k) * spcext(i,k)
            end do
         end do

!
!+WDC 8/30/2K 
! Hydrophilic organic carbon and sea-salt both have hydroscopic growth factors
!  Sea-salt moved here out of getsslt for simplicity
!
      else if ( trim(name).eq.'CPHIL' .or. trim(name).eq.'BBOCD' .or.
     $          trim(name).eq.'FFOCD' .or. trim(name).eq.'NOCD'  .or.
     $          trim(name).eq.'OCPHI' ) then

         call ocphi_extinction( rh, spcext)
         do k = 1, plev
            do i = 1, plon
               x(i,k) = mr(i,k) * spcext(i,k)
            end do
         end do

      else if ( trim(name) .eq. 'SSLT' ) then

         call sslt_extinction( rh, spcext )
         do k = 1, plev
            do i = 1, plon
               x(i,k) = mr(i,k) * spcext(i,k)
            end do
         end do

!
!-WDC
!
      else 

         aer_ext = ext_coeff( name )
         do k = 1, plev
            do i = 1, plon
               x(i,k) = mr(i,k) * aer_ext
            end do
         end do


      end if

      call cidry( x, pdel, q, od )


      return
      end

c##############################################################################

!
!+WDC 8/30/2K 
! Hydrophilic organic carbon has a hydroscopic growth factor
!
! mod by CTB 
      subroutine ocphi_extinction( rh, spcext )
C
c Compute hydrophilic organic carbon
C OPAC optics for water-soluble following Cooke et. al.(1999)
C 8/31/2000 Data from load_opac_aeronet and fit_opac_aeronet
C source = waso** optics 
C
C
      use precision
      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input
      real(r8) rh(plond,plev)

c Output
      real(r8) spcext(plond,plev)

c Local variables
      integer i, k, krh
      real(r8) aext(5)
      real(r8) ext_dry
      data aext/ 2.12077E+01, -8.52928E-02, -1.03633E+00, -5.26563E+02,  2.47152E+01 /
     $, ext_dry/ 2.3581343e3  / ! m^2/kg

      real(r8) aext_532(5)
      data aext_532/ 15.472636, -0.083201175, -1.0369290, -284.76629, 18.293469 /
      real(r8) ext_dry_532
      data ext_dry_532/ 3.06216e3 /

      real(r8) aext_550(5)
      data aext_550/ 15.548171, -0.083667377, -1.0368196, -286.36906, 18.306359 /
      real(r8) ext_dry_550
      data ext_dry_550/ 2.91374e3 /

      real(r8) aext_558(5)
      data aext_558/ 15.607089, -0.083882124, -1.0367792, -287.98471, 18.340030 /
      real(r8) ext_dry_558
      data ext_dry_558/ 2.85100e3 /

      real(r8) aext_630(5)
      data aext_630/ 15.803516, -0.085665864, -1.0364468, -290.29989, 18.254570 /
      real(r8) ext_dry_630
      data ext_dry_630/ 2.35950e3 /

c-----------------------------------------------------------------------

      spcext = 0.   ! initialize plon:plond 

      do k = 1, plev
         do i = 1, plon
#ifdef LITE
            spcext(i,k) = exp( aext_532(1) + aext_532(2)/(rh(i,k)+aext_532(3))
     $                       + aext_532(4)/(rh(i,k)+aext_532(5)) ) * ext_dry_532
#endif
#ifdef MODIS
            spcext(i,k) = exp( aext_550(1) + aext_550(2)/(rh(i,k)+aext_550(3))
     $                       + aext_550(4)/(rh(i,k)+aext_550(5)) ) * ext_dry_550
#endif
#ifdef MISR
            spcext(i,k) = exp( aext_558(1) + aext_558(2)/(rh(i,k)+aext_558(3))
     $                       + aext_558(4)/(rh(i,k)+aext_558(5)) ) * ext_dry_558
#endif
#ifdef AVHRR
            spcext(i,k) = exp( aext_630(1) + aext_630(2)/(rh(i,k)+aext_630(3))
     $                       + aext_630(4)/(rh(i,k)+aext_630(5)) ) * ext_dry_630
#endif
#ifdef OLDOPTICS
            spcext(i,k) = exp( aext(1) + aext(2)/(rh(i,k)+aext(3))
     $                       + aext(4)/(rh(i,k)+aext(5)) ) * ext_dry
#endif
         end do
      end do

      return
      end
!
!-WDC
!

c##############################################################################

      subroutine so4_extinction( rh, spcext)

c Compute SO4 specific extinction following Wong et al, 1998, eq. 11
c    Aerosol model: H2SO4 with LN dry PDF r_g = 0.05 microns, 
c                                         sigma_g = 2.03 microns
c    See Table 5, second line (r_g = 0.05, ln sigma_g = 0.7)

      use precision
      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input
      real(r8) rh(plond,plev)

c Output
      real(r8) spcext(plond,plev)

c Local variables
      integer i, k, krh
      real(r8) aext(5)
      real(r8) ext_dry
      data aext/ 1.124e1, -3.04e-1, -1.088, -1.776e2, 1.537e1 /
     $, ext_dry/ 3.203e3 / ! m^2/kg

      real(r8) aext_532(5)
      data aext_532/ 16.978767, -0.10997786, -1.0362092, -321.35314, 18.702310 /
      real(r8) ext_dry_532
      data ext_dry_532/ 5.14471e3 /

      real(r8) aext_550(5)
      data aext_550/ 17.134010, -0.11127778, -1.0360842, -324.17691, 18.692229 /
      real(r8) ext_dry_550
      data ext_dry_550/ 4.91311e3 /

      real(r8) aext_558(5)
      data aext_558/ 17.255250, -0.11174821, -1.0359958, -326.99307, 18.720622 /
      real(r8) ext_dry_558
      data ext_dry_558/ 4.80441e3 /

      real(r8) aext_630(5)
      data aext_630/ 17.612570, -0.11532377, -1.0354502, -331.02030, 18.554828 /
      real(r8) ext_dry_630
      data ext_dry_630/ 4.10697e3 /

c-----------------------------------------------------------------------


      spcext = 0.   ! initialize plon:plond 

      do k = 1, plev
         do i = 1, plon
#ifdef LITE
            spcext(i,k) = exp( aext_532(1) + aext_532(2)/(rh(i,k)+aext_532(3))
     $                       + aext_532(4)/(rh(i,k)+aext_532(5)) ) * ext_dry_532
#endif
#ifdef MODIS
            spcext(i,k) = exp( aext_550(1) + aext_550(2)/(rh(i,k)+aext_550(3))
     $                       + aext_550(4)/(rh(i,k)+aext_550(5)) ) * ext_dry_550
#endif
#ifdef MISR
            spcext(i,k) = exp( aext_558(1) + aext_558(2)/(rh(i,k)+aext_558(3))
     $                       + aext_558(4)/(rh(i,k)+aext_558(5)) ) * ext_dry_558
#endif
#ifdef AVHRR
            spcext(i,k) = exp( aext_630(1) + aext_630(2)/(rh(i,k)+aext_630(3))
     $                       + aext_630(4)/(rh(i,k)+aext_630(5)) ) * ext_dry_630
#endif
#ifdef OLDOPTICS
            spcext(i,k) = exp( aext(1) + aext(2)/(rh(i,k)+aext(3))
     $                       + aext(4)/(rh(i,k)+aext(5)) ) * ext_dry
#endif
         end do
      end do

      return
      end

c##############################################################################

!
!+WDC 8/30/2K 
! Sea-salt moved here out of getsslt for simplicity
!
      subroutine sslt_extinction( rh, spcext)
C  
C  OPAC optics for Coarse & accumulation mode sea-salt (6:1 mass ratio)
C    Data from load_opac_aeronet.pro and fit_opac_sslt_aeronet.pro for 630 nm
C    (8/31/2000)
C
      use precision
      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input

      real(r8) rh(plond,plev)

c Output
      real(r8) spcext(plond,plev)

c Local variables
      integer i, k, krh
      real(r8) aext(5)
      real(r8) ext_dry
      data aext/  1.50759E+01, -5.09172E-02, -1.02066E+00, -1.53318E+02,  1.01958E+01 /, 
     $      ext_dry  /   0.40192e3  / ! m^2/kg

      real(r8) aext_532(5)
      data aext_532/ 15.198924, -0.051143153, -1.0205953, -157.39744, 10.381994 /
      real(r8) ext_dry_532
      data ext_dry_532/ 0.399002e3 /

      real(r8) aext_550(5)
      data aext_550/ 15.189486, -0.051078460, -1.0205965, -157.00176, 10.362088 /
      real(r8) ext_dry_550
      data ext_dry_550/ 0.399965e3 /

      real(r8) aext_558(5)
      data aext_558/ 15.179063, -0.051048787, -1.0205986, -156.60625, 10.343104 /
      real(r8) ext_dry_558
      data ext_dry_558/ 0.400114e3 /

      real(r8) aext_630(5)
      data aext_630/ 15.208123, -0.050906078, -1.0206593, -156.17225, 10.295080 /
      real(r8) ext_dry_630
      data ext_dry_630/ 0.401900e3 /

c-----------------------------------------------------------------------

      spcext = 0.   ! initialize plon:plond 

      do k = 1, plev
         do i = 1, plon
#ifdef LITE
            spcext(i,k) = exp( aext_532(1) + aext_532(2)/(rh(i,k)+aext_532(3))
     $                       + aext_532(4)/(rh(i,k)+aext_532(5)) ) * ext_dry_532
#endif
#ifdef MODIS
            spcext(i,k) = exp( aext_550(1) + aext_550(2)/(rh(i,k)+aext_550(3))
     $                       + aext_550(4)/(rh(i,k)+aext_550(5)) ) * ext_dry_550
#endif
#ifdef MISR
            spcext(i,k) = exp( aext_558(1) + aext_558(2)/(rh(i,k)+aext_558(3))
     $                       + aext_558(4)/(rh(i,k)+aext_558(5)) ) * ext_dry_558
#endif
#ifdef AVHRR
            spcext(i,k) = exp( aext_630(1) + aext_630(2)/(rh(i,k)+aext_630(3))
     $                       + aext_630(4)/(rh(i,k)+aext_630(5)) ) * ext_dry_630
#endif
#ifdef OLDOPTICS
            spcext(i,k) = exp( aext(1) + aext(2)/(rh(i,k)+aext(3))
     $                       + aext(4)/(rh(i,k)+aext(5)) ) * ext_dry
#endif
         end do
      end do

      return
      end

!
!-WDC
!
c##############################################################################

      function ext_coeff( name )

c Return requested extinction coefficient.

      use precision

      implicit none

      real(r8) ext_coeff
c-----------------------------------------------------------------------
#include <dstgrd.h>
#include <dstodx.h>
c-----------------------------------------------------------------------

c Input argument:
      character*(*) name   ! keyword for aerosol

c Local variables:
c Extinction for CZ's dust type 1 (0.01 - 1 microns), type 2 (1 - 10 microns),
c type 3 (10 - 20 microns), type 4 (20 - 50 microns) from dstbdrad.f
c-drb      real(r8)
c     $  dst1, dst2, dst3, dst4        ! m^2/kg
c      data dst1/ 1302.82d0 /, dst2/ 124.263d0 /, dst3/ 39.8701d0 /
c     $     dst4/ 17.4849d0 /
c      save dst1, dst2, dst3, dst4

!
!+WDC
! Extinction for hydrophobic and hydrophilic carbon aerosol.
!   8/31/2000 Data taken from load_opac_aeronet and fit_opac_aeronet
!   ocpho = waso00 and bc = data from soot
!
      real(r8)
     $  bc, ocpho, bc_532, bc_550, bc_558, bc_630, ocpho_532, ocpho_550, ocpho_558, ocpho_630
      data bc/ 7.6339169e3 /, ocpho/ 2.3581343e3 /,
     $     bc_532/ 9.7282e3 /, bc_550/ 9.24736e3 /, bc_558/ 9.05065e3/, bc_630/ 7.63767e3 /,
     $     ocpho_532/ 3.06216e3 /, ocpho_550/ 2.91374e3 /, ocpho_558/ 2.851e3 /, ocpho_630/ 2.35950e3 /
      save bc, ocpho, bc_532, bc_550, bc_558, bc_630, ocpho_532, ocpho_550, ocpho_558, ocpho_630
!
!-WDC
!
c-----------------------------------------------------------------------
c-drb
c      if ( trim(name) .eq. 'DSTQ01' ) then
c         ext_coeff = dst1
c      else if ( trim(name) .eq. 'DSTQ02' ) then
c         ext_coeff = dst2
c      else if ( trim(name) .eq. 'DSTQ03' ) then
c         ext_coeff = dst3
c      else if ( trim(name) .eq. 'DSTQ04' ) then
c         ext_coeff = dst4


      if ( trim(name) .eq. 'DSTQ01' ) then
         ext_coeff = ext_cff_mss_dst_dgn(1)
      else if ( trim(name) .eq. 'DSTQ02' ) then
         ext_coeff = ext_cff_mss_dst_dgn(2)
      else if ( trim(name) .eq. 'DSTQ03' ) then
         ext_coeff = ext_cff_mss_dst_dgn(3)
      else if ( trim(name) .eq. 'DSTQ04' ) then
         ext_coeff = ext_cff_mss_dst_dgn(4)


!
!+WDC 8/30/2K
! all types of black carbon are treated as soot, and
! hydrophobic organic carbon is treated at unwettable "mixed" aerosol
!
      else if ( trim(name).eq.'BBBC'  .or. trim(name).eq.'BBBCD' .or.
     $          trim(name).eq.'FFBC'  .or. trim(name).eq.'FFBCD' .or. 
     $        trim(name).eq.'BCPHO' .or. trim(name).eq.'BCPHI') then
#ifdef LITE
         ext_coeff = bc_532
#endif
#ifdef MODIS
         ext_coeff = bc_550
#endif
#ifdef MISR
         ext_coeff = bc_558
#endif
#ifdef AVHRR
         ext_coeff = bc_630
#endif
#ifdef OLDOPTICS
         ext_coeff = bc
#endif
      else if ( trim(name).eq.'CPHOB' .or. trim(name).eq.'BBOC'  .or.
     $          trim(name).eq.'FFOC'  .or. trim(name).eq.'NOC'   .or.
     $          trim(name).eq.'OCPHO' ) then
#ifdef LITE
         ext_coeff = ocpho_532
#endif
#ifdef MODIS
         ext_coeff = ocpho_550
#endif
#ifdef MISR
         ext_coeff = ocpho_558
#endif
#ifdef AVHRR
         ext_coeff = ocpho_630
#endif
#ifdef OLDOPTICS
         ext_coeff = ocpho
#endif
!
!-WDC
!
      else
         write(*,*)'ext_coeff: aerosol keyword ',name,' not recognized'
         stop
      end if

      return
      end
