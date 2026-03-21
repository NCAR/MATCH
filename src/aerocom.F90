
   module aerocom
      use netcdf
      use pmgrid
      use precision
      implicit none
      save

      integer, parameter :: ndays = 366
      integer, parameter :: ndst = 4

      real(r8), dimension(plon, plat, ndays) :: so2_flux
      real(r8), dimension(plon, plat, ndays) :: dms_flux
      real(r8), dimension(plon, plat, ndays) :: pom_flux
      real(r8), dimension(plon, plat, ndays) :: bc_flux
      real(r8), dimension(ndst, plon, plat, ndays) :: dust_flux

      real(r8), dimension(plon, plat) :: so2_flux_interp
      real(r8), dimension(plon, plat) :: dms_flux_interp
      real(r8), dimension(plon, plat) :: pom_flux_interp
      real(r8), dimension(plon, plat) :: bc_flux_interp
      real(r8), dimension(ndst, plon, plat) :: dust_flux_interp

   contains

      subroutine aerocom_ini

         implicit none

#include <netcdf.inc>

         integer :: nc_id
         integer :: so2_flux_id, dms_flux_id, pom_flux_id, bc_flux_id, dust_flux_id

         write (*, *) 'aerocom_ini'
         call handle_ncerr( nf90_open('aerocom_emissions_2000_T62.nc', NF_NOWRITE, nc_id), &
            'aerocom_ini: error opening file '//'aerocom_emissions_2000_T62.nc' )

         call handle_ncerr( nf90_inq_varid(nc_id, 'so2_total', so2_flux_id), 'aerocom_ini: so2 id' )
         call handle_ncerr( nf90_inq_varid(nc_id, 'dms_total', dms_flux_id), 'aerocom_ini: dms id' )
         call handle_ncerr( nf90_inq_varid(nc_id, 'pom_total', pom_flux_id), 'aerocom_ini: pom id' )
         call handle_ncerr( nf90_inq_varid(nc_id, 'bc_total', bc_flux_id), 'aerocom_ini: bc id' )
         call handle_ncerr( nf90_inq_varid(nc_id, 'dust', dust_flux_id), 'aerocom_ini: dust id' )

         call handle_ncerr( nf90_get_var(nc_id, so2_flux_id, so2_flux), 'aerocom_ini: so2' )
         call handle_ncerr( nf90_get_var(nc_id, dms_flux_id, dms_flux), 'aerocom_ini: dms' )
         call handle_ncerr( nf90_get_var(nc_id, pom_flux_id, pom_flux), 'aerocom_ini: pom' )
         call handle_ncerr( nf90_get_var(nc_id, bc_flux_id, bc_flux), 'aerocom_ini: bc' )
         call handle_ncerr( nf90_get_var(nc_id, dust_flux_id, dust_flux), 'aerocom_ini: dust' )

      end subroutine aerocom_ini

      subroutine aerocom_interp(calday)

         use precision
         implicit none

         real(r8), intent(in) :: calday
         integer :: ddd

         write (*, *) 'aerocom_interp: ', calday
         ddd = floor(calday) + 1

         so2_flux_interp(:, :) = so2_flux(:, :, ddd)
         dms_flux_interp(:, :) = dms_flux(:, :, ddd)
         pom_flux_interp(:, :) = pom_flux(:, :, ddd)
         bc_flux_interp(:, :) = bc_flux(:, :, ddd)
         dust_flux_interp(:, :, :) = dust_flux(:, :, :, ddd)

       ! write (*, *) pom_flux_interp

      end subroutine aerocom_interp

      subroutine aerocom_carbon_get(klat, nlon, xx)

         implicit none

         integer, intent(in) :: klat
         integer, intent(in) :: nlon

         real(r8), dimension(nlon, 5), intent(out) :: xx

         integer :: klon

         do klon = 1, plon
            xx(klon, 1) = pom_flux_interp(klon, klat)
            xx(klon, 2) = bc_flux_interp(klon, klat)
            xx(klon, 3) = 0.0
            xx(klon, 4) = 0.0
            xx(klon, 5) = 0.0
         end do

      end subroutine aerocom_carbon_get


      subroutine aerocom_sox_get(klat, nlon, xx)

         implicit none

         integer, intent(in) :: klat
         integer, intent(in) :: nlon

         real(r8), dimension(nlon, 2), intent(out) :: xx

         integer :: klon

         do klon = 1, plon
            xx(klon, 1) = so2_flux_interp(klon, klat)
            xx(klon, 2) = 0.0
         end do

      end subroutine aerocom_sox_get

      subroutine aerocom_dms_get(klat, nlon, xx)

         implicit none

         integer, intent(in) :: klat
         integer, intent(in) :: nlon

         real(r8), dimension(nlon), intent(out) :: xx

         integer :: klon

         do klon = 1, plon
            xx(klon) = dms_flux_interp(klon, klat)
         end do

      end subroutine aerocom_dms_get

      subroutine aerocom_dust_get(klat, nlon, xx)

         implicit none

         integer, intent(in) :: klat
         integer, intent(in) :: nlon

         real(r8), dimension(nlon, ndst), intent(out) :: xx

         integer :: klon, kdst

         do kdst = 1, ndst
         do klon = 1, plon
            xx(klon, kdst) = dust_flux_interp(kdst, klon, klat)
         end do
         end do

      end subroutine aerocom_dust_get

   end module aerocom

