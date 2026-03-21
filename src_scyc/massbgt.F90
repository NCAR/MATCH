module massbgt

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Mass budget diagnostics for sulfur cycle.
   !
   ! It is assumed that the model calling this interface has been
   ! compiled so that 8 byte real data are being used.  On non CRAY
   ! machines this implies compiling with a "-r8" flag.
   ! 
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 

   use precision
   use pmgrid

   implicit none
   save
   private
   public :: &
      inimassbgt,    &! initialize module
      gamdry,        &! latitude contribution to global average mass (dry basis input)
      gamwet,        &! latitude contribution to global average mass (wet basis input)
      aal,           &! latitude contribution to area average of field
      budget2nc,     &! Write budget terms to netCDF file.
      closebudget     ! Close mass budget netCDF file.

   integer, parameter ::&
      nidstr=20              ! number of identifier strings (One for each place
                             ! in the code where a mass calculation is done.)
   real(r8) ::&
      rgravit,              &! 1./gravit
      cwava,                &! 1./(gravit*plon)
      gw(plat)               ! Gauss weights (these sum to 2)

   real(r8), dimension(plat,nidstr) ::&
      massl = 0.0            ! latitude contributions for each identifier

   integer ::&
      nrec = 1,             &! current record number
      ncid,                 &! ID for netCDF file
      vid(3+nidstr)          ! IDs for output variables

   character(len=64), dimension(nidstr) ::&
      idstr = (/         &! identifier strings
         'dms_a_fix  ',  &! dms after advection mass fixer
         'dmssnk     ',  &! dms chemistry sink
         'dms_a_tfilt',  &! dms after time filter
         'dms_b_tphys',  &! dms before tphys
         'dmssf      ',  &! dms surface flux
         'so2_a_fix  ',  &! so2 after advection mass fixer
         'so2srcg    ',  &! so2 source from gas phase chem
         'so2wet     ',  &! so2 wet deposition
         'so2_a_tfilt',  &! so2 after time filter
         'so2_b_tphys',  &! so2 before tphys
         'so2sf      ',  &! so2 surface flux
         'so2dry     ',  &! so2 dry deposition
         'so4_a_fix  ',  &! so4 after advection mass fixer
         'so4srca    ',  &! so4 source from aqueous phase chem
         'so4srcg    ',  &! so4 source from gas phase chem
         'so4wet     ',  &! so4 wet deposition
         'so4_a_tfilt',  &! so4 after time filter
         'so4_b_tphys',  &! so4 before tphys
         'so4sf      ',  &! so4 surface flux
         'so4dry     '/)  ! so4 dry deposition

!##############################################################################
contains
!##############################################################################

   subroutine inimassbgt( gravit, xgw )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Initialize module.  Open netCDF file 'mass_budget.nc' for diagnostic output.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: handle_ncerr

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         gravit,         &! Gravitational acceleration (m/s^2)
         xgw(plat)        ! Gauss weights

      ! Local variables:
      integer :: i, j, &
         recdid
      !-----------------------------------------------------------------------

      rgravit = 1./gravit
      cwava = 1./(plon*gravit)
      do j = 1, plat
         gw(j) = xgw(j)
      end do

      ! Create netCDF output file.
      call handle_ncerr( nf_create( 'mass_budget.nc', NF_CLOBBER, ncid ), &
         'inimassbgt: error opening file mass_budget.nc' )

      ! Define dimensions.
      call handle_ncerr( nf_def_dim( ncid, 'rec', NF_UNLIMITED, recdid ), &
         'inimassbgt:' )

      ! Define variables.
      call handle_ncerr( nf_def_var( ncid, 'nstep', NF_INT, 1, recdid, vid(1) ), &
         'inimassbgt:' )
      call handle_ncerr( nf_def_var( ncid, 'date', NF_INT, 1, recdid, vid(2) ), &
         'inimassbgt:' )
      call handle_ncerr( nf_def_var( ncid, 'datesec', NF_INT, 1, recdid, vid(3) ), &
         'inimassbgt:' )
      do i = 1, nidstr
         call handle_ncerr( nf_def_var( ncid, idstr(i), NF_DOUBLE, 1, recdid, vid(3+i) ), &
            'inimassbgt:' )
      end do

      ! End definitions.
      call handle_ncerr( nf_enddef( ncid ), 'inimassbgt:' )

   end subroutine inimassbgt

!##############################################################################

   subroutine gamdry( ident, lat, pdel, sh, as )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate latitude contribution to global average mass of tracer field.
      ! Input concentration units are dry mixing ratio (kg tracer / kg dry air).
      ! Output masses are in kg/m^2.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      character(len=*), intent(in) ::&
         ident                  ! identifier
      integer, intent(in) ::&
         lat                    ! latitude index
      real(r8), intent(in) ::&
         pdel(plond,plev),     &! pressure difference across layer (Pa)
         sh(plond,plev),       &! specific humidity (kg water/kg moist air)
         as(plond,plev)         ! advected species (dry mixing ratio)

      ! Local variables:
      integer :: i, k, id
      !-----------------------------------------------------------------------

      id = findindx( ident )

      massl(lat,id) = 0.0
      do k = 1, plev
         do i = 1, plon
            massl(lat,id) = massl(lat,id) + as(i,k)*(1.-sh(i,k))*pdel(i,k)
         end do
      end do
      massl(lat,id) =  massl(lat,id) * gw(lat)

   end subroutine gamdry

!##############################################################################

   subroutine gamwet( ident, lat, pdel, as )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate latitude contribution to global average mass of tracer field.
      ! Input concentration units are moist mixing ratio (kg tracer / kg moist air).
      ! Output masses are in kg/m^2.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      character(len=*), intent(in) ::&
         ident                  ! identifier
      integer, intent(in) ::&
         lat                    ! latitude index
      real(r8), intent(in) ::&
         pdel(plond,plev),     &! pressure difference across layer (Pa)
         as(plond,plev)         ! advected species (moist mixing ratio)

      ! Local variables:
      integer i, k, id
      !-----------------------------------------------------------------------

      id = findindx( ident )

      massl(lat,id) = 0.0
      do k = 1, plev
         do i = 1, plon
            massl(lat,id) = massl(lat,id) + as(i,k)*pdel(i,k)
         end do
      end do
      massl(lat,id) =  massl(lat,id) * gw(lat)

   end subroutine gamwet

!##############################################################################

   subroutine aal( ident, lat, x )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Calculate latitude contribution to area average of field.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      character(len=*), intent(in) ::&
         ident                ! identifier
      integer, intent(in) ::&
         lat                  ! latitude index
      real(r8), intent(in) ::&
         x(plond)             ! field to average

      ! Local variables:
      integer :: i, id
      !-----------------------------------------------------------------------

      id = findindx( ident )

      massl(lat,id) = 0.0
      do i = 1, plon
         massl(lat,id) = massl(lat,id) + x(i)
      end do
      massl(lat,id) =  massl(lat,id) * gw(lat)

   end subroutine aal

!##############################################################################

   real(r8) function gam( ident )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Return the global average mass (kg/m^2) of a tracer field.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      character(len=*), intent(in) ::&
         ident            ! identifier

      ! Local variables:
      integer :: j, id
      !-----------------------------------------------------------------------

      id = findindx( ident )

      gam = 0.0
      do j = 1, plat
         gam = gam + massl(j,id)
      end do
      gam =  gam * cwava * .5

   end function gam

!##############################################################################

   real(r8) function aa( ident )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Return the area average of a field.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      character(len=*), intent(in) ::&
         ident            ! identifier

      ! Local variables:
      integer j, id
      !-----------------------------------------------------------------------

      id = findindx( ident )

      aa = 0.0
      do j = 1, plat
         aa = aa + massl(j,id)
      end do
      aa =  aa * .5 / plon

   end function aa

!##############################################################################

   integer function findindx( ident )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Return index corresponding to ident string.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      character(len=*), intent(in) ::&
         ident            ! identifier

      ! Local variables:
      integer i
      !-----------------------------------------------------------------------

      findindx = 0
      do i = 1, nidstr
         if ( ident .eq. idstr(i) ) then
            findindx = i
            return
         end if
      end do

      write(6,*)'findindx: identifier not found :',ident,':'
      stop
   end function findindx

!##############################################################################

   subroutine budget2nc( nstep, date, sec )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Write budget terms to netCDF file.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: handle_ncerr

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      integer, intent(in) ::&
         nstep,   &! current timestep number
         date,    &! current date in yymmdd format
         sec       ! seconds past current date at 0Z

      ! Local variables.
      real(r8) area, dms2s, so22s, so42s
      character*10 mes
      !-----------------------------------------------------------------------

      mes = 'budget2nc:'
      area = 4._r8 * 3.14 * 6.37e6**2  ! global area in m^2
      dms2s = 1._r8/1.9375e9    ! convert kg DMS to Tg S
      so22s = 1._r8/2.e9        ! convert kg SO2 to Tg S
      so42s = 1._r8/3.e9        ! convert kg SO2 to Tg S

      call handle_ncerr( nf_put_var1_int( ncid, vid(1), nrec, nstep ),  mes )
      call handle_ncerr( nf_put_var1_int( ncid, vid(2), nrec, date ),  mes )
      call handle_ncerr( nf_put_var1_int( ncid, vid(3), nrec, sec ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(4), nrec, &
         gam( 'dms_a_fix' )*area*dms2s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(5), nrec, &
         gam( 'dmssnk' )*area*dms2s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(6), nrec, &
         gam( 'dms_a_tfilt' )*area*dms2s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(7), nrec, &
         gam( 'dms_b_tphys' )*area*dms2s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(8), nrec, &
         aa( 'dmssf' )*area*dms2s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(9), nrec, &
         gam( 'so2_a_fix' )*area*so22s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(10), nrec, &
         gam( 'so2srcg' )*area*so22s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(11), nrec, &
         gam( 'so2wet' )*area*so22s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(12), nrec, &
         gam( 'so2_a_tfilt' )*area*so22s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(13), nrec, &
         gam( 'so2_b_tphys' )*area*so22s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(14), nrec, &
         aa( 'so2sf' )*area*so22s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(15), nrec, &
         aa( 'so2dry' )*area*so22s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(16), nrec, &
         gam( 'so4_a_fix' )*area*so42s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(17), nrec, &
         gam( 'so4srca' )*area*so42s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(18), nrec, &
         gam( 'so4srcg' )*area*so42s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(19), nrec, &
         gam( 'so4wet' )*area*so42s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(20), nrec, &
         gam( 'so4_a_tfilt' )*area*so42s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(21), nrec, &
         gam( 'so4_b_tphys' )*area*so42s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(22), nrec, &
         aa( 'so4sf' )*area*so42s ),  mes )
      call handle_ncerr( nf_put_var1_double( ncid, vid(23), nrec, &
         aa( 'so4dry' )*area*so42s ),  mes )

      nrec = nrec + 1

   end subroutine budget2nc

!##############################################################################

   subroutine closebudget()

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Close mass budget netCDF file.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: handle_ncerr

      implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      call handle_ncerr( nf_close( ncid ), 'closebudget:' )

   end subroutine closebudget

!##############################################################################

end module massbgt
