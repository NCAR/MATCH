MODULE aod_assim_3d
!  
  use floats
  use analysis_3d
  use data_structures
  use pmgrid, only: plond, plon, plev, plat, pcnst, padv
  use error_messages, only: dealloc_err
  use precision
!
  IMPLICIT NONE
!
  PUBLIC :: Assimilate_AOD
  PUBLIC :: asmfluxout

  real(r8) :: tend(plon,plev,pcnst,plat) ! tendency
  real(r8) :: Ratio(plon,plev,plat)        ! ratio (in lon,lev,lat order)
  real(r8) :: nss_aod(plon, plat)     ! NOAA Pathfinder AOD - sea salt OD


CONTAINS

SUBROUTINE asmfluxout( lat, obuf )

! Write assimilation tendencies to history/output buffers.

  use histout, only: outfld
  use caer, only: caer_ncomp, caer_idx1


  IMPLICIT NONE


!-----------------------------------------------------------------------
#include <tracnam.h>
!-----------------------------------------------------------------------
#include <dstgrd.h>
!-----------------------------------------------------------------------

  
  INTEGER lat
  REAL(r8) obuf(*) 

  CHARACTER(len=3) trnum
  INTEGER m,mm
  INTEGER k

  m = 2 !SO4
  call outfld( trim(tracnam(m))//'ASM',tend(:,:,m,lat), plon, lat, obuf )

  do m = 1,dst_nbr
     mm = dst_idx_srt + m - 1
     write (trnum, '(i3)' ) m + 100
     call outfld( 'DST'//trnum(2:3)//'ASM',tend(:,:,mm,lat), plon, lat, obuf )
  end do 

  do m = 1, caer_ncomp()
     mm = caer_idx1() + m - 1
     call outfld( trim(tracnam(mm))//'ASM',tend(:,:,mm,lat), plon, lat, obuf )
  end do 

  call outfld( 'ASMRATIO',Ratio(:,:,lat),plon,lat,obuf )
  call outfld( 'ASMAOD',nss_aod(:,lat),plon,lat,obuf )

END SUBROUTINE asmfluxout

SUBROUTINE Assimilate_AOD( date, iSec, Ps, xt, T, as, delt )

!+
!
!
!-
   use sat_obs, only: satDataExists
   use caer, only: caer_ncomp, caer_idx1
   use precision

  IMPLICIT NONE

!-----------------------------------------------------------------------
#include <ssltcom.h>
!-----------------------------------------------------------------------



!
! INPUTS:
!
   INTEGER, INTENT(IN) :: date,   &! ncdate from MATCH (current date)
                          iSec     ! ncsec  from MATCH (current second)
   INTEGER, INTENT(IN) :: delt

   REAL(r8), INTENT(IN) :: &
        Ps(plond,plat), &          ! surface pressure field
        xt(plond,plev,padv,plat), &! misc advected vars. after slt on input.
                                   ! Contrib from physics added on output.     
        T(plond,plev,plat)         ! Atmospheric temperature

!                                       
! INPUTS/OUTPUTS:
!                                  
   REAL(r8), INTENT(INOUT) :: &
      as(plond,plev,pcnst,plat)    ! tracer mixing ratios
!
! PCNST indices for various mixing ratios
!
   integer, parameter :: I_SO4    = 2, &  ! SO4
                         I_DST1   = 5, &  ! dust class 1 (0.01-1 mu-m)
                         I_DST2   = 6, &  ! dust class 2 (1 -10  mu-m)
                         I_DST3   = 7, &  ! dust class 3 (10-20  mu-m)
                         I_DST4   = 8     ! dust class 4 (20-50  mu-m)
!
! INTERNALS:
!
   real(r8) ::    LatGrid(plat)   ! lat. of cell centers (deg)
   real(r8) ::    LonGrid(plon)   ! lon. of cell centers (deg)

   real(r8) :: &
        tadv(plond,plev) ,&      ! temperature corresponding to ptadv
        pint(plond,plev+1),&      ! pressure at layer interfaces        
        pmid(plond,plev),&       ! pressure at layer midpoints
        pdel(plond,plev)       ! pressure difference across layers

   INTEGER :: NLon, &              ! = plon
              NLat, &              ! = plat
              NLev                 ! = plev

   integer :: ilat                 ! latitude index for computing RH
   integer :: yyyymmdd             ! date with 4 digit year
   integer :: n_obs                ! = number of observations to assimilate

   TYPE(Orbits) :: Obs, &          ! structure to hold AOD data
                   Obs_G           ! Gridded obsservations?


    REAL(r8)  ::  p(plev,plat,plon), &            ! atmospheric pressure
                   dp(plev,plat,plon),  &         ! pressure differential
                   EXTzyx(plev,plat,plon), &       ! extinction (lev,lat,lon)
                   Ana(plev,plat,plon), &          ! analyzed extinction
                   Error(plev,plat,plon), &        ! relative (fractional) error in model flds
                   Ratiozyx(plev,plat,plon), &     ! ratio of Analyzed/Original extinction
                   RH(plond,plev,plat), &              ! relative humidity
                   EXT(plon,plev,plat), &          ! extinction (lon,lev,lat)
                   SPCEXT(plond,plev)

   real, parameter :: GRAVIT = 9.81! gravitational constant in MKS
   integer :: i, j, k, m, mm
   integer :: ier = 0

   integer :: imax,jmax,kmax
   real(r8) :: ratmax

! Externals:
   character(len=8) :: specie_name,name
   real(r8) :: ext_coeff
!-----------------------------------------------------------------------



   yyyymmdd = date
   if ( date < 2000101 ) yyyymmdd = 19000000 + date    ! add century to 2-digit dates
   if ( .not. satDataExists( yyyymmdd, iSec ) ) return

   CALL Get_AOD_Data( n_obs, Obs ) 

   write(*,*) "              # AOD obs     = ", n_obs
!
! do assimilation only if some observation are available
!
   if ( n_obs .gt. 0 ) then

      call getLatLon( 'degrees', LatGrid, LonGrid )
      
      NLon = plon
      NLat = plat
      NLev = plev
      
      CALL getPressure( ps, p, dp )
      
      p  = p /100.0 ! Pa -> mb
      dp = dp/100.0 ! Pa -> mb
!
! Calculate extinction EXTzyx here
!      
!   Compute RH
!

      do j = 1, plat
         call plevs0( ps(1,j), pint, pmid, pdel )
         call pt2t( xt(1,1,2,j), pmid, tadv )
         call relhum( ps(1,j), tadv, xt(1,1,1,j), RH(1,1,j) )
      end do


      ! Don't let relative humidity exceed .99
      ! must be consistent with physlic calculation

      do j = 1, plat
         do k = 1, plev
            do i = 1, plon
               if ( rh(i,k,j) > .99 ) rh(i,k,j) = .99
            end do
         end do
      end do

!-drb This works fine but is caught with debug flags on since plon:plond is undefined.
!      where ( rh > 0.99 )
!         rh = 0.99
!      endwhere


!
! Compute extinction (m^2/g): add in the various radiatively active tracers
!

      EXT = 0.0
      
!+drb  20020104 - implementing changes WBC 8/31/2K

      ! Sulfate
      do j = 1, nlat
         call so4_extinction( rh(1,1,j), SPCEXT )
         do k = 1, plev
            do i = 1, plon
               EXT(i,k,j) = EXT(i,k,j) + as(i,k,I_SO4,j) * SPCEXT(i,k)
            end do
         end do
     end do


     ! Sea Salt
      do j = 1, nlat
          call sslt_extinction(RH(1,1, j), SPCEXT)
          do k = 1, plev
             do i = 1, plon
                EXT(i,k,j) = EXT(i,k,j) +  sslt_profiles(i, k, j) * SPCEXT(i,k)
            end do
         end do
      end do      

          


     ! Dust
      do m = 1, 4
         mm = I_DST1 + m - 1
         EXT = EXT + as(1:plon,1:plev,mm,1:plat) * ext_coeff( specie_name(mm) )

              

      end do


      ! Carbon
      do m = 1, caer_ncomp()
         mm = caer_idx1() + m - 1
         name = specie_name(mm)
         if (   trim(name).eq.'BCPHI'  .or.   &
                trim(name).eq.'OCPHO' .or. trim(name).eq.'BCPHO'  ) then
         EXT = EXT + as(1:plon,1:plev,mm,1:plat) * ext_coeff( specie_name(mm) )
         
              

         else if (  trim(name).eq.'OCPHI'  ) then
             do j = 1, nlat
                call ocphi_extinction(rh(1,1, j), SPCEXT)
                do k = 1, plev
                   do i = 1, plon
                      EXT(i,k,j) = EXT(i,k,j) +  as(i, k, mm, j) * SPCEXT(i,k)
                   end do
                end do
             end do  
         endif


      end do
!-drb 20020104 - implementing changes WBC 8/31/2K

!   
!   Multiply by dp / g to get AOD per layer
!

      CALL xzy2zyx( EXT, EXTzyx )



      EXTzyx = EXTzyx * (dp * 100.0) / GRAVIT
!
! Perform analysis
!


!      Error(:,:,:) = 2.0
       Error(:,:,:) = 1.0 ! high value Ed4
!      Error(:,:,:) = 0.5 ! original value Ed2
!      Error(:,:,:) = 0.3
      
      CALL Analysis_XYZ_3d( NLev,NLat,NLon, Obs, P, dp, LatGrid, LonGrid, &
                            EXTzyx, Error, Ana, Obs_G)

! Hack from Boris to prevent creating negative mixing ratios.
      WHERE ( Ana  < 0. )
         Ana = 1.0E-5
      ENDWHERE


!
! Use the ratio of the analyzed extinction to original extinction
!     to modify mixing ratios
!


      Ratiozyx = Ana / EXTzyx
      CALL zyx2xzy( Ratiozyx, Ratio )

      print *,'range of ratio corrections',minval(Ratio(:plon,:plev,:plat)), &
           maxval(Ratio(:plon,:plev,:plat))
      imax = 1
      jmax = 1
      kmax = 1
      ratmax = 0.0
      do k = 1,plev
         do j = 1,plat
            do i = 1,plon
               if ( ratio(i,k,j).gt.ratmax ) then
                  imax = i
                  jmax = j
                  kmax = k
                  ratmax = ratio(i,k,j)
               endif
            end do  !i
         end do     !j
      end do        !k
      print *,'RATIO MAX at imax,kmax,jmax are ',ratmax,imax,kmax,jmax

      ! apply ratio from layer 11 down DWF
      do k = 1,10
         Ratio(:,k,:) = 1.0
      end do        !k
     
      !store to calculate tendency
      tend = as(1:plon,:,:,:)

      as(1:plon, 1:plev, I_SO4,    1:plat) = as(1:plon, 1:plev, I_SO4,    1:plat) * Ratio
      as(1:plon, 1:plev, I_DST1,   1:plat) = as(1:plon, 1:plev, I_DST1,   1:plat) * Ratio
      as(1:plon, 1:plev, I_DST2,   1:plat) = as(1:plon, 1:plev, I_DST2,   1:plat) * Ratio
      as(1:plon, 1:plev, I_DST3,   1:plat) = as(1:plon, 1:plev, I_DST3,   1:plat) * Ratio
      as(1:plon, 1:plev, I_DST4,   1:plat) = as(1:plon, 1:plev, I_DST4,   1:plat) * Ratio
      ! Carbon
      do m = 1, caer_ncomp()
         mm = caer_idx1() + m - 1
         as(1:plon,1:plev,mm,1:plat) = as(1:plon,1:plev,mm,1:plat) * Ratio
      end do

      !calculate tendency
      tend = ( as(1:plon,:,:,:) - tend ) / delt


!      print *,'as max after',as(imax,kmax,:,jmax)
!      print *,'extinctions at imax,jmax',Ana(:,jmax,imax)
      ratmax = 0.0
      do k = 1,plev
         ratmax = ratmax + Ana(k,jmax,imax)
      end do
!      print *,'exinctions at imax,jmax sum',ratmax

!
! Finish up by writing out old and new extinctions, deallocate storage
!
!      print *,'aod_assim_3d L234'
!      print *,'yyyymmdd',yyyymmdd
!      print *,'isec', isec
!-!      print *,'Obs',Obs
!-!      print *,'Obs_G',Obs_G
!-ok      print *, 'LatGrid',LatGrid
!-ok      print *,'LonGrid',LonGrid
!      print *,'P',size(P)
!      print *,'dp',size(dp)
!      print *,'EXTzyx',size(EXTzyx)
!      print *,'Ana',size(Ana)


!
       print *,'AOD_ASSIM output file is not being produced'
!     CALL NCDFOutput_3d( yyyymmdd, isec, Obs, Obs_G, LatGrid,LonGrid, &
!


      IF ( ASSOCIATED(Obs_G%Profile) ) DEALLOCATE(Obs_G%Profile, stat=ier)
      call dealloc_err( ier, 'Assimilate_AOD', 'Obs_G%Profile' )
      IF ( ASSOCIATED(Obs%Profile)   ) DEALLOCATE(Obs%Profile, stat=ier)
      call dealloc_err( ier, 'Assimilate_AOD', 'Obs%Profile' )

   endif
   return
!


END SUBROUTINE Assimilate_AOD


SUBROUTINE  xzy2zyx( A, B )
!+
!
!-
! INPUTS:
  REAL(r8), INTENT(IN) :: A(:,:,:)
! OUTPTS:
  REAL(r8), INTENT(OUT) :: B( SIZE(A(1,:,1)), SIZE(A(1,1,:)), SIZE(A(:,1,1)) )
!
  INTEGER :: i
!
  DO i=1,SIZE(A(:,1,1))
     B(:,:,i) = A(i,:,:)
  ENDDO   
!
RETURN
END SUBROUTINE  xzy2zyx




SUBROUTINE  zyx2xzy( A, B )
!+
!
!-
! INPUTS:
  REAL(r8), INTENT(IN) :: A(:,:,:)
! OUTPTS:
  REAL(r8), INTENT(OUT) :: B( SIZE(A(1,1,:)), SIZE(A(:,1,1)), SIZE(A(1,:,1)) )
!
  INTEGER :: i
!
  DO i=1,SIZE(A(1,1,:))
     B(i,:,:) = A(:,:,i)
  ENDDO   
!
RETURN
END SUBROUTINE  zyx2xzy




SUBROUTINE  xyz2zyx( A, B )
!+
!
!-
! INPUTS:
  REAL(r8), INTENT(IN) :: A(:,:,:)
! OUTPTS:
  REAL(r8), INTENT(OUT) :: B( SIZE(A(1,1,:)), SIZE(A(1,:,1)), SIZE(A(:,1,1)) )
!
  INTEGER :: i,j,k
!
  DO i=1,SIZE(A(1,1,:))
     B(i,:,:) = TRANSPOSE(A(:,:,i))
  ENDDO   

!
RETURN
END SUBROUTINE  xyz2zyx


SUBROUTINE  zyx2xyz( B, A )
!+
!
!-
! INPUTS:
  REAL(r8), INTENT(IN) :: B(:,:,:) 
! OUTPTS:
  REAL(r8), INTENT(OUT) :: A( SIZE(B(1,1,:)), SIZE(B(1,:,1)), SIZE(B(:,1,1)) )

!
  INTEGER :: i,j,k
!
  DO i=1,SIZE(B(:,1,1))
     A(:,:,i) = TRANSPOSE( B(i,:,:) )
  ENDDO   

!
RETURN
END SUBROUTINE  zyx2xyz


SUBROUTINE Get_AOD_Data( n_obs, Obs )
!
! This subroutine retrieves the AOD data from a specific
! file output by aerod.f. The data is filtered for cloudy regions.
! Each retrieved data is stored in the output array obs.
!
! JFL June 98
!
      use data_structures, only : orbits
      use floats
      use sat_obs, only: getNSSAOD
      use precision
!

      implicit none
!
! output
!
      type(orbits), intent(out) :: obs
      integer     , intent(out) :: n_obs
!
! local
!
      integer   :: n, &                    ! Index into observations
                   ilon, &                 ! Longitude index
                   ilat, &                 ! Latitude index
                   aodate, aosec           ! date and seconds for sat. obs

      real(r8) :: data_lat(plat), &       ! Latitudes of gridded sat. obs
                   data_lon(plon)       ! Longitude of gridded sat. obs

      integer :: ier = 0
!
! parameters
!
!      real(r8), parameter :: error   = 0.5   ! error in AOD is fixed at 50%
      real(r8), parameter :: error   = 0.2   ! error in AOD is fixed at 50%


      call getNSSAOD( aodate, aosec, data_lat, data_lon, nss_aod )

      n_obs = count(nss_aod .ge. 0.0)
!
! now, we know how many observations we have
! so that we can allocate the required memory
!
! CHANGE, removed following lines
!     if(associated(obs%profile)) deallocate(obs%profile, stat=ier )
!     call dealloc_err( ier, 'Get_AOD_Data', 'obs%profile' )
      if (n_obs .eq. 0) then 
         return
      endif
      allocate(obs%profile(n_obs))
!
! store data in orbit type structure
!
      obs%n            = n_obs
      obs%nodata       = -1.0
!
! the profile_grid is not used anywhere -- set to nonsense value
      obs%profile_grid = 0.340E+39
      n = 0
      do ilat = 1, plat
         do ilon = 1, plon
            if (nss_aod(ilon, ilat) .ge. 0.0) then
               n = n + 1
               obs%profile(n)%lat    = data_lat(ilat)
               obs%profile(n)%lon    = data_lon(ilon)
               obs%profile(n)%date   = aodate
               obs%profile(n)%secs   = aosec
               obs%profile(n)%v      = nss_aod(ilon,ilat)
               obs%profile(n)%sigmav = nss_aod(ilon,ilat) * error
!	       write(*,*) obs%profile(n)%lat, obs%profile(n)%lon, &
!			obs%profile(n)%date, obs%profile(n)%secs, &
!			obs%profile(n)%v(1)
            endif
         end do
      end do
      write (*, *) 'TEST - GET_AOD_DATA exit'
!
      return
      end subroutine GET_AOD_DATA
!

SUBROUTINE NCDFOutput_3d( yymmdd, isec, Obs, Obs_G, LatGrid,LonGrid, P, dp, EXT, EXT_New)
!
!
!
  use precision
  use netcdf 
  IMPLICIT NONE  
  INCLUDE 'netcdf.inc'
!
  INTEGER, INTENT(IN) :: yymmdd, isec
  TYPE(Orbits), INTENT(IN) :: Obs, Obs_G
  REAL(r8), INTENT(IN) :: LatGrid(:), LonGrid(:), P(:,:,:), EXT(:,:,:), EXT_New(:,:,:)
  REAL(r8), INTENT(IN) ::                        dP(:,:,:)
!
  CHARACTER*5 :: siSec
  CHARACTER*8 :: syyyymmdd
  INTEGER :: i
  REAL(r8), POINTER :: Obs_G_sigmaV(:,:), Obs_G_V(:,:)
  REAL(r8), POINTER :: Obs_sigmaV(:,:), Obs_V(:,:)
!
!
! INTERNALS:      
      INTEGER :: Status, ncid
      INTEGER :: Obs_N_DimID, Obs_G_N_DimID, N_MOP_Lev_DimID, NLat_DimID, NLon_DimID,  NLev_DimID
      INTEGER :: Date_ID, Sec_ID, NoData_ID, Lat_ID, Lon_ID, Lev_ID, Obs_Lat_ID, Obs_Lon_ID, Obs_Secs_ID, Obs_Date_ID
      INTEGER :: Obs_G_Lat_ID, Obs_G_Lon_ID, Obs_G_Secs_ID, Obs_G_Date_ID, Obs_P_ID, Obs_V_ID, Obs_sigmaV_ID
      INTEGER :: Obs_G_V_ID, Obs_G_sigmaV_ID, P_ID, EXT_ID, EXT_New_ID,dp_id
      INTEGER VDim1(1), VDim2(2), VDim3(3)
      CHARACTER*100 :: File

  character(len=15) :: routine = 'ncdfoutput_3d: '
  character(len=32) :: vname

!
  WRITE(siSec,'(I5.5)') iSec
!  WRITE(syyyymmdd,'(I8.8)') yymmdd + 19000000
  WRITE(syyyymmdd,'(I8.8)') yymmdd


  ALLOCATE( Obs_G_V( SIZE(Obs_G%Profile_Grid), Obs_G%N ) )
  ALLOCATE( Obs_G_sigmaV( SIZE(Obs_G%Profile_Grid), Obs_G%N ) )
  ALLOCATE( Obs_V( SIZE(Obs%Profile_Grid), Obs%N ) )
  ALLOCATE( Obs_sigmaV( SIZE(Obs%Profile_Grid), Obs%N ) )  

  DO i=1, SIZE(Obs_G%Profile_Grid)
     Obs_G_V(i,:) = Obs_G%Profile(:)%V(i)
     Obs_G_sigmaV(i,:) = Obs_G%Profile(:)%sigmaV(i)
  ENDDO   

  DO i=1, SIZE(Obs%Profile_Grid)
     Obs_V(i,:) = Obs%Profile(:)%V(i)
     Obs_sigmaV(i,:) = Obs%Profile(:)%sigmaV(i)
  ENDDO   


!
      File = 'AOD_ASSIM_'//syyyymmdd//'-'//siSec//'.nc'
      PRINT *, 'NCDFOutput% Creating a netCDF file ',TRIM(File)

      Status = NF_CREATE( TRIM(File), NF_CLOBBER, ncid )
      call handle_ncerr( status, routine//'error opening '//TRIM(File) )

!
! Define dimensions:
      Status = NF_DEF_DIM( ncid, 'Obs_N',  Obs%N,   Obs_N_DimID   )
      call handle_ncerr( status, routine//'error defining dimension Obs_N' )
      Status = NF_DEF_DIM( ncid, 'Obs_G_N',  Obs_G%N,   Obs_G_N_DimID   )
      call handle_ncerr( status, routine//'error defining dimension Obs_G_N' )

      Status = NF_DEF_DIM( ncid, 'N_MOP_Lev',  SIZE(Obs%Profile_Grid),   N_MOP_Lev_DimID   )
      call handle_ncerr( status, routine//'error defining dimension N_MOP_Lev' )

      Status = NF_DEF_DIM( ncid, 'lat',  SIZE(LatGrid),   NLat_DimID )
      call handle_ncerr( status, routine//'error defining dimension lat' )
      Status = NF_DEF_DIM( ncid, 'lon',  SIZE(LonGrid),   NLon_DimID )
      call handle_ncerr( status, routine//'error defining dimension lon' )
      Status = NF_DEF_DIM( ncid, 'lev',  SIZE(P(:,1,1)),   NLev_DimID )
      call handle_ncerr( status, routine//'error defining dimension lev' )

!
! Define scalar variables first:
      Status = NF_DEF_VAR( ncid, 'date', NF_INT, 0, VDim1, Date_ID )
      call handle_ncerr( status, routine//'error defining variable date' )

      Status = NF_DEF_VAR( ncid, 'time', NF_INT, 0, VDim1, Sec_ID )
      call handle_ncerr( status, routine//'error defining variable time' )

      Status = NF_DEF_VAR( ncid, 'no_data', NF_REAL, 0, VDim1, NoData_ID )
      call handle_ncerr( status, routine//'error defining variable no_data' )

!
! Define array variables:
      VDim1(1) = NLat_DimID
      Status = NF_DEF_VAR( ncid, 'lat', NF_FLOAT, 1, VDim1, Lat_ID )
      call handle_ncerr( status, routine//'error defining variable lat' )

      VDim1(1) = NLon_DimID
      Status = NF_DEF_VAR( ncid, 'lon', NF_FLOAT, 1, VDim1, Lon_ID )
      call handle_ncerr( status, routine//'error defining variable lon' )

      VDim1(1) = NLev_DimID
      Status = NF_DEF_VAR( ncid, 'lev', NF_FLOAT, 1, VDim1, Lev_ID )
      call handle_ncerr( status, routine//'error defining variable lev' )

      VDim1(1) = Obs_N_DimID
      Status = NF_DEF_VAR( ncid, 'obs_lat', NF_FLOAT, 1, VDim1, Obs_Lat_ID )
      call handle_ncerr( status, routine//'error defining variable obs_lat' )
      VDim1(1) = Obs_N_DimID
      Status = NF_DEF_VAR( ncid, 'obs_lon', NF_FLOAT, 1, VDim1, Obs_Lon_ID )
      call handle_ncerr( status, routine//'error defining variable obs_lon' )
      VDim1(1) = Obs_N_DimID
      Status = NF_DEF_VAR( ncid, 'obs_secs', NF_FLOAT, 1, VDim1, Obs_Secs_ID )
      call handle_ncerr( status, routine//'error defining variable obs_secs' )
      VDim1(1) = Obs_N_DimID
      Status = NF_DEF_VAR( ncid, 'obs_date', NF_INT, 1, VDim1, Obs_Date_ID )
      call handle_ncerr( status, routine//'error defining variable obs_date' )

      VDim1(1) = Obs_G_N_DimID
      Status = NF_DEF_VAR( ncid, 'obs_g_lat', NF_FLOAT, 1, VDim1, Obs_G_Lat_ID )
      call handle_ncerr( status, routine//'error defining variable obs_g_lat' )
      VDim1(1) = Obs_G_N_DimID
      Status = NF_DEF_VAR( ncid, 'obs_g_lon', NF_FLOAT, 1, VDim1, Obs_G_Lon_ID )
      call handle_ncerr( status, routine//'error defining variable obs_g_lon' )
      VDim1(1) = Obs_G_N_DimID
      Status = NF_DEF_VAR( ncid, 'obs_g_secs', NF_FLOAT, 1, VDim1, Obs_G_Secs_ID )
      call handle_ncerr( status, routine//'error defining variable obs_g_secs' )
      VDim1(1) = Obs_G_N_DimID
      Status = NF_DEF_VAR( ncid, 'obs_g_date', NF_INT, 1, VDim1, Obs_G_Date_ID )
      call handle_ncerr( status, routine//'error defining variable obs_g_date' )

      VDim1(1) = N_MOP_Lev_DimID
      Status = NF_DEF_VAR( ncid, 'obs_p', NF_FLOAT, 1, VDim1, Obs_P_ID )
      call handle_ncerr( status, routine//'error defining variable obs_p' )

      VDim2(1) = N_MOP_Lev_DimID
      VDim2(2) = Obs_N_DimID
      Status = NF_DEF_VAR(ncid, 'obs_v', NF_FLOAT, 2, VDim2, Obs_V_ID)
      call handle_ncerr( status, routine//'error defining variable obs_v' )
      Status = NF_DEF_VAR(ncid, 'obs_sigmav', NF_FLOAT, 2, VDim2, Obs_sigmaV_ID)
      call handle_ncerr( status, routine//'error defining variable obs_sigmav' )

      VDim2(1) = N_MOP_Lev_DimID
      VDim2(2) = Obs_G_N_DimID
      Status = NF_DEF_VAR(ncid, 'obs_g_v', NF_FLOAT, 2, VDim2, Obs_G_V_ID)
      call handle_ncerr( status, routine//'error defining variable obs_g_v' )
      Status = NF_DEF_VAR(ncid, 'obs_g_sigmav', NF_FLOAT, 2, VDim2, Obs_G_sigmaV_ID)
      call handle_ncerr( status, routine//'error defining variable obs_g_sigmav' )

      VDim3(1) = NLev_DimID
      VDim3(2) = NLat_DimID
      VDim3(3) = NLon_DimID
!      Status = NF_DEF_VAR(ncid, 'pressure', NF_FLOAT, 3, VDim3, P_ID)
!      CALL HANDLE_CDF_ERR( Status )
!      Status = NF_DEF_VAR(ncid, 'deltapressure', NF_FLOAT, 3, VDim3, dP_ID)
!      CALL HANDLE_CDF_ERR( Status )
      vname = 'EXT'
      Status = NF_DEF_VAR(ncid, trim(vname), NF_FLOAT, 3, VDim3,  EXT_ID)
      call handle_ncerr( status, routine//'error defining variable: '//trim(vname) )
      Status = NF_DEF_VAR(ncid, 'EXT_new', NF_FLOAT, 3, VDim3,  EXT_New_ID)
      call handle_ncerr( status, routine//'error defining variable EXT_new' )
!
! Leave definition mode:
      Status = NF_ENDDEF( ncid )
      call handle_ncerr( status, routine//'error leaving definition mode' )

!
! And write some variables:
 
      PRINT *, 'NCDFOutput% Dumping yymmdd ...'
      Status = NF_PUT_VAR_INT( ncid, Date_ID, yymmdd )
      call handle_ncerr( status, routine//'error writing date' )
      PRINT *, 'NCDFOutput% Dumping isec ...'
      Status = NF_PUT_VAR_INT( ncid, Sec_ID, isec )
      call handle_ncerr( status, routine//'error writing time' )

      PRINT *, 'NCDFOutput% Dumping LatGrid ...'
      status = NF90_PUT_VAR( ncid, Lat_ID, LatGrid )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: LatGrid' )

      PRINT *, 'NCDFOutput% Dumping LonGrid...'
      status = NF90_PUT_VAR( ncid,Lon_ID, LonGrid )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: LonGrid' )

      status = NF90_PUT_VAR( ncid,Lev_ID,P(:,1,1) )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: LevGrid' )

      PRINT *, 'NCDFOutput% Dumping  Obs%Profile_Grid...'
!      print *,'aod_assim_3d L613',Obs%Profile_Grid

      status = NF90_PUT_VAR( ncid,Obs_p_ID, Obs%Profile_Grid )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: Obs%Profile_Grid ' )

      PRINT *, 'NCDFOutput% Dumping Obs%NoData ...'
      status = NF90_PUT_VAR( ncid, NoData_ID, Obs%NoData )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: Obs%NoData' )

      PRINT *, 'NCDFOutput% Dumping Obs_Lon. MIN/MAX:', MINVAL(Obs%Profile%Lon), &
           MAXVAL(Obs%Profile%Lon)
      status = NF90_PUT_VAR( ncid,Obs_Lon_ID , Obs%Profile%Lon )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: Obs%Profile%Lon' )

      PRINT *, 'NCDFOutput% Dumping Obs_Lat. MIN/MAX:', MINVAL(Obs%Profile%Lat), &
           MAXVAL(Obs%Profile%Lat)
      status = NF90_PUT_VAR( ncid,Obs_Lat_ID, Obs%Profile%Lat)
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: Obs%Profile%Lat' )

!      Status = NF_PUT_VAR_INT( ncid, Obs_Date_ID, Obs%Profile%Date )
!      CALL HANDLE_CDF_ERR( Status )
!      call nfPutVarReal( ncid, Obs_Secs_ID, Obs%Profile%Secs )

      PRINT *, 'NCDFOutput% Dumping Obs_G_Lon. MIN/MAX:', MINVAL(Obs_G%Profile%Lon), &
           MAXVAL(Obs_G%Profile%Lon)
      status = NF90_PUT_VAR( ncid,Obs_G_Lon_ID , Obs_G%Profile%Lon)
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: Obs_G%Profile%Lon' )

      PRINT *, 'NCDFOutput% Dumping Obs_G_Lat. MIN/MAX:', MINVAL(Obs_G%Profile%Lat), &
           MAXVAL(Obs_G%Profile%Lat)
      status = NF90_PUT_VAR( ncid,  Obs_G_Lat_ID,Obs_G%Profile%Lat )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: Obs_G%Profile%Lat' )
!      Status = NF_PUT_VAR_INT( ncid, Obs_G_Date_ID, Obs_G%Profile%Date )
!      CALL HANDLE_CDF_ERR( Status )
!      call nfPutVarReal( ncid, Obs_G_Secs_ID, Obs_G%Profile%Secs )

      PRINT *, 'NCDFOutput% Dumping Obs_G_V. MIN/MAX:', MINVAL(Obs_G_V), MAXVAL(Obs_G_V)
      status = NF90_PUT_VAR( ncid, Obs_G_V_ID,Obs_G_V  )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: Obs_G_V' )

      PRINT *, 'NCDFOutput% Dumping Obs_G_sigmaV. MIN/MAX:', &
           MINVAL(Obs_G_sigmaV), MAXVAL(Obs_G_sigmaV)
      status = NF90_PUT_VAR( ncid, Obs_G_sigmaV_ID, Obs_G_sigmaV )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: Obs_G_sigmaV' )

      PRINT *, 'NCDFOutput% Dumping Obs_V. MIN/MAX:', MINVAL(Obs_V), MAXVAL(Obs_V)
      status = NF90_PUT_VAR( ncid,Obs_V_ID, Obs_V )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: Obs_V ' )

      PRINT *, 'NCDFOutput% Dumping Obs_sigmaV. MIN/MAX:', &
           MINVAL(Obs_sigmaV), MAXVAL(Obs_sigmaV)
      status = NF90_PUT_VAR( ncid,Obs_sigmaV_ID, Obs_sigmaV  )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: Obs_sigmaV ' )


!      PRINT *, 'NCDFOutput% Dumping P. MIN/MAX:', MINVAL(P), MAXVAL(P)
!      call nfPutVarReal( ncid, P_ID, P )
!      PRINT *, 'NCDFOutput% Dumping dP. MIN/MAX:', MINVAL(dP), MAXVAL(dP)
!      call nfPutVarReal( ncid, dP_ID, dP )
!      PRINT *, 'NCDFOutput% Dumping EXT. MIN/MAX:', MINVAL(EXT), MAXVAL(EXT)
      status = NF90_PUT_VAR( ncid,EXT_ID, EXT )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: EXT ' )

!     PRINT *, 'NCDFOutput% Dumping EXT_New. MIN/MAX:', MINVAL(EXT_New), MAXVAL(EXT_New)
      status = NF90_PUT_VAR( ncid,EXT_New_ID, EXT_New )
      call handle_ncerr( status, 'NF90_PUT_VAR: error return: EXT_new' )

! Close the file and leave:
      Status = NF_CLOSE( ncid )
      call handle_ncerr( status, routine//'error closing file' )

      PRINT *, 'NCDFOutput% Done. File: ', TRIM(File)


      DEALLOCATE( Obs_G_V )
      DEALLOCATE( Obs_G_sigmaV )
      DEALLOCATE( Obs_V )
      DEALLOCATE( Obs_sigmaV )  
!
      RETURN
      END SUBROUTINE NCDFOutput_3d

END MODULE aod_assim_3d
