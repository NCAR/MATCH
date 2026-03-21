
module covariances_3d

use floats
use data_structures
use geo
use f2v
use gaussian_weight
use obs_operator
use error_messages, only: dealloc_err
use precision

contains

subroutine obs_covariance_xyz(n_obs,nlev_obs,obs,latobs,lonobs &
                             ,npt_good,ind_good,o_xyz)
!
! for nadir-viewing instruments, independance in the horizontal is
! a good approximation.  So it is only necessary to calculate the
! values on the main diagonal.
!
! JFL Sep-98
!
      implicit none
! 
! inputs:
!
      integer     , intent(in) :: n_obs,nlev_obs
      type(orbits), intent(in) :: obs
      real(r8)   , intent(in) :: latobs(n_obs),lonobs(n_obs)
      integer     , intent(in) :: npt_good
      integer     , intent(in) :: ind_good(:,:)
!
! outputs:
!
      real(r8), intent(out) :: o_xyz(npt_good,npt_good)
! 
! internals:
!
      integer i,j,n
      real(r8) :: variance
      real(r8) :: dr        ! earthdist
      real(r8) :: l = 100.0 ! correlation length (km)
!      real(r8), parameter :: minvar = 0.0e-6 ! estimated accuracy or minimal variance
      real(r8), parameter :: minvar = 0.04 ! estimated accuracy or minimal variance
!
! reset covariance matrix
!
      o_xyz = 0.0
!
! define diagonal
!
      do n=1,npt_good
        i = ind_good(n,1)
        j = ind_good(n,2)
        variance = (obs%profile(i)%sigmav(j))**2 &
                 + minvar**2
        o_xyz(n,n) = variance
      enddo
!
!     print *,'done with obs_covariance_xyz'
!
      return
end subroutine obs_covariance_xyz

SUBROUTINE HxV( Nx, Ny, H, nH, V, HV )
  !+
  !     Product of the observation operator H and a vector V.
  !
!     H has 3 non-zero elements in each raw. 
!     H has Ny columns.
!
!     Indices if the non-zero elements are kept in nH(3,Ny). 
!     These indices range from 1 to Nx.
!-
      use precision
      IMPLICIT NONE
!
! INPUTS:
      INTEGER, INTENT(IN) :: Nx, Ny, nH(3,Ny)
      REAL(r8), INTENT(IN) :: V(Nx), H(3,Ny)
! OUTPUTS:
      REAL(r8), INTENT(OUT) :: HV(Ny)
!
! INTERNALS:
      INTEGER i,j,k
!
      DO i=1,Ny
       HV(i) = H(1,i)*V(nH(1,i)) + H(2,i)*V(nH(2,i)) + H(3,i)*V(nH(3,i))
      ENDDO
!
      RETURN
END SUBROUTINE HxV


subroutine interp_3d(nlon,nlat,nlev,n_g,nlev_obs    &
                    ,npt_good,ind_good,map_good     &
                    ,background,error_model         &
                    ,p_model,dp_model               &
                    ,hxy,nhxy,obs_profile_grid      &
                    ,background_int,error_model_int &
                    ,p_model_int,hz)
!
      use precision
!      implicit none
      implicit real(r8)(a-h,o-z)
!
! input
!
      integer  , intent(in) :: nlon,nlat,nlev,n_g,nlev_obs,nhxy &
                              ,npt_good
      integer  , intent(in) :: ind_good
      logical  , intent(in) :: map_good
      real(r8), intent(in) :: background,error_model,p_model,dp_model,hxy & 
                              ,obs_profile_grid

      dimension background      (nlev,nlat,nlon)   &
               ,error_model     (nlev,nlat,nlon)   &
               ,p_model         (nlev,nlat,nlon)   &
               ,dp_model        (nlev,nlat,nlon)   &
               ,hxy             (3,n_g)            &
               ,nhxy            (3,n_g)            &
               ,ind_good        (npt_good,2)       &
               ,map_good        (n_g,nlev_obs)     &
               ,obs_profile_grid(nlev_obs)
!
! output
!
      real(r8), intent(out) :: background_int,error_model_int,p_model_int &
                               ,hz
      dimension background_int (npt_good     )     &
               ,error_model_int(npt_good,nlev)     &
               ,p_model_int    (npt_good,nlev)     &
               ,hz             (npt_good,nlev)
!
! local
!
      real(r8) :: vector_1(nlon*nlat)
      real(r8) :: vector_2(nlev,n_g) &
                  ,vector_3(nlev,n_g) &
                  ,vector_4(nlev,n_g) &
                  ,vector_5(nlev,n_g)
      real(r8) :: wk_hz(nlev,nlev_obs)
      integer   :: i,j,k,npt,index
      logical   :: test = .false.
!
      npt = nlon * nlat
!
! interpolate various fields horizontally
!
!     write(113) n_g
!     write(113) hxy
!     write(113) nhxy
!
      do i=1,nlev ! loop over number of model levels
        call field2vector(nlat,nlon,error_model(i,:,:),vector_1)
        call hxv(npt,n_g,hxy,nhxy,vector_1,vector_2(i,:))
        call field2vector(nlat,nlon,p_model (i,:,:),vector_1)
        call hxv(npt,n_g,hxy,nhxy,vector_1,vector_3(i,:))
        call field2vector(nlat,nlon,dp_model(i,:,:),vector_1)
        call hxv(npt,n_g,hxy,nhxy,vector_1,vector_4(i,:))
        call field2vector(nlat,nlon,background(i,:,:),vector_1)
        call hxv(npt,n_g,hxy,nhxy,vector_1,vector_5(i,:))
      end do
!
      index = 0
      do i=1,n_g
        do j=1,nlev_obs
!
! test if this location has good observational data
! and map result into a vector
!
          test = map_good(i,j)
          if(.not.test) cycle
!
          index = index + 1
          if(index.gt.npt_good) then
            print *,'in interp_3d, index>npt_good'
            print *,index,npt_good
            call abort()
          endif
          error_model_int(index,:) = vector_2(:,i)
          p_model_int    (index,:) = vector_3(:,i)
!
! compute averaging kernel
!
          call obs_operator_z(obs_profile_grid,vector_3(:,i),vector_4(:,i),wk_hz)
          hz(index,:) = wk_hz(:,j)
!
! convolute with kernel
!
          background_int(index) = dot_product(vector_5(:,i),wk_hz(:,j))
!
        end do
!
      end do
!
! make sure all the points have been evaluated
!
      if(index.ne.npt_good) then
        print *,'in interp_3d,'
        print *,'index.ne.npt_good ',index,npt_good
        call abort()
      endif
!
!     print *,'done with interp_3d'
!
      return
end subroutine interp_3d

subroutine compute_hbht(n_g,nlev_obs       &
                       ,nlat,nlon,nlev     &
                       ,npt_good,ind_good  &
                       ,lats_g,lons_g      &
                       ,p_model_int        &
                       ,error_model_int,hz &
                       ,o_xyz)
!
      use precision
!      implicit none
      implicit real(r8)(a-h,o-z)
!
! input
!
      integer  ,intent(in)    :: n_g,nlev_obs,nlat,nlon,nlev,npt_good
      integer  ,intent(in)    :: ind_good
      real(r8),intent(in)    :: p_model_int,lats_g,lons_g,error_model_int
      real(r8),intent(inout) :: hz,o_xyz

      dimension lats_g         (n_g)               &
               ,lons_g         (n_g)               &
               ,ind_good       (npt_good,2)        &
               ,error_model_int(npt_good,nlev)     &
               ,p_model_int    (npt_good,nlev)     &
               ,hz             (npt_good,nlev)     &
               ,o_xyz          (npt_good,npt_good)
!
! local
!
      integer :: i,j,k,l,m
      integer :: ilat1,ilon1,ilev1,ilat2,ilon2,ilev2
      integer :: index1,index2
      real(r8) :: gauss_xy,gauss_z
      real(r8) :: dl
      real(r8) :: correl_length_xy = 100. ! km
      real(r8) :: correl_length_z  =    1.00
      real(r8), dimension(nlev) :: row_work
      real(r8), dimension(nlev) :: block_hbht
!
      LOOP_I : do i=1,npt_good
        index1 = ind_good(i,1)
        LOOP_J : do j=1,i
!
! horizontal correlation
!
          index2 = ind_good(j,1)
          dl = earthdist( lats_g(index1), lons_g(index1), &
                          lats_g(index2), lons_g(index2) )
          if(dl.ge.(2.*correl_length_xy)) cycle LOOP_J
          gauss_xy = gauss(dl,correl_length_xy)
!
          do k=1,nlev
            do l=1,nlev
!
! vertical correlation
!
              dl = abs(log(p_model_int(i,k)/p_model_int(j,l)))
              gauss_z = gauss(dl,correl_length_z)
!
              row_work(l) = sqrt(error_model_int(i,k)*error_model_int(j,l)) &
                          * gauss_xy * gauss_z
            end do
            block_hbht(k) = dot_product(row_work,hz(j,:))
          end do
          o_xyz(i,j) = o_xyz(i,j) + dot_product(block_hbht,hz(i,:))
        end do LOOP_J
      end do LOOP_I
!
! use symmetry to fill the rest of the array
!
      do i=1,npt_good
        do j=i+1,npt_good
          o_xyz(i,j) = o_xyz(j,i)
        end do
      end do
!
!     print *,'done with compute_hbht'
!
      return
end subroutine compute_hbht
!
subroutine compute_analysis(n_g,nlev_obs,nlev,nlat,nlon    &
                           ,npt_good,ind_good              &
                           ,lats_g,lons_g                  &
                           ,latgrd,longrd,error_model_int  &
                           ,background,error_model         &
                           ,p_model,p_model_int            &
                           ,hz,lambda,analysis)
!
! increment is equal to the scalar product of BHt and lambda
! (Eq(13) from Menard).
!
      use precision
!      implicit none
      implicit real(r8)(a-h,o-z)
!
! inputs
!
      integer  ,intent(in) :: n_g,nlev_obs,nlev,nlat,nlon &
                             ,npt_good,ind_good
      real(r8),intent(in) :: lats_g,lons_g           &
                             ,error_model_int &
                             ,background,error_model        &
                             ,p_model,p_model_int           &
                             ,hz,lambda

!-drb                             ,latgrd,longrd, &


      dimension ind_good(npt_good,2)
      dimension background (nlev,nlat,nlon) &
               ,error_model(nlev,nlat,nlon) &
               ,p_model    (nlev,nlat,nlon)
      dimension p_model_int(npt_good,nlev)
      dimension lats_g(n_g ),lons_g(n_g )
      real(r8), dimension(nlat) :: latgrd
      real(r8), dimension(nlon) :: longrd
      dimension lambda         (npt_good)      &
               ,hz             (npt_good,nlev) &
               ,error_model_int(npt_good,nlev)
!
! outputs
!
      real(r8),intent(out) :: analysis
      dimension analysis(nlev,nlat,nlon)
!
! local
!
      integer i,j,k,l,m
      integer ilat,ilon,ilev,index
      real(r8) :: increment (nlev,nlat,nlon)
      real(r8) :: row_work  (nlev    )
      real(r8) :: bht_lambda(npt_good)
      real(r8) :: gauss_xy  (npt_good)
      real(r8) :: ht_lambda (nlev    )
      real(r8) :: lon_model,lat_model,gauss_z
      real(r8) :: dl
      real(r8) :: correl_length_xy = 100.0 ! km
      real(r8) :: correl_length_z  =    1.00
!
! reset increment
!
      increment = 0.
!
      do k=1,nlon
        lon_model = longrd(k)
        do j=1,nlat
          lat_model = latgrd(j)
          gauss_xy = 0.
          do i=1,nlev
            bht_lambda = 0.
            LOOP_L : do l=1,npt_good
!
! horizontal correlation
!
              if(i.eq.1) then
                index = ind_good(l,1)
                dl = earthdist(lats_g(index),lons_g(index) &
                              ,lat_model    ,lon_model     )
                if(dl.ge.(2.*correl_length_xy)) then
                  gauss_xy(l) = 0.
                else
                  gauss_xy(l) = gauss(dl,correl_length_xy)
                endif
              endif
              if(gauss_xy(l).eq.0.) cycle LOOP_L
!
              do m=1,nlev
!
! vertical correlation
!
                dl = abs(log(p_model(i,j,k)/p_model_int(l,m)))
                gauss_z = gauss(dl,correl_length_z)
!
                row_work(m) = sqrt(error_model_int(l,m)*error_model(i,j,k)) &
                            * gauss_xy(l) * gauss_z
              end do
!
              ht_lambda(:)  = hz(l,:) * lambda(l)
              bht_lambda(l) = dot_product(row_work,ht_lambda)
!
            end do LOOP_L
            increment(i,j,k) = sum(bht_lambda)
!
          end do
        end do
      end do
!
! return analyzed field
!
      analysis = background + increment
!
      return
end subroutine compute_analysis
!
subroutine good_data_indx_3d(n_g,nlev_obs,obs_g,nodata &
                            ,npt_good,ind_good,map_good)
!
! Returns location of good data in the observations.
! In a single profile, one or more layers can be
! missing, but that does not mean that the whole
! profile is unusable.
!
      use precision
      implicit none
!
! inputs:
!
      integer     , intent(in) :: n_g,nlev_obs
      type(orbits), intent(in) :: obs_g
      real(r8)   , intent(in) :: nodata
!
! outputs:
!
      integer, pointer     :: ind_good(:,:)
      integer, intent(out) :: npt_good
      logical, intent(out) :: map_good(n_g,nlev_obs)
!
! internals:
!
      integer :: i,j
      integer :: counter
      integer, dimension(n_g*nlev_obs,2) :: work
      integer :: ier = 0
!
! reset map of observations
!
      map_good = .false.
!
      counter = 0
      do i=1,n_g
        do j=1,nlev_obs
          if(obs_g%profile(i)%v(j).eq.nodata) cycle
          counter = counter + 1
          work(counter,1) = i
          work(counter,2) = j
          map_good(i,j)   = .true.
        end do   
      end do   
!
! transfer to output arrays
!

      if(associated(ind_good)) deallocate(ind_good, stat=ier )
      call dealloc_err( ier, 'good_data_indx_3d', 'ind_good' )
      allocate(ind_good(counter,2))
      ind_good(:,1) = work(1:counter,1)
      ind_good(:,2) = work(1:counter,2)
      npt_good      = counter
!
      return
end subroutine good_data_indx_3d
!
end module covariances_3d
