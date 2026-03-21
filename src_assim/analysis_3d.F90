module analysis_3d
 use floats
 use data_structures
 use obs_operator
 use f2v
 use mathlib
 use covariances_3d
 use error_messages, only: dealloc_err
 use precision

 private 
 public :: orbits, std, analysis_xyz_3d

contains

subroutine analysis_xyz_3d(nlev_, nlat_, nlon_, &
                           obs, zgrid3_, dzgrid3_, latgrid_, longrid_, &
                           field_, sfield_, analysis_, obs_g)
!
! in this version the analysis is done in 3d, not split between
! the horizontal and the vertical as in the original version.
!
! JFL Sep 98
! 
     implicit none
!
! inputs:
!
     integer      , intent(in) :: nlev_, nlat_, nlon_
     type (orbits), intent(in) :: obs
     real(r8)    , intent(in) :: zgrid3_(nlev_,nlat_,nlon_)
     real(r8)    , intent(in) :: latgrid_(nlat_),longrid_(nlon_)
     real(r8)    , intent(in) :: dzgrid3_(nlev_,nlat_,nlon_)
     real(r8)    , intent(in) :: field_(nlev_,nlat_,nlon_),sfield_(nlev_,nlat_,nlon_)
!
! outputs:
!
     real(r8), intent(out) :: analysis_(nlev_,nlat_,nlon_) 
     type(orbits)           :: obs_g ! gridded observations (WDC: obs_g = obs!)
!
! internals:
!

!-drb2     real(r8), dimension(:), pointer :: lats_g, lons_g
     real(r8), allocatable  :: lats_g(:), lons_g(:)
     integer,   dimension(:,:), pointer :: ind_good
     logical, dimension(obs%n,size(obs_g%profile_grid)) :: map_good


     real(r8), allocatable :: lambda(:)

     integer  :: nhxy(3,obs%n)
     real(r8):: hxy(3,obs%n)
     real(r8), allocatable :: o_xyz(:,:)
     real(r8), allocatable :: field_int(:),sfield_int(:,:),p_int(:,:) &
                              ,hz(:,:)
!
     integer :: i, j, k, n_g, npt_good, index
     integer :: nlat, nlon, nlev, ni, nj, nk
     integer :: nlev_obs
     real(r8), allocatable :: latgrid(:), longrid(:)
     real(r8), allocatable :: zgrid3(:,:,:), dzgrid3(:,:,:)
     real(r8), allocatable :: field(:,:,:), sfield(:,:,:), analysis(:,:,:) 
     integer  , allocatable :: i_indx(:), j_indx(:)

     integer :: ii_indx(nlon_)
     real(r8)              :: latoff, lonoff, latmin, latmax, lonmin, lonmax
     real(r8)              :: minimum_variance_model = 0.

     integer :: alloc_error
     integer :: error = 0

!
! select only model grid points located close enough to the observations.
!
     latoff = 30.0
     lonoff = 30.0
!
     latmin =  minval( obs%profile%lat ) - latoff
     latmax =  maxval( obs%profile%lat ) + latoff
!
     print *, 'analysis_xyz% latmin, latmax: ', latmin, latmax
!
! select latitude indeces
!
     nj = 0
     do j=1,nlat_
        if ( (latgrid_(j) .le. latmax) .and. (latgrid_(j) .ge. latmin) ) nj = nj + 1
     enddo
     allocate ( j_indx(nj) )
     nj = 0
     do j=1,nlat_
        if ( (latgrid_(j) .le. latmax) .and. (latgrid_(j) .ge. latmin) ) then
           nj = nj + 1
           j_indx(nj) = j
        endif   
     enddo
!
! select longitude indeces
!
     ii_indx = 0
     ni = 0
     do i=1,nlon_
       do k=1,obs%n 
          lonmin= obs%profile(k)%lon - lonoff
          lonmax= obs%profile(k)%lon + lonoff
          if (lonmin .lt. 0  ) lonmin = 360 + lonmin
          if (lonmax .gt. 360) lonmax = lonmax - 360
          if ( insidelons(lonmin, lonmax, longrid_(i)) ) then
             ni = ni + 1
             ii_indx(ni) = i
             exit
          endif
       enddo
     enddo
!
     allocate(i_indx(ni))
     i_indx = ii_indx(1:ni)
!
     nk = nlev_
     nlon = ni
     nlat = nj
     nlev = nk
!    
     allocate(longrid           (nlon))
     allocate(latgrid      (nlat)     )
     allocate(zgrid3  (nlev,nlat,nlon))
     allocate(dzgrid3 (nlev,nlat,nlon))
     allocate(field   (nlev,nlat,nlon))
     allocate(sfield  (nlev,nlat,nlon))
     allocate(analysis(nlev,nlat,nlon))
!
! remap input data to reduced area
!
     longrid = longrid_( i_indx )
     latgrid = latgrid_(         j_indx )
     zgrid3  =  zgrid3_( 1:nlev, j_indx, i_indx )
     dzgrid3 = dzgrid3_( 1:nlev, j_indx, i_indx )
     field   = field_  ( 1:nlev, j_indx, i_indx )
     sfield  = sfield_ ( 1:nlev, j_indx, i_indx )
!
! put the variance in sfield
!
     do k=1,nlon
       do j=1,nlat
         do i=1,nlev
           sfield(i,j,k) = (sfield(i,j,k)*field(i,j,k))**2 &
                         + minimum_variance_model**2
         end do
       end do
     end do
!
! WDC -- replace logic for reducing resolution of observations
!        with identity operator -- the obs. are already on the model grid 
!
     n_g = obs%n
!-drb2     lats_g=>obs%profile%lat
!-drb2     lons_g=>obs%profile%lon
     allocate(lats_g(n_g), stat=error )
     call alloc_err( error, 'analysis_xyz_3d','lats_g',n_g)
     allocate(lons_g(n_g), stat=error )
     call alloc_err( error, 'analysis_xyz_3d','lons_g',n_g)
     
     lats_g(:) = obs%profile(:)%lat
     lons_g(:) = obs%profile(:)%lon

!
! prepare data structure (obs_g) to store the observations:
!
     obs_g%n = n_g
     obs_g%nodata = obs%nodata
     obs_g%profile_grid = obs%profile_grid
! CHANGE
!    if ( associated(obs_g%profile) ) deallocate(obs_g%profile, stat=error)
!    call dealloc_err( error, 'analysis_xyz_3d', 'obs_g%profile' )
     allocate( obs_g%profile(n_g) )
!-drb above logic is crap. quick change please

     obs_g%profile(:)%lat = lats_g
     obs_g%profile(:)%lon = lons_g
!-drb     obs_g%profile(:)%lat = obs%profile%lat
!-drb     obs_g%profile(:)%lon = obs%profile%lon

     obs_g%profile(:)%date   = obs%profile%date
     obs_g%profile(:)%secs   = obs%profile%secs
!++bug     do i=1,nlev
     do i=1,size(obs_g%profile_grid)
        obs_g%profile(:)%v(i)   = obs%profile(:)%v(i)
        obs_g%profile(:)%sigmav(i) = obs%profile(:)%sigmav(i)
     end do
!
! find the location of the good data
!
     nlev_obs = size(obs_g%profile_grid)


     call good_data_indx_3d(n_g,nlev_obs,obs_g,obs%nodata &
                           ,npt_good,ind_good,map_good)
!
! if no good observation is available, then return original field
!
     if(npt_good.eq.0) then
       print *,'no good data was available'
       print *,'return original field as analysis'
       analysis_ = field_

       ! deallocate these in case of returning early
       deallocate(lats_g, stat=error)
       call dealloc_err( error, 'analysis_xyz_3d','lats_g' )
       deallocate(lons_g, stat=error)
       call dealloc_err( error, 'analysis_xyz_3d','lons_g' )
       return
     endif
!
! compute observation covariance in 3d
!
     if(allocated(o_xyz)) deallocate(o_xyz, stat=error)
     call dealloc_err( error, 'analysis_xyz_3d', 'o_xyz' )
     allocate(o_xyz(npt_good,npt_good))
     call obs_covariance_xyz(n_g,nlev_obs,obs_g,lats_g,lons_g &
                            ,npt_good,ind_good,o_xyz)
!
! "interpolate" field and its variance to the location
! and vertical grid of the observations
!
! horizontal "interpolation" operator
!
     call obs_operator_xy(obs_g%profile(:)%lat &
                         ,obs_g%profile(:)%lon &
                         ,latgrid,longrid,     &
                          hxy,nhxy) 

     allocate( field_int(npt_good     ))
     allocate(sfield_int(npt_good,nlev))
     allocate(     p_int(npt_good,nlev))
     allocate(        hz(npt_good,nlev))
     call interp_3d(nlon,nlat,nlev,n_g,nlev_obs &
                   ,npt_good,ind_good,map_good  &
                   ,field,sfield,zgrid3,dzgrid3 &
                   ,hxy,nhxy                    &
                   ,obs_g%profile_grid(:)       &
                   ,field_int,sfield_int        &
                   ,p_int,hz)
!

!
! compute HBHt, point by point and store in o_xyz
!
     call compute_hbht(n_g,nlev_obs,nlat,nlon,nlev &
                      ,npt_good,ind_good           &
                      ,lats_g,lons_g,p_int         &
                      ,sfield_int,hz,o_xyz)
!
!    write(110) npt_good,nlev
!    write(110) o_xyz
!    write(110) field_int
!    write(110) sfield_int
!    write(110) lats_g,lons_g
!    write(110) ind_good
!    write(110) hz
!    write(110) p_int
!
! compute R.H.S. of PSAS equation (Eq(11) in Menard)
!
     do index=1,npt_good
       i = ind_good(index,1)
       j = ind_good(index,2)
       field_int(index) = obs_g%profile(i)%v(j) &
                        - field_int(index)
     end do
!
     allocate(lambda(size(field_int)))
!
! Solve linear system Ax=b (use LU-decomposition).
!
!                    A       b       x
!                    |       |       |
     call matsolve(o_xyz,field_int,lambda)
!
!    print *,'done with matsolve'
!    write(110) field_int
!    write(110) lambda
!
! free up some memory
!
     deallocate(o_xyz, stat=error )
     call dealloc_err( error, 'analysis_xyz_3d', 'o_xyz' )
     deallocate(field_int, stat=error)
     call dealloc_err( error, 'analysis_xyz_3d', 'field_int' )
!
! compute the analysis
!
     call compute_analysis(n_g,nlev_obs,nlev,nlat,nlon    &
                          ,npt_good,ind_good              &
                          ,lats_g,lons_g                  &
                          ,latgrid,longrid,sfield_int     &
                          ,field,sfield,zgrid3,p_int      &
                          ,hz,lambda,analysis)
!
! return analysis array for the whole domain
!
     analysis_                       = field_
     analysis_(1:nlev,j_indx,i_indx) = analysis
!
! free up memory
!
     deallocate(i_indx, stat=error    )
     call dealloc_err( error, 'analysis_xyz_3d', 'i_indx' )
     deallocate(j_indx, stat=error    )
     call dealloc_err( error, 'analysis_xyz_3d', 'j_indx' )
     deallocate(analysis, stat=error  )
     call dealloc_err( error, 'analysis_xyz_3d', 'analysis' )
     deallocate(lambda, stat=error    )
     call dealloc_err( error, 'analysis_xyz_3d', 'lambda' )

     deallocate(zgrid3, stat=error    )
     call dealloc_err( error, 'analysis_xyz_3d', 'zgrid3' )
     deallocate(dzgrid3, stat=error   )
     call dealloc_err( error, 'analysis_xyz_3d', 'dzgrid3' )
     deallocate(field, stat=error     )
     call dealloc_err( error, 'analysis_xyz_3d', 'field' )
     deallocate(sfield, stat=error    )
     call dealloc_err( error, 'analysis_xyz_3d', 'sfield' )
     deallocate(sfield_int, stat=error)
     call dealloc_err( error, 'analysis_xyz_3d', 'sfield_int' )
     deallocate(hz, stat=error        )
     call dealloc_err( error, 'analysis_xyz_3d', 'hz' )
     deallocate(p_int, stat=error     )
     call dealloc_err( error, 'analysis_xyz_3d', 'p_int' )


     deallocate(ind_good, stat=error  )
     call dealloc_err( error, 'analysis_xyz_3d', 'ind_good' )

     deallocate(lats_g, stat=error)
     call dealloc_err( error, 'analysis_xyz_3d','lats_g' )
     deallocate(lons_g, stat=error)
     call dealloc_err( error, 'analysis_xyz_3d','lons_g' )

!
     return
!
end subroutine analysis_xyz_3d

end module analysis_3d

