!
      module gaussian_weight
!
      use floats
      use precision
!
      implicit none
!
      real(r8), allocatable, dimension(:) :: gw_xy
      logical :: done_gaussian_weight = .false.
!
      contains
!
      subroutine gaussian_weight_xy(nx,latgrd,longrd,corrlength_xy)
!
      use geo
      use f2v
!
! compute gaussian weights for the calculation of the
! model covariance matrix (see model_covariance_xy)
!
! JFL : Aug 98
!
      use precision
      implicit none
!
! arguments
!
      integer  ,intent(in)  :: nx
      real(r8),intent(in)  :: corrlength_xy
      real(r8),intent(in)  :: latgrd(:),longrd(:)
!
! local
!
      integer   :: nlat,nlon,nw
      integer   :: i,j,index,ilat1,ilat2,ilon1,ilon2
      real(r8) :: dr
!
!     integer   :: npt_0 = 0
!
      nlat = size(latgrd)
      nlon = size(longrd)
!
! allocate memory
!
      nw = (nx*nx-nx)/2
!
! initialize gw_xy
!
      gw_xy = 0.
!
      do i=1,nx
        call indexv2f(nlat,nlon,i,ilat1,ilon1)
        do j=i+1,nx
          call indexv2f(nlat,nlon,j,ilat2,ilon2)
          dr = earthdist( latgrd(ilat1), longrd(ilon1), &
                          latgrd(ilat2), longrd(ilon2) )
          index = (i-1)*nx - (i*i+i+2)/2 + (j-1) + 2
          if(dr.lt.2*corrlength_xy) then
            gw_xy(index) = gauss(dr,corrlength_xy)
!         else
!           npt_0 = npt_0 + 1
          endif
        end do
      end do
!
!     print *,'number of elements with value 0 = ',npt_0,nx
!     stop
!
      done_gaussian_weight = .true.
!
      return
      end subroutine gaussian_weight_xy
!
! external function
!
      
      function gauss(dr,l)
!
      use precision
      implicit none
!
      real(r8) gauss
      real(r8), intent(in) :: dr, l
!
      real(r8) z
!
!     gauss = exp( - (dr/l)**2 )
!
! JFL Aug-98
!
! replace Gaussian by compactly-supported piecewise rational function
! (from Gaspari and Cohn 1998)
!
! sqrt (0.3) comes from Gaspari and Cohn
! sqrt (2.0) comes from the fact that the function approximates exp(-z^2/(2L^2))
! not exp(-z^2/L^2)
!
      z = abs(dr*sqrt(0.3*2.)/l)
!
      gauss = 0.
!
      if(z.le.1.) then
        gauss = (((((-0.25*z+0.5)*z)+0.625)*z)-5./3.) * z**2 &
              +  1.
      elseif(z.le.2.) then
        gauss = (((((((1./12.*z-0.5)*z)+0.625)*z)+5./3.)*z)-5.)*z &
              + 4. - 2./3.  / z
      endif
!
! make sure it is a positive value since the polynomial expression
! might become very slightly negative
!
      gauss = max(gauss,0._r8)
!
      return
      end function gauss
!
      end module gaussian_weight
