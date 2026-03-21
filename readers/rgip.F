c $Id$
! try to vectorize over z
      subroutine areaavg( nz, src, src_s, rst, rst_s, frac, idx)
      use precision

        implicit none

! input
        integer nz, src_s, rst_s
        real(r8) src(src_s, nz), frac(6)
        integer idx(6)
! output
        real(r8) rst(rst_s, nz)
! local
        integer z, i

        do z = 1, nz

            rst(1,z) = src(idx(1),z) * frac(1) +   
     $          src(idx(2),z) * frac(2) +  
     $          src(idx(3),z) * frac(3) +  
     $          src(idx(4),z) * frac(4) +  
     $          src(idx(5),z) * frac(5) +  
     $          src(idx(6),z) * frac(6)

        end do

        return
        end

    ! vectorize the time interpolation
      subroutine tinterp( dt, before, after, result, npts)

      use precision

      implicit none

    ! input
        integer npts
        real(r8) dt, before(npts), after(npts)
    ! output
        real(r8) result(npts)
    ! local
	integer i
        real(r8) dtm1

        dtm1 = 1 - dt

        do i=1,npts
            result(i) = before(i) * dt + after(i) * dtm1
        end do

        return
        end


