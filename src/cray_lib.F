c $Id$

      subroutine gather(n, a, b, index)
      use precision
      implicit real(r8) (a-h,o-z)
      dimension a(*), index(*), b(*)
      do 100 i=1,n
         a(i)=b(index(i))
  100 continue
      return
      end

c##############################################################################

      integer function ismax(n,sx,incx)
      use precision
      implicit real(r8) (a-h,o-z)
      real(r8) sx
      dimension sx(*)
      real(r8) max
      integer imax
      max = sx(1)
      imax = 1
      do 100 i=1+incx,n,incx
         if (max .lt. sx(i)) then
            max = sx(i)
            imax = i
         endif
  100 continue
      ismax = imax
      return
      end

c##############################################################################

      integer function ismin(n,sx,incx)
      use precision
c-drb      implicit real(r8) (a-h,o-z)
      implicit none

c Input
      integer n,incx
      real(r8) sx(n)

c Local
      real(r8) min
      integer imin,i
      min = sx(1)
      imin = 1
      do 100 i=1+incx,n,incx
         if (min .gt. sx(i)) then
            min = sx(i)
            imin = i
         endif
  100 continue
      ismin = imin
      return
      end

c##############################################################################

      function sasum( n, x, inc )
      use precision
      implicit real(r8) (a-h,o-z)
      integer n, inc, i
      real(r8) x(n)
      sasum = 0.
      if ( inc .eq. 1 ) then
         do i = 1, n
            sasum = sasum + abs( x(i) )
         end do
      else
         do i = 1, n
            sasum = sasum + abs( x(1+(i-1)*inc) )
         end do
      end if
      return
      end

c##############################################################################

      subroutine scatter(n, a, index, b)
      use precision
      implicit real(r8) (a-h,o-z)
      dimension a(*), index(*), b(*)
      do 100 i=1,n
         a(index(i))=b(i)
  100 continue
      return
      end

c##############################################################################

      function ssum( n, x, inc )
      use precision
      implicit real(r8) (a-h,o-z)
      integer n, inc, i
      real(r8) x(n)
      ssum = 0.
      if ( inc .eq. 1 ) then
         do i = 1, n
            ssum = ssum + x(i)
         end do
      else
         do i = 1, n
            ssum = ssum + x(1+(i-1)*inc)
         end do
      end if
      return
      end

c##############################################################################

      subroutine wheneq( n, array, inc, target, index, nval )

c Finds all array elements equal to the target.  The Cray library version
c says the array can be integer or real(r8).  The integer version is implemented
c here.  The algorithm is from the "Programmer's Library Reference Manual"
c SR-0113 D-01.

      use precision
      implicit none

c Input:

      integer
     $   n         ! Number of elements to be searched
     $,  array(n)  ! Array to be searched
     $,  inc       ! Increment between elements of the searched array
     $,  target    ! Value searched for in the array

c Output:

      integer
     $   index(n)  ! Array containing indices of array elements matching target
     $,  nval      ! Number of values put in index array

c Local:

      integer ina, i

      ina = 1
      nval = 0

      if ( inc .lt. 0 ) ina = -inc*(n-1) + 1

      do i = 1, n
         if ( array(ina) .eq. target ) then
            nval = nval + 1
            index(nval) = ina
         endif
         ina = ina + inc
      enddo

      return
      end

c##############################################################################

      subroutine whenflt( n, array, inc, target, index, nval )

c Finds all array elements less than the target.

      use precision
      implicit none

c Input:

      integer
     $   n         ! Number of elements to be searched
     $,  inc       ! Increment between elements of the searched array
      real(r8)
     $   array(n)  ! Array to be searched
     $,  target    ! Value searched for in the array

c Output:

      integer
     $   index(n)  ! Array containing indices of array elements matching target
     $,  nval      ! Number of values put in index array

c Local:

      integer ina, i

      ina = 1
      nval = 0

      if ( inc .lt. 0 ) ina = -inc*(n-1) + 1

      do i = 1, n
         if ( array(ina) .lt. target ) then
            nval = nval + 1
            index(nval) = ina
         endif
         ina = ina + inc
      enddo

      return
      end

c##############################################################################


      subroutine whenfgt(n,array,inc,target,index,nval)
      use precision
      implicit real(r8) (a-h,o-z)
      dimension index(*), array(*)
      ina=1
      nval=0
      if(inc .lt. 0) ina=(-inc)*(n-1)+1
      do 100 i=1,n
         if(array(ina) .gt. target) then
	    nval=nval+1
	    index(nval)=i
         end if
         ina=ina+inc
  100 continue
      return
      end

c##############################################################################
