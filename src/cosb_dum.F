c $Id$
      integer function crayopen(i1,i2,c1)
      integer i1,i2
      character*8 c1
      crayopen=0
      write(*,*)'crayopen: this routine should not have been called'
      stop
      end

      integer function crayread(i1,r1,i2,i3)
      use precision
      integer i1,i2,i3
      real(r8) r1
      crayread=0
      write(*,*)'crayread: this routine should not have been called'
      stop
      end

      integer function crayrew(i1)
      integer i1
      crayrew=0
      write(*,*)'crayrew: this routine should not have been called'
      stop
      end

      subroutine ctodpf(r1,r2,i1)
      use precision
      real(r8) r1,r2
      integer i1
      write(*,*)'ctodpf: this routine should not have been called'
      stop
      end

      integer function crayclose(i1)
      integer i1
      crayclose=0
      write(*,*)'crayclose: this routine should not have been called'
      stop
      end
