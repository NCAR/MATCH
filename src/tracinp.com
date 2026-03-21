c $Id$
c Reader for tracer initial conditions from a CCM format history file.

      integer, parameter ::
     $  mxfldsIn = 200    ! max number of fields in input history files

      character
     $  lpath*168    ! local path for tracer initial conditions file
     $, rpath*168    ! remote path for tracer initial conditions file
     $, mcflds(2,mxfldsIn)*8 ! field info from character header records

      common / tracinp_C / 
     $  lpath, rpath, mcflds

      logical
     $  cosbtra     ! .true. => file is COS blocked
     $, readNAS       ! t => read ICs for non-advected species
     $, readQIC     ! t => read initial conditions for the advected Q field
c                   !      rather than initializing from meteorology.

      common / tracinp_L /
     $  cosbtra, readNAS, readQIC

      integer
     $  luntra     ! Fortran unit number attached to file
     $, lenhdi     ! length of integer header record (number of integers)
     $, lenhdc     ! length of character header record (number of character*8)
     $, lenhdr     ! length of real header record (number of reals)
     $, maxsiz     ! length of data record (number of reals)
     $, nfld       ! number of fields in data records
     $, nlon       ! number of longitudes
     $, nlonw      ! number of longitudes written to history tape
     $, nlat       ! number of latitudes
     $, nlev       ! number of model layers
     $, tracdate   ! current date (in yymmdd format)
     $, tracsec    ! seconds relative to tracdate
     $, ndens      ! packing density for data records
     $, mflds(3,mxfldsIn) ! field info from integer header records

      common / tracinp_I /
     $  luntra, lenhdi, lenhdc, lenhdr, maxsiz, nfld, nlon, nlonw, nlat
     $, nlev, tracdate, tracsec, ndens, mflds
