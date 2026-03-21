c $Id$
c It is assumed in these routines that real(r8)s are 8-bytes.  On non-Cray
c machines this typically means compiling with a "-r8" flag.

      subroutine inimass2nc( etamid, tracnam )

c Open and initialize netCDF file for tracer mass quantities.

c N.B.  This routine must be called after the history module is initialized
c       since it queries that module for date information.

      use histout, only: qhsti
      use precision
      use netcdf

      implicit none
c---------------------------------------------------------------------
#include <pmgrid.h>
c---------------------------------------------------------------------
#include <netcdf.inc>
c---------------------------------------------------------------------
#include <mass2nc.com>
c---------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  etamid(plev)       ! hybrid coordinate levels
      character
     $  tracnam(pcnst)*8   ! specie names

c Local variables:
      integer
     $  k, m
     $, did(4), vid(7), tmp(2)
      integer
     $  ret
     $, nbdate, nbsec  ! base date and seconds for time units string
      real(r8)
     $  lev(plev)

c Externals
      character*10 idate2char
      character*8 isec2char
c-----------------------------------------------------------------------

      write (6,*) ' mass2nc disable till r8 issues addressed '
      return

#ifdef TURNOFF
c     Initial netCDF file record number (time sample number).
      nrec = 0

      do k = 1, plev
         lev(k) = 1000.*etamid(k)
      end do

      call handle_err( NF90_CREATE( 'gam.nc', NF90_CLOBBER, ncid ) )

c     Define dimensions.
      call handle_err( NF90_DEF_DIM( ncid, 'lev', plev, did(1) ) )
      call handle_err( NF90_DEF_DIM( ncid, 'trac', pcnst, did(2) ) )
      call handle_err( NF90_DEF_DIM( ncid, 'time', NF90_UNLIMITED, did(3) ))
      call handle_err( NF90_DEF_DIM( ncid, 'strlen', 8, did(4) ) )


c     Define variables.
      m = 1
      call handle_err( NF90_DEF_VAR( ncid, 'lev', NF90_FLOAT, did(1)
     $,                            vid(m) ) )
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'long_name'
     $,                'hybrid level at layer midpoints (1000*(A+B))' ))
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'units'
     $,                'hybrid_sigma_pressure' ) )
      m = 2
      tmp(1) = did(4)
      tmp(2) = did(2)
      call handle_err( NF90_DEF_VAR( ncid, 'tracnam', NF90_CHAR, tmp
     $,                            vid(m) ) )
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'long_name',
     $,                'tracer name' ))
      m = 3
      ret = qhsti( 'NBDATE', nbdate )
      ret = qhsti( 'NBSEC', nbsec )
      call handle_err( NF90_DEF_VAR( ncid, 'time', NF90_DOUBLE, did(3)
     $,                            vid(m) ) )
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'long_name', 4
     $,                'time' ))
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'units'
     $,   'days since '//idate2char(nbdate)//' '//isec2char(nbsec) ) )
      m = 4
      call handle_err( NF90_DEF_VAR( ncid, 'date', NF90_INT,  did(3)
     $,                            vid(m) ) )
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'long_name',
     $,                'date in YYYYMMDD format' ))
      m = 5
      call handle_err( NF90_DEF_VAR( ncid, 'datesec', NF90_INT, did(3)
     $,                            vid(m) ) )
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'long_name'
     $,                'seconds past date at 0Z' ))
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'units'
     $,                's' ) )
      m = 6
      call handle_err( NF90_DEF_VAR( ncid, 'gam', NF90_DOUBLE, did
     $,                            vid(m) ) )
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'long_name'
     $,                'global average mass per level' ))
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'units'
     $,                'kg/m2' ) )
      m = 7
      call handle_err( NF90_DEF_VAR( ncid, 'gammr', NF90_DOUBLE, did
     $,                            vid(m) ) )
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'long_name'
     $,                'global average mass mixing ratio per level' ))
      call handle_err( NF90_PUT_ATT( ncid, vid(m), 'units',
     $,                'kg/kg' ) )

c     Leave definition mode.
      call handle_err( NF90_ENDDEF( ncid ) )

c     Write non-record data to file.
      call handle_err( NF90_PUT_VAR_DOUBLE( ncid, vid(1), lev ) )
      call handle_err( NF90_PUT_VAR( ncid, vid(2), tracnam ) )
#endif

      return
      end

c#############################################################################

      subroutine wrmass2nc( trmk, mixk )

c Write record data to netCDF file.

      use precision
      use netcdf

      implicit none
c---------------------------------------------------------------------
#include <pmgrid.h>
c---------------------------------------------------------------------
#include <netcdf.inc>
c---------------------------------------------------------------------
#include <mass2nc.com>
c---------------------------------------------------------------------

c Input arguments:
      real(r8)
     $  trmk(plev,pcnst)          ! global average mass (kg/m2)
     $, mixk(plev,pcnst)          ! global average m.r. (kg/kg)

c Local variables:
      integer
     $  ret, ndcur, nscur, ncdate, ncsec
     $, vid(2), start(3), count(3)
      data
     $  start/ 1, 1, 1 /, count/ plev, pcnst, 1 /
      save start, count

      real(r8) time

c Externals:
      integer qctli
c---------------------------------------------------------------------

c     Query control module for current time/date info.
      ret = qctli( 'CURDAY', ndcur )
      ret = qctli( 'CURDAYSEC', nscur )
      ret = qctli( 'CURDATE', ncdate )
      ret = qctli( 'CURDATESEC', ncsec )

c-drb 23 Nov 2001 changed from time = ndcur + ncsec/86400. 
      time = ndcur + nscur/86400.

      nrec = nrec + 1
      start(3) = nrec

#ifdef TURNEDOFF
      call handle_err( NF90_INQ_VARID( ncid, 'time', vid ) )
      call handle_err( NF90_PUT_VAR( ncid, vid, nrec, time ) )
      call handle_err( NF90_INQ_VARID( ncid, 'date', vid ) )
      call handle_err( NF90_PUT_VAR( ncid, vid, nrec, ncdate ) )
      call handle_err( NF90_INQ_VARID( ncid, 'datesec', vid ) )
      call handle_err( NF90_PUT_VAR( ncid, vid, nrec, ncsec ) )
      call handle_err( NF90_INQ_VARID( ncid, 'gam', vid(1) ) )
      call handle_err( NF90_INQ_VARID( ncid, 'gammr', vid(2) ) )
      call handle_err( NF90_PUT_VAR( ncid, vid(1), start, count
     $,                                    trmk ) )
      call handle_err( NF90_PUT_VAR( ncid, vid(2), start, count
     $,                                    mixk ) )
#endif

      return
      end

c#############################################################################

      subroutine clmass2nc( )

c Close the netCDF file.

      use precision
      use netcdf

      implicit none
c---------------------------------------------------------------------
#include <pmgrid.h>
c---------------------------------------------------------------------
#include <netcdf.inc>
c---------------------------------------------------------------------
#include <mass2nc.com>
c---------------------------------------------------------------------

      call handle_err( NF90_CLOSE( ncid ) )
      return
      end
