c $Id$

c It is assumed that the model calling this interface has been
c compiled so that 8 byte real data are being used.  On non-CRAY
c machines this implies compiling with a "-r8" flag.

c###############################################################################

      subroutine volcemisini()

c Open netCDF file containing data for SO2 emissions from non-eruptive volcanoes.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <volcemis.h>
c-----------------------------------------------------------------------
#include <netcdf.inc>
c-----------------------------------------------------------------------

c Local variables:
      integer
     $  did
     $, nlon
     $, vid
c-----------------------------------------------------------------------

c     Open file.
      call handle_ncerr( 
     $  nf90_open( 'neVolc_emissions.nc', NF_NOWRITE, ncid )
     $, 'volcemisini: error opening file '//'neVolc_emissions.nc' )

c     Check that input data is a right resolution.
      call handle_ncerr( nf90_inq_dimid( ncid, 'lon', did )
     $,                  'volcemisini: ' )
      call handle_ncerr( nf90_inq_dimlen( ncid, did, nlon )
     $,                  'volcemisini: ' )
      if ( nlon .ne. plon ) then
         write(*,*)'volcemisini: change plon to ',nlon
         stop
      end if

      call handle_ncerr( nf90_inq_varid( ncid, 'SO2', vid )
     $,    'volcemisini: cannot find variable '//'SO2' )
      call handle_ncerr( nf90_get_var( ncid, vid, so2 )
     $,    'volcemisini: cannot read data for '//'SO2' )

      call handle_ncerr( nf90_inq_varid( ncid, 'MXVOLHGTIDX', vid )
     $,    'volcemisini: cannot find variable '//'MXVOLHGTIDX' )
      call handle_ncerr( nf90_get_var( ncid, vid, hgtidx )
     $,    'volcemisini: cannot read data for '//'MXVOLHGTIDX' )

      write(*,*)'volcemisini: read data.'

      return
      end

c#######################################################################

      subroutine volcemis( lat, gravit, rpdel, so2tend )

c Calculate SO2 tendency due to emissions from non-eruptive volcanoes.

c N.B. Time constraints dictated that this routine be hard-wired assuming
c      the vertical grid to be the 28 levels from the NCEP reanalysis.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <pmgrid.h>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  lat                  ! model latitude index
      real(r8)
     $  gravit
     $, rpdel(plond,plev)    ! reciprocal of pdel

c Output argument:
      real(r8)
     $  so2tend(plond,plev)  ! tracer tendency due to volcano emissions

c Local variables:
      integer i, k
     $, hgtidx(plond)
      real(r8)
     $  masfac               ! mass scale factor
     $, so2(plond)
     $, emis(plond,plev)     ! emissions distributed in vertical
c-----------------------------------------------------------------------

c     Input emissions are scaled from 1 Tg S/yr to 8 Tg S/yr.
c     The factor of 2 scales from kg S to kg SO2.
      masfac = 8. * 2.

c     Column emissions were divided as follows by Ulrike Lohmann
C   3 EMISSION LEVELS:
C    1. PRE-INTRA ERUPTION IN LEVEL IVOLC-HEIGHT (=TOP OF VOLCANO) 36%
C    2. POST-EXTRA ERUPTION IN LEVEL 15 -16 (CA 550-1736M)         36%
C    3. EXPLOSIVE ERUPTION IN LEVEL 10 - 11 (CA 5000-7900M)        28%

c     For the current model, the volcano height index is contained in the
c     data file, post-extra eruption levels are 21-25, and 
c     explosive eruption levels are 14-16

      call getvolcemis( lat, plond, so2, hgtidx )

      do k = 1, plev
         do i = 1, plon
            emis(i,k) = 0.0
         end do
      end do

c     Test
c      do i = 1, plon
c         emis(i,plev-1) = so2(i) * masfac*.5
c         emis(i,plev) = so2(i) * masfac*.5
c      end do

c     Add emissions at volcano height index.
      do i = 1, plon
         if ( hgtidx(i) .gt. 0 ) then
            emis(i,hgtidx(i)) = so2(i) * masfac * .36
         end if
      end do

c     Add emissions at post-extra eruption levels.  The factor of 5 is for
c     splitting the 36% up among 5 levels.
      do k = 21, 25
         do i = 1, plon
            emis(i,k) = emis(i,k) + so2(i) * masfac * .36 / 5.
         end do
      end do

c     Add emissions at explosive eruption levels.  The factor of 3 is for
c     splitting the 28% up among 3 levels.
      do k = 14, 16
         do i = 1, plon
            emis(i,k) = emis(i,k) + so2(i) * masfac * .28 / 3.
         end do
      end do

c     Scale emissions from kg/m2/s to kg/kg/s
      do k = 1, plev
         do i = 1, plon
            so2tend(i,k) = emis(i,k) * gravit*rpdel(i,k)
         end do
      end do

      return
      end

c#######################################################################

      subroutine getvolcemis( lat, nd1, so2lat, hgtidxlat )

c Return SO2 emission data for the requested latitude.

      use precision

      implicit none
c-----------------------------------------------------------------------
#include <volcemis.h>
c-----------------------------------------------------------------------

c Input arguments:
      integer
     $  lat        ! requested latitude index
     $, nd1        ! 1st dimension of output array

c Output argument:
      real(r8)
     $  so2lat(nd1)   ! SO2 emissions in kg S/m2/s
      integer
     $  hgtidxlat(nd1)! Index for layer that contains the top of 
c                     ! the volcano.

c Local variables:
      integer
     $  i
c-----------------------------------------------------------------------

      do i = 1, plon
         so2lat(i) = so2(i,lat)
         hgtidxlat(i) = hgtidx(i,lat)
      end do

      return
      end

c#######################################################################

