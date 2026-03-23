c $Id$
      module dyninp

CFortran public interface for the dynamics input data module:
C 9/18/96 Brian Eaton / John Caron
C
C 2/9/97 interface to nmc reanalysis grib files
C 3/21/97 remove linqin
C 4/11/97: add level parameter, add regridding capability
c 
c Notes:
c  not yet using dynnm

      use precision

      implicit none
      save
      private
      public ::
     $  inidyn    ! initialize the dynamics input module
     $, chkdyn    ! check that necessary input data is available
     $, interp1   ! Interpolate the fields whose values are required at the
                  ! begining of a timestep.
     $, interp2   ! Interpolate the fields whose values are required at the
                  ! timestep midpoints.
     $, qdyn      ! Query the dynamics module.

      integer :: nlat, nlon, nlev, nlond, nlatd, nlevp, nflds
      integer :: time_avg_offset

c -grb2d  offsets set to avg = 3 hr, fcst = 6 hr
      integer :: grb2d_avg_offset 
      integer :: grb2d_fcst_offset

      logical :: debug
      logical :: want_t, want_oro, want_phis, want_edot, want_omega,
     $	want_arch_vdiff, want_diag_vdiff, want_arch_conv,
     $	want_qrad, want_cldphys,
     $  xlimqin
      real(r8) :: latd(200)

      character(len=168) ::
     $  lpath    ! local path
     $, rpath    ! remote pathname

c##############################################################################
      contains
c##############################################################################

      subroutine inidyn( plon, plat, plev, plond, platd, plevp, 
     $  maxdynflds,
     $  lpth, rpth, date, datesec, offset, lun,
     $	lremove, rdt, rdoro, rdphis, calcedot, rdomega, 
     $  arvdiff, divdiff,
     $  arconvccm, addqrad, arcldphys,
     $  dynnm, dynsc, phi, latwts, hyam, hybm, hyai, hybi, phis, oro)

      implicit none

C input:
	integer plon, plat, plev, plond, platd, plevp, maxdynflds
	character lpth*(*)
	character rpth*(*)
	integer date, datesec, offset
	integer lun
	logical lremove
        logical rdt, rdoro, rdphis, calcedot, rdomega, arvdiff, divdiff
	logical arconvccm, addqrad, arcldphys
	character dynnm(maxdynflds)*32
	real(r8) dynsc(maxdynflds)

C output:
	real(r8) phi(plat)
	real(r8) latwts(plat)
	real(r8) hyam(plev)
	real(r8) hybm(plev)
	real(r8) hyai(plevp)
	real(r8) hybi(plevp)
	real(r8) phis(plond,plat)
	real(r8) oro(plond,plat)
C local
	integer i, ret, grib_open
        integer :: remove=0
c-grb2d  variables to hold offset times
        integer :: datesec_avg
        integer :: datesec_fcst
        integer :: date_avg
        integer :: date_fcst
        real*4 junk

	debug = .true.
	debug = .false.

	nlon = plon
	nlat = plat
	nlev = plev
	nlond = plond
	nlatd = platd
	nlevp = plevp
	nflds = maxdynflds
	time_avg_offset = offset

	want_t = rdt
	want_oro = rdoro
	want_phis = rdphis
	want_edot = .not. calcedot
	want_omega = rdomega
	want_arch_vdiff = arvdiff
	want_diag_vdiff = divdiff
	want_arch_conv = arconvccm
	want_qrad = addqrad
	want_cldphys = arcldphys

        lpath = lpth
        rpath = rpth
        if ( lremove ) remove=1
c-grb2d
        grb2d_avg_offset = 10800    !3 hours
        grb2d_fcst_offset = 21600   !6 hours

        if (time_avg_offset .ne. 0) then
           print *,'dyninp.F:INIDYN: warning: time_avg_offset != 0'
           print *,'can mess up grib file offset calculations'
        endif

	if (debug) then
	    print *, 'INIDYN'
 	    print *, 'dimensions ',plon,plat,plev,plond,platd,plevp, 
     $		maxdynflds
	    print *, 'local path = ', lpth
	    print *, 'remote path = ', rpth
	    print *, 'start = ', date, datesec, offset
	    print *, 'lun = ', lun
	    print *, 'remove = ', lremove
	    print *, 'read = ', 
     $		rdt,rdoro,rdphis,calcedot,rdomega,arvdiff, divdiff,
     $          arconvccm, addqrad, arcldphys
 
	    do i=1, maxdynflds
		print *, '    ', dynnm(i), dynsc(i)
	    enddo
	endif

        call modify_date_datesec(date,datesec, grb2d_fcst_offset,
     $       date_fcst,datesec_fcst)
        call modify_date_datesec(date,datesec, grb2d_avg_offset,
     $       date_avg,datesec_avg)
        if (debug) then
           print *,'inidyn:grb2d:modifying time for fcst from ',
     $          date,datesec,' to ',date_fcst,datesec_fcst
           print *,'inidyn: changed by:',
     $          grb2d_fcst_offset,'seconds'          
           print *,'inidyn:grb2d:modifying time for avg from: ',
     $          date,datesec,' to ',date_avg,datesec_avg
           print *,'inidyn: changed by:',grb2d_avg_offset,
     $          'seconds'
        end if

	ret = grib_open( lpth, rpth, 
     $       date_fcst,  datesec_fcst,
     $       plon, plat, plev,      
     $       remove, phi, latwts, hyam, hybm, hyai, hybi, 0)

	if (1 .ne. ret) then
	    print *, 'grib open failed'
	    call abort
	end if

c save lat in common block
	do i=1,nlat
	   latd(i) = phi(i)
	end do

 	call chkdyn( date, datesec, phis, oro)
      end subroutine inidyn

c##############################################################################

      subroutine chkdyn( date, datesec, phis, oro)
      implicit none

C input:
	integer date
  	integer datesec
C output:
	real(r8) phis(nlond,nlat)
	real(r8) oro(nlond,nlat)
c local
	integer grib_getdata
	integer i,j,ret
	integer from_2d_dataset, level
	real(r8) ice(nlon,nlat)


	integer time_phis, time_oro, time_ice
	data time_phis /-1/
	data time_oro /-1/
	data time_ice /-1/

	real(r8) latrad
c	-62.0 * 0.017453293 radians
	data latrad / -1.099 /

c-grb2d  variables to hold offset times
        integer :: datesec_avg
        integer :: datesec_fcst
        integer :: date_avg
        integer :: date_fcst

        if (debug) then
            print *, 'chkdyn date = ', date, datesec
	endif

       if (time_avg_offset .ne. 0) then
           print *,'dyninp.F:CHKDYN: warning: time_avg_offset != 0'
           print *,'can mess up grib file offset calculations'
        endif

	datesec = datesec - time_avg_offset

        call modify_date_datesec(date,datesec, grb2d_fcst_offset,
     $       date_fcst,datesec_fcst)
        call modify_date_datesec(date,datesec, grb2d_avg_offset,
     $       date_avg,datesec_avg)

      if (debug) then
           print *,'chkdyn:grb2d:modifying time for fcst from ',
     $          date,datesec,' to ',date_fcst,datesec_fcst
           print *,'chkdyn: changed by:',grb2d_fcst_offset,
     $          'seconds'          
           print *,'chkdyn:grb2d:modifying time for avg from: ',
     $          date,datesec,' to ',date_avg,datesec_avg
           print *,'chkdyn: changed by:',grb2d_avg_offset,
     $          'seconds'
        end if
    
c phis     
c Surface geopotential from the most recent time sample read
c   from the dynamics input file.  units: m^2/s^2

	if (want_phis) then
	    from_2d_dataset = 0
	    level = -1
	    ret = grib_getdata("HGT ", from_2d_dataset, level, date, 
     $		datesec, phis, 1, nlond, 0, time_phis)

	    if (ret .eq. 0) then
	      print *, 'ERROR : grib_getdata FAILED on phis'
	      call abort
	    end if

	    if (ret .eq. 1) then
c convert from geopotential = gravity * Geopotential height (HGT)
              do j=1,nlat
                do i=1,nlon
                  phis(i,j) = 9.80665 * phis(i,j)
                enddo
              enddo
	    end if

	endif

c oro
c Orography flag from the most recent time sample read
c   from the dynamics input file.  Has values: 0 over ocean,
c   1 over land, 2 over sea ice.

	if (want_oro) then
          from_2d_dataset = 1
          level = 0

	  ret = grib_getdata("LAND ", from_2d_dataset, level,
     $          date_fcst,datesec_fcst, oro, 1, nlond, 0, time_oro)

          if (ret .eq. 0) then
              print *, 'ERROR : grib_getdata FAILED on oro'
              call abort
          end if

	  if (ret.eq.1) then
             ret = 0
	    from_2d_dataset = 1
            level = 0
	    ret = grib_getdata("ICEC ", from_2d_dataset, level,  
     $           date_fcst,datesec_fcst, ice, 1, nlon, 0, time_ice )

            
            if (ret .eq. 0) then
              print *, 'ERROR : grib_getdata FAILED on ice'
              call abort
            end if


	    do j=1,nlat
	      do i=1,nlon
		if (ice(i,j) .gt. 0.5) then
		    oro(i,j) = 2.0
		else if (oro(i,j) .gt. .5) then
		    oro(i,j) = 1.0
		else
		    oro(i,j) = 0.0 
		end if
	      enddo
	    enddo
	  endif
	endif

      end subroutine chkdyn

c##############################################################################

      subroutine interp1(  date, datesec, 
     $	ps, t, q, zmu, zmd, zeu, hketa, hkbeta, cmfdqr,
     $  conicw, cldfr, concld, cme, nrain, nevapr, cwat)
      implicit none

c
c input:
   	integer date
	integer datesec
c output:
	real(r8) ps(nlond,nlat)
	real(r8) t(nlond,nlev,nlat)
	real(r8) q(nlond,nlev,nlat)
	real(r8) zmu(nlond,nlev,nlat)
	real(r8) zmd(nlond,nlev,nlat)
	real(r8) zeu(nlond,nlev,nlat)
        real(r8) hketa(nlond,nlev,nlat)
	real(r8) hkbeta(nlond,nlev,nlat)
   	real(r8) cmfdqr(nlond,nlev,nlat)
   	real(r8) conicw(nlond,nlev,nlat)
	real(r8) cldfr(nlond,nlevp,nlat)
	real(r8) concld(nlond,nlev,nlat)
	real(r8) cme(nlond,nlev,nlat)
	real(r8) nrain(nlond,nlev,nlat)
	real(r8) nevapr(nlond,nlev,nlat)
	real(r8) cwat(nlond,nlev,nlat)
c local
	integer grib_interp, ret
        integer from_2d_dataset, level
	integer i, j, k, lat

	real(r8) pint(nlond,nlevp), ! Pressure at model interfaces
     $     pmid(nlond,nlev),  ! Pressure at model levels
     $     pdel(nlond,nlev)   ! Layer thickness (pint(k+1) - pint(k))
	real(r8) es(nlond,nlev), qs(nlond,nlev)

        if (time_avg_offset .ne. 0) then
           print *,'dyninp.F:INTERP1: warning: time_avg_offset != 0'
           print *,'can mess up grib file offset calculations'
        endif

	datesec = datesec - time_avg_offset

c ps
c Surface pressure, units: Pa.
	from_2d_dataset = 0
	level = -1
        ret = grib_interp( "PRES ", from_2d_dataset, level, date, 
     $       datesec, ps, 1, nlond, 0)
  
        if (ret .eq. 0) then
          print *, 'ERROR : grib_interp FAILED on ps'
          call abort
        end if

c q         
c Specific humidity, units: kg water vapor/kg moist air.

        from_2d_dataset = 0
        level = -1

	ret = grib_interp( "SPFH ", from_2d_dataset, level, date, 
     $		datesec, q, 1, nlond*nlev, nlond)

        if (ret .eq. 0) then
          print *, 'ERROR : grib_interp FAILED on q'
          call abort
        end if

c prevent negetive Specific humidity
        do j=1,nlat
          do k=1,nlev
            do i=1,nlon
		if (q(i,k,j) .le. 2.e-12) then
		  q(i,k,j) = 2.e-12
		end if
	    end do
	  end do
	end do

c   Temperature, units: K.
	if (want_t) then
            from_2d_dataset = 0
            level = -1
	    ret = grib_interp( "TMP ", from_2d_dataset, level, date, 
     $		datesec, t, 1, nlond*nlev, nlond)
            if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on t'
              call abort
            end if

	endif

	if (want_arch_conv) then
	    print *, 'ERROR : convection not implemented '
	    call abort
	endif

	if (want_cldphys) then
	    print *, 'ERROR : cloud physics not implemented '
	    call abort
         endif
      end subroutine interp1

c##############################################################################

      subroutine interp2( date, datesec,
     $	uv, w, cgs, kvh, ts, taux, tauy, hflx, qflx,
c     $  fsds, fsus, flds, vwcsfc, snowh, qrad, omga)
     $  fsds, fsus, flds, vwcsfc, snowh, qrad)

      implicit none
c
c input:
   	integer date
	integer datesec
c output:
c	real(r8) uv(nlond,nlev,2,nlatd)
	real(r8) uv(nlond,nlev,2,nlat)
	real(r8) w(nlon,nlevp,nlat)
	real(r8) cgs(nlond,nlevp,nlat)
	real(r8) kvh(nlond,nlevp,nlat)
        real(r8) ts(nlond,nlat)
	real(r8) taux(nlond,nlat)
	real(r8) tauy(nlond,nlat)
	real(r8) hflx(nlond,nlat)
	real(r8) qflx(nlond,nlat)
	real(r8) fsds(nlond,nlat)      ! downward solar flux at surface (W/m2)
        real(r8) fsus(nlond,nlat)      ! upward solar flux at surface (W/m2)
        real(r8) flds(nlond,nlat)      ! downward longwave flux at surface (W/m2)
        real(r8) vwcsfc(nlond,nlat)    ! volumetric soil water content at surface (frac)
	real(r8) snowh(nlond,nlat)     ! water equivalent snow depth (m)
	real(r8) qrad(nlond,nlev,nlat)
        real(r8) omga(nlond,nlev,nlat)
c local
	integer i,j 
	integer grib_interp, ret
        integer from_2d_dataset, level
c-grb2d  variables to hold offset times
        integer :: datesec_avg
        integer :: datesec_fcst
        integer :: date_avg
        integer :: date_fcst



        if (time_avg_offset .ne. 0) then
           print *,'dyninp.F:INTERP2: warning: time_avg_offset != 0'
           print *,'can mess up grib file offset calculations'
        endif
	datesec = datesec - time_avg_offset

       call modify_date_datesec(date,datesec, grb2d_fcst_offset,
     $       date_fcst,datesec_fcst)
        call modify_date_datesec(date,datesec, grb2d_avg_offset,
     $       date_avg,datesec_avg)

        if (debug) then
           print *,'interp2:grb2d:modifying time for fcst from ',
     $          date,datesec,' to ',date_fcst,datesec_fcst
           print *,'interp2: changed by:',grb2d_fcst_offset,
     $          'seconds'          
           print *,'interp2:grb2d:modifying time for avg from: ',
     $          date,datesec,' to ',date_avg,datesec_avg
           print *,'interp2: changed by:',grb2d_avg_offset,
     $          'seconds'
        end if

c uv        
c uv(*,*,1,*) is the zonal wind,
c uv(*,*,2,*) is the meridional wind, units: m/s.

	from_2d_dataset = 0
	level = -1
	ret = grib_interp( "UGRD ", from_2d_dataset, level, date, 
     $		datesec, uv, 1, 2*nlond*nlev, nlond)
        if (ret .eq. 0) then
          print *, 'ERROR : grib_interp FAILED on u'
          call abort
        end if

        from_2d_dataset = 0
        level = -1
	ret = grib_interp( "VGRD ", from_2d_dataset, level, date, 
     $		datesec, uv(1,1,2,1), 1, 2*nlond*nlev, nlond)
        if (ret .eq. 0) then
          print *, 'ERROR : grib_interp FAILED on v'
          call abort
        end if

	if (want_edot) then
	    print *, 'ERROR : vertical winds (edot) not implemented '
	    call abort
	endif

	if (want_arch_vdiff) then
	    print *, 'ERROR : vert diffusion not implemented '
	    call abort
	endif

	if (want_diag_vdiff) then

c    surface temperature, units = K
c    note: level 0 = surface (kpds6=1)
          from_2d_dataset = 1
          level = 0
	  ret = grib_interp( "TMP ", from_2d_dataset, level,
     $         date_fcst,datesec_fcst, 
     $         ts, 1, nlond, 0)

          if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on surface temp'
              call abort
          end if

c    zonal surface stress, N/m^2
	  from_2d_dataset = 1
          level = 0
	  ret = grib_interp( "UFLX ", from_2d_dataset, level, 
     $         date_avg, datesec_avg,
     $         taux, 1, nlond, 0)
          if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on taux'
              call abort
          end if

c    meridional surface stress, N/m^2
          from_2d_dataset = 1
          level = 0
	  ret = grib_interp( "VFLX ", from_2d_dataset, level, 
     $         date_avg, datesec_avg,
     $         tauy, 1, nlond, 0)
          if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on tauy'
              call abort
          end if

c    surface sensible heat flux (positive up) units: w/m^2
          from_2d_dataset = 1
          level = 0
	  ret = grib_interp( "SHTFL ", from_2d_dataset, level, 
     $         date_avg, datesec_avg,
     $         hflx, 1, nlond, 0)
          if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on hflx'
              call abort
          end if

c    surface water vapor flux (positive up) kg/m^2/s
          from_2d_dataset = 1
          level = 0
          ret = grib_interp( "LHTFL ", from_2d_dataset, level, 
     $         date_avg, datesec_avg,
     $         qflx, 1, nlond, 0)
          if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on qflx'
              call abort
          end if

c    convert w/m^2 / J/kg == kg/m^2/s
          do j=1,nlat
              do i=1,nlon
                qflx(i,j) = qflx(i,j) / 2.31e6
              enddo
          enddo

c    downward solar flux at surface (W/m2)
          from_2d_dataset = 1
          level = 1
          ret = grib_interp( "DSWRF ", from_2d_dataset, level, 
     $         date_avg, datesec_avg,
     $         fsds, 1, nlond, 0)
           if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on fsds'
              call abort
          end if

c    upward solar flux at surface (W/m2)
          from_2d_dataset = 1
          level = 1
          ret = grib_interp( "USWRF ", from_2d_dataset, level,
     $         date_avg, datesec_avg, 
     $         fsus, 1, nlond, 0)
          if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on fsus'
              call abort
          end if

c    downward longwave flux at surface (W/m2)
          from_2d_dataset = 1
          level = 0
          ret = grib_interp( "DLWRF ", from_2d_dataset, level, 
     $         date_avg, datesec_avg, 
     $         flds, 1, nlond, 0)
          if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on flds'
              call abort
          end if

c    volumetric soil water content at surface (frac)
          from_2d_dataset = 1
          level = 1
          ret = grib_interp( "SOILW ", from_2d_dataset, level, 
     $         date_fcst, datesec_fcst,
     $         vwcsfc, 1, nlond, 0)
          if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on vwcsfc'
              call abort
          end if
c         Set missing values to 0.0
          do j=1,nlat
              do i=1,nlon
                vwcsfc(i,j) = max( vwcsfc(i,j), 0.0_r8 )
              enddo
          enddo

c    water equivalent accumulated snow depth (m)
          from_2d_dataset = 1
          level = 0
          ret = grib_interp( "WEASD ", from_2d_dataset, level, 
     $         date_fcst, datesec_fcst,
     $         snowh, 1, nlond, 0)
          if (ret .eq. 0) then
              print *, 'ERROR : grib_interp FAILED on snowh'
              call abort
          end if
c         Convert kg/m2 (=mm) --> m
          do j=1,nlat
              do i=1,nlon
                snowh(i,j) = snowh(i,j) * .001
              enddo
          enddo

	endif

	if (want_omega) then
	    print *, 'ERROR : vertical winds (omega) not implemented '
	    call abort
	endif

	if (want_qrad) then
	    print *, 'ERROR : qrad not implemented '
	    call abort
	endif

c	call grib_time
      end subroutine interp2

c##############################################################################

      subroutine qdyn( lpth, rpth, ntim)
      implicit none
c
c output:
	character lpth*(*)
	character rpth*(*)
	integer ntim

        lpth = lpath
        rpth = rpath
        ntim = 0

      return
      end subroutine qdyn

c##############################################################################

      end module dyninp
