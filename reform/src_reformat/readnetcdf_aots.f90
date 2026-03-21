	MODULE READNETCDF_AOTS
	use netcdf
	implicit none 
	integer ios
	integer ,PARAMETER :: nlon=192 ,nlat=94,ntime=1
        integer,PARAMETER ::  NREC=nlon*nlat,nf=4,np=11
	real, PARAMETER :: omit=-9999. 
	INTEGER*4 RECDIM,RCODE,START(10),COUNT(10),VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
	include 'netcdf.inc'

	character*2 collection
	CHARACTER*31 DUMMY,cname(5),getvarl(np)
	character*500 ifile1,ifile2
        INTEGER*4   base_time,IDATE(ntime)  
        character*8 CDATE(ntime)
  	REAL*4 v4(nlon,nlat,ntime)
	Integer*2 ipout(nlon,nlat,ntime,np)
	Integer*2 ipout1(nlon,nlat,np)


	data getvarl /&
      'AEROD', & ! 1
      'DSTODX01',& ! 2
      'DSTODX02',& ! 3
      'DSTODX03',& ! 4 
      'DSTODX04',& ! 5
      'SO4OD',& ! 6
      'SSLTOD',& ! 7 old var name 'TAUSSLT      '
      'BCPHIOD',& ! 8
      'BCPHOOD',& ! 9 
      'OCPHIOD',& ! 10
      'OCPHOOD'/ ! 11


!*****************************************************************
	real rlat(nlat),rlon(nlon)
	data rlat /                                                    &
      -88.54195,-86.65317,-84.75323, -82.85078,-80.94736,-79.0434,     &
      -77.13935, -75.23505, -73.33066, -71.42619,-69.52167,-67.617,    &
      -65.71251, -63.8079, -61.90326, -59.99861,-58.09395,-56.1892,    &
      -54.2846, -52.37991, -50.47522, -48.57052, -46.66582,-44.7611,   &
      -42.8564, -40.95169, -39.04697, -37.14225, -35.23753,-33.3328,   &
      -31.42808, -29.52336, -27.61863, -25.7139, -23.80917,-21.9044,   &
      -19.99971,-18.09498,-16.19024,-14.28551, -12.38078,-10.4760,     &
      -8.571308, -6.666573, -4.761838, -2.857103, -0.952368,0.95236,   &
      2.857103, 4.761838, 6.666573,8.571308,10.47604,12.38078,14.285,  &
      16.19024, 18.09498, 19.99971,21.90444,23.80917,25.7139, 27.618,  &
      29.52336, 31.42808, 33.33281,35.23753,37.14225,39.04697,40.951,  &
      42.8564, 44.76111, 46.66582,48.57052,50.47522,52.37991, 54.284,  &
      56.18928, 58.09395, 59.99861,61.90326,63.8079,65.71251, 67.617,  &
      69.52167, 71.42619, 73.33066,75.23505,77.13935,79.04349,80.947,  &
       82.85078, 84.75323, 86.65317,88.54195 /
	data rlon /                                                     &
      0,1.875,3.75,5.625,7.5,9.375,11.25,13.125,15,16.875,18.75,	&
      20.625,22.5,24.375,26.25,28.125,30,31.875,33.75,35.625,37.5,	&
      39.375,41.25,43.125,45,46.875,48.75,50.625,52.5,54.375,56.25,	&
      58.125,60,61.875,63.75,65.625,67.5,69.375,71.25,73.125,75,	&
      76.875,78.75,80.625,82.5,84.375,86.25,88.125,90,91.875,93.75,	&
      95.625,97.5,99.375,101.25,103.125,105,106.875,108.75,110.625,	&
      112.5,114.375,116.25,118.125,120,121.875,123.75,125.625,127.5,	&
      129.375,131.25,133.125,135,136.875,138.75,140.625,142.5,144.375,  &
      146.25,148.125,150,151.875,153.75,155.625,157.5,159.375,161.25,	&
      163.125,165,166.875,168.75,170.625,172.5,174.375,176.25,178.125,  &
      180,181.875,183.75,185.625,187.5,189.375,191.25,193.125,195,	&
      196.875,198.75,200.625,202.5,204.375,206.25,208.125,210,211.875,  &
      213.75,215.625,217.5,219.375,221.25,223.125,225,226.875,228.75,	&
      230.625,232.5,234.375,236.25,238.125,240,241.875,243.75,245.625,  &
      247.5,249.375,251.25,253.125,255,256.875,258.75,260.625,262.5,	&
      264.375,266.25,268.125,270,271.875,273.75,275.625,277.5,279.375,  &
      281.25,283.125,285,286.875,288.75,290.625,292.5,294.375,296.25,	&
      298.125,300,301.875,303.75,305.625,307.5,309.375,311.25,313.125,  &
      315,316.875,318.75,320.625,322.5,324.375,326.25,328.125,330,	&
      331.875,333.75,335.625,337.5,339.375,341.25,343.125,345,346.875,  &
      348.75,350.625,352.5,354.375,356.25,358.125/



	INTEGER:: i,j,ida,iyr,imo,imop,idap , idate1,iyrp,lenstr,itime,k,kkk
	INTEGER:: idate2,ip,nvars,ngatts,idx,ncid,nvs,ndsize,nvdim,ivarid,ntp,ndims
	

	contains
!========================================================================
	Subroutine readit
!input raw file
	call getarg(1,ifile1) !NETCDF MATCH ED2
	call getarg(2,ifile2) !BIN AOT
	
!	if ( collection == 'C4' .or. collection == 'C5' .or. collection == 'XX') then
!	print*,'NOTE: Will label output file as collection =',collection
!	else
!	stop ' What Collection is this C4 or C5'
!	endif
	
!	idx = index(ifile1,' ')-1
!output processed file
!	ifile2='./processed/aots/MATCH_TERRA_AOTS_MODIS.'//collection//'.'
	print*,trim(ifile2)
	
	 ios = nf90_open(trim(ifile1), nf90_nowrite,NCID)
	if ( ios .ne. 0) stop ' Bad nf90_open'
!      NCID=NCOPN(trim(ifile1),NCNOWRIT,RCODE)
	print*, ncid
      CALL NCINQ(NCID,NDIMS,NVARS,NGATTS,RECDIM,RCODE)
	print*, NCID,NDIMS,NVARS,NGATTS,RECDIM,RCODE

!----Date
      Ios = nf90_inq_varid(NCID,'date          ',ivarid )
!	 ivarid = ncvid(ncid,'date          ',rcode)

	print*, ivarid
	CALL NCVINQ(NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
	print*,NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE
	  do  J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  	enddo

      CALL NCVGT(NCID,ivarid,START,COUNT,IDATE,RCODE)

!EARLIER VERSIONS The information contained in IDATE represents the avg of
!the previous 24 hours ie (20010315) = 0z march 14th thru 0z march15
!!
	print*,'DATE:',idate
	do itime =1,ntime
	idate1= idate(itime)
	call idate_chop(idate1,iyr,imo,ida)
	print*,'YMD',iyr,imo,ida

!	call  previous_day(iyr,imo,ida, iyrp,imop,idap )
	iyrp = iyr ; imop = imo ; idap = ida ! SAME DAY USED FOR MODIS 

	idate2 = iyrp*10000+imop*100+idap
	write (cdate(itime),'(I8)')idate2
	print*,'PREVIOUS',itime,idate1,idate2,cdate(itime)
	
! At This point we give the file name (CDATE) that of the
!! PREVIOUS DAY so
! the (20010314) represents 0z march 14th thru 0z march15
! centered arounf 20000314 (12z)
! 

	enddo
	
!---- Contituents
	CONSTITUENT : do ip=1,np
	print'(i3,a1,a31)',ip,'-',getvarl(ip)
 !     ivarid = ncvid(ncid,getvarl(ip),rcode)
      Ios = nf90_inq_varid(NCID,getvarl(ip) ,ivarid )      
      CALL NCVINQ(NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
!  	print*,NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE
      LENSTR=1

      do  J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  	enddo

      CALL NCVGT(NCID,ivarid,START,COUNT,v4(1,1,1),RCODE)

!!! PACK & OUTPUT------------------------------------------------------


	TIME : do k=1, ntime
	do j=1,nlat
	do i=1,nlon

	if (v4 (i,j,k) <= 3.0 )then
	 ipout(i,j,k,ip) = v4(i,j,k)*10000
	else
	 ipout(i,j,k,ip) = 3.0 * 10000
	endif

	enddo;enddo

	enddo TIME 
	enddo CONSTITUENT
!---------------------------------------------------------------------

	TIME2 : do k  =1,ntime
	open(51,file=trim(ifile2),form='Unformatted', &
       access='direct',recl=nlon*nlat*np*2)
	ipout1= ipout(1:nlon,1:nlat,k,1:np)
	write(51,rec=1) ipout1
	close(51)
	print*,k,IDATE(k)
	enddo TIME2
	

	print*,' Complete READIT'
!-------------------------------------------------------
	

	return

	end subroutine readit
!---------------------------------
	subroutine idate_chop(idate,iyr,imo,ida)
	integer idate,iyr,imo,ida
	integer kkk
	kkk=idate
	iyr = kkk/10000
	kkk = kkk-iyr*10000
	imo = kkk/100
	kkk = kkk-imo*100
	ida = kkk
	if ( ida == 42 ) ida = 31
	if ( ida == 41 .and. imo /= 2) ida = 30
	if ( ida == 41 .and. imo == 2) ida = 29
!	print*,iyr,imo,ida
	end subroutine idate_chop
!============================================================
	subroutine previous_day(iyr,imo,ida, iyrp,imop,idap )
	integer iyr,imo,ida,iyrp,imop,idap
	integer ideom(12)
	data ideom/31,28,31,30,31,30,31,31,30,31,30,31/


	  iyrp= iyr
	  imop= imo
	  idap= ida
	if ( ida >1 ) then
	  idap= ida-1
	  return
	 else

	if ( imo > 1 ) then
	 imop= imo-1
	 idap = ideom(imop)
	if ( imo == 3 .and. ida==1 ) then
	 if (mod(iyr,4  ) == 0 ) idap = 29
	 if (mod(iyr,100) == 0 ) idap = 28
	 if (mod(iyr,400) == 0 ) idap = 29
	endif

	else
	 iyrp=iyr-1
	 imop=12
	 idap=ideom(imop)
	endif

	endif


	end subroutine previous_day
!================================================================
	END MODULE READNETCDF_AOTS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!================================================================
subroutine process_aots
 USE READNETCDF_AOTS
 call readit
end subroutine process_aots
!------------------------------------
subroutine process_vertical
 USE READNETCDF_VERTICAL
 call process_hdf_match_profile
end subroutine process_vertical
!------------------------------------
 call process_aots
 call process_vertical
 call process_sage_aot
 call process_sage_profile
	stop
	end
