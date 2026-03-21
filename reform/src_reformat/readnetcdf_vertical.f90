	MODULE READNETCDF_VERTICAL
	use netcdf
	implicit none
	integer ios
	integer ,parameter :: nlon=192 ,nlat=94,nlev=28,ntime=1
        integer, parameter :: NREC=nlon*nlat, nf=4,np=10
	real ,parameter    :: omit=-9999.
	INTEGER*4 RECDIM,RCODE,START(10),COUNT(10),VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
	
	character*31 DUMMY,cname(5)
	character*14 getvarl(0:np)
	character*500 ifile1,ifile2
        character*2 collection
        integer*4   base_time  ,IDATE(ntime)  
        character*8 CDATE(ntime)
  	real*4 v4(nlon,nlat,nlev,ntime)
	real*8 v8(nlon,nlat,nlev,ntime)
	real*4 ps(nlon,nlat,ntime)

        real*8 ps8(nlon,nlat,ntime)


	integer*1  ipout(nlon,nlat,nlev,np)
		

	integer ip,i,j,k,l
	integer idate2,iyr,imo,ida

	real hybi(nlev+1),psfc,hybm(nlev)

!	include '/CERES/sarb/home/rose/netcdf/netcdf.inc.f90'
	include 'netcdf.inc'
	
	data getvarl / &
      'PS',   &
      'DSTQ01',    & !1
      'DSTQ02',    & !2
      'DSTQ03',    & !3
      'DSTQ04',    & !4
      'SO4' ,   & !5
      'SSLT',    & !6
      'BCPHI',    & !7
      'BCPHO',    & !8
      'OCPHI',    & !9
      'OCPHO'/      !10

!	data getvarl /
!     & 'AEROD       ',  ! 1
!     & 'DSTODX01     ', ! 2
!     & 'DSTODX02     ', ! 3
!     & 'DSTODX03     ', ! 4 
!     & 'DSTODX04     ', ! 5
!     & 'SO4OD        ', ! 6
!     & 'SSLTOD       ', ! 7 old var name 'TAUSSLT      '
!     & 'BCPHIOD      ', ! 8
!     & 'BCPHOOD      ', ! 9 
!     & 'OCPHIOD      ', ! 10
!     & 'OCPHOOD      '/ ! 11


!*****************************************************************
	real rlat(nlat),rlon(nlon)
	data rlat /                                                   &
      -88.54195,-86.65317,-84.75323, -82.85078,-80.94736,-79.0434,    &
      -77.13935, -75.23505, -73.33066, -71.42619,-69.52167,-67.617,   &
      -65.71251, -63.8079, -61.90326, -59.99861,-58.09395,-56.1892,   &
      -54.2846, -52.37991, -50.47522, -48.57052, -46.66582,-44.7611,  &
      -42.8564, -40.95169, -39.04697, -37.14225, -35.23753,-33.3328,  &
      -31.42808, -29.52336, -27.61863, -25.7139, -23.80917,-21.9044,  &
      -19.99971,-18.09498,-16.19024,-14.28551, -12.38078,-10.4760,     &
      -8.571308, -6.666573, -4.761838, -2.857103, -0.952368,0.95236,  &
      2.857103, 4.761838, 6.666573,8.571308,10.47604,12.38078,14.285, &
      16.19024, 18.09498, 19.99971,21.90444,23.80917,25.7139, 27.618, &
      29.52336, 31.42808, 33.33281,35.23753,37.14225,39.04697,40.951, &
      42.8564, 44.76111, 46.66582,48.57052,50.47522,52.37991, 54.284, &
      56.18928, 58.09395, 59.99861,61.90326,63.8079,65.71251, 67.617, &
      69.52167, 71.42619, 73.33066,75.23505,77.13935,79.04349,80.947, &
       82.85078, 84.75323, 86.65317,88.54195 /
	data rlon /                                                   &
      0,1.875,3.75,5.625,7.5,9.375,11.25,13.125,15,16.875,18.75,       &
      20.625,22.5,24.375,26.25,28.125,30,31.875,33.75,35.625,37.5,     &
      39.375,41.25,43.125,45,46.875,48.75,50.625,52.5,54.375,56.25,    &
      58.125,60,61.875,63.75,65.625,67.5,69.375,71.25,73.125,75,       &
      76.875,78.75,80.625,82.5,84.375,86.25,88.125,90,91.875,93.75,    &
      95.625,97.5,99.375,101.25,103.125,105,106.875,108.75,110.625,    &
      112.5,114.375,116.25,118.125,120,121.875,123.75,125.625,127.5,   &
      129.375,131.25,133.125,135,136.875,138.75,140.625,142.5,144.375, &
      146.25,148.125,150,151.875,153.75,155.625,157.5,159.375,161.25,  &
      163.125,165,166.875,168.75,170.625,172.5,174.375,176.25,178.125, &
      180,181.875,183.75,185.625,187.5,189.375,191.25,193.125,195,     &
      196.875,198.75,200.625,202.5,204.375,206.25,208.125,210,211.875, &
      213.75,215.625,217.5,219.375,221.25,223.125,225,226.875,228.75,  &
      230.625,232.5,234.375,236.25,238.125,240,241.875,243.75,245.625, &
      247.5,249.375,251.25,253.125,255,256.875,258.75,260.625,262.5,   &
      264.375,266.25,268.125,270,271.875,273.75,275.625,277.5,279.375, &
      281.25,283.125,285,286.875,288.75,290.625,292.5,294.375,296.25,  &
      298.125,300,301.875,303.75,305.625,307.5,309.375,311.25,313.125, &
      315,316.875,318.75,320.625,322.5,324.375,326.25,328.125,330,     &
      331.875,333.75,335.625,337.5,339.375,341.25,343.125,345,346.875, &
      348.75,350.625,352.5,354.375,356.25,358.125/

	data hybi / &
      0.00100,0.00657,0.01386,0.02309,0.03469,0.0492, &
      0.06723,0.08945,0.11654,0.14916,0.18783,0.23286, &
      0.28421,0.34137,0.40334,0.4686,0.53529,0.60135, &
      0.66482,0.72401,0.77773,0.82527,0.86642,0.90135, &
      0.93054,0.95459,0.97418,0.99,1.0/

data hybm /0.0027,0.0101,0.0183,0.0288,0.0418,0.058,&
           0.0782,0.1028,0.1326,0.1682,0.2101,0.2582,&
           0.3125,0.372, 0.4357,0.5017,0.5681,0.6329,&
           0.6943,0.7508,0.8014,0.8458,0.8838,0.9159,&
	   0.9425,0.9644,0.9821,0.995/

	contains
!========================================================================
	Subroutine process_hdf_match_profile

	integer NCID,NDIMS,NVARS,NGATTS
	integer LENSTR,NDSIZE,NTP,NVDIM,NVS
	integer iyrp,imop,idap
	integer idate1,idx,itime,ivarid
	integer RECDIM,RCODE

!input raw file
	call getarg(1,ifile1) ! MATCH NETCDF ED2	
	call getarg(3,ifile2) !!BIN Profile
	
!	if ( collection == 'C4' .or. collection == 'C5' .or. collection == 'XX') then
!	print*,'NOTE: Will label output file as collection =',collection
!	else
!	stop ' What Collection is this C4 or C5 or XX'
!	endif
	
!	idx = index(ifile1,' ')-1
!output processed file
!	ifile2='./processed/vertical/MATCH_TERRA_VERTICAL_MODIS.'//collection//'.'
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
	print'(a4,3i5)','YMD',iyr,imo,ida

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
	CONSTITUENT : do ip=0,np
	print*,ip,'-',getvarl(ip)
 !     ivarid = ncvid(ncid,trim(getvarl(ip)),rcode)
      Ios = nf90_inq_varid(NCID,getvarl(ip) ,ivarid ) 
      CALL NCVINQ(NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
  	print*,NCID,ivarid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE
      LENSTR=1

      do  J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
	print*,'XXXX',j,NDSIZE
  	enddo

	if (ip == 0 ) then
	 
	 if ( collection == 'XX' ) then
	  CALL NCVGT(NCID,ivarid,START,COUNT,ps(1,1,1),RCODE)
!	  ps=ps8
	 else
	  CALL NCVGT(NCID,ivarid,START,COUNT,ps(1,1,1),RCODE)
	 endif
	 
	 
	else
	  if ( collection == 'XX' ) then
	    CALL NCVGT(NCID,ivarid,START,COUNT,v4(1,1,1,1),RCODE)
!	    v4=v8
	  else
	    CALL NCVGT(NCID,ivarid,START,COUNT,v4(1,1,1,1),RCODE)
	  endif
	endif

	
!!! PACK & OUTPUT------------------------------------------------------

	
	if (ip >0 ) call cvt(ip)


 
	enddo CONSTITUENT
!---------------------------------------------------------------------

!	stop ' PREMATURE'

	TIME2 : do k  =1,ntime
	open(51,file=trim(ifile2),form='Unformatted', &
        access='direct',recl=nlon*nlat*nlev*np)
	write(51,rec=1) ipout
	close(51)
	print*,k,IDATE(k)
	enddo TIME2
	

	print*,' Complete READIT'	
	return

	end subroutine process_hdf_match_profile
!========================================================================

	subroutine idate_chop(idate,iyr,imo,ida)
	integer kkk,idate,iyr,imo,ida
	
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
!========================================================================

	subroutine cvt(ip)
	
	real pmid(nlev),pctn(nlev),amt(nlev)
	real s1 
	integer ip,i,j,k,l

         
	TIME : do k=1, ntime
	do j=1,nlat
	do i=1,nlon

	s1 =0
	do l=1,nlev
        psfc = ps(i,j,k)*0.01
	pmid(l)   = psfc*hybm(l)
	amt(l)    = v4(i,j,l,k) 
	enddo
!---------------------
	s1 = maxval(amt)

	if ( s1 > 1E-35) then
	do l=1,nlev
	pctn(l) = 100.0* amt(l)/s1
!	print'(2f8.2,f10.3)',pbnd(l),pbnd(l+1),pct(l)
	enddo
	else
	 pctn(1:nlev-1) =0.
	 pctn(nlev) =100.
	endif

!PACK
	do l=1,nlev
	if ( pctn(l) < 0.0 .or. pctn(l) > 100.001) STOP ' OUTOF RANGE PCT' 
	 ipout(i,j,l,ip) = pbv( pctn(l) )
	enddo

!	print'(28i4,I3,f6.1)',nint(pctn),ip,psfc
	
	enddo;enddo

	enddo TIME 
        
	end subroutine cvt

!-------------------------------
integer*1 function pbv(x)
real x
pbv = nint(x *2.50)  !!! PACK (0:100)  to Byte nearest 0.40
end function pbv
!-------------------------------
real function upbv(b1)
! unpack a PBV packed byte a 0 to 100 value
real x
integer*1 b1

if  ( b1 < 0 ) then
 x= b1+256
else
 x=b1
endif

upbv = x*0.4000
end function upbv
!================================================================

	END MODULE READNETCDF_VERTICAL
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!================================================================

!	USE READNETCDF_VERTICAL
	
!	call process_hdf_match_profile

!!!!	call read_processed

!	stop
!	end
