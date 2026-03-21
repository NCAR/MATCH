Module match_sageaot
  
USE SAGE_CLIM , only: rd_sage_clim,ilat_sage_q,sage
implicit none
 character*300 filein,fileout

 integer ,parameter :: nlon=192 ,nlat=94,ntyp=11
 real ,parameter    :: omit=-9999.

 integer*2  ipaer0(nlon,nlat,ntyp) 
 integer*2  ipaer1(nlon,nlat,ntyp)
 
 integer imonth 
 real, parameter :: tiny = 0.0001
CONTAINS
!------------------------------------------------------------------------
subroutine IN_FILE_aot

 open(51,file=trim(filein),form='Unformatted',action='read',status='old', &
 access='direct',recl=nlon*nlat*ntyp*2)
 read(51,rec=1) ipaer0
 close(51)

 
end subroutine IN_FILE_aot
!------------------------------------------------------------------------
subroutine OUT_FILE_aot

 open(52,file=trim(fileout),form='Unformatted', &
 access='direct',recl=nlon*nlat*ntyp*2)
 write(52,rec=1) ipaer1
 close(52)
 
 
! print '(2(10i6,3x))',nbad
end subroutine OUT_FILE_aot
!------------------------------------------------------------------------

!======================================================
subroutine doit_aot
integer ilat,ilon,ip
real sageaot
real aots0(ntyp),aots1(ntyp)

LAT:do ilat=1,nlat

sageaot    = sage%aot(ilat_sage_q(ilat) ,imonth)

 LON: do ilon=1,nlon

  aots0(1:ntyp) = ipaer0(ilon,ilat,1:ntyp) / 10000.0
 
  aots1(1)  = aots0(1) ! TOTAL AOT
  aots1(2)  = aots0(2)  ! ' DSTQ01'
  aots1(3)  = aots0(3)  ! ' DSTQ02' 
  aots1(4)  = aots0(4)  ! ' DSTQ03'
  aots1(5)  = aots0(5)  ! ' DSTQ04'
  aots1(6)  = sageaot   ! '    SO4' --> SAGE
  aots1(7)  = aots0(7)  ! '   SSLT'
  aots1(8)  = aots0(8)  ! '  BCPHI'
  aots1(9)  = aots0(9)  ! '  BCPHO'
  aots1(10) = (aots0(10) + aots0(6) )- sageaot ! '  OCPHI' --> (OCPHI+SO4 - SAGE)
  aots1(11) = aots0(11) ! '  OCPHO'


  if ( aots1(10) < tiny)  aots1(10) = tiny
  
  if ( sum( aots1(2:11)) < sageaot+tiny .or. aots1(1) < sageaot+tiny ) aots1(1) = sageaot+tiny
 
 
 
  ipaer1(ilon,ilat,1:ntyp) = aots1(1:ntyp) * 10000.0

 enddo LON
enddo LAT
end subroutine doit_aot



end MODULE match_sageaot

!=======================================================
subroutine process_sage_aot
USE match_sageaot
implicit none
 character*2 MM

 call getarg(2,filein) !BIN AOT
 
 call getarg(5,fileout) !Bin Saged AOT
 
 call getarg(4,MM)  ! MONTH 01-12

  read(MM,'(i2.2)') imonth
!  print*,'MM=',MM ,'  imonth=',imonth

  call rd_sage_clim 

!  print'(i4,32f8.5)',imonth,sage%aot(1:32,imonth)

 
! print*,'IN:',trim(filein)
 call IN_FILE_aot
  
 call doit_aot
 
!  print*,'OUT:',trim(fileout)

   call OUT_FILE_aot
end  subroutine process_sage_aot
 
