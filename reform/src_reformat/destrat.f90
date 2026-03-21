Module match_destrat
  
USE SAGE_CLIM , only: rd_sage_clim,ilat_sage_q,sage
implicit none
 character*300 filein,fileout

 integer ,parameter :: nlon=192 ,nlat=94,nlev=28,np=10
 real ,parameter    :: omit=-9999.

 integer*1  bdat(nlon,nlat,nlev,np)
 
 integer nbad(np,2)
 integer imonth 
CONTAINS
!------------------------------------------------------------------------
subroutine IN_FILE_destrat

 open(51,file=trim(filein),form='Unformatted',action='read',status='old', &
 access='direct',recl=nlon*nlat*nlev*np)
 read(51,rec=1) bdat
 close(51)

 nbad =0
 
end subroutine IN_FILE_destrat
!------------------------------------------------------------------------
subroutine OUT_FILE_destrat

 open(52,file=trim(fileout),form='Unformatted', &
 access='direct',recl=nlon*nlat*nlev*np)
 write(52,rec=1) bdat
 close(52)
 
 
 print '(2(10i6,3x))',nbad
end subroutine OUT_FILE_destrat
!------------------------------------------------------------------------


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

!======================================================
subroutine doit
integer ilat,ilon,ilev,ip
real s1,s2,s3
real prf1(nlev),prf0(nlev),sage0(nlev)


LAT:do ilat=1,nlat
sage0(1:nlev) = 0.0
sage0(1:10)    = sage%prf(1:10, ilat_sage_q(ilat) ,imonth)

LON: do ilon=1,nlon

PARM: do ip =1,np



if ( ip == 5) then !SO4 SAGE
prf1(1:nlev) = sage0(1:nlev)
else
 do ilev = 1,nlev  
  prf0(ilev) = upbv( bdat(ilon,ilat,ilev,ip) ) 
  prf1(ilev) =prf0(ilev)
 enddo

 s1 = maxval(prf1)
 s2 = maxval(prf1(11:nlev))
 prf1(1:10) =0.0 ! Zero out Stratosphere below 0.18783Level , 0.2101Layer

if ( s2 < 100 .and. s2 > 0.0 ) prf1(11:nlev) = 100.0* prf1(11:nlev)/s2  ! Reset max value to 100%

s3 = maxval(prf1) !Check
if (s3 == 0.0 ) then ! In case of no AOT at all just put AOT in trop below .187
 prf1(1:10) =0.0
 prf1(11:nlev) =100.
endif







   if ( s2 < 10 .and. s3 .ne. 0 ) then
    nbad(ip,1)=nbad(ip,1) +1 
!    print'(3i4,3f6.1,3x,10f6.1,2x,18f6.1)',ilon,ilat,ip ,s1,s2,s3,prf0
   endif
   if (s3 == 0 ) then
       nbad(ip,2)=nbad(ip,2) +1 
!      print*,'s3 =0 '

   endif
   
endif ! elseSO4 SAGE


!print'(I4,28f6.1)',ip,prf1



 do ilev = 1,nlev  
  bdat(ilon,ilat,ilev,ip)  =  pbv( prf1(ilev) )   
 enddo


enddo PARM

enddo LON
enddo LAT
end subroutine doit



end MODULE match_destrat

!=======================================================
subroutine process_sage_profile
USE match_destrat
implicit none
 character*2 MM
integer im,ilat

 call getarg(3,filein) ! Bin Profile

 call getarg(6,fileout) !Bin Saged Profile
 
 call getarg(4,MM) !Month 01-12

  read(MM,'(i2.2)') imonth
!  print*,'MM=',MM ,'  imonth=',imonth

  call rd_sage_clim 

!  do ilat=1,32 ; print'(f8.1,f8.4,10f6.1)',sage%lat(ilat),sage%aot(ilat,imonth),sage%prf(1:10,ilat,imonth); enddo

 
! print*,'IN:',trim(filein)
 call IN_FILE_destrat
 
 
 
 call doit
 
 
 
 
!  print*,'OUT:',trim(fileout)

   call OUT_FILE_destrat
end  subroutine process_sage_profile
 
