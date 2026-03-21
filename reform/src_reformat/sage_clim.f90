Module SAGE_Clim
implicit none
integer,parameter :: nlev=10,nlat=32,nmth=12
type SAGE_TYP
sequence
real prf(nlev,nlat,nmth)
real aot(nlat,nmth)
real lat(nlat)
end type SAGE_TYP

TYPE (SAGE_TYP) sage

 contains
!=======================================================
subroutine rd_sage_clim
 character*80 line
integer im ,ilat

open(21,file='/homedir/dfillmor/MATCH/reform/sage.2000_2004.matchprf.2.txt',action='read',status='old')


do im=1,12
 read(21,*) line
 do ilat=1,32
 read(21,*) sage%lat(ilat), sage%aot(ilat,im),sage%prf(1:10,ilat,im)
 enddo
enddo

 close(21)

end subroutine rd_sage_clim
!======================================================
integer function ilat_sage_q(ilatm)
integer ilatm
integer icvt(94)

data icvt( 1: 8) / 8*1 / 
data icvt( 9:10) / 2*2 /
data icvt(11:13) / 3*3 /
data icvt(14:15) / 2*4 /
data icvt(16:18) / 3*5 /
data icvt(19:21) / 3*6 /
data icvt(22:23) / 2*7 /
data icvt(24:26) / 3*8 /
data icvt(27:29) / 3*9 /
data icvt(30:31) / 2*10/
data icvt(32:34) / 3*11/
data icvt(35:36) / 2*12/
data icvt(37:39) / 3*13/
data icvt(40:42) / 3*14/
data icvt(43:44) / 2*15/
data icvt(45:47) / 3*16/
data icvt(48:50) / 3*17/
data icvt(51:52) / 2*18/
data icvt(53:55) / 3*19/
data icvt(56:58) / 3*20/
data icvt(59:60) / 2*21/
data icvt(61:63) / 3*22/
data icvt(64:65) / 2*23/
data icvt(66:68) / 3*24/
data icvt(69:71) / 3*25/
data icvt(72:73) / 2*26/
data icvt(74:76) / 3*27/
data icvt(77:79) / 3*28/
data icvt(80:81) / 2*29/
data icvt(82:84) / 3*30/
data icvt(85:86) / 2*31/
data icvt(87:94) / 8*32/

 ilat_sage_q = icvt(ilatm)

end function ilat_sage_q
!======================================================
integer function ilat_sage(ilatm)
integer is,ilatm

   REAL,      PARAMETER :: rlatca ( 94) =                             &
      (/ -88.54195,-86.65317,-84.75323, -82.85078,-80.94736,-79.0434,      &
         -77.13935, -75.23505, -73.33066, -71.42619,-69.52167,-67.617,     &
         -65.71251, -63.8079, -61.90326, -59.99861,-58.09395,-56.1892,     &
         -54.2846, -52.37991, -50.47522, -48.57052, -46.66582,-44.7611,    &
         -42.8564, -40.95169, -39.04697, -37.14225, -35.23753,-33.3328,    &
         -31.42808, -29.52336, -27.61863, -25.7139, -23.80917,-21.9044,    &
         -19.99971,-18.09498,-16.19024,-14.28551, -12.38078,-10.4760,      &
          -8.571308, -6.666573, -4.761838, -2.857103, -0.952368,0.95236,    &
           2.857103, 4.761838, 6.666573,8.571308,10.47604,12.38078,14.285,  &
          16.19024, 18.09498, 19.99971,21.90444,23.80917,25.7139, 27.618,   &
          29.52336, 31.42808, 33.33281,35.23753,37.14225,39.04697,40.951,   &
          42.8564, 44.76111, 46.66582,48.57052,50.47522,52.37991, 54.284,   &  
          56.18928, 58.09395, 59.99861,61.90326,63.8079,65.71251, 67.617,   &
          69.52167, 71.42619, 73.33066,75.23505,77.13935,79.04349,80.947,   &
          82.85078, 84.75323, 86.65317,88.54195 /)

ilat_sage =-1

 if( rlatca(ilatm) < sage%lat(1)-2.5 ) then
  ilat_sage = 1
 elseif  ( rlatca(ilatm) > sage%lat(nlat)+2.5 ) then
  ilat_sage = nlat
 else
  do is=1,nlat 
   if ( rlatca(ilatm) >= sage%lat(is)-2.5 .and. rlatca(ilatm) <= sage%lat(is)+2.5 ) then
    ilat_sage = is
    endif
  enddo
endif

end function ilat_sage

end Module SAGE_Clim
!========================================================================
subroutine test_Sage

USE SAGE_CLIM , only: rd_sage_clim, ilat_sage,ilat_sage_q,sage
implicit none
integer ilatm,ilat ,im 
 call rd_sage_clim


do ilatm=1,94
print*,ilatm,ilat_sage_q(ilatm), ilat_sage(ilatm)-ilat_sage_q(ilatm)
enddo

 do im=1,12
 do ilat=1,32
  print'(f8.1,f8.4,10f6.1)',sage%lat(ilat),sage%aot(ilat,im),sage%prf(1:10,ilat,im)
 enddo
 enddo
end subroutine test_Sage
!===================================================================
! call test_sage ; end
