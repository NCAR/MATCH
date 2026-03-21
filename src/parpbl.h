c
c $Id$
c $Author$
c
C
C Constants for the surface temperature calculation
C
      real(r8) zzocen   ! Ocean aerodynamic roughness length
      real(r8) zzsice   ! Sea ice aerodynamic roughness length
      real(r8) xkar     ! Von Karman constant
      real(r8) hmixmn   ! Minimum boundary layer height
      real(r8) ric      ! Critical Richardson number
      real(r8) rimax    ! Maximum Richardson number
      real(r8) epsi     ! Small factor to prevent exact zeros
      real(r8) zref     ! 10m reference height 
      real(r8) umin     ! Minimum wind speed at bottom level
C
      parameter ( zzocen =    0.0001  ,
     $            zzsice =    0.0400   ,
     $              xkar =    0.4      ,
     $            hmixmn =  500.0      ,
     $               ric =    3.05     ,
     $             rimax =    0.50*ric ,
     $              epsi =    1.0e-12  ,
     $              zref =   10.0      ,
     $              umin =    1.0      )  
