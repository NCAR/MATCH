c $Id$
C
C Common block and statement functions for saturation vapor pressure
C look-up procedure, J. J. Hack, February 1990
C
      integer plenest  ! length of saturation vapor pressure table
      parameter (plenest=250)
C
C Table of saturation vapor pressure values es from tmin degrees
C to tmax+1 degrees k in one degree increments.  ttrice defines the
C transition region where es is a combination of ice & water values
C
      common/comes0/estbl(plenest),tmin,tmax,ttrice,pcf(6),
     $             epsqs,rgasv,hlatf,hlatv,cp,icephs
C
      real(r8) estbl,  ! table values of saturation vapor pressure
     $     tmin,   ! min temperature (K) for table
     $     tmax,   ! max temperature (K) for table
     $     ttrice, ! transition range from es over water to es over ice
     $     pcf,    ! polynomial coeffs -> es transition water to ice
     $     epsqs,  ! Ratio of h2o to dry air molecular weights 
     $     rgasv,  ! Gas constant for water vapor
     $     hlatf,  ! Latent heat of vaporization
     $     hlatv,  ! Latent heat of fusion
     $     cp      ! specific heat of dry air
      logical
     $     icephs  ! false => saturation vapor pressure over water only
C
C Dummy variables for statement functions
C
      real(r8) td,     ! dummy variable for function evaluation
     $     tlim,   ! intermediate variable for es look-up with estbl4 
     $     estblf, ! statement function es look-up
     $     estbl4  ! statement function es look-up
C
C Statement functions used in saturation vapor pressure table lookup
C there are two ways to use these three statement functions.
C For compilers that do a simple in-line expansion:
C => ttemp = tlim(t)
C    es    = estbl4(ttemp)
C
C For compilers that provide real optimization:
C => es    = estblf(t)
C
      tlim(td) = max(min(td,tmax),tmin)
C
      estblf(td) =  (tmin + int(tlim(td)-tmin) - tlim(td) + 1.0)
     $            *estbl(int(tlim(td)-tmin)+1)
     $            -(tmin + int(tlim(td)-tmin) - tlim(td)      )
     $            *estbl(int(tlim(td)-tmin)+2)
C
      estbl4(td) =  (tmin+int(td-tmin)+1.0-td)*estbl(int(td-tmin)+1)
     $            + ( td-(tmin+int(td-tmin)) )*estbl(int(td-tmin)+2)
C
