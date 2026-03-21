c     $Header$
      
c     Purpose: dstrad.F contains radiation subroutines for the dust model
      
      subroutine dst_trn(
     $     mpp_dst,             ! I [kg m-2] Dust path above kth interface level
     $     dst_trn_ttl)         ! O [frc] Total transmission through all size classes between interfaces k1 and k2
c     Purpose: Compute dust transmissions needed in absorptivity/emissivity calculations
c     dst_trn() is called by CCM:physics/radclw() when doabsems is .true.

      use precision

      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstrad.h> /* Dust radiative properties */
c     Parameters
c     Input
      real(r8) mpp_dst(plond,plevp,dst_nbr) ! I [kg m-2] Dust path above kth interface level
c     Output
      real(r8) dst_trn_ttl(plond,plevp,plevp,bnd_nbr_LW) ! O [frc] Total transmission through all size classes between interfaces k1 and k2
c     Local
      integer bnd_idx           ! [idx] LW band index
      integer i                 ! [idx] lon index
      integer k1                ! [idx] lev index
      integer k2                ! [idx] lev index
      integer m                 ! [idx] cst index
      real(r8) dst_pth_dlt          ! [kg m-2] Dust path between interface levels k1 and k2
      real(r8) odap_dst             ! [frc] Path absorption optical depth of this size class
      real(r8) odap_dst_ttl         ! [frc] Total path absorption optical depth of all size classes
      
      do bnd_idx=1,bnd_nbr_LW
         do k1=1,plevp          ! NB: end at plevp
            do k2=1,plevp       ! NB: end at plevp
               do i=1,plon
                  odap_dst_ttl=0.0 ! [frc]
                  do m=1,dst_nbr
                     dst_pth_dlt=abs(mpp_dst(i,k1,m)-mpp_dst(i,k2,m)) ! [kg m-2] Dust path between interface levels k1 and k2
                     odap_dst=abs_cff_mss_dst(bnd_idx,m)*dst_pth_dlt ! [frc] Absorption optical depth of this size class between interfaces k1 and k2
                     odap_dst_ttl=odap_dst_ttl+odap_dst ! [frc] Total absorption optical depth of all size classes between interfaces k1 and k2
                  end do        ! end loop over cst
c     fxm: Scale aerosol extinction to account for scattering as in CLT99
                  dst_trn_ttl(i,k1,k2,bnd_idx)=exp(-1.66*odap_dst_ttl) ! [frc] Total transmission through all size classes between interfaces k1 and k2
               end do           ! end loop over lon
            end do              ! end loop over lev (k2)
         end do                 ! end loop over lev (k2)
      end do                    ! end loop over bnd
      
      return
      end                       ! end dst_trn()
      
      
