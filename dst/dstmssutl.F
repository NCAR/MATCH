c     $Header$
      
c     Purpose: Mass budget utilities for the dust model
      
      subroutine dst_zero(q_dst)
c     Purpose: Set dust tracers to zero
c     dst_zero() is currently not called (since inidat() does this by default)

      use precision

      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
c     Input/Output
      real(r8) q_dst(plond,plev,dst_nbr) ! kg dust/kg dry air
c     Local
      integer i                 ! lon index
      integer k                 ! lev index
      integer m                 ! cst index 
c     Initialize dust tracer to zero
      do m=1,dst_nbr
         do k=1,plev
            do i=1,plon
               q_dst(i,k,m)=0.0 ! [kg kg-1]
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst
      return
      end                       ! end dst_zero()
      
      subroutine mmrd2mpc(q_H2O_vpr,q_dry,pdel,mpc_dst,mpc_dst_ttl)
c     Purpose: Given a dry mass mixing ratio q_dry, and a specific humidity q_H2O_vpr 
c     (i.e., moist mass mixing ratio of water vapor), mmrd2mpc computes column mass path.
c     mmrd2mpc returns mass path for each constituent separately, in array mpc_dst, 
c     and total mass path of all constituents in array mpc_dst_ttl.
c     NB: q_dry, but not q_H2O_vpr, is a dry mass mixing ratio.
c     mmrd2mpc() is called by CCM:dynamics/eul/linemsbc(), MATCH:src/physlic()

      use precision

      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
c     Input
      real(r8) pdel(plond,plev)     ! [Pa] pdel(k)=pint(k+1)-pint(k)
      real(r8) q_H2O_vpr(plond,plev)    ! [kg H2O vapor kg-1 moist air] Specific humidity
      real(r8) q_dry(plond,plev,dst_nbr) ! [kg tracer kg-1 dry air] Tracer dry mass mixing ratio
c     Output
      real(r8) mpc_dst_ttl(plond)   ! [kg m-2] Total column mass path of dust
      real(r8) mpc_dst(plond,dst_nbr) ! [kg m-2] Column mass path of dust
c     Local
      integer i                 ! lon index
      integer k                 ! lev index
      integer m                 ! cst index 
c     Initialize column mass to zero
      call vec_set(mpc_dst,plond*dst_nbr,0.0_r8) ! [kg m-2]
      call vec_set(mpc_dst_ttl,plond,0.0_r8) ! [kg m-2]
      do m=1,dst_nbr
         do k=1,plev
            do i=1,plon
               mpc_dst(i,m)=mpc_dst(i,m)+q_dry(i,k,m)*(1.0-q_H2O_vpr(i,k))*pdel(i,k)*grv_sfc_rcp
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst
c     Add all dust constituents to obtain total
      call dst_add_lon(mpc_dst,mpc_dst_ttl)
      return
      end                       ! end mmrd2mpc()

      subroutine dst_add_lon_lev(q,q_ttl)
c     Purpose: Given a 3-D array of an additive property (e.g., mixing ratio, flux)
c     which is dimensioned (plond,plev,dst_nbr), dst_add_lon_lev() computes 
c     and returns the total property (e.g., mixing ratio, flux), obtained 
c     by simply adding along the (dust) constituent dimension.

      use precision

      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
c     Input
      real(r8) q(plond,plev,dst_nbr) ! Total property
c     Output
      real(r8) q_ttl(plond,plev)    ! Property for each size class
c     Local
      integer i                 ! lon index
      integer k                 ! lev index
      integer m                 ! cst index 
c     Initialize total property to zero then integrate
      call vec_set(q_ttl,plond*plev,0.0_r8)
      do m=1,dst_nbr
         do k=1,plev
            do i=1,plon
               q_ttl(i,k)=q_ttl(i,k)+q(i,k,m)
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst
      return
      end                       ! end dst_add_lon_lev()

      subroutine dst_add_lon(q,q_ttl)
c     Purpose: Given a 3-D array of an additive property (e.g., mixing ratio, flux)
c     which is dimensioned (plond,dst_nbr), dst_add_lon() computes 
c     and returns the total property (e.g., mixing ratio, flux), obtained 
c     by simply adding along the (dust) constituent dimension.

      use precision

      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
c     Input
      real(r8) q(plond,dst_nbr)     ! Total property
c     Output
      real(r8) q_ttl(plond)         ! Property for each size class
c     Local
      integer i                 ! lon index
      integer m                 ! cst index 
c     Initialize total property to zero then integrate
      call vec_set(q_ttl,plond,0.0_r8)
      do m=1,dst_nbr
         do i=1,plon
            q_ttl(i)=q_ttl(i)+q(i,m)
         end do                 ! end loop over lon
      end do                    ! end loop over cst
      return
      end                       ! end dst_add_lon()

      subroutine dst_pth(
     $     prs_dlt,             ! I [Pa] Pressure thickness
     $     q_dst,               ! I [kg kg-1] Dust mixing ratio
     $     mpc_dst,             ! O [kg m-2] Column mass path of dust
     $     mpc_dst_ttl,         ! O [kg m-2] Total column mass path of dust
     $     mpl_dst,             ! O [kg m-2] Layer dust amount
     $     mpp_dst,             ! O [kg m-2] Dust path above kth interface level
     $     odxc_dst,            ! O [frc] Column dust optical depth
     $     odxc_dst_ttl)        ! O [frc] Total column dust optical depth
c     Purpose: Compute dust paths needed in radiation routines and for outfld()
c     dst_pth() is called by CCM:physics/tphysbc()
c     fxm: mpp_dst requires a lot of memory. Perhaps mpp_dst should get its own routine called only when doabsems == true.?

      use precision

      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstcst.h> /* Physical constants for dust routines */
#include <dstodx.h> /* Optical depth information */
c     Input
      real(r8) prs_dlt(plond,plev)     ! I [Pa] Pressure thickness
      real(r8) q_dst(plond,plev,dst_nbr) ! I [kg kg-1] Dust mixing ratio
c     Output
      real(r8) mpc_dst(plond,dst_nbr) ! O [kg m-2] Column mass path of dust
      real(r8) mpc_dst_ttl(plond)   ! O [kg m-2] Total column mass path of dust
      real(r8) mpl_dst(plond,plev,dst_nbr) ! O [kg m-2] Layer dust amount
      real(r8) mpp_dst(plond,plevp,dst_nbr) ! O [kg m-2] Dust path above kth interface level
      real(r8) odxc_dst(plond,dst_nbr) ! O [frc] Column dust optical depth
      real(r8) odxc_dst_ttl(plond)  ! O [frc] Total column dust optical depth
c     Local
      integer i                 ! lon index
      integer k                 ! lev index
      integer m                 ! cst index 
      real(r8) mpl_air(plond,plev)  ! [kg m-2] Air mass path in layer
      
c     Compute necessary derived fields
      do k=1,plev
         do i=1,plon
c     Compute mass of air currently in gridbox
            mpl_air(i,k)=prs_dlt(i,k)*grv_sfc_rcp ! [kg m-2]
         end do                 ! end loop over lon
      end do                    ! end loop over lev
      
c     Initialize dust path to zero
      call vec_set(mpc_dst,plond*dst_nbr,0.0_r8) ! [kg m-2]
      call vec_set(mpp_dst,plond*plevp*dst_nbr,0.0_r8) ! [kg m-2]

c     Compute local and column integrated dust paths
      do m=1,dst_nbr
         do k=1,plev
            do i=1,plon
               mpl_dst(i,k,m)=q_dst(i,k,m)*prs_dlt(i,k)*grv_sfc_rcp ! [kg kg-1] --> [kg m-2]
               mpc_dst(i,m)=mpc_dst(i,m)+mpl_dst(i,k,m) ! [kg m-2]
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst

c     Compute partial integrated dust paths above interfaces
      do m=1,dst_nbr
         do k=1,plev            
            do i=1,plon
               mpp_dst(i,k+1,m)=mpp_dst(i,k,m)+mpl_dst(i,k,m) ! [kg m-2] Dust path above kth interface level
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst
      
c     Initialize integrated column totals
      call vec_set(odxc_dst_ttl,plond,0.0_r8) ! [frc]
      call vec_set(mpc_dst_ttl,plond,0.0_r8) ! [kg m-2]
c     Compute column optical depths
      do m=1,dst_nbr
         do i=1,plon
            mpc_dst_ttl(i)=mpc_dst_ttl(i)+mpc_dst(i,m)
            odxc_dst(i,m)=mpc_dst(i,m)*ext_cff_mss_dst_dgn(m)
            odxc_dst_ttl(i)=odxc_dst_ttl(i)+odxc_dst(i,m)
         end do                 ! end loop over lon
      end do                    ! end loop over cst
      return
      end                       ! end dst_pth()

      subroutine mpc2odxc(mpc_dst,odxc_dst,odxc_dst_ttl)
c     Given a column mass path mpc_dst [kg m-2], mpc2odxc() computes and returns
c     column optical depths at the specified channel separately, in array odxc_dst,
c     and the total column optical depth from all constituents in array odxc_dst_ttl.
c     mpc2odxc() is called by CCM:dynamics/eul/linemsbc(), MATCH:src/physlic()

      use precision

      implicit none
#include <pmgrid.h>
#include <dstgrd.h> /* Dust grid sizes */
#include <dstodx.h> /* Optical depth information */
c     Input
      real(r8) mpc_dst(plond,dst_nbr) ! [kg m-2] Column mass path of dust
c     Output
      real(r8) odxc_dst(plond,dst_nbr) ! [frc] Column dust optical depth
      real(r8) odxc_dst_ttl(plond) ! [frc] Total column dust optical depth
c     Local
      integer i                 ! [idx] lon index
      integer m                 ! [idx] cst index 

c     Initialize total column optical depth to zero
      call vec_set(odxc_dst_ttl,plond,0.0_r8) ! [frc]
c     Compute individual optical depths of each size class
c     Sum them to obtain total optical depth
      do m=1,dst_nbr
         do i=1,plon
            odxc_dst(i,m)=mpc_dst(i,m)*ext_cff_mss_dst_dgn(m) ! [frc] 
            odxc_dst_ttl(i)=odxc_dst_ttl(i)+odxc_dst(i,m) ! [frc] 
         end do                 ! end loop over lon
      end do                    ! end loop over cst
      return
      end                       ! end mpc2odxc()

      subroutine mmrm2mmrd(
     $     plev,                ! I [nbr] Number of levels
     $     plon,                ! I [nbr] Number of longitudes
     $     plond,               ! I [nbr] First dimension of arrays
     $     q_H2O_vpr,           ! I [kg H2O vapor kg-1 moist air] Specific humidity
     $     q_mst,               ! I [kg tracer kg-1 moist air] Tracer moist mass mixing ratio
     $     trc_nbr,             ! I [nbr] Number of tracers
     $     q_dry)               ! O [kg tracer kg-1 dry air] Tracer dry mass mixing ratio
c     Purpose: Convert moist mass mixing ratios to dry mass mixing ratios
c     Given moist mass mixing ratio q_mst [kg tracer kg-1 moist air], 
c     and a specific humidity q_H2O_vpr [kg H2O vapor kg-1 moist air], 
c     mmrm2mmrd() converts q_mst to dry mass mixing ratio q_dry [kg tracer kg-1 dry air]. 
c     mmrm2mmrd() is called by CCM:dynamics/eul/linemsbc(), MATCH:src/physlic()
c     mmrm2mmrd() is inverse of mmrd2mmrm()

      use precision

      implicit none
c     Input
      integer,intent(in)::plond ! I [nbr] First dimension of arrays
      integer,intent(in)::plon  ! I [nbr] Number of longitudes
      integer,intent(in)::plev  ! I [nbr] Number of levels
      integer,intent(in)::trc_nbr ! I [nbr] Number of tracers
      real(r8),intent(in)::q_H2O_vpr(plond,plev) ! I [kg H2O vapor kg-1 moist air] Specific humidity
      real(r8),intent(in)::q_mst(plond,plev,trc_nbr) ! I [kg tracer kg-1 moist air] Tracer moist mass mixing ratio
c     Output
      real(r8),intent(out)::q_dry(plond,plev,trc_nbr) ! O [kg tracer kg-1 dry air] Tracer dry mass mixing ratio
c     Local
      integer lon_idx           ! [idx] Longitude index
      integer lev_idx           ! [idx] Level index
      integer trc_idx           ! [idx] Tracer index 
      do trc_idx=1,trc_nbr
         do lev_idx=1,plev
            do lon_idx=1,plon
               q_dry(lon_idx,lev_idx,trc_idx)=q_mst(lon_idx,lev_idx,trc_idx)/(1.0-q_H2O_vpr(lon_idx,lev_idx))
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst
      return
      end                       ! end mmrm2mmrd()

      subroutine mmrd2mmrm(
     $     plev,                ! I [nbr] Number of levels
     $     plon,                ! I [nbr] Number of longitudes
     $     plond,               ! I [nbr] First dimension of arrays
     $     q_H2O_vpr,           ! I [kg H2O vapor kg-1 moist air] Specific humidity
     $     q_dry,               ! I [kg tracer kg-1 dry air] Tracer dry mass mixing ratio
     $     trc_nbr,             ! I [nbr] Number of tracers
     $     q_mst)               ! O [kg tracer kg-1 moist air] Tracer moist mass mixing ratio
c     Purpose: Convert dry mass mixing ratios to moist mass mixing ratios
c     Given dry mass mixing ratio q_dry [kg tracer kg-1 dry air], 
c     and a specific humidity q_H2O_vpr [kg H2O vapor kg-1 moist air], 
c     mmrd2mmrm() converts q_dry to moist mass mixing ratio q_mst [kg tracer kg-1 moist air]. 
c     mmrd2mmrm() is called by CCM:dynamics/eul/linemsbc(), MATCH:src/physlic()
c     mmrd2mmrm() is inverse of mmrm2mmrd()

      use precision

      implicit none
c     Input
      integer,intent(in)::plond ! I [nbr] First dimension of arrays
      integer,intent(in)::plon  ! I [nbr] Number of longitudes
      integer,intent(in)::plev  ! I [nbr] Number of levels
      integer,intent(in)::trc_nbr ! I [nbr] Number of tracers
      real(r8),intent(in)::q_H2O_vpr(plond,plev) ! I [kg H2O vapor kg-1 moist air] Specific humidity
      real(r8),intent(in)::q_dry(plond,plev,trc_nbr) ! I [kg tracer kg-1 dry air] Tracer dry mass mixing ratio
c     Output
      real(r8),intent(out)::q_mst(plond,plev,trc_nbr) ! O [kg tracer kg-1 moist air] Tracer moist mass mixing ratio
c     Local
      integer lon_idx           ! [idx] Longitude index
      integer lev_idx           ! [idx] Level index
      integer trc_idx           ! [idx] Tracer index 
      do trc_idx=1,trc_nbr
         do lev_idx=1,plev
            do lon_idx=1,plon
               q_mst(lon_idx,lev_idx,trc_idx)=q_dry(lon_idx,lev_idx,trc_idx)*(1.0-q_H2O_vpr(lon_idx,lev_idx))
            end do              ! end loop over lon
         end do                 ! end loop over lev
      end do                    ! end loop over cst
      return
      end                       ! end mmrd2mmrm()

