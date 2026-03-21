c     $Header$
      
c     Purpose: dstszdstbd.F contains block data to initialize mineral dust 
c     size distributions
c     dstszdstbd is currently initialized in CCM:control/ccm3(), MATCH:src/main()
      
#include <params.h>

      block data dstszdstbd
c     Purpose: dstszdstbd initializes aerosol microphysical properties known at compile time
      use precision

      implicit none
#include <dstgrd.h> /* Dust grid sizes */
#include <dstszdst.h> /* Diagnostic aerosol microphysical properties */
      
c     Set properties of source soil distribution
c     Following parameters are from She84 p. 75 Table 1 with gsd = 2.0 modification from SBG98 p. 10581 Table 1
c     These are Shettle's "Background Desert Model" modes
      data dmt_vma_src / 0.0111e-6,  2.524e-6, 42.10e-6 / ! [m] Mass median diameter She84 p. 75 Table 1 
      data gsd_anl_src / 1.89     ,  2.0     ,  2.13    / ! [frc] Geometric standard deviation She84 p. 75 Table 1 
      data mss_frc_src / 2.6e-4   ,  0.781   ,  0.219   / ! [frc] Mass fraction She84 p. 75 Table 1 
      
c     Following parameters are from BSM96 p. 73 Table 2 
c     These are D'Almeida's (1987) "Background" modes
c      data dmt_vma_src / 0.832e-6 ,  4.82e-6 , 19.38e-6 / ! [m] Mass median diameter BSM96 p. 73 Table 2
c      data gsd_anl_src / 2.10     ,  1.90    ,  1.60    / ! [frc] Geometric standard deviation BSM96 p. 73 Table 2 
c      data mss_frc_src / 0.036    ,  0.957   ,  0.007   / ! [frc] Mass fraction BSM96 p. 73 Table 2 
      
c     Following parameters are from PaG77 p. 2080 Table 1 (converted from radius to diameter)
c     Source soils assumed to have three modes dubbed, in increasing size order, modes C, A, and B
c     Mode A is mineral dust transport mode, seen in source regions and downwind
c     Mode B is seen in the source soil itself, and in the atmosphere during dust events
c     Mode C is seen most everywhere, but does not usually correlate with local dust amount
c     Mode C is usually a global, aged, background, anthropogenic aerosol, typically rich in sulfate and black carbon
c     Sometimes, however, mode C has a mineral dust component
c     Modes C and B are averages from PaG77 Table 1 p. 2080
c     Mode A is based on summary recommendation: rds_mdn_sfc = 1.5 and gsd = 2.2
c     mss_frc_src is unavailable in PaG77 because N_0 is not given
c     Use mss_frc_src from BSM96 p. 73 Table 2 instead
c     PaG77 nomenclature:  / Mode C   ,  Mode A  ,  Mode B  /
c     data dmt_vma_src / 0.27e-6  ,  5.6e-6  ,  57.6e-6 / ! [m] Mass median diameter PaG77 p. 2080 Table 1 
c     data gsd_anl_src / 1.88     ,  2.2     ,  1.62    / ! [frc] Geometric standard deviation PaG77 p. 2080 Table 1 
c     data mss_frc_src / 0.036    ,  0.957   ,  0.007   / ! [frc] Mass fraction BSM96 p. 73 Table 2 
      
      end                       ! end dstszdstbd
      

