c $Id$
      common /mass2nc0/ ncid, nrec

      integer
     $  ncid       ! netCDF file ID
     $, nrec       ! current record being written
