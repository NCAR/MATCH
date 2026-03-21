#ifndef PARAMS_H
#define PARAMS_H

#define PLON  192
#define PLAT  94
#define PLEV  28


#define CONVTRAN

#define DI_VDIFF 
#define DI_CONV_CCM
#define DI_CLOUD_PHYS

#define SCYC_NBR     4
#define SCYC_IDX_SRT 1
#define DST_NBR 4 
#define DST_IDX_SRT 5
#define CAER_NBR     4
#define CAER_IDX_SRT 9
#define PCNST  SCYC_NBR+DST_NBR+CAER_NBR+1
#define MODIS




#ifdef SCYC
#ifdef DST
#ifdef CAER
#define AER
#endif   
#endif   
#endif


#define PNATS  0

#endif
