/*$Id: ccmhdr.h,v 1.1.1.1 1998/02/17 23:44:00 eaton Exp $*/

#ifndef CCMHDR_H
#define CCMHDR_H

#include "version.h"
#include "fortio.h"

typedef struct Field {
  char name[9];   /* Name of field (8 chars + a null byte) */
  Hint lev;       /* Type of levels, 0 = single level field, 1 = multilevel
                     field at layer interfaces, 2 = multilevel fields at
		     layer midpoints */
  Hint tim;       /* Time representation, 0 = instantaneous, 1 = averaged,
		     2 = minimum value at each gridpoint since last write,
		     3 = maximum value at each gridpoint since last write */
  Hint off;       /* Offset for field in data record */
  Hint pdens;     /* Packing density (1, 2, or 4) */
  char units[9];  /* Physical units of field */
  char* long_name; /* descriptive name for field */
  char* udunit;   /* units conforming to UNIDATA's udunits specification */
  Hint spval_flg; /* Special value flag.  0 if not used, 1 otherwise. */
} Field;

typedef struct CCM_HDR {
  Hint *ih;       /* Storage for internal representation of integer header */
  char *ch;       /* Storage for internal representation of char header */
  Hfloat *fh;     /* Storage for internal representation of float header */
  Hint lenhdi;    /* Length of header record 1 */
  Hint mftyp;     /* Format code for this history tape */
  Hint mfilh;     /* Time sample relative to beginning of this tape */
  Hint mfilth;    /* Maximum number of time samples on a tape */
  Hint nrbd;      /* Number of records before data records */
  Hint maxsiz;    /* Length of data record */
  Hint ndavu;     /* Length of data record after unpacking */
  Hint nspher;    /* Horizonal domain flag: 1 for global, 2 for hemispheric */
  Hint nlon;      /* Number of longitude points in model grid */
  Hint nlonw;     /* Number of longitude data values written in data records */
  Hint norec;     /* Number of latitudes.  Also number of data records */
  Hint nlev;      /* Number of vertical levels */
  Hint ntrm;      /* Spectral truncation parameter, M */
  Hint ntrn;      /* Spectral truncation parameter, N */
  Hint ntrk;      /* Spectral truncation parameter, K */
  Hint nfldh;     /* Number of fields */
  Hint nsteph;    /* Time step number for this time sample */
  Hint nstprh;    /* Time step number at start of run (restarts are new runs) */
  Hint nitslf;    /* Number of time steps since previous time sample written */
  Hint ndbase;    /* Base day number */
  Hint nsbase;    /* Base number of seconds relative to ndbase */
  Hint ndcur;     /* Current day number (corresponds with nsteph) */
  Hint nscur;     /* Current number of seconds relative to ndcur */
  Hint nbdate;    /* Base date */
  Hint nbsec;     /* Base number of seconds relative to nbdate */
  Hint ncdate;    /* Current date (corresponds with nsteph) */
  Hint ncsec;     /* Current number of seconds relative to ncdate */
  Hint mdt;       /* Model timestep in seconds */
  Hint mhisf;     /* Number of timesteps between time sample writes */
  Hint mfstrt;    /* Type of start: 0=initial run, 1=restart, 
                     2=regeneration, 3=branch run */
  Hint lenhdc;    /* Length of character header in 8 byte words */
  Hint lenhdr;    /* Length of real header (number of floating point values) */
  Hfloat *mphyab; /* Pointer to hybrid A+B coords (2*nlev+1 values) */
  Hfloat *mphya;  /* Pointer to hybrid A coords (2*nlev+1 values) */
  Hfloat *mphyb;  /* Pointer to hybrid B coords (2*nlev+1 values) */
  Hfloat *mplat;  /* Pointer to latitude values (norec values) */
  Hfloat *mpwts;  /* Pointer to Gaussian weights (norec values) */

  Hfloat* etamid; /* hybrid A+B at layer midpoints */
  Hfloat* lon;    /* longitude values */

  char mcase[9];
  char mcstit[81];
  char lnhstc[81];
  char ldhstc[9];
  char lthstc[9];
  char lshstc[9];
  char lnhstp[81];
  char ldhstp[9];
  char lthstp[9];
  char lshstp[9];
  char lnhstf[81];
  char ldhstf[9];
  char lthstf[9];
  char lshstf[9];
  char lnhsti[81];
  char ldhsti[9];
  char lthsti[9];
  char lshsti[9];
  char lnhsta[81];
  char ldhsta[9];
  char lthsta[9];
  char lshsta[9];
  char lnhstt[81];
  char ldhstt[9];
  char lthstt[9];
  char lshstt[9];
  char lnhsts[81];
  char ldhsts[9];
  char lthsts[9];
  char lshsts[9];
  char lnhsto[81];
  char ldhsto[9];
  char lthsto[9];
  char lshsto[9];
  Field *fld;
  int ts_count;    /* Count the number of time samples read.  This counter is */
                   /* continuous across files.                                */
} CCM_HDR;

/* long_name and units for standard CCM history fields */

static struct ccm_names {
  char* std_name;
  char* std_unit;
  char* long_name;
  char* udunit;
} ccm_name_tab[] = {
  "PHIS",     "M2/S2",    "surface geopotential",                           "m2/s2",
  "PS",       "PA",       "surface pressure",                               "Pa",
  "T",        "K",        "temperature",                                    "K",
  "U",        "M/S",      "zonal wind component",                           "m/s",
  "V",        "M/S",      "meridional wind component",                      "m/s",
  "ETADOT",   "1/S",      "vertical (eta) velocity",                        "s-1",
  "Q",        "KG/KG",    "specific humidity",                              "kg/kg",
  "HA01",     "KG/KGS",   "horizontal advective tendency of water vapor",   "kg/kg/s",
  "VA01",     "KG/KGS",   "vertical advective tendency of water vapor",     "kg/kg/s",
  "DF01",     "KG/KGS",   "SLT tendency of water vapor",                    "kg/kg/s",
  "TA01",     "KG/KGS",   "total advective tendency of water vapor",        "kg/kg/s",
  "VD01",     "KG/KGS",   "vertical diffusion tendency of water vapor",     "kg/kg/s",
  "DC01",     "KG/KGS",   "convective adjustment tendency of water vapor",  "kg/kg/s",
  "TE01",     "KG/KGS",   "total tendency of water vapor",                  "kg/kg/s",
  "SS01",     "KG/KGS",   "source/sink tendency of water vapor",            "kg/kg/s",
  "SF01",     "KG/M2/S",  "surface flux of water vapor",                    "kg/m2/s",
  "DUH",      "K/S",      "U horizontal diffusive heating",                 "K/s",
  "DVH",      "K/S",      "V horizontal diffusive heating",                 "K/s",
  "DTH",      "K/S",      "T horizontal diffusive heating",                 "K/s",
  "NDQH",     "KG/KGS",   "Q horizontal diffusion",                         "kg/kg/s",
  "ORO",      "FLAG",     "ocean (0), land (1), sea ice (2) flag",          NULL,
  "WET",      "M",        "Soil moisture",                                  "m",
  "SNOWH",    "M",        "Water equivalent snow depth",                    "m",
  "PRECL",    "M/S",      "Large-scale (stable) precipitation rate",        "m/s",
  "PRECC",    "M/S",      "Convective precipitation rate",                  "m/s",
  "PRECSL",   "M/S",    "Large-scale (stable) snow rate (water equivalent)","m/s",
  "PRECSC",   "M/S",      "Convective snow rate (water equivalent)",        "m/s",
  "SHFLX",    "W/M2",     "Surface sensible heat flux",                     "W/m2",
  "LHFLX",    "W/M2",     "Surface latent heat flux",                       "W/m2",
  "QFLX",     "KG/M2/S",  "Surface water flux",                             "kg/m2/s",
  "PBLH",     "M",        "PBL height",                                     "m",
  "USTAR",    "M/S",      "Surface friction velocity",                      "m/s",
  "TREFHT",   "K",        "Reference height temperature",                   "K",
  "Q10",      "KG/KG",    "10-meter specific humidity",                     "kg/kg",
  "CGH",      "K/M",      "Counter-gradient term for heat in PBL",          "K/m",
  "CGQ",      "1/M",      "Counter-gradient term for moisture in PBL",      "m-1",
  "CGS",      "S/M2", "Counter-gradient coefficient on surface kinematic fluxes","s/m2",
  "TPERT",    "K",        "Perturbation temperature (eddies in PBL)",       "K",
  "QPERT",    "KG/KG",    "Perturbation specific humidity (eddies in PBL)", "kg/kg",
  "KVH",      "M2/S",     "Vertical diffusion diffusivities (heat/moisture)","m2/s",
  "KVM",      "M2/S",     "Vertical diffusion diffusivities (momentum)",    "m2/s",
  "DUV",      "M/S2",     "U vertical diffusion",                           "m/s2",
  "DVV",      "M/S2",     "V vertical diffusion",                           "m/s2",
  "DTV",      "K/S",      "T vertical diffusion",                           "K/s",
  "FSNS",     "W/M2",     "Net solar flux at surface",                      "W/m2",
  "FLNS",     "W/M2",     "Net longwave flux at surface",                   "W/m2",
  "FLNT",     "W/M2",     "Net longwave flux at top",                       "W/m2",
  "FSNT",     "W/M2",     "Net solar flux at top",                          "W/m2",
  "CLOUD",    "FRACTION", "Cloud fraction",                                 NULL,
  "SETLWP",   "G/M2",     "Specified liquid water path lengths",            "gram/m2",
"CLDLWP","G/M2","Actual cloud liquid water path length (account for cloud fraction)","gram/m2",
  "EFFCLD",   "FRACTION", "Effective cloud fraction",                       NULL,
  "FLNTC",    "W/M2",     "Clearsky net longwave flux at top",              "W/m2",
  "FSNTC",    "W/M2",     "Clearsky net solar flux at top",                 "W/m2",
  "FLNSC",    "W/M2",     "Clearsky net longwave flux at surface",          "W/m2",
  "FSNSC",    "W/M2",     "Clearsky net solar flux at surface",             "W/m2",
  "OMEGA",    "PA/S",     "Vertical pressure velocity",                     "Pa/s",
  "DQP",      "KG/KGS",   "Specific humidity tendency due to precipitation","kg/kg/s",
  "TAUX",     "N/M2",     "X-component (east-west) of surface stress",      "N/m2",
  "TAUY",     "N/M2",     "Y-component (north-south) of surface stress",    "N/m2",
  "SRFRAD",   "W/M2",     "Net radiative flux at surface",                  "W/m2",
  "QRS",      "K/S",      "Solar heating rate",                             "K/s",
  "QRL",      "K/S",      "Longwave heating rate",                          "K/s",
  "CLDTOT",   "FRACTION", "Vertically-integrated, random overlap, total cloud amount",    NULL,
  "CLDLOW",   "FRACTION", "Vertically-integrated, random overlap, low cloud amount",      NULL,
  "CLDMED",   "FRACTION", "Vertically-integrated, random overlap, mid-level cloud amount",NULL,
  "CLDHGH",   "FRACTION", "Vertically-integrated, random overlap, high cloud amount",     NULL,
  "TOTLWP",   "FRACTION", "Vertically-integrated actual liquid water path length",        NULL,
  "LWSH",     "M",        "Liquid water scale height",                      "m",
  "TS1",      "K",        "Surface temperature (level 1)",                  "K",
  "TS2",      "K",        "Subsurface temperature (level 2)",               "K",
  "TS3",      "K",        "Subsurface temperature (level 3)",               "K",
  "TS4",      "K",        "Subsurface temperature (level 4)",               "K",
  "TS",       "K",        "Surface temperature",                            "K",
  "SOLIN",    "W/M2",     "Solar insolation",                               "W/m2",
  "UTEND",    "M/S2",     "U tendency",                                     "m/s2",
  "VTEND",    "M/S2",     "V tendency",                                     "m/s2",
  "TTEND",    "K/S",      "T tendency",                                     "K/s",
  "LPSTEN",   "PA/S",     "Surface pressure tendency",                      "Pa/s",
  "UTGW",     "M/S2",     "U tendency - gravity wave drag",                 "m/s2",
  "VTGW",     "M/S2",     "V tendency - gravity wave drag",                 "m/s2",
  "TAUGWX",   "N/M2",     "East-west gravity wave drag surface stress",     "N/m2",
  "TAUGWY",   "N/M2",     "North-south gravity wave drag surface stress",   "N/m2",
  "DTCOND",   "K/S",      "T tendency - convective adjustment",             "K/s",
  "CMFDT",    "K/S",      "T tendency - moist convetion",                   "K/s",
  "CMFDQ",    "KG/KGS",   "Q tendency - moist convection",                  "kg/kg/s",
  "CMFDQF",   "KG/KGS",   "Q tendency - moist convective rainout",          "kg/kg/s",
  "CMFMC",    "KG/M2S",   "Moist convection mass flux",                     "kg/m2/s",
  "CMFSL",    "W/M2",     "Moist convection liquid water static energy flux","W/m2",
  "CMFLQ",    "W/M2",     "Moist convection total water flux",              "W/m2",
  "CNVCLD",   "FRACTION", "Random overlap total convective cloud amount",   NULL,
  "VT",       "KM/S",     "Meridional heat transport",                      "K m/s",
  "VZ",       "M3/S3",    "Meridional transport",                           "m3/s3",
  "VQ",       "M/S",      "Meridional water transport",                     "m/s",
  "VVPUU",    "M2/S2",    "Kinetic Energy",                                 "m2/s2",
  "ALB",      "FRACTION", "Albedo",                                         NULL,
  "ALBCLR",   "FRACTION", "Clear sky albedo",                               NULL,
  "RELHUM",   "FRACTION", "Relative humidity",                              NULL,
  "Z2",       "M",        "Geopotential Height (above sea level)",          "m",
  "MQ",       "KG/M2",    "water mass",                                     "kg/m2",
  "PSL",      "PASCALS",  "Sea level pressure",                             "Pa",
  "OMEGAUP",  "PA/S",     "Average of Omega when omega is < 0 (up)",        "Pa/s",
  "NUMOMGUP", "FRACTION", "Percentage of time omega is < 0 (up)",           NULL,
  "CLOUDUP",  "FRACTION", "Average of Cloud when omega is < 0 (up)",        NULL,
  "DPSLON",   "M/S2",     "Longitude Pressure Gradient",                    "m/s2",
  "DPSLAT",   "M/S2",     "Latitude Pressure Gradient",                     "m/s2",
  NULL,NULL,NULL,NULL
};


/*
Function declarations:
*/

extern int init_hdr( FFILE ff, CCM_HDR *hdr );
extern int next_hdr( FFILE ff, CCM_HDR *hdr );

#endif
