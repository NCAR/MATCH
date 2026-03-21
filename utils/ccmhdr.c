/*$Id: ccmhdr.c,v 1.2 1998/03/30 19:59:41 eaton Exp $*/

/* Initialize and update the ccm header structure from ccm history tape
   format files.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "ccmhdr.h"
#include "fortio.h"

#define NINTS 32    /* no. of ints read from header to get size info */

void charr2str( char* in, int nc1 );

/* Offsets for ccm1 and ccm2 integer header */
static const int lenhdi=0, mftyp=1, mfilh=2, mfilth=3, nrbd=4, maxsiz=5,
            ndavu=6,
            nspher=7, nlon=8, nlonw=9, norec=10, nlev=11, ntrm=12, ntrn=13,
            ntrk=14, nfldh=15, nsteph=16, nstprh=17, nitslf=18, ndbase=19,
            nsbase=20, ndcur=21, nscur=22, nbdate=23, nbsec=24, ncdate=25,
            ncsec=26, mdt=27, mhisf=28, mfstrt=29, lenhdc=30, lenhdr=31,
            mpsig=32, mplat=33, mpwts=34, mpflds=35, mpcfld=36;

/*===========================================================================*/


int init_hdr( FFILE ff, CCM_HDR *hdr )
{
  /*
    Initialize the header structure.  The information from the first 3 header
    records of a CCM1 or CCM2 history tape are stored in the header structure.
    The input file is left positioned at the end of the 3rd record.

    Input arguments:
  
    ff      Input stream.

    Output arguments:

    hdr     Initialized header structure.

    Return values:

     0  Successful return.
    -1  Error reading first 32 words of integer header.
    -2  Error allocating memory for integer header record.
    -3  Error allocating memory for character header record.
    -4  Error allocating memory for real header record.
    -5  Error allocating memory for data requiring conversion.
    -6  Error reading integer header.
    -7  Error reading character header.
    -8  Error reading float header.
    -11 Error allocating memory for field structures.
  */

  /* Should get rid of this declaration in all cases */
#ifndef RS6000
  /* Offsets for ccm1 and ccm2 integer header */
  extern const int lenhdi, mftyp, mfilh, mfilth, nrbd, maxsiz, ndavu,
                   nspher, nlon, nlonw, norec, nlev, ntrm, ntrn,
                   ntrk, nfldh, nsteph, nstprh, nitslf, ndbase,
                   nsbase, ndcur, nscur, nbdate, nbsec, ncdate,
                   ncsec, mdt, mhisf, mfstrt, lenhdc, lenhdr,
                   mpsig, mplat, mpwts, mpflds, mpcfld;
#endif

  /* Character header info */
  int mcase, mcstit,
      lnhstc, ldhstc, lthstc, lshstc,  /* Current history tape */
      lnhstp, ldhstp, lthstp, lshstp,  /* Previous history tape */
      lnhstf, ldhstf, lthstf, lshstf,  /* First history tape of case */
      lnhsti, ldhsti, lthsti, lshsti,  /* Initial history tape */
      lnhsta, ldhsta, lthsta, lshsta,  /* Analysed history tape */
      lnhstt, ldhstt, lthstt, lshstt,  /* Time-invariant boundary dataset */
      lnhsts, ldhsts, lthsts, lshsts,  /* Time-variant SST dataset */
      lnhsto, ldhsto, lthsto, lshsto,  /* Ozone dataset */
      mcflds;

  Hint ihbuf[NINTS*2];   /* storage for NINTS 8-byte integers */
  Hint ihtmp[NINTS];     /* storage for integer conversion if necessary */

  int ihlen,   /* length of integer header (no. of int) */
      chlen,   /* length of character header (no. of char[8]) */
      fhlen,   /* length of float header (no. of floats) */
      nbytes;
  void *buf;   /* storage for data that requires conversion */
  int i, n, st, itmp;
  int fo;
  char *cptr;
  Hfloat dlon;
/*---------------------------------------------------------------------------*/

  /* Read enough of the first integer header record to determine the lengths
     of the headers. */

  if ( need_icv( ff ) ) {
    if ( rdcv_irec( ff, NINTS, ihbuf, ihtmp ) != NINTS ) return -1;
  } else {
    if ( read_rec( ff, ihtmp, NINTS*ff.isz ) != NINTS*ff.isz ) return -1;
  }

  rewind_file( ff );

  ihlen = ihtmp[lenhdi];
  chlen = ihtmp[lenhdc];
  fhlen = ihtmp[lenhdr];

  /* Allocate space for header records. */

  if ( (hdr->ih = (Hint *)malloc( ihlen*sizeof(Hint) )) == NULL ) return -2;
  if ( (hdr->ch = (char *)malloc( chlen*8 )) == NULL ) return -3;
  if ( (hdr->fh = (Hfloat *)malloc( fhlen*sizeof(Hfloat) )) == NULL ) return -4;

  /* Allocate space for data conversion. */

  if ( need_icv( ff ) || need_fcv( ff ) ) {
    nbytes = (ihlen*ff.isz > fhlen*ff.fsz) ? ihlen*ff.isz : fhlen*ff.fsz;
    if ( (buf = malloc( nbytes )) == NULL ) return -5;
  }

  /* Read header records. */

  if ( need_icv( ff ) ) {
    if ( rdcv_irec( ff, ihlen, buf, hdr->ih ) != ihlen ) return -6;
  } else {
    if ( read_rec( ff, hdr->ih, ihlen*ff.isz ) != ihlen*ff.isz ) return -6;
  }

  if ( read_rec( ff, hdr->ch, chlen*8 ) != chlen*8 ) return -7;

  if ( need_fcv( ff ) ) {
    if ( (itmp=rdcv_frec( ff, fhlen, buf, hdr->fh )) != fhlen ) {
      if ( itmp == -3 ) {
	printf("\ninit_hdr: ***warning*** out of range floats in header\n");
      } else {
	printf("init_hdr: %d floats requested, %d floats returned\n", fhlen,
	       itmp );
	return -8;
      }
    }
  } else {
    if ( read_rec( ff, hdr->fh, fhlen*ff.fsz ) != fhlen*ff.fsz ) return -8;
  }

  /* Initialize time sample counter */

  hdr->ts_count = 1;

  /* Initialize data from integer header record */

  hdr->lenhdi  = *(hdr->ih+lenhdi);
  hdr->mftyp   = *(hdr->ih+mftyp);
  hdr->mfilh   = *(hdr->ih+mfilh);
  hdr->mfilth  = *(hdr->ih+mfilth);
  hdr->nrbd    = *(hdr->ih+nrbd);
  hdr->maxsiz  = *(hdr->ih+maxsiz);
  hdr->ndavu   = *(hdr->ih+ndavu);
  hdr->nspher  = *(hdr->ih+nspher);
  hdr->nlon    = *(hdr->ih+nlon);
  hdr->nlonw   = *(hdr->ih+nlonw);
  hdr->norec   = *(hdr->ih+norec);
  hdr->nlev    = *(hdr->ih+nlev);
  hdr->ntrm    = *(hdr->ih+ntrm);
  hdr->ntrn    = *(hdr->ih+ntrn);
  hdr->ntrk    = *(hdr->ih+ntrk);
  hdr->nfldh   = *(hdr->ih+nfldh);
  hdr->nsteph  = *(hdr->ih+nsteph);
  hdr->nstprh  = *(hdr->ih+nstprh);
  hdr->nitslf  = *(hdr->ih+nitslf);
  hdr->ndbase  = *(hdr->ih+ndbase);
  hdr->nsbase  = *(hdr->ih+nsbase);
  hdr->ndcur   = *(hdr->ih+ndcur);
  hdr->nscur   = *(hdr->ih+nscur);
  hdr->nbdate  = *(hdr->ih+nbdate);
  hdr->nbsec   = *(hdr->ih+nbsec);
  hdr->ncdate  = *(hdr->ih+ncdate);
  hdr->ncsec   = *(hdr->ih+ncsec);
  hdr->mdt     = *(hdr->ih+mdt);
  hdr->mhisf   = *(hdr->ih+mhisf);
  hdr->mfstrt  = *(hdr->ih+mfstrt);
  hdr->lenhdc  = *(hdr->ih+lenhdc);
  hdr->lenhdr  = *(hdr->ih+lenhdr);

  /* Initialize pointers to data in real header record */

  hdr->mphyab  = hdr->fh + *(hdr->ih+mpsig) - 1;

  if ( *(hdr->ih+mplat) - *(hdr->ih+mpsig)  ==  3*(2*hdr->nlev+1) ){/* hybrid coordinates */
    hdr->mphya = hdr->mphyab + 2*hdr->nlev + 1;
    hdr->mphyb = hdr->mphya  + 2*hdr->nlev + 1;

    /* If top of model is at zero pressure, modify so that it isn't */
    if ( *hdr->mphyab == 0. ) {
      *hdr->mphyab = *(hdr->mphyab+1) / 10.;
      *hdr->mphya = *(hdr->mphya+1) / 10.;
      *hdr->mphyb = *(hdr->mphyb+1) / 10.;
    }

  } else {                                             /* sigma coordinates */
    hdr->mphya = (Hfloat *)NULL;
    hdr->mphyb = hdr->mphyab;
  }
  hdr->mplat   = hdr->fh + *(hdr->ih+mplat) - 1;
  hdr->mpwts   = hdr->fh + *(hdr->ih+mpwts) - 1;


  /* Allocate space for coordinate arrays */
  hdr->etamid = (Hfloat*)malloc( sizeof(Hfloat)*hdr->nlev );
  hdr->lon = (Hfloat*)malloc( sizeof(Hfloat)*hdr->nlon );

  /* Build coordinate arrays */
  for ( i = 0; i < hdr->nlev; ++i ) hdr->etamid[i] = *(hdr->mphyab + 2*i + 1);
  dlon = 360. / hdr->nlon;
  for ( i=0; i < hdr->nlon; i++ ) hdr->lon[i] = i * dlon;


  /* Initialize the offsets into the character header that depend on the
     history tape type. */

  if ( hdr->mftyp%10 == 2 ) {

    mcase=0, mcstit=8,
    lnhstc=88, ldhstc=168, lthstc=176, lshstc=184,
    lnhstp=192, ldhstp=272, lthstp=280, lshstp=288,
    lnhstf=296, ldhstf=376, lthstf=384, lshstf=392,
    lnhsti=400, ldhsti=480, lthsti=488, lshsti=496,
    lnhsta=504, ldhsta=584, lthsta=592, lshsta=600,
    mcflds=608;

  } else {

    if ( hdr->mftyp%10 != 3 ) {
      fprintf( stderr,
      "init_hdr: history tape type %d set in header\n"
      " continue processing assuming type 43\n", hdr->mftyp );
      hdr->mftyp = 43;
    }

    mcase=0, mcstit=8,
    lnhstc=88, ldhstc=168, lthstc=176, lshstc=184,
    lnhstf=192, ldhstf=272, lthstf=280, lshstf=288,
    lnhsti=296, ldhsti=376, lthsti=384, lshsti=392,
    lnhstt=400, ldhstt=480, lthstt=488, lshstt=496,
    lnhsts=504, ldhsts=584, lthsts=592, lshsts=600,
    lnhsto=608, ldhsto=688, lthsto=696, lshsto=704,
    mcflds=712;

  }

  /* Initialize data from character header record */

  memcpy( hdr->mcase, hdr->ch+mcase, 8 );
  memcpy( hdr->mcstit, hdr->ch+mcstit, 80 );
  memcpy( hdr->lnhstc, hdr->ch+lnhstc, 80 );
  memcpy( hdr->ldhstc, hdr->ch+ldhstc, 8 );
  memcpy( hdr->lthstc, hdr->ch+lthstc, 8 );
  memcpy( hdr->lshstc, hdr->ch+lshstc, 8 );
  memcpy( hdr->lnhstf, hdr->ch+lnhstf, 80 );
  memcpy( hdr->ldhstf, hdr->ch+ldhstf, 8 );
  memcpy( hdr->lthstf, hdr->ch+lthstf, 8 );
  memcpy( hdr->lshstf, hdr->ch+lshstf, 8 );
  memcpy( hdr->lnhsti, hdr->ch+lnhsti, 80 );
  memcpy( hdr->ldhsti, hdr->ch+ldhsti, 8 );
  memcpy( hdr->lthsti, hdr->ch+lthsti, 8 );
  memcpy( hdr->lshsti, hdr->ch+lshsti, 8 );
  charr2str( hdr->mcase, 8 );
  charr2str( hdr->mcstit, 80 );
  charr2str( hdr->lnhstc, 80 );
  charr2str( hdr->ldhstc, 8 );
  charr2str( hdr->lthstc, 8 );
  charr2str( hdr->lshstc, 8 );
  charr2str( hdr->lnhstf, 80 );
  charr2str( hdr->ldhstf, 8 );
  charr2str( hdr->lthstf, 8 );
  charr2str( hdr->lshstf, 8 );
  charr2str( hdr->lnhsti, 80 );
  charr2str( hdr->ldhsti, 8 );
  charr2str( hdr->lthsti, 8 );
  charr2str( hdr->lshsti, 8 );

  if ( hdr->mftyp%10 == 2 ) {

    memcpy( hdr->lnhstp, hdr->ch+lnhstp, 80 );
    memcpy( hdr->ldhstp, hdr->ch+ldhstp, 8 );
    memcpy( hdr->lthstp, hdr->ch+lthstp, 8 );
    memcpy( hdr->lshstp, hdr->ch+lshstp, 8 );
    memcpy( hdr->lnhsta, hdr->ch+lnhsta, 80 );
    memcpy( hdr->ldhsta, hdr->ch+ldhsta, 8 );
    memcpy( hdr->lthsta, hdr->ch+lthsta, 8 );
    memcpy( hdr->lshsta, hdr->ch+lshsta, 8 );
    charr2str( hdr->lnhstp, 80 );
    charr2str( hdr->ldhstp, 8 );
    charr2str( hdr->lthstp, 8 );
    charr2str( hdr->lshstp, 8 );
    charr2str( hdr->lnhsta, 80 );
    charr2str( hdr->ldhsta, 8 );
    charr2str( hdr->lthsta, 8 );
    charr2str( hdr->lshsta, 8 );

  } else if ( hdr->mftyp%10 == 3 ) {

    memcpy( hdr->lnhstt, hdr->ch+lnhstt, 80 );
    memcpy( hdr->ldhstt, hdr->ch+ldhstt, 8 );
    memcpy( hdr->lthstt, hdr->ch+lthstt, 8 );
    memcpy( hdr->lshstt, hdr->ch+lshstt, 8 );
    memcpy( hdr->lnhsts, hdr->ch+lnhsts, 80 );
    memcpy( hdr->ldhsts, hdr->ch+ldhsts, 8 );
    memcpy( hdr->lthsts, hdr->ch+lthsts, 8 );
    memcpy( hdr->lshsts, hdr->ch+lshsts, 8 );
    memcpy( hdr->lnhsto, hdr->ch+lnhsto, 80 );
    memcpy( hdr->ldhsto, hdr->ch+ldhsto, 8 );
    memcpy( hdr->lthsto, hdr->ch+lthsto, 8 );
    memcpy( hdr->lshsto, hdr->ch+lshsto, 8 );
    charr2str( hdr->lnhstt, 80 );
    charr2str( hdr->ldhstt, 8 );
    charr2str( hdr->lthstt, 8 );
    charr2str( hdr->lshstt, 8 );
    charr2str( hdr->lnhsts, 80 );
    charr2str( hdr->ldhsts, 8 );
    charr2str( hdr->lthsts, 8 );
    charr2str( hdr->lshsts, 8 );
    charr2str( hdr->lnhsto, 80 );
    charr2str( hdr->ldhsto, 8 );
    charr2str( hdr->lthsto, 8 );
    charr2str( hdr->lshsto, 8 );

  }

  /* Allocate memory and store field info in array of Field structures. */

  hdr->fld = (Field *)malloc( sizeof(Field)*hdr->nfldh );
  if ( hdr->fld  == NULL ) return -11;

  fo = *(hdr->ih+mpflds) - 1;     /* offset for field info in integer header */

  for ( n=0; n < hdr->nfldh; ++n ) {
    
    /* add null bytes to field name and units name */

    memcpy( (hdr->fld+n)->name, hdr->ch+mcflds+16*n, 8 );
    cptr = (hdr->fld+n)->name;
    for ( i=0; i <= 8; ++i, ++cptr ) {
      if ( *cptr == ' ' || i == 8 ) { *cptr = '\0'; break; }
    }

    memcpy( (hdr->fld+n)->units, hdr->ch+mcflds+16*n+8, 8 );
    cptr = (hdr->fld+n)->units;
    for ( i=0; i <= 8; ++i, ++cptr ) {
      if ( *cptr == ' ' || i == 8 ) { *cptr = '\0'; break; }
    }

    (hdr->fld+n)->lev = *(hdr->ih+fo+3*n) % 10;
    (hdr->fld+n)->tim = *(hdr->ih+fo+3*n) / 10;
/* for type 2, should calculate pointer into packed record for packed data */
    (hdr->fld+n)->off = *(hdr->ih+fo+3*n+1) -1;/* convert Fortran index to C index*/
    (hdr->fld+n)->pdens = *(hdr->ih+fo+3*n+2);
    (hdr->fld+n)->spval_flg = 0;

    (hdr->fld+n)->long_name = (hdr->fld+n)->name;
    (hdr->fld+n)->udunit = (hdr->fld+n)->units;

    /* search table and add long_name and units */

    for ( i=0; ccm_name_tab[i].std_name != NULL; i++ ) {

      if ( strcmp( ccm_name_tab[i].std_name, (hdr->fld+n)->name ) == 0 ) {

        (hdr->fld+n)->long_name = ccm_name_tab[i].long_name;

        if ( strcmp( ccm_name_tab[i].std_unit, (hdr->fld+n)->units ) == 0 ) {
          (hdr->fld+n)->udunit = ccm_name_tab[i].udunit;
        }
	break;
      }
    }
  }

  /* free any memory that was required for data conversions */
  if ( need_icv( ff ) || need_fcv( ff ) ) {
    free( buf );
  }

  return 0;
}                                                                /* init_hdr */
/*===========================================================================*/

int next_hdr( FFILE ff, CCM_HDR *hdr )
{
  /*
    Read the next header (3 records) in a CCM history tape and update the
    hdr structure

    Input arguments:
  
    ff      Input stream.

    Output arguments:

    hdr     Updated header structure.

    Return values:

     0  Successful return.
    -5  Error allocating memory for data requiring conversion.
    -6  Error reading integer header.
    -7  Error reading character header.
    -8  Error reading float header.
  */

#ifndef RS6000
  /* Offsets for ccm1 and ccm2 integer header */
  extern const int mfilh, nsteph, nitslf, ndcur, nscur, ncdate, ncsec;
#endif

  int ihlen, chlen, fhlen, itmp;
  Hint *ih;
  char *ch;
  Hfloat *fh;
  int nbytes;
  void *buf;   /* storage for data that requires conversion */
/*---------------------------------------------------------------------------*/

  ihlen = hdr->lenhdi;
  chlen = hdr->lenhdc;
  fhlen = hdr->lenhdr;
  ih = hdr->ih;
  ch = hdr->ch;
  fh = hdr->fh;

  /* Allocate space for data conversion. */

  if ( need_icv( ff ) || need_fcv( ff ) ) {
    nbytes = (ihlen*ff.isz > fhlen*ff.fsz) ? ihlen*ff.isz : fhlen*ff.fsz;
    if ( (buf = malloc( nbytes )) == NULL ) return -5;
  }

  /* Read header records. */

  if ( need_icv( ff ) ) {
    if ( rdcv_irec( ff, ihlen, buf, ih ) != ihlen ) return -6;
  } else {
    if ( read_rec( ff, ih, ihlen*ff.isz ) != ihlen*ff.isz ) return -6;
  }

  if ( read_rec( ff, ch, chlen*8 ) != chlen*8 ) return -7;

  if ( need_fcv( ff ) ) {
    if ( (itmp=rdcv_frec( ff, fhlen, buf, fh )) != fhlen ) {
      if ( itmp == -3 ) {
/*
	printf("\nnext_hdr: ***warning*** out of range floats in header\n");
*/
      } else {
	printf("next_hdr: %d floats requested, %d floats returned\n", fhlen,
	       itmp );
	return -8;
      }
    }
  } else {
    if ( read_rec( ff, fh, fhlen*ff.fsz ) != fhlen*ff.fsz ) return -8;
  }

  /* If top of model is at zero pressure, modify so that it isn't */
  if ( *hdr->mphyab == 0. ) {
    *hdr->mphyab = *(hdr->mphyab+1) / 10.;
    if ( hdr->mphya != NULL ) *hdr->mphya = *(hdr->mphya+1) / 10.;
    if ( hdr->mphya != NULL ) *hdr->mphyb = *(hdr->mphyb+1) / 10.;
  }

  /* Increment time sample counter */

  hdr->ts_count += 1;

  /* Update data from integer header record */

  hdr->mfilh   = *(ih+mfilh);
  hdr->nsteph  = *(ih+nsteph);
  hdr->nitslf  = *(ih+nitslf);
  hdr->ndcur   = *(ih+ndcur);
  hdr->nscur   = *(ih+nscur);
  hdr->ncdate  = *(ih+ncdate);
  hdr->ncsec   = *(ih+ncsec);

  return 0;
}                                                                /* next_hdr */
/*===========================================================================*/

void charr2str( char* in, int nc1 )
{
  /* Convert an input FORTRAN character array to a string.  The input array
     is assumed to be long enough that the added null byte won't overwrite
     any input data.

     Input argument:

     in    input array
     nc1   maximum number of characters in the input array
  */

  char* end;
  int len;

  /* Start at end of input FORTRAN character array, and determine position of
     last significant (non-blank or non-null) character.
  */
  end = in + nc1 - 1;
  for ( len = nc1; len > 0; --len, --end ) {
    if ( *end != ' ' && *end != '\0' ) break;
    if ( len == 1 ) {
      len = 0;
      break;
    }
  }

  /* Add null byte after last significant character */

  *(in+len) = '\0';
}
