/*$Id: iscosb.h,v 1.1.1.1 1998/02/17 23:44:02 eaton Exp $*/

#ifndef ISCOSB_H
#define ISCOSB_H

typedef struct {
  unsigned int type  :  4;
  unsigned int junk  : 27;
  unsigned int bnhi  :  1;
  unsigned int bnlo  : 23;
  unsigned int fwi   :  9;
} BCW;

/*
  The BCW bitfield is used to extract the type of control word and the
  block number from a block control word.  In an 8 byte Cray control word
  the block number would be read as a single 24 bit piece of data.  In
  BCW the block number is split into two pieces, bnhi and bnlo, so that
  the block number field does not span a word boundary on a machine with
  4-byte words.

  See the Cray manpage BLOCKED(4F) for more information.
*/


extern int iscosb( const char *filepath );

#endif
