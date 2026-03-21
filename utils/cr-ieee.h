/*$Id: cr-ieee.h,v 1.2 1998/03/30 19:59:43 eaton Exp $*/

#ifndef CR_IEEE_H
#define CR_IEEE_H

int
ctospi( void *in, void *out, int n );

void
sptoci( void *in, void *out, int n );

int
ctospf( void *in, void *out, int n );

int
sptocf( void *in, void *out, int n );

int
ctodpf( void *in, void *out, int n );

int
dptocf( void *in, void *out, int n );

#endif
