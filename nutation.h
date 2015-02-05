//
//  nutation.h
//  tstastro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//

#ifndef __nutation__
#define __nutation__



#include "common.h"
#include "astrocst.h"
#include "cnvfunc.h"

/*
 * Terms definition for 'nutation & obliquity'.
 */

typedef struct NUT_TERMS_T
{
   
   int i1, i2, i3, i4, i5;
   double a,b,c,d;
   
} nut_terms_t;



/*
 * Terms periodique pour nutation & obliquity.
 */

typedef struct NUT_PERIOTERMS_T
{
   
   //__ Number of terms in serie.
   long long unsigned  n;
   
   //__ Terms.
   nut_terms_t       * t;
   
} nut_perioterms_t;


/*
 *  Function prototypes.
 */

nut_perioterms_t **  _nut_perioterms(void);
int                  loadNutation(char * f, nut_perioterms_t * pt);
int                  nutation(nut_perioterms_t * pt, double jd, double * d_psi, double * d_epsilon, double * epsilon0);

#endif /* defined(__nutation__) */
