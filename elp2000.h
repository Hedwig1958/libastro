//
//  elp2000.h
//  tstastro
//
//  Created by Corto Maltese on 30/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//

#ifndef __elp2000__
#define __elp2000__

#include "common.h"
#include "nutation.h"
#include "cnvfunc.h"
#include "coordinate.h"


#define  NUM_ELP_PERIOTERMS           36
#define  C1                         60.0
#define  C2                       3600.0


/*
 * Terms definition for 'ELP2000 theory'.
 */
typedef struct ELP_TERMS_T1
{
   
   int i1, i2, i3, i4;
   double a,b1,b2,b3,b4,b5,b6;
   
} elp_terms_t1;


/*
 * Terms definition for 'ELP2000 theory'.
 */
typedef struct ELP_TERMS_T2
{
   
   int i1, i2, i3, i4, i5;
   double phi, a, p;
   
} elp_terms_t2;


/*
 * Terms definition for 'ELP2000 theory'.
 */
typedef struct ELP_TERMS_T3
{
   
   int i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11;
   double phi, a, p;
   
} elp_terms_t3;


/*--- ---*/
typedef struct ELP_PERIOTERMS_T
{
   
   /*--- Number of terms in serie ---*/
   int       n[NUM_ELP_PERIOTERMS];
   /*--- Terms ---*/
   void    * e[NUM_ELP_PERIOTERMS];
   
} elp_perioterms_t;


/*
 * Structure d'accueil pour le calcul des positions geocentriques.
 * Structure pour la lune.
 */
typedef struct LUNAR_DATA_T
{
   
   //__ Nutation in longitude & obliquity.
   nut_perioterms_t   * nutation;
   
   //__ Elp2000.
   elp_perioterms_t   * moon;
   
} lunar_data_t;



/*
 * Locale circonstance for one eclipse.
 */
typedef struct LUNAR_SHIFT_T
{
   
   double lambda;
   double beta;
   
} lunar_shift_t;



/*
 *  Function prototypes.
 */


lunar_data_t   ** _lunar_data(void);
int               loadElp2000(elp_perioterms_t * pt, char * work_directory);
int               fileElp2000(char * f, void ** t);
sph_circ_t        elp2000(elp_perioterms_t * pt, double jd);
sph_circ_t        elp2000_bis(elp_perioterms_t * pt, double jd);



#endif /* defined(__elp2000__) */
