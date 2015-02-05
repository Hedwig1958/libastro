//
//  vsop87.h
//  tstastro
//
//  Created by Corto Maltese on 30/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//

#ifndef __vsop87__
#define __vsop87__


#include "cnvfunc.h"
#include "common.h"
#include "nutation.h"
#include "coordinate.h"



#define  NUM_VSOP_PERIOTERMS          18


/*
 * Terms definition for 'VSOP87 theory'.
 */
typedef struct VSOP_TERMS_T
{
   int idx;
   double a,b,c;
   
} vsop_terms_t;


/*
 * Tableau pour les termes periodiques.
 */
typedef struct VSOP_PERIOTERMS_T
{
   
   //__Number of terms in serie.
   long long unsigned    n;
   
   //__ Terms.
   vsop_terms_t        * t;
   
} vsop_perioterms_t;


/*
 * Structure d'accueil pour le calcul des positions geocentriques.
 * Structure pour les planetes.
 */
typedef struct PLANETS_DATA_T
{
   
   //__ Nutation in longitude & obliquity.
   nut_perioterms_t    * nutation;
   
   //__ vsop87 terms.
   vsop_perioterms_t   * earth;
   vsop_perioterms_t   * planet;
   
} planets_data_t;




/*
 *  Function prototypes.
 */


planets_data_t ** _planets_data(void);
int               loadVsop87(char * f, vsop_perioterms_t * pt);
sph_circ_t        vsop87(vsop_perioterms_t * pt, double jd);


#endif /* defined(__vsop87__) */
