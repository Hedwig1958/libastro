//
//  stars.h
//  tstastro
//
//  Created by Corto Maltese on 30/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//

#ifndef __stars__
#define __stars__

#include "common.h"
#include "nutation.h"
#include "cnvfunc.h"


/*
 * Definition d'une structure de type 'STAR_REC_T'.
 */
typedef struct STAR_REC_T
{
   
   //__ Skymap Identifier.
   long   skymap;
   
   //__ RA & DEC in radians, proper motion for RA in second and for DEC in arcsecond.
   double ra,dec,pm_ra,pm_dec,vmag;
   
   //__ Henry Draper number.
   char   hd[7];
   char   hd_m;
   char   hd_u;
   
   //__ SAO number.
   char   sao[7];
   char   sao_m;
   
   //__ Durchmusterung number.
   char   dm[12];
   char   dm_m;
   char   dm_u;
   
   //__ Washington double stars number.
   char   wds[11];
   char   wds_m;
   char   wds_u;
   
   //__ Star Name.
   char   star_name[11];
   char   var_name[11];
   
}  star_rec_t;



/*
 * Constantes pour le traitement du Catalogue of d'etoiles.
 */
typedef struct STARS_DATA_T
{
   
   //__ Number of record in file.
   long long unsigned   n;
   
   //__ Data records.
   star_rec_t         * s;
   
} stars_data_t;


/*
 *  Function prototypes.
 */

int    dumpStarsData(stars_data_t * stars);
int    loadStars(char * f, stars_data_t * sd);


#endif /* defined(__stars__) */
