//
//  coordinate.h
//  astro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//




#ifndef __coordinate__
#define __coordinate__

#include "astrocst.h"
#include "cnvfunc.h"


//__ Definition d'une structure de type 'REC_COORD_T'.
typedef struct REC_COORD_T
{
   
   double x;
   double y;
   double z;
   
} rec_coord_t;


//__ Definition d'une structure de type 'SPH_COORD_T'.
typedef struct SPH_COORD_T
{
   
   double l;
   double b;
   double radius;
   
} sph_coord_t;


//__ Definition d'une structure de type 'EQU_COORD_T'.
typedef struct EQU_COORD_T
{
   
   double alpha;   //__ right ascension
   double delta;   //__ declinaison
   
} equ_coord_t;


//__ Definition d'une structure de type 'GEO_COORD_T'.
typedef struct GEO_COORD_T
{
   
   double longitude;
   double latitude;
   double altitude;
   
} geo_coord_t;


//__ Definition d'une structure de type 'HORIZ_COORD_T'.
typedef struct HORIZ_COORD_T
{
   
   double azimuth;
   double altitude;
   
} horiz_coord_t;


//__ Definition d'une structure de type 'HELIO_COORD_T'.
typedef struct HELIO_COORD_T
{
   
   double latitude;
   double longitude;
   double radius;
   
} helio_coord_t;


//__ Definition d'une structure de type 'REC_CIRC_T'.
typedef struct REC_CIRC_T
{
   
   double       jd;
   int          type;
   
   rec_coord_t  pos;
   
} rec_circ_t;


//__ Definition d'une structure de type 'SPH_CIRC_T'.
typedef struct SPH_CIRC_T
{
   
   double       jd;
   int          type;
   
   sph_coord_t  pos;
   
} sph_circ_t;


//__ Definition d'une structure de type 'GEO_CIRC_T'.
typedef struct GEO_CIRC_T
{
   
   char *       name;   
   geo_coord_t  pos;
   
}  geo_circ_t;


//__ Definition d'une structure de type 'HOIZ_CIRC_T'.
typedef struct HORIZ_CIRC_T
{
   
   double       jd;
   int          type;
   double       theta0;   //__ Temps Siderale Apparent a Greenwich.
   double       theta;    //__ Temps Siderale Apparent locale.
   double       eta;      //__ Angle Horaire.
   double       q;        //__ Angle Parallactique.
   
   horiz_coord_t  pos;
   
} horiz_circ_t;


//__ Definition d'une structure de type 'EQU_CIRC_T'.
typedef struct EQU_CIRC_T
{
   
   double      jd;
   int         type;
   double      parallax;
   double      semidiam;
   
   sph_coord_t sph;
   equ_coord_t equ;
   
} equ_circ_t;


//__ Definition d'une structure de type 'HELIO_CIRC_T'.
typedef struct HELIO_CIRC_T
{
   
   double        jd;
   int           type;
   
   helio_coord_t pos;
   
} helio_circ_t;



/*
 *  Function prototypes.
 */

equ_coord_t       eprec(double jd0, double jd, equ_coord_t * p0);
equ_coord_t       eprec2000(double jd, equ_coord_t * p0);
sph_circ_t        sprec(sph_circ_t * p0, double jd0);
sph_circ_t        sprec2000(sph_circ_t * p0);
sph_circ_t        rectosph(rec_circ_t *r);



#endif /* defined(__coordinate__) */

