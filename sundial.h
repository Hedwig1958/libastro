//
//  sundial.h
//  astro
//
//  Created by Corto Maltese on 30/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//

#ifndef __sundial__
#define __sundial__



#include "common.h"
#include "astro.h"



/*
 * Sundial variables container.
 */
typedef struct SUN_DATA_T
{
   
   double jd;
   double eta;
   double theta0;
   double event_time;
   equ_coord_t pos;
   
} sun_data_t;



/*
 * Shadow coordinate.
 */
typedef struct SHADOW_COORD_T
{
   
   fbool  illum;;
   double x;
   double y;

} shadow_coord_t;


/*
 * Stylus size & coordinate.
 */
typedef struct STYLUS_PARAM_T
{
   double length;
   double x0;
   double y0;
   double u;
   double psi;
   
} stylus_param_t;


/*
 *  Function prototypes.
 */

sun_data_t        sunPositionFromHourAngle(nut_perioterms_t  * nut, planets_data_t * earth, geo_circ_t * geo, double  jd, double eta);
void              planarSundial(nut_perioterms_t  * nut, planets_data_t * earth, geo_circ_t * g, double declinaison, double zenithal, double length, double year);
void              sunRiseTransitSetting(nut_perioterms_t * nut, planets_data_t * earth, geo_circ_t * geo, double jd);
stylus_param_t    stylus(double latitude, double declinaison, double zenithal, double length);
shadow_coord_t    shadow(double latitude, double delta, double declinaison, double zenithal, double eta, double length);



#endif /* defined(__astro__sundial__) */
