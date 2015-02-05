//
//  astro.h
//  astro
//
//  Created by Corto Maltese on 30/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//



#ifndef __astro__
#define __astro__

#include "linear.h"
#include "sideraltime.h"
#include "stars.h"
#include "vsop87.h"
#include "elp2000.h"
#include "sundial.h"



/*
 * Cantantes pour usage generale.
 */

#define  OPTIQUE                    true
#define  BARYCENTRE                false
#define  RON_VONDRAK_ELEM             36
#define  ITERATION_LIGHT_TIME_MOON     2
#define  ITERATION_LIGHT_TIME_SUN      2
#define  ITERATION_LIGHT_TIME_PLANETS  2
#define  LOW_ACCURATE_FOR_SEASONS      1
#define  HIGH_ACCURATE_FOR_SEASONS     3


/*
 * Constantes pour les saisons.
 */
#define SPRING                         0
#define SUMMER                         1
#define AUTUMN                         2
#define WINTER                         3


/*
 * Constantes pour les phases lunaire.
 */
#define NEW_MOON                     0.0
#define FIRST_QUARTER                0.25
#define FULL_MOON                    0.5
#define LAST_QUARTER                 0.75



/*
 * Constante for Rising, Transit and Setting.
 */

#define EPHEM_NUM_ITERATION            5
#define EPHEM_DAY                      1
#define NOT_CIRCUMPOLAR            false
#define CIRCUMPOLAR                 true
#define RISING                         1
#define SETTING                        2



/*
 * Ephemerides data for 'STARS', 'SUN', 'PLANET', 'MOON'.
 */
typedef struct EPHEM_DATA_T
{
   
   double      jd;
   int         type;
   fbool       circumpolar;
   double      h0;
   double      theta0;
   double      eta;
   double      azimuth;
   
   equ_coord_t pos;
   
} ephem_data_t;


/*
 * Structure for seasons elements.
 */
typedef struct SEASONS_DATA_T
{
   
   double     jd;
   int        type;
   
   equ_circ_t circ;
   
} seasons_data_t;


/*
 * functions prototypes.
 */

int            *  _seasons(void);
double         *  _phases(void);
double            _seasons_longitude(double jd);
double            _lunarphases_longitude(double jd);
double            appSideraltime(nut_perioterms_t  * nut, double jd);
horiz_circ_t      equtohoriz(nut_perioterms_t  * nut, equ_coord_t * equ, geo_coord_t * geo, double jd, int type);
double            sunAberration(double radius, double jd);
double            parallacticAngle(double eta, double delta, double phi);
double            moonIlluminatedFraction(equ_circ_t sun, equ_circ_t moon);
double            moonPhaseAngle(equ_circ_t sun, equ_circ_t moon);
double            moonPositionAngle(equ_circ_t sun, equ_circ_t moon);
helio_circ_t      sunEphemeris(planets_data_t * data, double jd);
equ_coord_t       aberrRonVondrak(equ_coord_t * pos, double jd);
equ_circ_t        appGeoPosPlanets(planets_data_t * data, double jd);
equ_circ_t        appGeoPosSun(planets_data_t * data, double jd);
equ_circ_t        appGeoPosMoon(lunar_data_t * data, double jd, lunar_shift_t * lf);
equ_circ_t        appGeoPosStars(nut_perioterms_t  * nut, star_rec_t * rec, double jd);
double            altitude0(int type,double parallax, double refraction, double semi_diam, double altitude_obs, double distance_obs, double altitude_hor);
ephem_data_t      transit(nut_perioterms_t  * nut, equ_circ_t * equ, geo_circ_t * geo, double eta);
ephem_data_t      rising(nut_perioterms_t  * nut, equ_circ_t * equ, geo_circ_t * geo, double h0);
ephem_data_t      setting(nut_perioterms_t  * nut, equ_circ_t * equ, geo_circ_t * geo, double h0);
seasons_data_t    seasons(planets_data_t * eart, double year, int type, int accurate);
double            lunarphases(planets_data_t * d1, lunar_data_t * d2, double phase, double year, int accurate);


#endif /* defined(__astro__) */






