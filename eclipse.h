//
//  eclipse.h
//  tstastro
//
//  Created by Corto Maltese on 30/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//

#ifndef __eclipse__
#define __eclipse__

#include "elp2000.h"
#include "vsop87.h"
#include "astro.h"
/*
 * Constants for eclipses.
 */
#define LUNAR_ECLIPSE                0.5
#define SOLAR_ECLIPSE                0.0

#define NO_ECLIPSE                     0
#define PARTIAL                        1
#define CENTRAL_ANNULAR                2
#define CENTRAL_TOTAL                  3
#define CENTRAL_ANNULAR_TOTAL          4
#define NON_CENTRAL_ANNULAR            5
#define NON_CENTRAL_TOTAL              6


#define  HOURS_IN_DAY               24.0
#define  BESSEL_SHIFT_TIME     (3.0/HOURS_IN_DAY)
#define  BESSEL_HOURS_IN_DAY   (1.0/HOURS_IN_DAY)
#define  BESSEL_MIDDLE_ELEM            7

#define  NUM_SAMPLES                   7

/*
 * Stucture d'accueil pour les donnÈes relatives
 * aux eclipses.
 */
typedef struct ECLIPSES_T
{
   
   int    type;
   double jd;
   double gamma;
   double u;
   double magnitude;
   
} eclipses_t;


/*
 * Structure d'accueil les les donnÈes
 * de type 'Besselian elements for eclipses'.
 */
typedef struct BESSEL_ELEM_T
{
   
   double  jd;
   char    type[2];
   double  gamma;
   double  k;
   double  saros;
   double  t0;
   double  x[4];
   double  y[4];
   double  z[4];
   double  d[3];
   double  m[5];
   double  f1;
   double  f2;
   double  u1[3];
   double  u2[3];
   
} bessel_elem_t;


/*
 * Locale circonstance for one eclipse.
 */
typedef struct LOC_CIRCONST_T
{
   
   double kcnt;
   double jd;
   double h;
   double longitude;
   double latitude;
   double width;
   double duration;
   double moon_ratio;
   
} loc_circonst_t;



/*
 *  Function prototypes.
 */

double            _eclipses_ra(double jd);
eclipses_t        eclipses(double type, double year);
bessel_elem_t     besselianElements(planets_data_t * d1, lunar_data_t * d2, lunar_shift_t * lf, double k0, double jd, double dt);
int               eclipseCircunsTime(bessel_elem_t * bess, loc_circonst_t * loc, double jd, double dt);
int               eclipseCircunsLongitude(bessel_elem_t * bess, loc_circonst_t * loc, double lambda, double tol
                                          , double i, double g, double dt);
int               eclipseExtremePoints(bessel_elem_t * bess, loc_circonst_t * ext, double dt);
int               eclipseMaximumPoint(bessel_elem_t * bess, loc_circonst_t * ext, double dt);
int               eclipseNoonPoint(bessel_elem_t * bess, loc_circonst_t * ext, double dt);



#endif /* defined(__eclipse__) */
