//
//  sideraltime.c
//  tstastro
//
//  Created by Corto Maltese on 7/12/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//


#include <stdio.h>
#include <math.h>

#include "sideraltime.h"


/*
 * Sideral time for Greenwish at 0h.
 */

double sideraltime0(double jd)
{
   
   double t, theta0;
   
   //__ Time in Julian centuries from epoch 2000 - Meeus 11.1
   t=(floor(jd)+0.5-2451545.0)/36525.0;
   
   //__ Mean Sideral Time at Greenwich at 0h UT - Meeus 11.3
   theta0=100.46061837+t*(36000.770053608+t*(0.000387933-t*(1/38710000)));
   
   return fmod(dgtord(theta0),(M_PI*2));
}



/*
 * Sideral time for Greenwish at any instant.
 */

double sideraltime(double jd)
{
   
   double t, theta0;
   
   //__ Time in Julian centuries from epoch 2000 - Meeus 11.1
   t=(jd-EPOCH_2000)/JULIAN_CENTURIES;
   
   //__ Mean Sideral Time at Greenwich at 0h UT - Meeus 11.4
   theta0=280.46061837+360.98564736629*(jd-2451545.0)+t*t*(0.000387933-t*(1/38710000));
   
   return fmod(dgtord(theta0),(M_PI*2));
}




