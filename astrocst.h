//
//  astrocst.h
//  astro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//


#ifndef __astrocst__
#define __astrocst__

# include "fbool.h"

/*
 * Constantes pour le calcul 'h0'.
 */

#define  HIGH                             1
#define  LOW                              0

#define PLANET                            1
#define STARS                             2
#define SUN                               3
#define MOON                              4


#define  UA                     149597870.0              //__ Distance moyenne terre/soleil.
#define  EARTH_RADIUS               6378.14
#define  EARTH_RADIUS_METER         6378140              //__ Earth radius in meter.
#define  PI_0                      8.794148              //__ Parallaxe horizontale du soleil a UA en seconde ARC.
#define  S_0                         959.63              //__ Demi-diametre solaire en seconde ARC.
#define  LIGHT_TIME            0.0057755183              //__ Light time constant.
#define  EPOCH_2000               2451545.0              //__ Epoch 2000.
#define  JULIAN_YEARS                365.25              //__ Julian years.
#define  JULIAN_CENTURIES           36525.0              //__ Julian centuries.
#define  JULIAN_MILLENIA           365250.0              //__ Julian millenia.

#define  K_0_UMBRA_IAU            0.2725076
#define  K_0_UMBRA_MEEUS          0.2722740
#define  K_0_PENUMBRA_MEEUS       0.2724810
#define  K_0_UMBRA_NASA_T         0.2725076
#define  K_0_UMBRA_NASA_A         0.2722810

#define  SIDERAL_CST_MTTOST       1.0027379
#define  SIDERAL_CST_STTOMT       0.9972696

#define  REFRACTION_RADAU        (0.61*M_PI/180.0)
#define  REFRACTION_NAUTICAL     (0.56666667*M_PI/180.0)

#define  SUN_SEMI_DIAMETER       (0.26666667*M_PI/180.0)
#define  MOON_SEMI_DIAMETER      (0.26666667*M_PI/180.0)




#endif /* defined(__astrocst__) */


