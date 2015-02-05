//
//  date.h
//  astro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//

#ifndef __date__
#define __date__


#include "astrocst.h"

#ifndef   __DATE_T__
#define   __DATE_T__


//__ Definition d'une structure de type 'date_t'.

typedef struct DATE_T
{
   short   year, month, day, hour, min, sec;
   double  sec100;
   
} date_t;

#endif /* defined(__DATE_T__) */



/*
 *  Function prototypes.
 */


double     deltatime(double jd);
int        jdtojul(date_t * date, double jd);
int        jdtogreg(date_t * date, double jd);
double     gregtojd(date_t * d);
double     jultojd(date_t * d);
char *     ptime(date_t * date);
char *     pdate(date_t * date);
int        setdate(date_t * date, short y, short m, short d, short h, short mi, short s, double s100);


#endif /* defined(__date__) */
