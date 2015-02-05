//
//  cnvfunc.h
//  tstastro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//



#ifndef __cnvfunc__
#define __cnvfunc__


#include "date.h"
#include "fbool.h"


#define  DG_0_2PI                      1
#define  DG_PI_PI                      2


#define  ANGLE_2PI                 false
#define  ANGLE_PI                   true


/*--- Prototype des fonctions de convertion ---*/

double         rdtodg(double rd);
double         dgtord(double dg);
double         rdtosec(double rd);
double         sectord(double sec);
double         hmstord(double hms);
double         rdtohms(double rd);
char *         daystrhms(double day, char * buff);
char *         hstrhms(double hour, char * buff);
char *         rdstrhms(double dg, char * buff);
char *         rdstrdms(double rd, char * buff, int type);
double         rdreduce(double rd, fbool f_decli);
double         gregtojd(date_t * date);
double         jultojd(date_t * date);
int            jdtogreg(date_t * date, double jd);
int            jdtojul(date_t * date, double jd);
double         interpolation(double y1, double y2, double y3, double n);
double         interpolationBessel(double y0, double y1, double y2, double y3, double n);



#endif /* defined(__cnvfunc__) */
