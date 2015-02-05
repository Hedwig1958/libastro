//
//  nutation.c
//  tstastro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>


#include "nutation.h"



/*
 * Nutation periods terms.
 */

nut_perioterms_t ** _nut_perioterms(void)
{
   static nut_perioterms_t * n;
   return &n;
}



/*
 * Loading terms for nutation in longitude & obliquity.
 */

int loadNutation(char * f, nut_perioterms_t * pt)
{
   
   FILE           * data;
   long unsigned    size;
   struct stat      s;
   
   
   //__ Ouverture du fichier contenant des data.
   if(!(data = fopen(f,"rb")))
     {
      fprintf(__ERROR_STREAM__,"This file '%s' cannot be correctly opened.\n",f);
      exit(EXIT_FAILURE);
     }
   
   //__ Search size of file.
   if(stat(f,&s))
     {
      fprintf(__ERROR_STREAM__,"This file '%s' cannot be correctly opened for size information.\n",f);
      exit(EXIT_FAILURE);
     }
   
   //__ Dynamic allocation of memory table ___*/
   if(!(pt->t = (nut_terms_t *) malloc(s.st_size)))
     {
      fprintf(__ERROR_STREAM__,"Memory size cannot be correctly allocated in function loadNutation\n");
      exit(EXIT_FAILURE);
     }
   
   /*___ Update structure 'perioterms' ___*/
   pt->n = s.st_size/sizeof(nut_terms_t);
   size = fread(pt->t,1,s.st_size,data);
   
   //__ Clossing 'data file'.
   fclose(data);
   
   return EXIT_SUCCESS;
}




/*
 * Nutation in longitude & Obliquity.
 */

int nutation(nut_perioterms_t * pt, double jd, double * d_psi, double * d_epsilon, double * epsilon0)
{
   double         tau,tau100;
   double         d,m,mp,f,omega;
   double         x;
   nut_terms_t  * t;
   
   //__  Time in Julian century from epoch 2000.
   tau=(jd-EPOCH_2000)/JULIAN_CENTURIES;
   
   tau100=tau/100;
   
   //__ Calculate terms 'd,m,mp,f,omega'.
   d = 297.85036+tau*445267.111480-tau*tau*0.0019142+tau*tau*tau/189474;
   m = 357.52772+tau*35999.050340-tau*tau*0.0001603-tau*tau*tau/300000;
   mp = 134.96298+tau*477198.867398+tau*tau*0.0086972+tau*tau*tau/56250;
   f = 93.27191+tau*483202.017538-tau*tau*0.0036825+tau*tau*tau/327270;
   omega = 125.04452-tau*1934.136261+tau*tau*0.0020708+tau*tau*tau/450000;
   
   //__ Init sone arguments.
   *d_psi = *d_epsilon=*epsilon0=0.0;
   
   //__ Calculate nutation in longitude & obliquity.
   for(int i=0;i<pt->n;i++)
     {
      t = &pt->t[i];
      x = t->i1*d+t->i2*m+t->i3*mp+t->i4*f+t->i5*omega;
      *d_psi += t->a*sin(dgtord(x))+t->b*sin(dgtord(x))*tau;
      *d_epsilon += t->c*cos(dgtord(x))+t->d*cos(dgtord(x))*tau;
     }
   
   //__ Convert in radians.
   *d_psi = dgtord(*d_psi/36000000);
   *d_epsilon = dgtord(*d_epsilon/36000000);
   
   //__ The mean obliquity of the ecliptic.
   *epsilon0 = 84381.448-tau100*(4680.93-tau100*(1.55+tau100*(1999.25-tau100*
      (51.38-tau100*(249.67-tau100*(39.05+tau100*(7.12+tau100*(27.87+tau100*(5.79+tau100*2.45)))))))));
   
   //__ Convert in radians.
   *epsilon0 = sectord(*epsilon0);
   
   return EXIT_SUCCESS;
}




