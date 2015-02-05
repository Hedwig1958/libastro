//
//  vsop87.c
//  tstastro
//
//  Created by Corto Maltese on 30/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>


#include "vsop87.h"


/*
 * Planets data type.
 */

planets_data_t ** _planets_data(void)
{
   static planets_data_t * d;
   return &d;
}


/*
 * Loading terms for planets heliocentrics positions.
 */

int loadVsop87(char * f, vsop_perioterms_t * pt)
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
   
   //__ Search size of file ___*/
   if(stat(f,&s))
     {
      fprintf(__ERROR_STREAM__,"This file '%s' cannot be correctly opened for size information.\n",f);
      exit(EXIT_FAILURE);
     }
   
   //__ Dynamic allocation of memory table.
   if(!(pt->t = (vsop_terms_t *) malloc(s.st_size)))
     {
      fprintf(__ERROR_STREAM__,"Memory size cannot be correctly allocated in function load_vsop87\n");
      exit(EXIT_FAILURE);
     }
   
   //__ Update structure 'perioterms'.
   pt->n = s.st_size/sizeof(vsop_terms_t);
   size = fread(pt->t,1,s.st_size,data);
   
   //__ Clossing 'data file'.
   fclose(data);
   
   return EXIT_SUCCESS;
}



/*
 * Variations Seculaire of Orbitals Planets.
 *
 * Ref : Astronomical algorithms, Jean Meeus.
 *
 */

sph_circ_t vsop87(vsop_perioterms_t * pt, double jd)
{
   vsop_terms_t  * p;
   double          tau, latitude, longitude, radius, t[NUM_VSOP_PERIOTERMS];
   static sph_circ_t event;
   
   //__ Initialisation of 't[]'.
   memset(t,0,sizeof(t));
   
   //__ Time in Julian millenia from epoch 2000.
   tau=(jd-EPOCH_2000)/JULIAN_MILLENIA;
   
   //__ Compute terms l0->l5, b0->b5, r0->r5.
   for(int i=0;i<pt->n;i++)
     {
      p = &pt->t[i];
      t[p->idx] += p->a*cos(p->b+p->c*tau);
     }
   
   //__ Compute longitude in radians t[0->5] equivalent at l0 -> l5.
   longitude = t[0]+tau*(t[1]+tau*(t[2]+tau*(t[3]+tau*(t[4]+tau*t[5]))));
   
   //__ Compute latitude  in radians t[6->11] equivalent at b0 -> b5.
   latitude = t[6]+tau*(t[7]+tau*(t[8]+tau*(t[9]+tau*(t[10]+tau*t[11]))));
   
   //__ Compute radius in UA t[12->17] equivalent at r0 -> r5.
   radius = t[12]+tau*(t[13]+tau*(t[14]+tau*(t[15]+tau*(t[16]+tau*t[17]))));
   
   event.pos.l = longitude;
   event.pos.b = latitude;
   event.pos.radius = radius;
   
   return event;
}

