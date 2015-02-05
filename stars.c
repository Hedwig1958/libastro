//
//  stars.c
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


#include "stars.h"

/*
 * Dump data for stars
 */

int dumpStarsData(stars_data_t * stars)
{
   int  i;
   char buff1[20];
   char buff2[20];
   
   for(i=0;i<stars->n;i++)
      printf("%d %s %s %s\n",i ,rdstrhms(stars->s[i].ra,buff1), rdstrdms(stars->s[i].dec,buff2,DG_PI_PI), stars->s[i].star_name);
   
   return EXIT_SUCCESS;
}



/*
 * Loading stars data, positions, magnitude.
 */

int loadStars(char * f, stars_data_t * sd)
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
   
   //__ Dynamic allocation of memory table.
   if(!(sd->s=(star_rec_t*)malloc(s.st_size)))
     {
      fprintf(__ERROR_STREAM__,"Memory size cannot be correctly allocated in function loadStars\n");
      exit(EXIT_FAILURE);
     }
   
   //__ Update.
   sd->n = s.st_size/sizeof(star_rec_t);
   size = fread(sd->s,1,s.st_size,data);
   
   //__ Clossing 'data file'.
   fclose(data);
   
   return EXIT_SUCCESS;
}

