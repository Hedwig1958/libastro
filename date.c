//
//  date.c
//  tstastro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "date.h"


/*
 * Delta Time.
 */

double deltatime(double jd)
{
   
   double dt = 0.0,t;
   
   //__  Time in Julian centuries from epoch 2000.
   t=(jd-EPOCH_2000)/JULIAN_CENTURIES;
   
   //__ Pour periode avant 948 AD.
   if(t<9.48-20.0)
      dt=2715.6+t*(573.36+t*46.5);
   
   //__  Pour periode apres 948 AD et avant 1600 AD.
   if(t>=9.48-20.0 && t<16-20)
      dt=50.6+t*(67.5+t*22.5);
   
   if(t>=16-20)
      dt=-15+pow(jd-2382148,2)/41048480;
   
   return dt;
}



/*
 * Convert 'julian day' into Julian date in format 'date_t'
 */

int jdtojul(date_t * date, double jd)
{
   
   double a,b,c,d,e,f,z;
   
   memset(date,0,sizeof(date_t));
   
   //__ 'z' integer part of 'jd + 0.5'.
   z=floor(jd+0.5);
   
   //__ 'f' decimal part of 'jd + 0.5'.
   f=jd+0.5-z;
   
   //__ Pour le calendrier Julien 'a' vaut 'z'.
   a=z;
   b=a+1524;
   c=floor((b-122.1)/365.25);
   d=floor(365.25*c);
   e=floor((b-d)/30.6001);
   
   //__ Compute date.
   
   date->day=(short)(b-d-floor(30.6001*e));
   date->month=e<14 ? (short)e-1 : (short)e-13;
   date->year=date->month>2.0 ? (short)c-4716 : (short)c-4715;
   
   //__ Compute hour.
   if(f>0.0)
     {
      date->hour=(short)floor(f*24.0);
      f-=date->hour/24.0;
      date->min=(short)floor(f*1440.0);
      f-=date->min/1440.0;
      date->sec=(short)floor(f*86400.0);
      f-=date->sec/86400.0;
      date->sec100=f*8640000.0;
     }
   
   //__ the julian day '2299161' correspond at Gregorian date 15 octber 1582.
   return EXIT_SUCCESS;
}



/*
 * Convert 'julian day' into Gregorian date in format 'date_t'
 */

int jdtogreg(date_t * date, double jd)
{
   
   double alpha, a,b,c,d,e,f,z;
   
   memset(date,0,sizeof(date_t));
   
   //__ 'z' partie entiere de 'jd + 0.5'.
   z=floor(jd+0.5);
   
   //__ 'f' partie decimale de 'jd + 0.5'.
   f=jd+0.5-z;
   
   //__ Nombre de siecles avec date pivot ''.
   alpha=floor((z-1867216.25)/36524.25);
   
   a=z+1+alpha-floor(alpha/4);
   b=a+1524;
   c=floor((b-122.1)/365.25);
   d=floor(365.25*c);
   e=floor((b-d)/30.6001);
   
   //__ Mise a jour de la date.
   date->day=(short)(b-d-floor(30.6001*e));
   date->month = e<14 ? (short)e-1 : (short)e-13 ;
   date->year = date->month>2.0 ? (short)c-4716 : (short)c-4715 ;
   
   //__ Mise a jour de l'heure.
   if(f>0.0)
     {
      date->hour=(short)floor(f*24.0);
      f-=date->hour/24.0;
      date->min=(short)floor(f*1440.0);
      f-=date->min/1440.0;
      date->sec=(short)floor(f*86400.0);
      f-=date->sec/86400.0;
      date->sec100=f*8640000.0;
      
      //__ Arrondi au 1/1000 sec.
      if(date->sec100>=99.9)
        {
         date->sec100=0.0;
         date->sec++;
         if(date->sec==60)
           {
            date->sec=0;
            date->min++;
            if(date->min==60)
              {
               date->min=0;
               date->hour++;
               if(date->hour==24)
                 {
                  date->hour=0;
                 }
              }
           }
        }
     }
   
   return EXIT_SUCCESS;
}



/*
 * Convert Gregorian date of type 'date_t' to 'julian day'
 */

double gregtojd(date_t * d)
{
   short   century,year=d->year,month=d->month,delay;
   
   //__ Correction for january and february.
   if(d->month<=2)
     {
      year--;
      month+=12;
     }
   
   //__ Compute the century.
   century=year/100;
   
   //__ Gregorian Delay.
   delay=2-century+century/4;
   
   //__ Compute the julian day.
   return floor((365.25*(4716+year)))+floor((30.6001*(month+1)))+d->day+delay-1524.5+
   ((double)d->hour/24)+((double)d->min)/1440+((double)d->sec)/86400+((double)d->sec100)/8640000;
}



/*
 * Convert Julian date of type 'date_t' to 'julian day'.
 */

double jultojd(date_t * d)
{
   short   century,year=d->year,month=d->month;
   
   //__ Correction for january and february.
   if(d->month<=2)
     {
      year--;
      month+=12;
     }
   
   //__ Compute the century.
   century=year/100;
   
   //__ Compute the julian day.
   return floor((365.25*(year+4716)))+floor((30.6001*(month+1)))+d->day-1524.5+
   ((double)d->hour)/24+((double)d->min)/1440+((double)d->sec)/86400+((double)d->sec100)/8640000;
}



/*
 * Impression de la partie 'date' d'une date de type 'date_t'.
 */

char * pdate(date_t * date)
{
   static char buff[11];
   sprintf(buff,"%02d/%02d/%04d",date->day,date->month,date->year);
   return buff;
}



/*
 * Impression de la partie 'heure' d'une date de type 'date_t'.
 */

char * ptime(date_t * date)
{
   static char buff[14];
   sprintf(buff,"%02d:%02d:%02d,%04.1f",date->hour,date->min,date->sec,date->sec100);
   return buff;
}



/*
 * Positionne la date avec certaines verification sur la validite.
 */

int setdate(date_t * date, short y, short m, short d, short h, short mi, short s, double s100)
{
   
   /*--- Ajouter des tests de validites des dates ---*/
   date->year=y;
   date->month=m;
   date->day=d;
   date->hour=h;
   date->min=mi;
   date->sec=s;
   date->sec100=s100;
   
   return EXIT_SUCCESS;
}




