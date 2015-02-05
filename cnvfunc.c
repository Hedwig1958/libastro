//
//  cnvfunc.c
//  tstastro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//


#include <stdio.h>
#include <math.h>

#include "cnvfunc.h"


/*
 * transform radian to degrees.
 */

double rdtodg(double rd)
{
   return rd / M_PI * 180.0;
}


/*
 * transform degrees to radian.
 */

double dgtord(double dg)
{
   return dg  * M_PI / 180.0;
}

/*
 * transform radian to second.
 */

double rdtosec(double rd)
{
   return rd / M_PI * 648000.0;
}


/*
 * transform second to radian.
 */

double sectord(double dg)
{
   return dg  * M_PI / 648000.0;
}

/*
 * transform h.m.s to radian.
 */

double hmstord(double hms)
{
   return hms  * M_PI / 12.0;
}


/*
 * transform radian to h.m.s.
 */

double rdtohms(double rd)
{
   return rd / M_PI * 12.0;
}

/*
 * Convert 'D.dec' to 'H,M,S.dec'.
 */

char * daystrhms(double day, char * buff)
{
   double x,h,m;
   
   if(day<0.0)
      day+=1;
   
   x=day*24;
   h=floor(x);
   x-=h;
   x*=60.0;
   m=floor(x);
   x-=m;
   x*=60.0;
   
   /*--- ---*/
   sprintf(buff,"%02.0fh%02.0fm%06.3fs",h,m,x);
   return buff;
}


/*
 * Convert 'H.dec' to 'H,M,S.dec'.
 */

char * hstrhms(double hour, char * buff)
{
   
   double x,h,m;
   
   if(hour<0.0)
      hour+=24.0;
   
   x=hour;
   h=floor(x);
   x-=h;
   x*=60.0;
   m=floor(x);
   x-=m;
   x*=60.0;
   
   //__ result output.
   sprintf(buff,"%02.0fh%02.0fm%06.3fs",h,m,x);
   return buff;
}



/*
 * Convert angle in radians to string H,M,S.dec.
 */

char * rdstrhms(double rd, char * buff)
{
   
   double x,h,m;
   
   x=atan2(sin(rd),cos(rd));
   
   if(x<0.0)
      x+=(M_PI*2);
   
   //__ convert to degrees.
   x=rdtodg(x);
   
   x/=15.0;
   h=floor(x);
   x-=h;
   x*=60.0;
   m=floor(x);
   x-=m;
   x*=60.0;
   
   //__result output.
   sprintf(buff,"%02.0fh%02.0fm%06.3fs",h,m,x);
   
   return buff;
}


/*
 * Convert angle in radians to D,M,S.dec.
 */

char * rdstrdms(double rd, char * buff, int type)
{
   double x,d,m,r;
   
   /*
    ---------------------------------------------
    arg : type
    
    DG_0_2PI    :    0    <= angle <= + 2*PI
    DG_PI_PI    :  - PI   <= angle <= + PI
    
    ---------------------------------------------
    */
   
   x=atan2(sin(rd),cos(rd));
   
   if(type == DG_PI_PI)
     {
      //__ -PI <= angle <= +PI
      if(x > M_PI)
         x -= (M_PI*2);
     }
   
   //__ convert to degrees.
   x=rdtodg(x);
   
   r=fabs(x);
   d=floor(r);
   r-=d;
   r*=60.0;
   m=floor(r);
   r-=m;
   r*=60.0;
   
   //__ transformation to string value.
   sprintf(buff,"%c%02.0fd%02.0f'%05.2f\"",x >= 0.0 ?' ':'-',d,m,r);
   
   return buff;
   
}


/*
 * Reduction d'un angle en radian.
 */

double rdreduce(double rd, fbool f_decli)
{
   
   double x;
   
   x=fmod(rd,(M_PI*2));
   
   if(!f_decli && x<0.0)
     {
      x+=(M_PI*2);
     }
   else
     {
      if(f_decli && x<-M_PI)
        {
         x+=(M_PI*2);
        }
      else
        {
         if(f_decli && x>M_PI)
           {
            x-=(M_PI*2);
           }
        }
     }
   
   return x;
}


/*
 * Interpolation - Meeus 3.2
 */

double interpolation(double y1, double y2, double y3, double n)
{
   double a,b,c;
   
   a=y2-y1;
   b=y3-y2;
   c=b-a;
   return y2+n/2*(a+b+n*c);
}

/*
 * Interpolation, formule de Bessel.
 */

double interpolationBessel(double y0, double y1, double y2, double y3, double n)
{
   double d11, d12, d13;
   double d21, d22;
   
   //__ Differences premieres.
   d11=y1-y0;
   d12=y2-y1;
   d13=y3-y2;
   
   //__ Differences secondes.
   d21=d12-d11;
   d22=d13-d12;
   
   return y1+(n*d12)-(0.25*n*(1-n)*(d21+d22));
}

