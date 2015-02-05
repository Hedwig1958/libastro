//
//  coordinate.c
//  astro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//


#include <stdio.h>
#include <math.h>

#include "coordinate.h"


/*
 * Transformation des coordonnees rectangulaire en coordonnees spheriques.
 */

sph_circ_t rectosph(rec_circ_t * r)
{
   sph_circ_t s;
   
   s.jd=r->jd;
   s.type=r->type;
   s.pos.radius=sqrt(r->pos.x*r->pos.x+r->pos.y*r->pos.y+r->pos.z*r->pos.z);
   s.pos.l=atan2(r->pos.y,r->pos.x);
   s.pos.b=asin(r->pos.z/s.pos.radius);
   
   return s;
}


/*
 * Transformation epoque to another.
 */

equ_coord_t eprec(double jd0, double jd, equ_coord_t * p0)
{
   
   double       t,dt;
   double       dzeta,theta,z;
   double       a,b,c;
   equ_coord_t  p;
   
   //__ Time in Julian century from epoch 2000.
   t=(jd-EPOCH_2000)/JULIAN_CENTURIES;
   dt=(jd-jd0)/JULIAN_CENTURIES;
   
   //__ Accurate reduction of positions from one equinox to another - Meeus 20.2
   dzeta=(2306.2181+1.39656*t-0.000139*t*t)*dt
   +(0.30188-0.000344*t)*dt*dt
   +0.017998*dt*dt*dt;
   dzeta=dzeta/dgtord(3600);
   
   
   z=(2306.2181+1.39656*t-0.000139*t*t)*dt
   +(1.09468+0.000066*t)*dt*dt
   +0.018203*dt*dt*dt;
   z=z/dgtord(3600);
   
   theta=(2004.3109+0.85330*t-0.000217*t*t)*dt
   -(0.42665+0.000217*t)*dt*dt
   -0.041833*dt*dt*dt;
   theta=theta/dgtord(3600);
   
   //__ See Meeus 20.3
   a=cos(p0->delta)*sin(p0->alpha+dzeta);
   b=cos(theta)*cos(p0->delta)*cos(p0->alpha+dzeta)-sin(theta)*sin(p0->delta);
   c=sin(theta)*cos(p0->delta)*cos(p0->alpha+dzeta)+cos(theta)*sin(p0->delta);
   
   //__ Update output values - Meeus 20.4
   p.alpha=atan2(a,b)+z;
   p.delta=asin(c);
   
   return p;
}



/*
 * Transformation coordonnees equatoriales de l'epoque 2000 vers une autre.
 */

equ_coord_t eprec2000(double jd, equ_coord_t * p0)
{
   
   double       t;
   double       dzeta,theta,z;
   double       a,b,c;
   equ_coord_t  p;
   
   //__  Time in Julian century from epoch 2000
   t=(jd-EPOCH_2000)/JULIAN_CENTURIES;
   
   //__ Accurate reduction of positions from one equinox to another - Meeus 20.2
   dzeta=2306.2181*t+0.30188*t*t+0.017998*t*t*t;
   dzeta=dzeta/dgtord(3600);
   
   z=2306.2181*t+1.09468*t*t+0.018203*t*t*t;
   z=z/dgtord(3600);
   
   theta=2004.3109*t-0.42665*t*t-0.041833*t*t*t;
   theta=theta/dgtord(3600);
   
   //__ See Meeus 20.3
   a=cos(p0->delta)*sin(p0->alpha+dzeta);
   b=cos(theta)*cos(p0->delta)*cos(p0->alpha+dzeta)-sin(theta)*sin(p0->delta);
   c=sin(theta)*cos(p0->delta)*cos(p0->alpha+dzeta)+cos(theta)*sin(p0->delta);
   
   //__ Update output values - Meeus 20.4
   p.alpha=atan2(a,b)+z;
   p.delta=asin(c);
   
   return p;
}


/*
 * Transformation d'une epoque vers une autre.
 */

sph_circ_t sprec2000(sph_circ_t * p0)
{
   
   double       t;
   double       aeta,api,ap;
   double       a,b,c;
   sph_circ_t  p;
   
   //__ Time in Julian century from epoch 2000.
   t=(p0->jd-EPOCH_2000)/JULIAN_CENTURIES;
   
   
   //__ Accurate reduction of positions from one equinox to another - Meeus 20.6
   aeta=t*(47.0029-t*(0.03302+t*0.000060));
   aeta=aeta/dgtord(3600);
   
   api=174.876384*60-t*(869.8089+t*0.03536);
   api=api/dgtord(3600);
   
   ap=t*(5029.0966+t*(1.11113-t*0.000006));
   ap=ap/dgtord(3600);
   
   //__ See Meeus 20.7
   a=cos(aeta)*cos(p0->pos.b)*sin(api-p0->pos.l)-sin(aeta)*sin(p0->pos.b);
   b=cos(p0->pos.b)*cos(api-p0->pos.l);
   c=cos(aeta)*sin(p0->pos.b)+sin(aeta)*cos(p0->pos.b)*sin(api-p0->pos.l);
   
   //__ Update output values - Meeus 20.4
   p.type=p0->type;
   p.jd=p0->jd;
   p.pos.l=-atan2(a,b)+ap+api;
   p.pos.b=asin(c);
   p.pos.radius=p0->pos.radius;
   
   return p;
}


