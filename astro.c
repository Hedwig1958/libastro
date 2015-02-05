//
//  astro.c
//  astro
//
//  Created by Corto Maltese on 30/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "astro.h"

#define VERBOSE_MOON_PHASES
#define VERBOSE_SEASONS
#define VERBOSE_RON_VONDRAK
#define VERBOSE_SUN_PHYS_OBS_EPHEMERIS
#define VERBOSE_SUN_ABERRATION


/*
 * Apparent Sideral time for Greenwich at any instant.
 */

double appSideraltime(nut_perioterms_t  * nut, double jd)
{
   double t, theta0;
   double d_psi,d_epsilon,epsilon0;
   
   //__ Time in Julian centuries from epoch 2000 - Meeus 11.1
   t=(jd-EPOCH_2000)/JULIAN_CENTURIES;
   
   //__ Mean Sideral Time at Greenwich at 0h UT - Meeus 11.4
   theta0=280.46061837+360.98564736629*(jd-2451545.0)+t*t*(0.000387933-t*(1/38710000));
   
   //__ Corrections for nutation in longitude & obliquity - Meeus 21.1, 21.2
   nutation(nut,jd,&d_psi,&d_epsilon,&epsilon0);
   theta0=fmod(dgtord(theta0),(M_PI*2))+d_psi*cos(epsilon0);
   
   return theta0;
}


/*
 * Transform equatorial coordinate in horizontal coordinate.
 */

horiz_circ_t equtohoriz(nut_perioterms_t  * nut, equ_coord_t * equ, geo_coord_t * geo, double jd, int type)
{
   
   horiz_circ_t   horiz;
   
   //__ Temps Sideral apparent a Greenwich.
   horiz.theta0=appSideraltime(nut,jd);
   
   //__ Angle Horaire.
   horiz.theta=horiz.theta0-geo->longitude;
   horiz.eta=horiz.theta0-geo->longitude-equ->alpha;
   
   //__ Angle Parallactique.
   horiz.q=atan2(sin(horiz.eta),tan(geo->latitude)*cos(equ->delta)-sin(equ->delta)*cos(horiz.eta));
   
   //__ Coordonnees locales.
   horiz.type = type;
   horiz.pos.azimuth = atan2(sin(horiz.eta),cos(horiz.eta)*sin(geo->latitude)-tan(equ->delta)*cos(geo->latitude));
   horiz.pos.altitude = asin(sin(geo->latitude)*sin(equ->delta)+cos(geo->latitude)*cos(equ->delta)*cos(horiz.eta));
   
   return horiz;
   
}


/*
 * Aberration en longitude pour le soleil.
 */

double sunAberration(double radius, double jd)
{
   double tau = (jd-EPOCH_2000)/JULIAN_MILLENIA;
   double dlambda = 3548.193
   + 118.568                   * sin(dgtord(  87.5287 +  359993.7286 * tau ))
   +   2.476                   * sin(dgtord(  85.0561 +  719987.4571 * tau ))
   +   1.376                   * sin(dgtord(  27.8502 + 4452671.1152 * tau ))
   +   0.119                   * sin(dgtord(  73.1375 +  450368.8564 * tau ))
   +   0.114                   * sin(dgtord( 337.2264 +  329644.6718 * tau ))
   +   0.086                   * sin(dgtord( 222.5400 +  659289.3436 * tau ))
   +   0.078                   * sin(dgtord( 162.8136 + 9224659.7915 * tau ))
   +   0.054                   * sin(dgtord(  82.5823 + 1079981.1857 * tau ))
   +   0.052                   * sin(dgtord( 171.5189 +  225184.4282 * tau ))
   +   0.034                   * sin(dgtord(  30.3214 + 4092677.3866 * tau ))
   +   0.033                   * sin(dgtord( 119.8105 +  337181.4711 * tau ))
   +   0.023                   * sin(dgtord( 247.5418 +  299295.6151 * tau ))
   +   0.023                   * sin(dgtord( 325.1526 +  315559.5560 * tau ))
   +   0.021                   * sin(dgtord( 155.1241 +  675553.2846 * tau ))
   +   7.311 * tau             * sin(dgtord( 333.4515 +  359993.7286 * tau ))
   +   0.305 * tau             * sin(dgtord( 330.9814 +  719987.4571 * tau ))
   +   0.010 * tau             * sin(dgtord( 328.5170 + 1079981.1857 * tau ))
   +   0.309 * tau * tau       * sin(dgtord( 241.4518 +  359993.7286 * tau ))
   +   0.021 * tau * tau       * sin(dgtord( 205.0482 +  719987.4571 * tau ))
   +   0.004 * tau * tau       * sin(dgtord( 297.8610 + 4452671.1152 * tau ))
   +   0.010 * tau * tau * tau * sin(dgtord( 154.7066 +  359993.7286 * tau ));
   
   //__ Sun aberratgion.
   double aberr = -LIGHT_TIME*radius*dlambda;
   return sectord(aberr);
   
}


/*
 * Parallactic Angle.
 * eta   : angle Horaire.
 * delta : declinaison.
 * phi   : latitude.
 */

double parallacticAngle(double eta, double delta, double phi)
{
   return sin(eta)/(tan(phi)*cos(delta)-sin(delta)*cos(eta));
}


/*
 * Ephemeris for Physical Observation of the Sun.
 */

helio_circ_t sunEphemeris(planets_data_t * data, double jd)
{
   int    i;
   double x,y;
   double theta;   // Sideral period for daily rotation.
   double iota;    // Inclinaison de l'equateur solaire.
   double kappa;   // longitude du noeud ascendant de l'equateur solaire.
   double aberr;
   
   double tau,eta;
   double d_psi,d_epsilon,epsilon0,epsilon;
   double lambda,lambda_p;
   
   sph_circ_t     e = vsop87(data->earth,jd);
   helio_circ_t   helio;
   
   //__ Iteration for light time.
   for(i=0;i<ITERATION_LIGHT_TIME_SUN;i++)
     {
      tau = LIGHT_TIME*e.pos.radius;
      e = vsop87(data->earth,jd-tau);
     }
   
   //__ Sun aberration.
   aberr = sunAberration(e.pos.radius,jd);
   
   //__ Corrections for nutation in longitude & obliquity - Meeus 21.1, 21.2
   nutation(data->nutation,jd,&d_psi,&d_epsilon,&epsilon0);
   epsilon = epsilon0+d_epsilon;
   lambda = e.pos.l+M_PI+aberr;
   lambda_p = lambda+d_psi;
   
   //__ Calcule theta, iota & kappa.
   theta = (jd-2398220)*((M_PI*2)/25.38);
   iota = dgtord(7.25);
   kappa = dgtord(73.6667)+dgtord(1.3958333)*((jd-2396758)/36525);
   
   //__ Calcule l'angle x & y.
   x = atan(-cos(lambda_p)*tan(epsilon));
   y = atan(-cos(lambda-kappa)*tan(iota));
   
   helio.jd = jd;
   helio.pos.radius = x+y;
   helio.pos.latitude = asin(sin(lambda-kappa)*sin(iota));
   eta = atan2(-sin(lambda-kappa)*cos(iota),-cos(lambda-kappa));
   helio.pos.longitude = eta-theta;
   
#ifdef VERBOSE_SUN_PHYS_OBS_EPHEMERIS
  {
   char               buff[20];
   printf("EPHEMERIS FOR PHYSICAL OBSERVATION OF THE SUN\n");
   printf("Theta       : %s\n",rdstrdms(theta,buff,DG_PI_PI));
   printf("Iota        : %s\n",rdstrdms(iota,buff,DG_PI_PI));
   printf("Kappa       : %s\n",rdstrdms(kappa,buff,DG_PI_PI));
   printf("Radius      : %f\n",e.pos.radius);
   printf("Sun Aberr.  : %s\n",rdstrdms(aberr,buff,DG_PI_PI));
   printf("Epsilon     : %s\n",rdstrdms(epsilon,buff,DG_PI_PI));
   printf("Lambda      : %s\n",rdstrdms(lambda,buff,DG_0_2PI));
   printf("Lambda'     : %s\n",rdstrdms(lambda_p,buff,DG_0_2PI));
   printf("position    : %s\n",rdstrdms(helio.pos.radius,buff,DG_PI_PI));
   printf("latitude    : %s\n",rdstrdms(helio.pos.latitude,buff,DG_PI_PI));
   printf("longitude   : %s\n",rdstrdms(helio.pos.longitude,buff,DG_PI_PI));
  }
#endif
   
   return helio;
   
}


/*
 * Illuminated Fraction of the Moon's disk.
 */

double moonIlluminatedFraction(equ_circ_t sun, equ_circ_t moon)
{
   double psi = acos(sin(sun.equ.delta)*sin(moon.equ.delta)+cos(sun.equ.delta)*cos(moon.equ.delta)*cos(sun.equ.alpha-moon.equ.alpha));
   double iota = atan2(sun.sph.radius*UA*sin(psi),moon.sph.radius-sun.sph.radius*UA*cos(psi));
   
   return (1+cos(iota))/2;
}


/*
 * Phase Angle of the Moon's.
 */

double moonPhaseAngle(equ_circ_t sun, equ_circ_t moon)
{
   double d = cos(moon.equ.delta)*sin(sun.equ.delta)-sin(moon.equ.delta)*cos(sun.equ.delta)* cos(moon.equ.alpha-sun.equ.alpha);
   double n = cos(sun.equ.delta)*sin(moon.equ.alpha-sun.equ.alpha);
   
   return atan2(d,n);
}



/*
 * Position angle of the Moon's bright limb.
 */

double moonPositionAngle(equ_circ_t sun, equ_circ_t moon)
{
   return atan2(cos(sun.equ.delta)*sin(sun.equ.alpha-moon.equ.alpha),sin(sun.equ.delta)*
                cos(moon.equ.delta)-cos(sun.equ.delta)*sin(moon.equ.delta)*cos(sun.equ.alpha-moon.equ.alpha));
}




/*
 * The Ron-Vondrak expression for aberration.
 */

equ_coord_t aberrRonVondrak(equ_coord_t * pos, double jd)
{
   
   equ_coord_t  p;
   double       d_alpha, d_delta;
   double       t=(jd-EPOCH_2000)/JULIAN_CENTURIES;
   double       l_2,l_3,l_4,l_5,l_6,l_7,l_8,lp,d,mp,f;
   double       args[RON_VONDRAK_ELEM];
   double       xp,yp,zp;
   double       c=17314463350.0;
   
   l_2=3.1761467+1021.3285546*t;
   l_3=1.7534703+628.3075849*t;
   l_4=6.2034809+334.0612431*t;
   l_5=0.5995465+52.9690965*t;
   l_6=0.8740168+21.3299095*t;
   l_7=5.4812939+7.4781599*t;
   l_8=5.3118863+3.8133036*t;
   lp=3.8103444+8399.6847337*t;
   d=5.1984667+7771.3771486*t;
   mp=2.3555559+8328.6914289*t;
   f=1.6279052+8433.4661601*t;
   
   //__ Calculs pour les arguments.
   args[0]=l_3;
   args[1]=2*l_3;
   args[2]=l_5;
   args[3]=lp;
   args[4]=3*l_3;
   args[5]=l_6;
   args[6]=f;
   args[7]=lp+mp;
   args[8]=2*l_5;
   args[9]=2*l_3-l_5;
   args[10]=3*l_3-8*l_4+3*l_5;
   args[11]=5*l_3-8*l_4+3*l_5;
   args[12]=2*l_2-l_3;
   args[13]=l_2;
   args[14]=l_7;
   args[15]=l_3-2*l_5;
   args[16]=l_8;
   args[17]=l_3+l_5;
   args[18]=2*l_2-2*l_3;
   args[19]=l_3-l_5;
   args[20]=4*l_3;
   args[21]=3*l_3-2*l_5;
   args[22]=l_2-2*l_3;
   args[23]=2*l_2-3*l_3;
   args[24]=2*l_6;
   args[25]=2*l_2-4*l_3;
   args[26]=3*l_3-2*l_4;
   args[27]=lp+2*d-mp;
   args[28]=8*l_2-12*l_3;
   args[29]=8*l_2-14*l_3;
   args[30]=2*l_4;
   args[31]=3*l_2-4*l_3;
   args[32]=2*l_3-2*l_5;
   args[33]=3*l_2-3*l_3;
   args[34]=2*l_3-2*l_4;
   args[35]=lp-2*d;
   
   //__ Calculs X'.
   xp=((-1719914-2*t)*sin(args[0])-25*cos(args[0]));
   xp+=((6434+141*t)*sin(args[1])+(28007-107*t)*cos(args[1]));
   xp+=(715*sin(args[2]));
   xp+=(715*sin(args[3]));
   xp+=((486-5*t)*sin(args[4])-(236-4*t)*cos(args[4]));
   xp+=(159*sin(args[5]));
   xp+=(0);
   xp+=(39*sin(args[7]));
   xp+=(33*sin(args[8])-10*cos(args[8]));
   xp+=(31*sin(args[9])+cos(args[9]));
   xp+=(8*sin(args[10])-28*cos(args[10]));
   xp+=(8*sin(args[11])-28*cos(args[11]));
   xp+=(21*sin(args[12]));
   xp+=(-19*sin(args[13]));
   xp+=(17*sin(args[14]));
   xp+=(16*sin(args[15]));
   xp+=(16*sin(args[16]));
   xp+=(11*sin(args[17])-cos(args[17]));
   xp+=(-11*cos(args[18]));
   xp+=(-11*sin(args[19])-2*cos(args[19]));
   xp+=(-7*sin(args[20])-8*cos(args[20]));
   xp+=(-10*sin(args[21]));
   xp+=(-9*sin(args[22]));
   xp+=(-9*sin(args[23]));
   xp+=(-9*cos(args[24]));
   xp+=(-9*cos(args[25]));
   xp+=(8*sin(args[26]));
   xp+=(8*sin(args[27]));
   xp+=(-4*sin(args[28])-7*cos(args[28]));
   xp+=(-4*sin(args[29])-7*cos(args[29]));
   xp+=(-6*sin(args[30])-5*cos(args[30]));
   xp+=(-sin(args[31])-cos(args[31]));
   xp+=(4*sin(args[32])-6*cos(args[32]));
   xp+=(-7*cos(args[33]));
   xp+=(5*sin(args[34])-5*cos(args[34]));
   xp+=(5*sin(args[35]));
   
   //__ Calculs Y'.
   yp=((25-13*t)*sin(args[0])+(1578089+156*t)*cos(args[0]));
   yp+=((25697-95*t)*sin(args[1])-(5904-130*t)*cos(args[1]));
   yp+=(6*sin(args[2])-657*cos(args[2]));
   yp+=(-656*cos(args[3]));
   yp+=((-216-4*t)*sin(args[4])-(446+5*t)*cos(args[4]));
   yp+=(2*sin(args[5])-147*cos(args[5]));
   yp+=(26*cos(args[6]));
   yp+=(-36*cos(args[7]));
   yp+=(-9*sin(args[8])-30*cos(args[8]));
   yp+=(sin(args[9])-28*cos(args[9]));
   yp+=(25*sin(args[10])+8*cos(args[10]));
   yp+=(-25*sin(args[11])-8*cos(args[11]));
   yp+=(-19*cos(args[12]));
   yp+=(17*cos(args[13]));
   yp+=(-16*cos(args[14]));
   yp+=(15*cos(args[15]));
   yp+=(sin(args[16])-15*cos(args[16]));
   yp+=(-sin(args[17])-10*cos(args[17]));
   yp+=(-10*sin(args[18]));
   yp+=(-2*sin(args[19])+9*cos(args[19]));
   yp+=(-8*sin(args[20])+6*cos(args[20]));
   yp+=(9*cos(args[21]));
   yp+=(-9*cos(args[22]));
   yp+=(-8*cos(args[23]));
   yp+=(-8*sin(args[24]));
   yp+=(8*sin(args[25]));
   yp+=(-8*cos(args[26]));
   yp+=(-7*cos(args[27]));
   yp+=(-6*sin(args[28])+4*cos(args[28]));
   yp+=(6*sin(args[29])-4*cos(args[29]));
   yp+=(-4*sin(args[30])+5*cos(args[30]));
   yp+=(-2*sin(args[31])-7*cos(args[31]));
   yp+=(-5*sin(args[32])-4*cos(args[32]));
   yp+=(-6*sin(args[33]));
   yp+=(-4*sin(args[34])-5*cos(args[34]));
   yp+=(-5*cos(args[35]));
   
   //__ Calculs Z'.
   zp=((10+32*t)*sin(args[0])+(684185-358*t)*cos(args[0]));
   zp+=((11141-48*t)*sin(args[1])-(2559-55*t)*cos(args[1]));
   zp+=(-15*sin(args[2])-282*cos(args[2]));
   zp+=(-285*cos(args[3]));
   zp+=(-94*sin(args[4])-193*cos(args[4]));
   zp+=(-6*sin(args[5])-61*cos(args[5]));
   zp+=(-59*cos(args[6]));
   zp+=(-16*cos(args[7]));
   zp+=(-5*sin(args[8])-13*cos(args[8]));
   zp+=(-12*cos(args[9]));
   zp+=(11*sin(args[10])+3*cos(args[10]));
   zp+=(-11*sin(args[11])-3*cos(args[11]));
   zp+=(-8*cos(args[12]));
   zp+=(8*cos(args[13]));
   zp+=(-7*cos(args[14]));
   zp+=(sin(args[15])+7*cos(args[15]));
   zp+=(-3*sin(args[16])-6*cos(args[16]));
   zp+=(-sin(args[17])-5*cos(args[17]));
   zp+=(-4*sin(args[18]));
   zp+=(-sin(args[19])+4*cos(args[19]));
   zp+=(-3*sin(args[20])+3*cos(args[20]));
   zp+=(4*cos(args[21]));
   zp+=(-4*cos(args[22]));
   zp+=(-4*cos(args[23]));
   zp+=(-3*sin(args[24]));
   zp+=(3*sin(args[25]));
   zp+=(-3*cos(args[26]));
   zp+=(-3*cos(args[27]));
   zp+=(-3*sin(args[28])+2*cos(args[28]));
   zp+=(3*sin(args[29])-2*cos(args[29]));
   zp+=(-2*sin(args[30])+2*cos(args[30]));
   zp+=(sin(args[31])-4*cos(args[31]));
   zp+=(-2*sin(args[32])-2*cos(args[32]));
   zp+=(-3*sin(args[33]));
   zp+=(-2*sin(args[34])-2*cos(args[34]));
   zp+=(-2*cos(args[35]));
   
#ifdef VERBOSE_RON_VONDRAK
   printf("\n*** RON VONDRAK ***\n\n");
   printf("T: %.12f\n\n",t);
   printf("L2: %.12f\n",l_2);
   printf("L3: %.12f\n",l_3);
   printf("L4: %.12f\n",l_4);
   printf("L5: %.12f\n",l_5);
   printf("L6: %.12f\n",l_6);
   printf("L7: %.12f\n",l_7);
   printf("L8: %.12f\n",l_8);
   printf("lp: %.12f\n",lp);
   printf("d: %.12f\n",d);
   printf("mp: %.12f\n",mp);
   printf("f: %.12f\n\n",f);
   printf("X': %.12f\n",xp);
   printf("Y': %.12f\n",yp);
   printf("Z': %.12f\n",zp);
#endif
   
   d_alpha = (yp*cos(pos->alpha)-xp*sin(pos->alpha))/(c*cos(pos->delta));
   d_delta = -((xp*cos(pos->alpha)+yp*sin(pos->alpha))*sin(pos->delta)-zp*cos(pos->delta))/c;
   
   p.alpha = pos->alpha+d_alpha;
   p.delta = pos->delta+d_delta;
   
   return p;
}



/*
 * Calculation for apparent geocentric position of the planets.
 *
 * Ref : Astronomical algorithms, Jean Meeus.
 *
 */

equ_circ_t appGeoPosPlanets(planets_data_t * data, double jd)
{
   double       x = 0.0,y = 0.0,z = 0.0;
   double       sun_theta,sun_beta;
   double       aberr_e,aberr_pi;
   double       d_lambda,d_beta;
   double       radius = 0.0;
   double       lambda,lambda_p,beta;
   double       kappa = 20.49552;
   double       d_psi,d_epsilon,epsilon0,epsilon;
   double       t,tau = 0.0;
   
   sph_circ_t   sphEarth = vsop87(data->earth,jd);
   sph_circ_t   sphPlanet;
   equ_circ_t   equPLanet;

   
   //__ Time in Julian centuries from epoch 2000.
   t = (jd-EPOCH_2000)/JULIAN_CENTURIES;
   
   //__ Iteration for light-time.
   for(int i=0;i<ITERATION_LIGHT_TIME_PLANETS;i++)
     {
      sphPlanet = vsop87(data->planet,jd-tau);
      
      //__ Coordinate x,y & z.
      x = sphPlanet.pos.radius*cos(sphPlanet.pos.b)*cos(sphPlanet.pos.l)-sphEarth.pos.radius*cos(sphEarth.pos.b)*cos(sphEarth.pos.l);
      y = sphPlanet.pos.radius*cos(sphPlanet.pos.b)*sin(sphPlanet.pos.l)-sphEarth.pos.radius*cos(sphEarth.pos.b)*sin(sphEarth.pos.l);
      z = sphPlanet.pos.radius*sin(sphPlanet.pos.b)                     -sphEarth.pos.radius*sin(sphEarth.pos.b);
      
      //__ Radius & light-time correction.
      radius = sqrt(x*x+y*y+z*z);
      tau = LIGHT_TIME*radius;
      
      //__ Save first radius.
      if(i==0)
         equPLanet.sph.radius = radius;
     }

   //__ Sun longitude 'theta' & latitude 'beta'.
   sun_theta = sphEarth.pos.l+M_PI;
   sun_beta = -sphEarth.pos.b;
   
   //__ Planet lambda & beta.
   lambda = atan2(y,x);
   beta = atan(z/sqrt(x*x+y*y));
   
   //__ The effect of aberration - Meeus 22.2.
   aberr_e = 0.016708617-t*(0.000042037-0.0000001236*t);
   aberr_pi = dgtord(102.93735+t*(0.71953+0.00046*t));
   d_lambda = -kappa*cos(sun_theta-lambda)+aberr_e*kappa*cos(aberr_pi-lambda)/cos(beta);
   d_beta = -kappa*sin(beta)*(sin(sun_theta-lambda)-aberr_e*sin(aberr_pi-lambda));
   lambda += sectord(d_lambda);
   beta += sectord(d_beta);
   
   //__ Corrections for reduction to the FK5 system - Meeus 31.3
   lambda_p = dgtord(rdtodg(lambda)-1.397*t-0.00031*t*t);
   d_lambda = -0.09033+0.03916*(cos(lambda_p)+sin(lambda_p))*tan(beta);
   d_beta = 0.03916*(cos(lambda_p)-sin(lambda_p));
   lambda += sectord(d_lambda);
   beta += sectord(d_beta);
   
   //__ Corrections for nutation in longitude & obliquity - Meeus 21.1, 21.2
   nutation(data->nutation,jd,&d_psi,&d_epsilon,&epsilon0);
   epsilon = epsilon0+d_epsilon;
   lambda += d_psi;
   
   //__ Transformations from ecliptical into equatorial coordinates - Meeus 12.3, 12.4
   equPLanet.type = PLANET;
   equPLanet.jd = jd;
   equPLanet.equ.alpha = rdreduce(atan2((sin(lambda)*cos(epsilon)-tan(beta)*sin(epsilon)),cos(lambda)),ANGLE_2PI);
   equPLanet.equ.delta = rdreduce(asin(sin(beta)*cos(epsilon)+cos(beta)*sin(epsilon)*sin(lambda)),ANGLE_PI);
   equPLanet.sph.l = rdreduce(lambda,ANGLE_2PI);
   equPLanet.sph.b = rdreduce(beta,ANGLE_PI);
   equPLanet.parallax = asin(EARTH_RADIUS/(equPLanet.sph.radius*UA));
   equPLanet.semidiam = 0.004652418/equPLanet.sph.radius;
   
   return equPLanet;
   
}


/*
 * Calculation for apparent geocentric position of the sun.
 *
 * Ref : Astronomical algorithms, Jean Meeus.
 *
 */

equ_circ_t appGeoPosSun(planets_data_t * data, double jd)
{
   double       t;
   double       lambda,lambda_p,theta,beta;
   double       d_theta,d_beta;
   double       d_psi,d_epsilon,epsilon0,epsilon;
   double       aberr;
   
   sph_circ_t   sphEarth = vsop87(data->earth,jd);
   equ_circ_t   equSun;
   
   //__ Time in Julian centuries from epoch 2000.
   t = (jd-EPOCH_2000)/JULIAN_CENTURIES;
   
   //__ Sun Radius.
   equSun.sph.radius = sphEarth.pos.radius;
   
   //__ Corrections for sun aberration.
   aberr = sunAberration(sphEarth.pos.radius,jd);

   //__ Sun longitude 'theta' & latitude 'beta'.
   theta = sphEarth.pos.l+M_PI+aberr;
   beta = -sphEarth.pos.b;
   
   //__ Corrections for reduction to the FK5 system - Meeus 31.3.
   lambda_p = dgtord(rdtodg(theta)-1.397*t-0.00031*t*t);
   d_theta = sectord(-0.09033);
   d_beta = sectord(0.03916*(cos(lambda_p)-sin(lambda_p)));
   theta += d_theta;
   beta += d_beta;
   
   
   //__ Corrections for nutation in longitude & obliquity - Meeus 21.1, 21.2.
   nutation(data->nutation,jd,&d_psi,&d_epsilon,&epsilon0);
   epsilon = epsilon0+d_epsilon;
   lambda = theta+d_psi;
   
   //__ Transformations from ecliptical into equatorial coordinates - Meeus 12.3, 12.4.
   equSun.type = SUN;
   equSun.jd = jd;
   equSun.equ.alpha = rdreduce(atan2((sin(lambda)*cos(epsilon)-tan(beta)*sin(epsilon)),cos(lambda)),ANGLE_2PI);
   equSun.equ.delta = rdreduce(asin(sin(beta)*cos(epsilon)+cos(beta)*sin(epsilon)*sin(lambda)),ANGLE_PI);
   equSun.sph.l = rdreduce(lambda,ANGLE_2PI);
   equSun.sph.b = rdreduce(beta,ANGLE_PI);
   equSun.parallax = asin(EARTH_RADIUS/(equSun.sph.radius*UA));
   equSun.semidiam = 0.004652418/equSun.sph.radius;
   
   return equSun;
   
}



/*
 * Calculation for apparent geocentric position of the moon.
 *
 * Ref : Astronomical algorithms, Jean Meeus.
 *
 */

equ_circ_t appGeoPosMoon(lunar_data_t * data, double jd, lunar_shift_t * ls)
{
   double       t,tau;
   double       lambda,lambda_p,beta;
   double       d_theta,d_beta;
   double       d_psi,d_epsilon,epsilon0,epsilon;
   
   sph_circ_t   sphMoon = elp2000(data->moon,jd);
   equ_circ_t   equMoon;
   
   
   //__ Time in Julian centuries from epoch 2000.
   t = (jd-EPOCH_2000)/JULIAN_CENTURIES;
   
   //__ Iteration for light-time.
   equMoon.sph.radius = sphMoon.pos.radius;
   
   for(int i=0;i<ITERATION_LIGHT_TIME_MOON;i++)
     {
      tau = LIGHT_TIME*sphMoon.pos.radius/UA;
      sphMoon = elp2000(data->moon,jd-tau);
     }
   
   //__ Moon longitude 'lambda' & latitude 'beta'.
   lambda = sphMoon.pos.l;
   beta = sphMoon.pos.b;
   
   //__ Corrections for reduction to the FK5 system - Meeus 31.3.
   lambda_p = dgtord((rdtodg(lambda))-t*(1.397-0.00031*t));
   d_theta = sectord(-0.09033);
   d_beta = sectord(0.03916*(cos(lambda_p)-sin(lambda_p)));
   lambda += d_theta;
   beta += d_beta;
   
   //__ Corrections for nutation in longitude & obliquity - Meeus 21.1, 21.2.
   nutation(data->nutation,jd,&d_psi,&d_epsilon,&epsilon0);
   lambda += d_psi;
   epsilon = epsilon0+d_epsilon;
   
   //__ Lunar shift.
   lambda += sectord(ls->lambda*cos(epsilon));
   beta += sectord(ls->beta);
   
   //__ Transformations from ecliptical into equatorial coordinates - Meeus 12.3, 12.4.
   equMoon.type = MOON;
   equMoon.jd = jd;
   equMoon.equ.alpha = rdreduce(atan2((sin(lambda)*cos(epsilon)-tan(beta)*sin(epsilon)),cos(lambda)),ANGLE_2PI);
   equMoon.equ.delta = rdreduce(asin(sin(beta)*cos(epsilon)+cos(beta)*sin(epsilon)*sin(lambda)),ANGLE_PI);
   equMoon.sph.l = rdreduce(lambda,ANGLE_2PI);
   equMoon.sph.b = rdreduce(beta,ANGLE_PI);
   equMoon.parallax = asin(EARTH_RADIUS/equMoon.sph.radius);
   equMoon.semidiam = asin(K_0_UMBRA_IAU*EARTH_RADIUS/equMoon.sph.radius);
   
   return equMoon;
   
}


/*
 * Calculation for apparent geocentric position of the stars.
 *
 * Ref : Astronomical algorithms, Jean Meeus.
 *
 */

equ_circ_t appGeoPosStars(nut_perioterms_t  * nut, star_rec_t * rec, double jd)
{
   double      tau;
   double      d_psi, d_epsilon, epsilon0;
   double      d_alpha, d_delta;
   equ_circ_t  equStar;
   
   //__ Time in julian year from epoch 2000.
   tau = (jd-EPOCH_2000)/JULIAN_YEARS;
   
   //__ Stars position and Proper motion from epoch 2000.
   equStar.type = STARS;
   equStar.jd = jd;
   equStar.equ.alpha = rec->ra+(hmstord(rec->pm_ra*tau/(3600*15)));
   equStar.equ.delta = rec->dec+(dgtord(rec->pm_dec*tau/3600));
   
   //__ Aberration, 'Ron-Vondrak'.
   equStar.equ = aberrRonVondrak(&equStar.equ,jd);
   
   //__ Precession.
   equStar.equ = eprec2000(jd,&equStar.equ);
   
   //__ Corrections for nutation in longitude & obliquity - Meeus 21.1, 21.2.
   nutation(nut,jd,&d_psi,&d_epsilon,&epsilon0);
   d_alpha = (cos(epsilon0)+sin(epsilon0)*sin(equStar.equ.alpha)*tan(equStar.equ.delta))*d_psi-(cos(equStar.equ.alpha)*tan(equStar.equ.delta))*d_epsilon;
   d_delta = sin(epsilon0)*cos(equStar.equ.alpha)*d_psi+sin(equStar.equ.alpha)*d_epsilon;
   
   equStar.equ.alpha += d_alpha;
   equStar.equ.delta += d_delta;
   
   return equStar;
   
}


/*
 * The geometric altitude of the center of the body at the time of
 * apparent rising or setting.
 *
 * ref : Guide des donnŽes astronomiques.
 *       EDP sciences.
 */

double altitude0(int type, double parallax, double refraction, double semi_diam, double altitude_obs, double distance_obs, double altitude_hor)
{
   double n1 = 0.0;
   double n2 = 0.0;
   double h0 = 0.0;
   
   //__ observer altitude.
   n1=acos(EARTH_RADIUS_METER/(EARTH_RADIUS_METER+altitude_obs));
   
   //__ horizon altitude.
   if (altitude_hor!=0.0 && distance_obs!=0.0)
      n2 = atan(altitude_hor/distance_obs);
   
   //__ Determination of 'h0'.
   switch(type)
  {
     case PLANET:
     case STARS:
      h0 = -refraction-n1+n2;
      break;
     case SUN:
      h0 = -refraction-semi_diam-n1+n2;
      break;
     case MOON:
      h0 = parallax-refraction-semi_diam-n1+n2;
      break;
  }
   
   //__ Transfert result.
   return h0;
   
}


/*
 * Calculation for sun or moon transit.
 *
 * ref : Guide des donnŽes astronomiques.
 *       EDP sciences.
 *
 * hour angle : eta
 */

ephem_data_t transit(nut_perioterms_t  * nut, equ_circ_t * equ, geo_circ_t * geo, double eta)
{
   double         n = 0.0;
   double         m = 0.0;
   
   equ_circ_t     tmp;
   ephem_data_t   ephem;
   
   //__ Apparent Sideral time for Greenwich at any instant.
   double theta0 = appSideraltime(nut,equ[EPHEM_DAY].jd);
   tmp.equ.alpha = equ[EPHEM_DAY].equ.alpha;
   tmp.equ.delta = equ[EPHEM_DAY].equ.delta;
   
   //__ check if body is circumpolar.
   fbool circumpolar = fabs(cos(geo->pos.latitude)*cos(tmp.equ.delta))>1.0;
   
   for(int i=0;i<EPHEM_NUM_ITERATION;i++)
     {
      
      //__ Raw transit time.
      m = tmp.equ.alpha+geo->pos.longitude-theta0+eta;
      n = (SIDERAL_CST_STTOMT*m)/(M_PI*2);
      
      //__ Interpolation for alpha & delta.
      tmp.equ.alpha = interpolationBessel(equ[0].equ.alpha,equ[1].equ.alpha,equ[2].equ.alpha,equ[3].equ.alpha,n);
      tmp.equ.delta = interpolationBessel(equ[0].equ.delta,equ[1].equ.delta,equ[2].equ.delta,equ[3].equ.delta,n);
      
     }
   
   //__ Transfert result.
   ephem.jd = equ[EPHEM_DAY].jd + n;
   ephem.type = equ[EPHEM_DAY].type;
   ephem.circumpolar =  circumpolar;
   ephem.eta = eta;
   ephem.theta0 = theta0;
   ephem.azimuth = 0.0;
   ephem.pos.alpha = tmp.equ.alpha;
   ephem.pos.delta = tmp.equ.delta;
   
   return ephem;
   
}


/*
 * Calculation for sun or moon rising.
 *
 * ref : Guide des donnŽes astronomiques
 *       EDP sciences.
 */

ephem_data_t rising(nut_perioterms_t  * nut, equ_circ_t * equ, geo_circ_t * geo, double h0)
{
   double         n = 0.0;
   double         m = 0.0;
   double         eta0 = 0.0;
   double         azimuth = 0.0;
   
   equ_circ_t     tmp;
   ephem_data_t   ephem;
   
   
   //__ Apparent Sideral time for Greenwich at any instant.
   double theta0 = appSideraltime(nut,equ[EPHEM_DAY].jd);
   tmp.equ.alpha = equ[EPHEM_DAY].equ.alpha;
   tmp.equ.delta = equ[EPHEM_DAY].equ.delta;
   
   //__ check if body is circumpolar.
   fbool circumpolar = fabs(cos(geo->pos.latitude)*cos(tmp.equ.delta))>1.0;
   
   if (!circumpolar)
     {
      for(int i=0;i<EPHEM_NUM_ITERATION;i++)
        {
         
         //__ Approximate rising time.
         eta0 = acos((sin(h0)-sin(geo->pos.latitude)*sin(tmp.equ.delta))/(cos(geo->pos.latitude)*cos(tmp.equ.delta)));
         m = tmp.equ.alpha-eta0;
         n = (m+geo->pos.longitude-theta0)/(M_PI*2);
         n = n*SIDERAL_CST_STTOMT;
         
         //__ Interpolation for alpha & delta.
         tmp.equ.alpha = interpolationBessel(equ[0].equ.alpha,equ[1].equ.alpha,equ[2].equ.alpha,equ[3].equ.alpha,n);
         tmp.equ.delta = interpolationBessel(equ[0].equ.delta,equ[1].equ.delta,equ[2].equ.delta,equ[3].equ.delta,n);
        }
      
      //__ Compute the azimuth.
      azimuth = -acos((sin(h0*tan(geo->pos.latitude)-(sin(tmp.equ.delta/cos(geo->pos.latitude)))))/cos(h0));
     }
   
   //__ Transfert result.
   ephem.jd = equ[EPHEM_DAY].jd + n;
   ephem.type = equ[EPHEM_DAY].type;
   ephem.circumpolar =  circumpolar;
   ephem.h0 = h0;
   ephem.eta = eta0;
   ephem.theta0 = theta0;
   ephem.azimuth = azimuth;
   ephem.pos.alpha = tmp.equ.alpha;
   ephem.pos.delta = tmp.equ.delta;
   
   return ephem;
}


/*
 * Calculation for sun or moon setting.
 *
 * ref : Guide des donnŽes astronomiques
 *       EDP sciences
 */

ephem_data_t setting(nut_perioterms_t  * nut, equ_circ_t * equ, geo_circ_t * geo, double h0)
{
   double         n = 0.0;
   double         m = 0.0;
   double         eta0 = 0.0;
   double         azimuth = 0.0;
   
   equ_circ_t     tmp;
   ephem_data_t   ephem;

   //__ Apparent Sideral time for Greenwich at any instant.
   double theta0 = appSideraltime(nut,equ[EPHEM_DAY].jd);
   tmp.equ.alpha = equ[EPHEM_DAY].equ.alpha;
   tmp.equ.delta = equ[EPHEM_DAY].equ.delta;
   
   //__ check if body is circumpolar.
   fbool circumpolar = fabs(cos(geo->pos.latitude)*cos(tmp.equ.delta))>1.0;
   
   if (!circumpolar)
     {
      for(int i=0;i<EPHEM_NUM_ITERATION;i++)
        {
         
         //__ Approximate rising time.
         eta0 = acos((sin(h0)-sin(geo->pos.latitude)*sin(tmp.equ.delta))/(cos(geo->pos.latitude)*cos(tmp.equ.delta)));
         m = tmp.equ.alpha+eta0;
         n = (m+geo->pos.longitude-theta0)/(M_PI*2);
         n = n*SIDERAL_CST_STTOMT;
         
         //__ Interpolation for alpha & delta.
         tmp.equ.alpha = interpolationBessel(equ[0].equ.alpha,equ[1].equ.alpha,equ[2].equ.alpha,equ[3].equ.alpha,n);
         tmp.equ.delta = interpolationBessel(equ[0].equ.delta,equ[1].equ.delta,equ[2].equ.delta,equ[3].equ.delta,n);
        }
      
      //__ Compute the azimuth.
      azimuth = acos((sin(h0*tan(geo->pos.latitude)-(sin(tmp.equ.delta/cos(geo->pos.latitude)))))/cos(h0));
     }
   
   //__ Transfert result.
   ephem.jd = equ[EPHEM_DAY].jd + n;
   ephem.type = equ[EPHEM_DAY].type;
   ephem.circumpolar =  circumpolar;
   ephem.h0 = h0;
   ephem.eta = eta0;
   ephem.theta0 = theta0;
   ephem.azimuth = azimuth;
   ephem.pos.alpha = tmp.equ.alpha;
   ephem.pos.delta = tmp.equ.delta;

   return ephem;
}



/*
 * Seasons.
 */
int * _seasons(void)
{
   static int k;
   return &k;
}



/*
 * Seasons longitude.
 */

double _seasons_longitude(double jd)
{
   double l = appGeoPosSun(*_planets_data(),jd).sph.l;
   double longitude = 0.0;
   
   switch(*_seasons())
  {
     case SPRING:
      longitude = -cos(l);
     case SUMMER:
      longitude = -sin(l);
     case AUTUMN:
      longitude = cos(l);
     case WINTER:
      longitude = sin(l);
  }
   
   return longitude;
}


/*
 * Equinoxes & Solstices.
 */

seasons_data_t seasons(planets_data_t * eart, double year, int type, int accurate)
{
   double          tau = 0.0;
   double          jd0 = 0.0;
   double          jd = 0.0;
   double          y,t,w,dl,s;
   equ_circ_t      pos;
   seasons_data_t  elem;
   
   //__ Init type.
   elem.type = type;
   
   if(year>=1000)
      y = (year-2000)/1000;
   else
      y = year/1000;
   
   /*
    * Select season, Spring, Summer, Autumn, winter.
    */
   
   switch(type)
  {
      
     case SPRING:
      
      //__ March equinox, Jean Meeus 26.[A,B].1
      if(year>=1000)
         jd0 = 2451623.80984+y*(365242.37404+y*(0.05169-y*(0.00411-y*0.00057)));
      else
         jd0 = 1721139.29189+y*(365242.13740+y*(0.06134+y*(0.00111-y*0.00071)));
      
      break;
      
     case SUMMER:
      
      //__ June solstice, Jean Meeus 26.[A,B].2
      if(year>=1000)
         jd0 = 2451716.56767+y*(365241.62603+y*(0.00325+y*(0.00888-y*0.00030)));
      else
         jd0 = 1721233.25401+y*(365241.72562-y*(0.05323+y*(0.00907-y*0.00025)));
      
      break;
      
     case AUTUMN:
      
      //__ September equinox, Jean Meeus 26.[A,B].3
      if(year>=1000)
         jd0 = 2451810.21715+y*(365242.01767-y*(0.11575+y*(0.00337+y*0.00078)));
      else
         jd0 = 1721325.70455+y*(365242.49558-y*(0.11677-y*(0.00297-y*0.00074)));
      
      break;
      
     case WINTER:
      
      //__ December solstice, Jean Meeus 26.[A,B].4
      if(year>=1000)
         jd0 = 2451900.05952+y*(365242.74049-y*(0.06223-y*(0.00823+y*0.00032)));
      else
         jd0 = 1721414.39987+y*(365242.88257-y*(0.00769-y*(0.00933-y*0.00006)));
      
      break;
      
  }
   
   //__ Jean Meeus 26.C
   t = (jd0-2451545.0)/36525;
   w = 35999.373*t-2.47;
   dl = 1+0.0334*cos(dgtord(w))+0.0007*cos(dgtord(2*w));
   
   s = 485*cos(dgtord(324.96+  1934.136*t))
   +203*cos(dgtord(337.23+ 32964.467*t))
   +199*cos(dgtord(342.08+    20.186*t))
   +182*cos(dgtord( 27.85+445267.112*t))
   +156*cos(dgtord( 73.14+ 45036.886*t))
   +136*cos(dgtord(171.52+ 22518.443*t))
   + 77*cos(dgtord(222.54+ 65928.934*t))
   + 74*cos(dgtord(296.72+  3034.906*t))
   + 70*cos(dgtord(243.58+  9037.513*t))
   + 58*cos(dgtord(119.81+ 33718.147*t))
   + 52*cos(dgtord(297.17+   150.678*t))
   + 50*cos(dgtord( 21.02+  2281.226*t))
   + 45*cos(dgtord(247.54+ 29929.562*t))
   + 44*cos(dgtord(325.15+ 31555.956*t))
   + 29*cos(dgtord( 60.93+  4443.417*t))
   + 18*cos(dgtord(155.12+ 67555.328*t))
   + 17*cos(dgtord(288.79+  4562.452*t))
   + 16*cos(dgtord(198.04+ 62894.029*t))
   + 14*cos(dgtord(199.76+ 31436.921*t))
   + 12*cos(dgtord( 95.39+ 14577.848*t))
   + 12*cos(dgtord(287.11+ 31931.756*t))
   + 12*cos(dgtord(320.81+ 34777.259*t))
   +  9*cos(dgtord(227.73+  1222.114*t))
   +  8*cos(dgtord( 15.45+ 16859.074*t));
   
   jd=jd0+(0.00001*s)/dl;
   
   //__ Higher accurancy.
   for (int i=0;i<=(accurate ? HIGH_ACCURATE_FOR_SEASONS : LOW_ACCURATE_FOR_SEASONS);i++)
     {
      jd+=tau;
      pos = appGeoPosSun(eart,jd);
      tau = 58*sin(M_PI/2.0*type - pos.sph.l);
     }
   
   //__ transfert result value.
   elem.jd = jd;
   elem.circ = pos;
   
   return elem;
}


/*
 * Use for the compute for lunar phases.
 */

double * _phases(void)
{
   static double k;
   return &k;
}


/*
 * for lunar phases longitude.
 */

double _lunarphases_longitude(double jd)
{
   
   lunar_shift_t lf = {0.0, 0.0};
   
   double ls = appGeoPosSun(*_planets_data(),jd).sph.l;
   double lm = appGeoPosMoon(*_lunar_data(),jd,&lf).sph.l;
   double p = *_phases();
   
   return pow(sin(ls)-sin(lm+p*(M_PI*2)),2);
}


/*
 * Calcul des phases lunaires.
 */

double lunarphases(planets_data_t * d1, lunar_data_t * d2, double phase, double year, int accurate)
{
   
   double       k;
   double       jd,t,m,mp,f,e,w,omega;
   double       corr = 0.0;
   double       c_add;
   date_t       date;
   
   /*___ Jean Meeus 47.2 ___*/
   k = floor((year-2000)*12.3685);
   k += phase;
   
#ifdef VERBOSE_MOON_PHASES
   printf("k:%f\n",k);
#endif
   
   /*___ Jean Meeus 47.3 ___*/
   t = k/1236.85;
   
   /*___ Calcul approximatif du jd pour les phases de la lune - Jean Meeus 47.1 ___*/
   jd = 2451550.09765+29.530588853*k+t*t*(0.0001337-t*(0.000000150+t*0.00000000073));
   
   /*___ Sun's mean anomaly - Jean Meeus 47.4 ___*/
   m = 2.5534+29.10535669*k-t*t*(0.0000218+0.00000011*t);
   
   /*___ Moon's mean anomaly - Jean Meeus 47.5 ___*/
   mp = 201.5643+385.81693528*k+t*t*(0.0107438+0.00001239*t-0.000000058*t*t);
   
   /*___ Moon's arument of latitude - Jean Meeus 47.6---*/
   f = 160.7108+390.67050274*k-t*t*(0.0016341+0.00000227*t-0.000000011*t*t);
   
   /*___ Longitude of the ascending node of the lunar orbit - Jean Meeus 47.7 ___*/
   omega = 124.7746-1.56375580*k+t*t*(0.0020691+0.00000215*t);
   
   /*___ Eccentricity of the Earth - Jean Meeus 45.6 ___*/
   e = 1-t*(0.002516+0.0000074*t);
   
   /*___ Additional corrections ___*/
   w = 0.00306-0.00038*e*cos(dgtord(m))+0.00026*cos(dgtord(mp))-0.00002*cos(dgtord(mp-m))
   +0.00002*cos(dgtord(mp+m))+0.00002*cos(dgtord(2*f));
   
   /*
    _________________________________
    Planetary arguments & Additionnal
    corrections for all phases
    _________________________________
    */
   
   c_add = 0.000325*sin(dgtord(299.77+0.107408*k-0.009173*t*t))
   +0.000165*sin(dgtord(251.88+ 0.016321*k))
   +0.000164*sin(dgtord(251.83+26.651886*k))
   +0.000126*sin(dgtord(349.42+36.412478*k))
   +0.000110*sin(dgtord( 84.66+18.206239*k))
   +0.000062*sin(dgtord(141.74+53.303771*k))
   +0.000060*sin(dgtord(207.14+ 2.453732*k))
   +0.000056*sin(dgtord(154.84+ 7.306860*k))
   +0.000047*sin(dgtord( 34.52+27.261239*k))
   +0.000042*sin(dgtord(207.19+ 0.121824*k))
   +0.000040*sin(dgtord(291.34+ 1.844379*k))
   +0.000037*sin(dgtord(161.72+24.198154*k))
   +0.000035*sin(dgtord(239.56+25.513099*k))
   +0.000023*sin(dgtord(331.55+ 3.592518*k));
   
   
   /*___ New Moon ___*/
   if(phase==NEW_MOON)
      corr = -0.40720    *sin(dgtord(mp))
      +0.17241*e  *sin(dgtord(m))
      +0.01608    *sin(dgtord(2*mp))
      +0.01039    *sin(dgtord(2*f))
      +0.00739*e  *sin(dgtord(mp-m))
      -0.00514*e  *sin(dgtord(mp+m))
      +0.00208*e*e*sin(dgtord(2*m))
      -0.00111    *sin(dgtord(mp-2*f))
      -0.00057    *sin(dgtord(mp+2*f))
      +0.00056*e  *sin(dgtord(2*mp+m))
      -0.00042    *sin(dgtord(3*mp))
      +0.00042*e  *sin(dgtord(m+2*f))
      +0.00038*e  *sin(dgtord(m-2*f))
      -0.00024*e  *sin(dgtord(2*mp-m))
      -0.00017    *sin(dgtord(omega))
      -0.00007    *sin(dgtord(mp+2*m))
      +0.00004    *sin(dgtord(2*mp-2*f))
      +0.00004    *sin(dgtord(3*m))
      +0.00003    *sin(dgtord(mp+m-2*f))
      +0.00003    *sin(dgtord(2*mp+2*f))
      -0.00003    *sin(dgtord(mp+m+2*f))
      +0.00003    *sin(dgtord(mp-m+2*f))
      -0.00002    *sin(dgtord(mp-m-2*f))
      -0.00002    *sin(dgtord(3*mp+m))
      +0.00002    *sin(dgtord(4*mp));
   
   /*___ Full Moon ___*/
   if(phase==FULL_MOON)
      corr = -0.40614    *sin(dgtord(mp))
      +0.17302*e  *sin(dgtord(m))
      +0.01614    *sin(dgtord(2*mp))
      +0.01043    *sin(dgtord(2*f))
      +0.00734*e  *sin(dgtord(mp-m))
      -0.00515*e  *sin(dgtord(mp+m))
      +0.00209*e*e*sin(dgtord(2*m))
      -0.00111    *sin(dgtord(mp-2*f))
      -0.00057    *sin(dgtord(mp+2*f))
      +0.00056*e  *sin(dgtord(2*mp+m))
      -0.00042    *sin(dgtord(3*mp))
      +0.00042*e  *sin(dgtord(m+2*f))
      +0.00038*e  *sin(dgtord(m-2*f))
      -0.00024*e  *sin(dgtord(2*mp-m))
      -0.00017    *sin(dgtord(omega))
      -0.00007    *sin(dgtord(mp+2*m))
      +0.00004    *sin(dgtord(2*mp-2*f))
      +0.00004    *sin(dgtord(3*m))
      +0.00003    *sin(dgtord(mp+m-2*f))
      +0.00003    *sin(dgtord(2*mp+2*f))
      -0.00003    *sin(dgtord(mp+m+2*f))
      +0.00003    *sin(dgtord(mp-m+2*f))
      -0.00002    *sin(dgtord(mp-m-2*f))
      -0.00002    *sin(dgtord(3*mp+m))
      +0.00002    *sin(dgtord(4*mp));
   
   
   /*___ First & Last quarters ___*/
   if(phase==FIRST_QUARTER||phase==LAST_QUARTER)
      corr = -0.62801    *sin(dgtord(mp))
      +0.17172*e  *sin(dgtord(m))
      -0.01183*e  *sin(dgtord(mp+m))
      +0.00862    *sin(dgtord(2*mp))
      +0.00804    *sin(dgtord(2*f))
      +0.00454*e  *sin(dgtord(mp-m))
      +0.00204*e*e*sin(dgtord(2*m))
      -0.00180    *sin(dgtord(mp-2*f))
      -0.00070    *sin(dgtord(mp+2*f))
      -0.00040    *sin(dgtord(3*mp))
      -0.00034*e  *sin(dgtord(2*mp-m))
      +0.00032*e  *sin(dgtord(m+2*f))
      +0.00032*e  *sin(dgtord(m-2*f))
      -0.00028*e  *sin(dgtord(mp+2*m))
      +0.00027*e  *sin(dgtord(2*mp+m))
      -0.00017    *sin(dgtord(omega))
      -0.00005    *sin(dgtord(mp-m-2*f))
      +0.00004    *sin(dgtord(2*mp+2*f))
      -0.00004    *sin(dgtord(mp+m+2*f))
      +0.00004    *sin(dgtord(mp-2*m))
      +0.00003    *sin(dgtord(mp+m-2*f))
      +0.00003    *sin(dgtord(3*m))
      +0.00002    *sin(dgtord(2*mp-2*f))
      +0.00002    *sin(dgtord(mp-m+2*f))
      -0.00002    *sin(dgtord(3*mp+m));
   
   
   jd+=(corr+c_add);
   
   if(phase==FIRST_QUARTER)
      jd += w;
   
   if(phase==LAST_QUARTER)
      jd -= w;
   
   jdtogreg(&date,jd);
   
   return jd;
}




