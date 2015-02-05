//
//  sundial.c
//  tstastro
//
//  Created by Corto Maltese on 30/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//



#include <stdio.h>
#include <math.h>

#include "sundial.h"


/*
 * Solar declinaison for a given hour angle.
 *
 * julian day : jd
 * hour angle : eta
 */

sun_data_t sunPositionFromHourAngle(nut_perioterms_t  * nut, planets_data_t * earth, geo_circ_t * geo, double  jd, double eta)
{
   date_t             d;
   equ_circ_t         tpos[4];
   equ_circ_t         sunpos;
   sun_data_t         sun;
   
   //__ Calculate 4 sun positions around the giving jd.
   tpos[0] = appGeoPosSun(earth,jd-1);
   tpos[1] = appGeoPosSun(earth,jd);
   tpos[2] = appGeoPosSun(earth,jd+1);
   tpos[3] = appGeoPosSun(earth,jd+2);
   
   //__ Calculate the event time for the giving eta.
   ephem_data_t e = transit(nut,tpos,geo,eta);
   jdtogreg(&d,e.jd);
   
   //__ Calculate de sun parameters for this time.
   sunpos = appGeoPosSun(earth,e.jd);
   
   sun.jd = jd;
   sun.eta = eta;
   sun.event_time = e.jd;
   sun.theta0 = e.theta0;
   sun.pos.alpha = sunpos.equ.alpha;
   sun.pos.delta = sunpos.equ.delta;
   
   return sun;
}

/*
 * Planar Sundial.
 *
 * gnomonic declinaison    : declinaison
 * zenithal distance       : zenithal
 * length of polar stylus  : stylus
 *
 */



void planarSundial(nut_perioterms_t * nut, planets_data_t * earth, geo_circ_t * geo, double declinaison, double zenithal, double length, double year)
{
   char            buff[20];
   double          hourAngle=dgtord(0.0);
   int             accurate = true;

   seasons_data_t  seasonsArray[4];
   date_t          d;
   sun_data_t      sun;
   stylus_param_t  sty;

   
   //__ Calculation for size & coordinate of the stylus.
   sty = stylus(geo->pos.latitude,declinaison,zenithal,length);
   printf("Stylus size & coordinate, length: %.8f, x0: %.8f, y0: %.8f, u: %.8f, psi: %s\n",sty.length,sty.x0,sty.y0,sty.u,rdstrdms(sty.psi,buff,DG_PI_PI ));

   //__ Calculation for the seasons dates, spring, summer, autumn and winter.
   seasonsArray[SPRING] = seasons(earth,year,SPRING,accurate);
   jdtogreg(&d,seasonsArray[SPRING].jd);
   printf("SPRING DATE : %s %s\n",pdate(&d),ptime(&d));
   printf("Delta : %s, longitude %s\n",rdstrdms(seasonsArray[SPRING].circ.equ.delta,buff,DG_PI_PI),rdstrdms(seasonsArray[SPRING].circ.sph.l,buff,DG_0_2PI));
   sun = sunPositionFromHourAngle(nut,earth,geo,seasonsArray[SPRING].jd,hourAngle);
   printf("Delta : %s for hour angle: %f\n",rdstrdms(sun.pos.delta,buff,DG_PI_PI),hourAngle);
   
   //__ Summer.
   seasonsArray[SUMMER] = seasons(earth,year,SUMMER,accurate);
   jdtogreg(&d,seasonsArray[SUMMER].jd);
   printf("SUMMER DATE : %s %s\n",pdate(&d),ptime(&d));
   printf("Delta : %s, longitude %s\n",rdstrdms(seasonsArray[SUMMER].circ.equ.delta,buff,DG_PI_PI),rdstrdms(seasonsArray[SUMMER].circ.sph.l,buff,DG_0_2PI));
   sun = sunPositionFromHourAngle(nut,earth,geo,seasonsArray[SUMMER].jd,hourAngle);
   printf("Delta : %s for hour angle: %f\n",rdstrdms(sun.pos.delta,buff,DG_PI_PI),hourAngle);
   
   //__ Calculation the coordinates of the shadow for each hour angle.
   shadow_coord_t sha = shadow(geo->pos.latitude,sun.pos.delta,declinaison,zenithal,dgtord(0.0),length);
   printf("Shadow coordinate, below plane: %d x: %.8f, y: %.8f\n", sha.illum, sha.x, sha.y);

   //__ Autumn.
   seasonsArray[AUTUMN] = seasons(earth,year,AUTUMN,accurate);
   jdtogreg(&d,seasonsArray[AUTUMN].jd);
   printf("AUTUMN DATE : %s %s\n",pdate(&d),ptime(&d));
   printf("Delta : %s, longitude %s\n",rdstrdms(seasonsArray[AUTUMN].circ.equ.delta,buff,DG_PI_PI),rdstrdms(seasonsArray[AUTUMN].circ.sph.l,buff,DG_0_2PI));
   sun = sunPositionFromHourAngle(nut,earth,geo,seasonsArray[AUTUMN].jd,hourAngle);
   printf("Delta : %s for hour angle: %f\n",rdstrdms(sun.pos.delta,buff,DG_PI_PI),hourAngle);
   
   //__ Winter.
   seasonsArray[WINTER] = seasons(earth,year,WINTER,accurate);
   jdtogreg(&d,seasonsArray[WINTER].jd);
   printf("WINTER DATE : %s %s\n",pdate(&d),ptime(&d));
   printf("Delta : %s, longitude %s\n",rdstrdms(seasonsArray[WINTER].circ.equ.delta,buff,DG_PI_PI),rdstrdms(seasonsArray[WINTER].circ.sph.l,buff,DG_0_2PI));
   sun = sunPositionFromHourAngle(nut,earth,geo,seasonsArray[WINTER].jd,hourAngle);
   printf("Delta : %s for hour angle: %f\n",rdstrdms(sun.pos.delta,buff,DG_PI_PI),hourAngle);
   
   
   //__ Calculation for the rise and set for each dates.
   sunRiseTransitSetting(nut,earth,geo,seasonsArray[SPRING].jd);
   sunRiseTransitSetting(nut,earth,geo,seasonsArray[SUMMER].jd);
   sunRiseTransitSetting(nut,earth,geo,seasonsArray[AUTUMN].jd);
   sunRiseTransitSetting(nut,earth,geo,seasonsArray[WINTER].jd);
   
}


/*
 * Sun Rise Transit Setting.
 */
void sunRiseTransitSetting(nut_perioterms_t * nut, planets_data_t * earth, geo_circ_t * geo, double jd)
{
   char            buff[20];
   equ_circ_t      sun[4];
   ephem_data_t    ephem;
   date_t          date;
   double          h0;
   
   sun[0] = appGeoPosSun(earth,jd);
   sun[1] = appGeoPosSun(earth,jd+1);
   sun[2] = appGeoPosSun(earth,jd+2);
   sun[3] = appGeoPosSun(earth,jd+3);
   
   //__ h0.
   h0 = altitude0(sun[0].type,0.0,REFRACTION_RADAU,0.0,0.0,0.0,0.0);
   printf("Ephemris for the sun at %s, h0 : %s\n",geo->name,rdstrdms(ephem.h0,buff,DG_PI_PI));
   
   //__ Sun transit.
   ephem = transit(nut,sun,geo,0.0);
   jdtogreg(&date,ephem.jd);
   printf("Transit : %s %s\n",pdate(&date),ptime(&date));
   printf("hour angle: %s\n",rdstrdms(ephem.eta,buff,DG_0_2PI));

   //__ Sun rising.
   ephem = rising(nut,sun,geo,h0);
   jdtogreg(&date,ephem.jd);
   printf("Rise    : %s %s, azimuth : %s\n",pdate(&date),ptime(&date),rdstrdms(ephem.azimuth,buff,DG_PI_PI));
   printf("hour angle: %s\n",rdstrdms(ephem.eta,buff,DG_0_2PI));
   
   //__ Sun setting.
   ephem = setting(nut,sun,geo,h0);
   jdtogreg(&date,ephem.jd);
   printf("Setting : %s %s, azimuth : %s\n",pdate(&date),ptime(&date),rdstrdms(ephem.azimuth,buff,DG_PI_PI));
   printf("hour angle: %s\n",rdstrdms(ephem.eta,buff,DG_0_2PI));
   
}


/*
 * Calculate shadow coordinate.
 *
 * geographical latitude  : latitude
 * sun declinaison        : delta
 * gnomonic declinaison   : declinaison
 * zenithal distance      : zenithal
 * hour angle             : eta
 * length of polar stylus : length
 *
 */
shadow_coord_t shadow(double latitude, double delta, double declinaison, double zenithal, double eta, double length)
{
   double p;
   double q;
   double nx;
   double ny;
   
   shadow_coord_t sha = {false,0.0,0.0};
   
   //__ calulation for p & q.
   p = sin(latitude)*cos(zenithal)-cos(latitude)*sin(zenithal)*cos(declinaison);
   q = sin(declinaison)*sin(zenithal)*sin(eta)+(cos(latitude)*cos(zenithal)+sin(latitude)*sin(zenithal)*cos(declinaison))*cos(eta)+p*tan(delta);
   
   sha.illum = q>=0.0;
   
   if (sha.illum)
     {
      
      //__ nx & ny.
      nx = cos(declinaison)*sin(eta)-sin(declinaison)*(sin(latitude)*cos(eta)-cos(latitude)*tan(delta));
      ny = cos(zenithal)*sin(declinaison)*sin(eta)-(cos(latitude)*sin(zenithal)-sin(latitude)*cos(zenithal)*cos(declinaison))*cos(eta)
                           -(sin(latitude)*sin(zenithal)+cos(latitude)*cos(zenithal)*cos(declinaison))*tan(delta);
      
      //__ shadow coordinate, x & y.
      sha.x = length*(nx/q);
      sha.y = length*(ny/q);
      
     }
   
   return sha;
}


/*
 * Calculate stylus size & coordinate.
 *
 * geographical latitude  : latitude
 * gnomonic declinaison   : declinaison
 * zenithal distance      : zenithal
 * length of polar stylus : length
 *
 */
stylus_param_t stylus(double latitude, double declinaison, double zenithal, double length)
{
   double p;
   double q;
   double nx;
   double ny;
   
   stylus_param_t sty;
   
   //__ calulation for p,q,nx,ny.
   p = sin(latitude)*cos(zenithal)-cos(latitude)*sin(zenithal)*cos(declinaison);
   q = (cos(latitude)*cos(zenithal)+sin(latitude)*sin(zenithal)*cos(declinaison));
   nx = -sin(declinaison)*(sin(latitude));
   ny = -(cos(latitude)*sin(zenithal)-sin(latitude)*cos(zenithal)*cos(declinaison))-(sin(latitude)*sin(zenithal));

   //__ stylus coordinates & size.
   sty.length = length;
   sty.x0 = (length/p)*cos(latitude)*sin(declinaison);
   sty.y0 = -(length/p)*(sin(latitude)*sin(zenithal)+cos(latitude)*cos(zenithal)*cos(declinaison));
   sty.u = length/fabs(p);
   sty.psi = asin(fabs(p));
   
   return sty;
}

