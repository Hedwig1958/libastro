//
//  tstastro.c
//  astro
//
//  Created by Corto Maltese on 30/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "tstastro.h"

#define PRINT             1

#define RUN_MOON          1
#define RUN_SUN           2
#define RUN_MARS          3
#define RUN_VENUS         4
#define RUN_STARS         5
#define RUN_SEASONS       6
#define RUN_LUNAR_PHASES  7
#define RUN_NUTATION      8
#define RUN_SUNDIAL       9



int main(int argc, char * argv[])
{
   
   char               file[128];
   char               buff[20];
   char               buff1[20];
   char               buff2[20];
   double             jd;
   double             h0;
   planets_data_t     planet_1;
   planets_data_t     planet_2;
   planets_data_t     planet_3;
   lunar_data_t       satel_1;
   stars_data_t       data_s18;
   stars_data_t       data_s14;
   stars_data_t       data_s4;
   stars_data_t       data_s5;
   stars_data_t       data_s6;
   stars_data_t       data_s2;
   nut_perioterms_t   nut;
   vsop_perioterms_t  earth;
   vsop_perioterms_t  venus;
   vsop_perioterms_t  mars;
   elp_perioterms_t   moon;
   date_t             d1,d2;
   sph_circ_t         sph_moon;
   sph_circ_t         sph_earth;
   equ_circ_t         epos;
   horiz_circ_t       hpos;
   equ_circ_t         star, alpheratz, vega, sirius, capella, thuban, aldebaran, northpole;
   equ_circ_t         tpos[4];
   equ_circ_t         tpos1[4];
   geo_circ_t *       tmp;
   geo_circ_t         paris={"Paris, FR",{dgtord(-2-20.2/60),dgtord(48+50.0/60+11/3600),0.0}};
   geo_circ_t         stonehenge={"Stonehenge, GB",{dgtord(1.75),dgtord(51.15),0.0}};
   geo_circ_t         oslo={"Oslo, S",{dgtord(-10.0),dgtord(58.5),0.0}};
   geo_circ_t         le_caire={"Le Caire, Egypte",{dgtord(-31.2333),dgtord(30.0166),0.0}};
   geo_circ_t         boston={"Boston, US",{dgtord(71.0833),dgtord(42.3333),0.0}};
   geo_circ_t         washington={"Washington, US",{dgtord(77.06555),dgtord(38.9213),0.0}};
   geo_circ_t         uccle={"Uccle, BE",{dgtord(-4.35),dgtord(50.80),0.0}};
   geo_circ_t         oisquercq={"Oisquercq, BE",{dgtord(-4.21),dgtord(50.66),0.0}};
   geo_circ_t         bologne={"Bologne, Italie",{dgtord(-11.33),dgtord(44.50),0.0}};
   geo_circ_t         zurich={"Zurich, Suisse",{dgtord(-8.55),dgtord(47.38),0.0}};
   geo_circ_t         meridian={"Meridian",{0.0,dgtord(52),0.0}};
   geo_circ_t         greenwich_0={"Greenwich 0 lat",0,0.0,{0.0,dgtord(52),0.0}};
   geo_circ_t         greenwich={"Greenwich, GB",{0.0,dgtord(51+28/60+22/3600),0.0}};
   geo_circ_t         greenwich_47={"Greenwich 47 lat",{0.0,dgtord(47.38),0.0}};
   geo_circ_t         saopaulo={"Sao Paulo, Bresil",{dgtord(46.6236111111),dgtord(-23.651166666),0.0}};
   geo_circ_t         argeles={"Argeles-Sur-Mer, FR",{dgtord(-3.033333),dgtord(42.65),0.0}};
   geo_circ_t         nord52={"52° nord, BE",{dgtord(0.0),dgtord(52.0),0.0}};
   geo_circ_t         nord30={"30° nord, BE",{dgtord(0.0),dgtord(30.0),0.0}};
   ephem_data_t       ephem;
   star_rec_t         pole={0L,0.0,dgtord(90.0),0.0,0.0,0.0};
   helio_circ_t       sun;
   
   char *work_directory=getenv("WORK_DIRECTORY");
   
   
   /*--- Loading nutation & obliquity terms ---*/
   sprintf(file,"%s%s",work_directory?work_directory:"","nutation.dat");
   loadNutation(file,&nut);
   planet_1.nutation=&nut;
   planet_2.nutation=&nut;
   planet_3.nutation=&nut;
   satel_1.nutation=&nut;
   
   /*--- Test for the earth ---*/
   sprintf(file,"%s%s",work_directory?work_directory:"","earth.dat");
   loadVsop87(file,&earth);
   planet_1.earth=&earth;
   //jd=2451545.0;
   //sph_earth=vsop87(planet_1.earth,jd);
   
   
   /*--- Test for venus ---*/
   sprintf(file,"%s%s",work_directory?work_directory:"","venus.dat");
   loadVsop87(file,&venus);
   planet_2.earth=&earth;
   planet_2.planet=&venus;
   //jd=2451545.0; //___ 01/01/2000 ___
   //jd=2122820.0; //___ 19/12/1099 ___
   //jd=2268920.0; //___ 19/12/1499 ___
   //sph_earth=vsop87(planet_2.planet,jd);
   
   /*--- Test for mars ---*/
   sprintf(file,"%s%s",work_directory?work_directory:"","mars.dat");
   loadVsop87(file,&mars);
   planet_3.earth=&earth;
   planet_3.planet=&mars;
   //jd=2451545.0; //___ 01/01/2000 ___
   //jd=2122820.0; //___ 19/12/1099 ___
   //jd=2268920.0; //___ 19/12/1499 ___
   //sph_earth=vsop87(planet_3.planet,jd);
   
   /*--- Test for the moon ---*/
   loadElp2000(&moon,work_directory);
   satel_1.moon=&moon;
   //setdate(&d1,1828,10,5,0,0,0,0);
   //jd=gregtojd(&d1);
   //sph_moon=elp2000(satel_1.moon,jd);
   
   
   
   switch(RUN_SUN)
  {
     case RUN_SUN:
    {
     int i;
     //double dt=59.0;
     double dt=0.0;
     
     
     /*--- Ephemeris for the SUN ---*/
     printf("\n/*--- Ephemeris for the SUN ---*/\n\n");
     //setdate(&d1,2005,4,8,20,15,38,20.1);
     //setdate(&d1,2005,10,3,10,10,40,84.6);
     //setdate(&d1,2010,1,1,0,0,0,0);
     //setdate(&d1,2100,1,1,0,0,0,0);
     //setdate(&d1,1990,1,1,0,0,0,0);
     //setdate(&d1,2000,1,1,0,0,0,0);
     //setdate(&d1,2007,5,11,0,0,0,0);
     //setdate(&d1,2006,1,1,0,0,0,0);
     //setdate(&d1,1998,3,21,0,0,0,0);
     setdate(&d1,2001,5,10,0,0,0,0);
     //setdate(&d1,2015,2,4,0,0,0,0);
     jd=gregtojd(&d1);
     jd+=(dt/86400);
     
     tmp=&greenwich;
     
     epos=appGeoPosSun(&planet_1,jd);
     printf("\nDate: %s a %s\n",pdate(&d1),ptime(&d1));
     printf("SUN\nLongitude: %s\n",rdstrdms(epos.sph.l,buff,DG_0_2PI));
     printf("Latitude: %s\n",rdstrdms(epos.sph.b,buff,DG_PI_PI));
     printf("... Alpha    : %s\n",rdstrhms(epos.equ.alpha,buff));
     printf("... Delta    : %s\n",rdstrdms(epos.equ.delta,buff,DG_PI_PI));
     printf("... Radius   : %.9f UA\n",epos.sph.radius);
     printf("... Parallax : %s\n",rdstrdms(epos.parallax,buff,DG_PI_PI));
     printf("... Semi Diam: %f\n",epos.semidiam);
     hpos=equtohoriz(&nut,&epos.equ,&greenwich_47.pos,epos.jd,epos.type);
     printf("Azinuth: %s\n",rdstrdms(hpos.pos.azimuth,buff,DG_0_2PI));
     printf("Altitude: %s\n",rdstrdms(hpos.pos.altitude,buff,DG_PI_PI));
     printf("Temps Siderale apparent à Greenwich: %s\n",rdstrhms(hpos.theta0,buff));
     printf("Temps Siderale apparent a %s : %s\n",tmp->name,rdstrhms(hpos.theta,buff));
     printf("Angle Horaire: %s\n",rdstrhms(hpos.eta,buff));
     printf("Angle Parallactique: %s\n",rdstrdms(hpos.q,buff,DG_PI_PI));
     
     //setdate(&d1,2014,12,21,0,0,0,0);
     //setdate(&d1,2015,2,4,0,0,0,0);
     //setdate(&d1,2014,11,16,0,0,0,0);
     jd=gregtojd(&d1);
     jd+=(dt/86400);
     
     tpos[0]=appGeoPosSun(&planet_1,jd-1);
     tpos[1]=appGeoPosSun(&planet_1,jd);
     tpos[2]=appGeoPosSun(&planet_1,jd+1);
     tpos[3]=appGeoPosSun(&planet_1,jd+2);
     
     for(i=0;i<4;i++)
        printf("Alpha: %s, Delta: %s, Radius: %.9f \n"
               ,rdstrhms(tpos[i].equ.alpha,buff1)
               ,rdstrdms(tpos[i].equ.delta,buff2,DG_PI_PI)
               ,tpos[i].sph.radius);
     
     ephem=transit(&nut,tpos,tmp,0.0);
     jdtogreg(&d1,ephem.jd);
     printf("Temps Siderale apparent à Greenwich: %s\n",rdstrhms(ephem.theta0,buff));
     printf("Transit du SOLEIL a %s le %s a %s\n",tmp->name,pdate(&d1),ptime(&d1));
     h0=altitude0(tpos[0].type,0.0,REFRACTION_RADAU,0.0,0.0,0.0,0.0);
     ephem=rising(&nut,tpos,tmp,h0);
     jdtogreg(&d1,ephem.jd);
     printf("Lever du SOLEIL a %s le %s a %s, h0 : %s\n",tmp->name,pdate(&d1),ptime(&d1),rdstrdms(ephem.h0,buff,DG_PI_PI));
     printf("Azinuth: %s\n",rdstrdms(ephem.azimuth,buff,DG_PI_PI));
     ephem=setting(&nut,tpos,tmp,h0);
     jdtogreg(&d1,ephem.jd);
     printf("Coucher du SOLEIL a %s le %s a %s, h0 : %s\n",tmp->name,pdate(&d1),ptime(&d1),rdstrdms(ephem.h0,buff,DG_PI_PI));
     printf("Azinuth: %s\n",rdstrdms(ephem.azimuth,buff,DG_PI_PI));
     
     sun=sunEphemeris(&planet_1,jd);
     printf("SUN PHYS OBS EPHEMERIS, position    : %s\n",rdstrdms(sun.pos.radius,buff,DG_PI_PI));
     printf("SUN PHYS OBS EPHEMERIS, latitude    : %s\n",rdstrdms(sun.pos.latitude,buff,DG_PI_PI));
     printf("SUN PHYS OBS EPHEMERIS, longitude   : %s\n",rdstrdms(sun.pos.longitude,buff,DG_0_2PI));
     
     
     
    }
      break;
     case RUN_MOON:
    {
     int i;
     //double dt=65.0;
     double dt=0.0;
     //lunar_shift_t ls={0.5,-0.25};
     lunar_shift_t ls={0.0,0.0};
     
     printf("\n/*--- Ephemeris for the MOON ---*/\n\n");
     //setdate(&d1,2005,4,8,20,15,38,20.1);
     //setdate(&d1,2005,10,3,10,10,40,84.6);
     //setdate(&d1,2100,1,1,0,0,0,0);
     //setdate(&d1,2007,2,10,0,0,0,0);
     setdate(&d1,1993,1,7,0,0,0,0);
     //setdate(&d1,1990,1,1,0,0,0,0);
     //setdate(&d1,1900,1,1,0,0,0,0);
     //setdate(&d1,2005,5,21,0,0,0,0);
     //setdate(&d1,2047,10,17,0,0,0,0);
     //setdate(&d1,1828,10,5,0,0,0,0);
     //setdate(&d1,1992,4,12,0,0,0,0);
     jd=gregtojd(&d1);
     jd+=(dt/86400);
     
     
     epos=appGeoPosMoon(&satel_1,jd,&ls);
     printf("\nDates: %s a %s\n",pdate(&d1),ptime(&d1));
     
     printf("MOON\nLongitude: %s, %.15f\n",rdstrdms(epos.sph.l,buff,DG_0_2PI),epos.sph.l);
     printf("Latitude: %s, %.15f\n",rdstrdms(epos.sph.b,buff,DG_PI_PI),epos.sph.b);
     printf("... Alpha    : %s\n",rdstrhms(epos.equ.alpha,buff));
     printf("... Delta    : %s\n",rdstrdms(epos.equ.delta,buff,DG_PI_PI));
     printf("... Radius   : %.9f Km\n",epos.sph.radius);
     printf("... Parallax : %s\n",rdstrdms(epos.parallax,buff,DG_PI_PI));
     printf("... Semi Dia : %s\n",rdstrdms(epos.semidiam,buff,DG_PI_PI));
     printf("%f\n",epos.jd);
     
     hpos=equtohoriz(&nut,&epos.equ,&paris.pos,epos.jd,epos.type);
     printf("\nMOON\n\nAzinuth: %s\n",rdstrdms(hpos.pos.azimuth,buff,DG_0_2PI));
     printf("Altitude: %s\n",rdstrdms(hpos.pos.altitude,buff,DG_PI_PI));
     
     /*--- Ephemerides pour Moon ---*/
     if (true)
       {
        
        tpos[0]=appGeoPosMoon(&satel_1,jd,&ls);
        tpos[1]=appGeoPosMoon(&satel_1,jd+1,&ls);
        tpos[2]=appGeoPosMoon(&satel_1,jd+2,&ls);
        tpos[3]=appGeoPosMoon(&satel_1,jd+3,&ls);
        
        for(i=0;i<4;i++)
           printf("Alpha: %s, Delta: %s, Radius: %.9f \n"
                  ,rdstrhms(tpos[i].equ.alpha,buff1)
                  ,rdstrdms(tpos[i].equ.delta,buff2,DG_PI_PI)
                  ,tpos[i].sph.radius);
        
        tmp=&greenwich;
        ephem=transit(&nut,tpos,tmp,0.0);
        jdtogreg(&d1,ephem.jd);
        printf("Transit de la LUNE a %s le %s a %s\n",tmp->name,pdate(&d1),ptime(&d1));
        h0=altitude0(tpos[0].type,0.0,REFRACTION_RADAU,MOON_SEMI_DIAMETER,0.0,0.0,0.0);
        ephem=rising(&nut,tpos,tmp,h0);
        jdtogreg(&d1,ephem.jd);
        printf("Lever de la LUNE a %s le %s a %s\n",tmp->name,pdate(&d1),ptime(&d1));
        ephem=setting(&nut,tpos,tmp,h0);
        jdtogreg(&d1,ephem.jd);
        printf("Coucher de la LUNE a %s le %s a %s\n",tmp->name,pdate(&d1),ptime(&d1));
        
        tpos1[0]=appGeoPosSun(&planet_1,jd);
        printf("Illuminated Moon Disk : %.4f\n",moonIlluminatedFraction(tpos1[0],tpos[0]));
        printf("Phase angle Moon : %s\n",rdstrdms(moonPhaseAngle(tpos1[0],tpos[0]),buff,DG_0_2PI));
        printf("Position angle Moon limb : %s\n",rdstrdms(moonPositionAngle(tpos1[0],tpos[0]),buff,DG_0_2PI));
       }
     
    }
      break;
     case RUN_VENUS:
    {
     
     /*--- Teste for venus ---*/
     //setdate(&d1,1987,4,10,19,21,0,0);
     //setdate(&d1,2100,1,1,0,0,0,0);
     setdate(&d1,2010,3,20,0,0,0,0);
     //setdate(&d1,1997,4,10,0,0,0,0);
     //setdate(&d1,1992,12,20,0,0,0,0);
     jd=gregtojd(&d1);
     
     epos=appGeoPosPlanets(&planet_2,jd);
     printf("\nDate: %s a %s\n",pdate(&d1),ptime(&d1));
     printf("VENUS\n\nLongitude: %s\n",rdstrdms(epos.sph.l,buff,DG_0_2PI));
     printf("Latitude: %s\n",rdstrdms(epos.sph.b,buff,DG_PI_PI));
     printf("... Alpha    : %s\n",rdstrhms(epos.equ.alpha,buff));
     printf("... Delta    : %s\n",rdstrdms(epos.equ.delta,buff,DG_PI_PI));
     printf("... Radius   : %.9f UA\n",epos.sph.radius);
     printf("... Parallax : %f\n",rdtodg(epos.parallax)*3600);
     printf("... Semi Diam: %f\n",epos.semidiam);
     hpos=equtohoriz(&nut,&epos.equ,&washington.pos,epos.jd,epos.type);
     printf("Azinuth: %s\n",rdstrdms(hpos.pos.azimuth,buff,DG_0_2PI));
     printf("Altitude: %s\n",rdstrdms(hpos.pos.altitude,buff,DG_PI_PI));
     
     /*--- Ephemerides pour venus ---*/
     
     tpos[0]=appGeoPosPlanets(&planet_2,jd-1);
     tpos[1]=appGeoPosPlanets(&planet_2,jd);
     tpos[2]=appGeoPosPlanets(&planet_2,jd+1);
     tpos[3]=appGeoPosPlanets(&planet_2,jd+2);
     
     tmp=&greenwich_47;
     ephem=transit(&nut,tpos,tmp,0.0);
     jdtogreg(&d1,ephem.jd);
     printf("Transit de Venus a %s le %s a %s\n",tmp->name,pdate(&d1),ptime(&d1));

    }
      
      break;
     case RUN_MARS:
    {
     /*--- Test for mars ---*/
     //setdate(&d1,2004,6,8,8,43,4,96);
     //setdate(&d1,2100,1,1,0,0,0,0);
     //setdate(&d1,1900,1,1,0,0,0,0);
     //setdate(&d1,2005,6,10,0,0,0,0);
     setdate(&d1,2005,5,21,0,0,0,0);
     jd=gregtojd(&d1);
     
     epos=appGeoPosPlanets(&planet_3,jd);
     printf("\nDate: %s a %s\n",pdate(&d1),ptime(&d1));
     printf("MARS\n\nLongitude: %s\n",rdstrdms(epos.sph.l,buff,DG_0_2PI));
     printf("Latitude: %s\n",rdstrdms(epos.sph.b,buff,DG_PI_PI));
     printf("... Alpha    : %s\n",rdstrhms(epos.equ.alpha,buff));
     printf("... Delta    : %s\n",rdstrdms(epos.equ.delta,buff,DG_PI_PI));
     printf("... Radius   : %.9f UA\n",epos.sph.radius);
     printf("... Parallax : %f\n",rdtodg(epos.parallax)*3600);
     printf("... Semi Diam: %f\n",epos.semidiam);
     hpos=equtohoriz(&nut,&epos.equ,&uccle.pos,epos.jd,epos.type);
     printf("Azinuth: %s\n",rdstrdms(hpos.pos.azimuth,buff,DG_0_2PI));
     printf("Altitude: %s\n",rdstrdms(hpos.pos.altitude,buff,DG_PI_PI));
     
     /*--- Ephemerides pour mars ---*/
     tpos[0]=appGeoPosPlanets(&planet_3,jd-1);
     tpos[1]=appGeoPosPlanets(&planet_3,jd);
     tpos[2]=appGeoPosPlanets(&planet_3,jd+1);
     tpos[3]=appGeoPosPlanets(&planet_3,jd+2);
     
     tmp=&greenwich_47;
     ephem=transit(&nut,tpos,tmp,0.0);
     jdtogreg(&d1,ephem.jd);
     printf("Transit de Mars a %s le %s a %s\n",tmp->name,pdate(&d1),ptime(&d1));
    }
      
      break;
      
     case RUN_STARS:
    {
     /*
      --- DUMP STARS DATA ---
      alpha Ursae minoris, Polaris  : sky2.dat [133]
      beta Ursae minoris, Kochab    : sky14.dat [223]
      alpha Orionis, Betelgeuse     : sky5.dat [327]
      alpha Ceti, Menkar            : sky3.dat  [9]
      alpha 1 libra                 : sky14.dat [222]
      alpha Tauri, aldebaran        : sky4.dat [175]
      alpha Draconis, Thuban        : sky14.dat [14]
      alpha Lyrae, Vega             : sky18.dat [218]
      iota Cygni                    : sky19.dat [158]
      
      */
     
     sprintf(file,"%s%s",work_directory?work_directory:"","sky2.dat");
     loadStars(file,&data_s2);
     
     if(!PRINT)
       {
        dumpStarsData(&data_s2);
       }
     
     //setdate(&d1,2028,11,13,4,33,36,0);
     //setdate(&d1,2007,6,14,21,9,0,0);
     setdate(&d1,2001,10,10,1,20,0,0);
     //setdate(&d1,2030,1,1,0,0,0,0);
     //setdate(&d1,1950,1,1,0,0,0,0);
     //setdate(&d1,1994,7,11,19,56,0,0);
     //setdate(&d1,2006,7,3,11,0,0,0);
     jd=gregtojd(&d1);
     
     star=appGeoPosStars(&nut,&data_s2.s[133],jd);
     printf("Alpha: %s\n",rdstrhms(star.equ.alpha,buff));
     printf("Delta: %s\n",rdstrdms(star.equ.delta,buff,DG_PI_PI));
     
     tmp=&greenwich;
     
     hpos=equtohoriz(&nut,&star.equ,&greenwich.pos,star.jd,star.type);
     printf("Azinuth: %s\n",rdstrdms(hpos.pos.azimuth,buff,DG_0_2PI));
     printf("Altitude: %s\n",rdstrdms(hpos.pos.altitude,buff,DG_PI_PI));
     printf("Temps Siderale apparent a Greenwich : %s\n",rdstrhms(hpos.theta0,buff));
     printf("Temps Siderale apparent locale a %s: %s\n",tmp->name,rdstrhms(hpos.theta,buff));
     printf("Angle Horaire: %s\n",rdstrhms(hpos.eta,buff));
     
     if (!PRINT)
       {
        double date = -13000;
        int i;
        for(i=0;i<2700;i++)
          {
           setdate(&d1,date,1,1,0,0,0,0);
           jd=gregtojd(&d1);
           star=appGeoPosStars(&nut,&data_s2.s[133],jd);
           printf("%f\t%f\t%f\n",date,star.equ.alpha,star.equ.delta);
           date+=10;
          }
       }
    }
      break;
     case RUN_SEASONS:
      
    {
     
     date_t    d;
     double    year=2014.00;
     
     jdtogreg(&d,seasons(&planet_1,year,SPRING,true).jd);
     printf("SEASON SPRING: %s %s\n",pdate(&d),ptime(&d));
     
     jdtogreg(&d,seasons(&planet_1,year,SUMMER,true).jd);
     printf("SEASON SUMMER: %s %s\n",pdate(&d),ptime(&d));
     
     jdtogreg(&d,seasons(&planet_1,year,AUTUMN,true).jd);
     printf("SEASON AUTUMN: %s %s\n",pdate(&d),ptime(&d));
     
     jdtogreg(&d,seasons(&planet_1,year,WINTER,true).jd);
     printf("SEASON WINTER: %s %s\n\n",pdate(&d),ptime(&d));
     
     printf("----\n");
    }
      
      break;
      
      
     case RUN_LUNAR_PHASES:
      
    {
     
     date_t    d;
     int       accurate=(1?HIGH:LOW);
     double    lunaison=2012.60;
     
     jdtogreg(&d,lunarphases(&planet_1,&satel_1,NEW_MOON,lunaison,accurate));
     printf("NEW MOON      : %s %s\n\n",pdate(&d),ptime(&d));
     
     jdtogreg(&d,lunarphases(&planet_1,&satel_1,FIRST_QUARTER,lunaison,accurate));
     printf("FIRST QUARTER : %s %s\n\n",pdate(&d),ptime(&d));
     
     jdtogreg(&d,lunarphases(&planet_1,&satel_1,FULL_MOON,lunaison,accurate));
     printf("FULL MOON     : %s %s\n\n",pdate(&d),ptime(&d));
     
     jdtogreg(&d2,lunarphases(&planet_1,&satel_1,LAST_QUARTER,lunaison,accurate));
     printf("LAST QUARTER  : %s %s\n\n",pdate(&d),ptime(&d));
     printf("----\n");
     
    }
      break;
      
      
     case RUN_NUTATION:
    {
     /*___ Nutation ___*/
     double d_psi, d_epsilon, epsilon0;
     jd=2446895.5;
     nutation(&nut,jd,&d_psi,&d_epsilon,&epsilon0);
     printf("Nutation: %s, %.8f, %.8f", rdstrdms(epsilon0,buff,DG_PI_PI ),rdtosec(d_psi),rdtosec(d_epsilon));
     break;
    }
      
     case RUN_SUNDIAL:
    {
     
     /*___ Calculate solar declinaison for sundial ___*/
     //setdate(&d1,2028,11,13,4,33,36,0);
     //setdate(&d1,2007,6,14,21,9,0,0);
     setdate(&d1,2001,3,19,0,0,0,0);
     //setdate(&d1,2030,1,1,0,0,0,0);
     //setdate(&d1,1950,1,1,0,0,0,0);
     //setdate(&d1,1994,7,11,19,56,0,0);
     //setdate(&d1,2006,7,3,11,0,0,0);
     jd=gregtojd(&d1);
          
     planarSundial(&nut,&planet_1,&greenwich,dgtord(-15.0),dgtord(90.0),1.0,2001);
     
    }
     default:
      break;
  }
   
   return EXIT_SUCCESS;
   
}









