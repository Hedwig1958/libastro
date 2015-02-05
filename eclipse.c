//
//  eclipse.c
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


#include "eclipse.h"

/*
 * for eclipses.
 */
double _eclipses_ra(double jd)
{
   
   lunar_shift_t lf = {0.5, -0.25};
   
   double ls=appGeoPosSun(*_planets_data(),jd).equ.alpha;
   double lm=appGeoPosMoon(*_lunar_data(),jd,&lf).equ.alpha;
   double p=*_phases();
   
   return pow(sin(ls)-sin(lm+p*(M_PI*2)),2);
}



/*
 * Calcul des eclipses.
 */
eclipses_t  eclipses(double type, double year)
{
   
   double        jd,k,t,corr=0,phase;
   double        m,mp,f,f1,omega,e,a1;
   double        p,q;
   double        w;
   eclipses_t    ecl;
   
   //__ Meeus 47.3
   k=floor((year-2000)*12.3685);
   
   //__ Meeus 47.2
   k+=type;
   t=k/1236.85;
   
   //__ Calcul approximatif du jd pour les phases de la lune - Meeus 47.1
   jd=2451550.09765+29.530588853*k+t*t*(0.0001337-t*(0.000000150+t*0.00000000073));
   
   //__ Sun's mean anomaly - Jean Meeus 47.4
   m=2.5534+29.10535669*k-t*t*(0.0000218+0.00000011*t);
   
   //__ Moon's mean anomaly - Jean Meeus 47.5
   mp=201.5643+385.81693528*k+t*t*(0.0107438+0.00001239*t-0.000000058*t*t);
   
   //__ Moon's arument of latitude - Jean Meeus 47.6
   f=160.7108+390.67050274*k-t*t*(0.0016341+0.00000227*t-0.000000011*t*t);
   
   //__ Longitude of the ascending node of the lunar orbit - Jean Meeus 47.7
   omega=124.7746-1.56375580*k+t*t*(0.0020691+0.00000215*t);
   
   //__ Eccentricity of the Earth - Jean Meeus 45.6
   e=1-t*(0.002516+0.0000074*t);
   
   //__.
   f1=f-0.02665*sin(dgtord(omega));
   a1=299.77+0.107408*k-0.009173*t*t;
   
#ifdef VERBOSE_ECLIPSES
   printf("k: %f, jd: %f, m: %f, mp: %f, f: %f, omega: %f\n",k,jd,m,mp,f,omega);
   printf("f1: %f, a1: %f, e: %f\n",f1,a1,e);
#endif
   
   //__ Time of maximun of solar & lunar eclipses.
   if(type==SOLAR_ECLIPSE)
     {
      corr= -0.4075*sin(dgtord(mp))
      +0.1721*e  *sin(dgtord(m));
      phase=NEW_MOON;
     }
   
   if(type==LUNAR_ECLIPSE)
     {
      corr= -0.4065    *sin(dgtord(mp))
      +0.1727*e  *sin(dgtord(m));
      phase=FULL_MOON;
     }
   
   corr+=(0.0161    *sin(dgtord(2*mp))
          -0.0097    *sin(dgtord(2*f1))
          +0.0073*e  *sin(dgtord(mp-m))
          -0.0050*e  *sin(dgtord(mp+m))
          -0.0023*e*e*sin(dgtord(2*m))
          +0.0021    *sin(dgtord(mp-2*f1))
          +0.0012    *sin(dgtord(mp+2*f1))
          +0.0006*e  *sin(dgtord(2*mp+m))
          -0.0004    *sin(dgtord(3*mp))
          -0.0003*e  *sin(dgtord(m+2*f1))
          +0.0003    *sin(dgtord(a1))
          +0.0002*e  *sin(dgtord(m-2*f1))
          -0.0002*e  *sin(dgtord(2*mp-m))
          -0.0002    *sin(dgtord(omega)));
   
   ecl.jd=jd+corr;
   
#ifdef VERBOSE_ECLIPSES
   printf("jd (corrected): %f, correction: %f\n",ecl.jd,corr);
#endif
   
   //__ Calcul de p,q.
   p=0.2070*e*sin(dgtord(m))
   +0.0024*e*sin(dgtord(2*m))
   -0.0392*  sin(dgtord(mp))
   +0.0116*  sin(dgtord(2*mp))
   -0.0073*e*sin(dgtord(mp+m))
   +0.0067*e*sin(dgtord(mp-m))
   +0.0118*  sin(dgtord(2*f1));
   
   q=5.2207
   -0.0048*e*cos(dgtord(m))
   +0.0020*e*cos(dgtord(2*m))
   -0.3299*  cos(dgtord(mp))
   -0.0060*e*cos(dgtord(mp+m))
   +0.0041*e*cos(dgtord(mp-m));
   
   //__ Calcul de w, gamma, u.
   w=fabs(cos(dgtord(f1)));
   ecl.gamma=(p*cos(dgtord(f1))+q*sin(dgtord(f1)))*(1-0.0048*w);
   ecl.u=0.0059
   +0.0046*e*cos(dgtord(m))
   -0.0182*  cos(dgtord(mp))
   +0.0004*  cos(dgtord(2*mp))
   -0.0005*  cos(dgtord(mp+m));
   
#ifdef VERBOSE_ECLIPSES
   printf("p: %f, q: %f, gamma: %f, u: %f\n",p,q,ecl.gamma,ecl.u);
#endif
   
   //__ Initialisations.
   ecl.magnitude=0.0;
   
   //__ Determination du type de l'eclipse.
   if(fabs(ecl.gamma)>1.5433+ecl.u)
      ecl.type=NO_ECLIPSE;
   else
      if(fabs(ecl.gamma)<0.9972)
        {
         //__ Eclipse is 'CENTRAL'.
         if(ecl.u<0.0)
            ecl.type=CENTRAL_TOTAL;
         else
            //__ Annular.
            if(ecl.u>0.0047)
               ecl.type=CENTRAL_ANNULAR;
            else
              {
               //__ Annular or annular/total.
               double omega=0.00464*sqrt(1-ecl.gamma*ecl.gamma);
               if(omega>0.0)
                 {
                  if(ecl.u<omega)
                     ecl.type=CENTRAL_ANNULAR_TOTAL;
                  else
                     ecl.type=CENTRAL_ANNULAR;
                 }
              }
        }
      else
         if(fabs(ecl.gamma)>0.9972 && fabs(ecl.gamma+ecl.u)<1.5433)
           {
            //__ Eclipse is 'NON_CENTRAL TOTAL OR ANNULAR'.
            if(fabs(ecl.gamma)>0.9972 && fabs(ecl.gamma)<0.9972+ecl.u)
              {
               //__ Eclipse is 'NON_CENTRAL TOTAL OR ANNULAR'.
               if(ecl.u<0.0)
                  ecl.type=NON_CENTRAL_ANNULAR;
               else
                  ecl.type=NON_CENTRAL_TOTAL;
              }
            else
              {
               //__ Eclipse is 'PARTIAL'.
               ecl.type=PARTIAL;
               ecl.magnitude=(1.5433+ecl.u-fabs(ecl.gamma))/(0.5461+2*ecl.u);
              }
           }
   return ecl;
}



/*
 * Besselian elements.
 */
bessel_elem_t besselianElements(planets_data_t * d1, lunar_data_t * d2, lunar_shift_t * ls, double k0, double jd, double dt)
{
   
   int            i;
   double         t;
   double         a,b,d,g;
   double         a1,a2,a3,a4;
   double         x,y,z;
   double         m;
   double         u1,u2;
   double         f1,f2;
   double         exp_x[NUM_SAMPLES*2];
   double         exp_y[NUM_SAMPLES*2];
   double         exp_z[NUM_SAMPLES*2];
   double         exp_d[NUM_SAMPLES*2];
   double         exp_m[NUM_SAMPLES*2];
   double         exp_u1[NUM_SAMPLES*2];
   double         exp_u2[NUM_SAMPLES*2];
   double         exp_f1[NUM_SAMPLES*2];
   double         exp_f2[NUM_SAMPLES*2];
   double         coef[4] = {1.0, 1.0, 1.0, 1.0};
   
   
   equ_circ_t     sun, moon;
   bessel_elem_t  elem;
   
   //__ Initialisation de elem & t0
   memset(&elem,0, sizeof(bessel_elem_t));
   elem.t0=jd;
   
   //__ t-n ... t-1, t0, t+1 ... t+n
   for(i=0;i<NUM_SAMPLES;i++)
     {
      t=-BESSEL_SHIFT_TIME+(((double)i)*BESSEL_HOURS_IN_DAY);
      
      //__ Calcul des positions equatoriale apparente du soleil et de la lune.
      sun=appGeoPosSun(d1,jd+t+(dt/86400));
      moon=appGeoPosMoon(d2,jd+t+(dt/86400),ls);
      
      //__ Declinaison (radian) et ascension droite (radian) de l'axe soleil/lune.
      b=moon.sph.radius/UA/sun.sph.radius;
      a1=cos(sun.equ.delta)*cos(sun.equ.alpha)-b*cos(moon.equ.delta)*cos(moon.equ.alpha);
      a2=cos(sun.equ.delta)*sin(sun.equ.alpha)-b*cos(moon.equ.delta)*sin(moon.equ.alpha);
      a3=sin(sun.equ.delta)-b*sin(moon.equ.delta);
      a4=cos(sun.equ.delta)-b*cos(moon.equ.delta);
      
      a=atan2(a2,a1);
      d=atan2(a3,a4);
      g=a3/sin(d);
      
      exp_d[i*2]=(double)i;
      exp_d[i*2+1]=d;
      
      //__ Position du x,y,z du centre de la lune.
      x=moon.sph.radius/EARTH_RADIUS*(cos(moon.equ.delta)*sin(moon.equ.alpha-a));
      y=moon.sph.radius/EARTH_RADIUS*(sin(moon.equ.delta)*cos(d)-cos(moon.equ.delta)*sin(d)*cos(moon.equ.alpha-a));
      z=moon.sph.radius/EARTH_RADIUS*(sin(moon.equ.delta)*sin(d)+cos(moon.equ.delta)*cos(d)*cos(moon.equ.alpha-a));
      
      //__ Memorisation des valeurs de x,y,z.
      exp_x[i*2]=(double)i;
      exp_x[i*2+1]=x;
      exp_y[i*2]=(double)i;
      exp_y[i*2+1]=y;
      exp_z[i*2]=(double)i;
      exp_z[i*2+1]=z;
      
      //__ Angle Horaire & Sideral time - Les eclipses de soleil, pg: 111.
      m=rdreduce(appSideraltime(d1->nutation,jd+t)-a,ANGLE_PI);
      exp_m[i*2]=(double)i;
      exp_m[i*2+1]=m;
      
      //__ les demi-angles au sommet (radian).
      f1=asin((sin(sectord(S_0))+k0*sin(sectord(PI_0)))/(g*sun.sph.radius));
      f2=asin(-(sin(sectord(S_0))-k0*sin(sectord(PI_0)))/(g*sun.sph.radius));
      
      exp_f1[i*2]=(double)i;
      exp_f1[i*2+1]=f1;
      exp_f2[i*2]=(double)i;
      exp_f2[i*2+1]=f2;
      
      //__ Rayons des sections circulaires des cones.
      u1=z*tan(f1)+k0/cos(f1);
      u2=z*tan(f2)+k0/cos(f2);
      
      //__ Memorisation des valeurs de u1, u2.
      exp_u1[i*2]=(double)i;
      exp_u1[i*2+1]=u1;
      exp_u2[i*2]=(double)i;
      exp_u2[i*2+1]=u2;
     }
   
   //__ Developpements polynomiaux.
   leastsqr(coef,4,exp_x,NUM_SAMPLES,elem.x);
   leastsqr(coef,4,exp_y,NUM_SAMPLES,elem.y);
   leastsqr(coef,4,exp_z,NUM_SAMPLES,elem.z);
   leastsqr(coef,3,exp_d,NUM_SAMPLES,elem.d);
   leastsqr(coef,4,exp_m,NUM_SAMPLES,elem.m);
   
   //__ Set the fourth element of 'm'.
   elem.m[4]=dgtord(-0.00417807);
   leastsqr(coef,3,exp_u1,NUM_SAMPLES,elem.u1);
   leastsqr(coef,3,exp_u2,NUM_SAMPLES,elem.u2);
   elem.f1=exp_f1[BESSEL_MIDDLE_ELEM];
   elem.f2=exp_f2[BESSEL_MIDDLE_ELEM];
   
   
#ifdef VERBOSE_BESSEL
  {
   date_t date;
   jdtogreg(&date,elem.t0);
   printf("\n\n");
   printf("******************************\n");
   printf("****  ELEMENTS OF BESSEL  ****\n");
   printf("******************************\n");
   printf("t0: %0.7f, %s %s\n",elem.t0,pdate(&date),ptime(&date));
   printf("X: %15.8f %15.8f %15.8f %15.8f\n",elem.x[0],elem.x[1],elem.x[2],elem.x[3]);
   printf("Y: %15.8f %15.8f %15.8f %15.8f\n",elem.y[0],elem.y[1],elem.y[2],elem.y[3]);
   printf("Z: %15.8f %15.8f %15.8f %15.8f\n",elem.z[0],elem.z[1],elem.z[2],elem.z[3]);
   printf("D: %15.8f %15.8f %15.8f\n",elem.d[0],elem.d[1],elem.d[2]);
   printf("M: %15.8f %15.8f %15.8f %15.8f\n",elem.m[0]/RDTODG,elem.m[1]/RDTODG,elem.m[2]/RDTODG,elem.m[3]/RDTODG);
   printf("tan f1      : %16.8f\n",tan(elem.f1));
   printf("tan f2      : %16.8f\n",tan(elem.f2));
   printf("u1: %16.8f %16.8f %16.8f\n",elem.u1[0],elem.u1[1],elem.u1[2]);
   printf("u2: %16.8f %16.8f %16.8f\n\n",elem.u2[0],elem.u2[1],elem.u2[2]);
   printf("******************************\n");
  }
#endif
   
   return elem;
}



/*
 * Circunstances of eclipse for given times.
 */
int eclipseCircunsTime(bessel_elem_t * bess, loc_circonst_t * loc, double jd, double dt)
{
   
   double tmp,t;
   double kcnt;
   double x,xp,y,yp,y1,d,m,n;
   double u1,u1p,u2,u2p;
   double h,p,k;
   double omega,phi,sin_phi_1,lambda,eta;
   double a,moon_ratio,b,b1,b2,c;
   
   //__ Initialisation de la structure d'accueil.
   memset(loc,0,sizeof(loc_circonst_t));
   
   //__ Time 't' in hour from reference hour 't0'.
   loc->jd=jd;
   
   t=jd-bess->t0+BESSEL_SHIFT_TIME;
   t*=HOURS_IN_DAY;
   
   //__ Besselian elements 'x','y','d' and 'm'.
   x=bess->x[0]+t*(bess->x[1]+t*(bess->x[2]+t*bess->x[3]));
   y=bess->y[0]+t*(bess->y[1]+t*(bess->y[2]+t*bess->y[3]));
   d=bess->d[0]+t*(bess->d[1]+t*bess->d[2]);
   m=bess->m[0]+t*(bess->m[1]+t*(bess->m[2]+t*bess->m[3]))+bess->m[4]*dt;
   u1=bess->u1[0]+t*(bess->u1[1]+t*bess->u1[2]);
   u2=bess->u2[0]+t*(bess->u2[1]+t*bess->u2[2]);
   
   //__ Hourly variations.
   xp=bess->x[1]+t*(2*bess->x[2]+t*3*bess->x[3]);
   yp=bess->y[1]+t*(2*bess->y[2]+t*3*bess->y[3]);
   
   /*
    -------------------------------------------------
    Evaluation du coefficient de centralite.
    Si 'kcnt' n'existe pas (tmp <= 0.0), il n'existe
    pas d'eclipse centrale pour le temps donne.
    -------------------------------------------------
    */
   
   omega=1/sqrt(1-0.006694385*cos(d)*cos(d));
   y1=omega*y;
   kcnt=(tmp=1-x*x-y1*y1)>0.0 ? sqrt(tmp) : 0.0;
   
   if(kcnt>0.0)
     {
      loc->kcnt=kcnt;
      
      //__ Longitude, latitude and Sun altitude of the eclipse.
      b1=omega*sin(d);
      b2=0.99664719*omega*cos(d);
      p=rdtodg(bess->m[1])/57.2957795;
      b=yp-p*x*sin(d);
      c=xp+p*y*sin(d);
      a=c-p*kcnt*cos(d);
      n=sqrt(a*a+b*b);
      
      eta=atan2(x,(kcnt*b2-y1*b1));
      sin_phi_1=kcnt*b1+y1*b2;
      phi=atan(1.00336409*tan(asin(sin_phi_1)));
      lambda=m-eta+1.0027379*(M_PI*2)*dt/86400;
      h=sin(d)*sin(phi)+cos(d)*cos(phi)*cos(eta);
      
      //__ Calculate width of the path.
      k=kcnt*kcnt+pow((x*a+y*b),2)/(n*n);
      
      //__ Moon ratio.
      u1p=u1-kcnt*tan(bess->f1);
      u2p=u2-kcnt*tan(bess->f2);
      moon_ratio=(u1p+u2p)/(u1p-u2p);
      
      //__ MaJ de la structure de sortie.
      loc->h=asin(h);
      loc->longitude=lambda;
      loc->latitude=phi;
      loc->duration=-7200.0*u2p/n;
      loc->width=12756.0*fabs(u2p)/sqrt(k);
      loc->moon_ratio=moon_ratio;
      
     }
   
#ifdef VERBOSE_CIRCUNS_TIME
  {
   char buff[128];
   date_t date;
   jdtogreg(&date,loc->jd);
   printf("\n\n");
   printf("***********************************************************\n");
   printf("****  CIRCUNSTANCE OF CENTRAL ECLIPSE FOR GIVEN TIMES  ****\n");
   printf("***********************************************************\n");
   printf("JD: %0.8f, %s %s\n",loc->jd,pdate(&date),ptime(&date));
   printf("t     : %15.8f\n",t);
   printf("X     : %15.8f\n",x);
   printf("Y     : %15.8f\n",y);
   printf("D     : %15.8f\n",d/RDTODG);
   printf("M     : %15.8f\n",m/RDTODG);
   printf("u1    : %15.8f\n",u1);
   printf("u2    : %15.8f\n",u2);
   printf("Xp    : %15.8f\n",xp);
   printf("Yp    : %15.8f\n",yp);
   printf("OMEGA : %15.8f\n",omega);
   printf("B     : %15.8f\n",kcnt);
   printf("\n\n");
   if(kcnt>0.0)
     {
      printf("p     : %15.8f\n",p);
      printf("b     : %15.8f\n",b);
      printf("c     : %15.8f\n",c);
      printf("u1p   : %15.8f\n",u1p);
      printf("u2p   : %15.8f\n",u2p);
      printf("a     : %15.8f\n",a);
      printf("n     : %15.8f\n",n);
      printf("ETA   : %15.8f\n",eta/RDTODG);
      printf("PHI   : %15.8f\n",phi/RDTODG);
     }
  }
#endif
   
   return EXIT_SUCCESS;
}



/*
 * Extreme limits of the central line.
 */
int eclipseExtremePoints(bessel_elem_t * bess, loc_circonst_t * ext, double dt)
{
   
   int      i,j;
   double   vr,tmp;
   double   offset,t[2]=
  {
   0.0,0.0
  }
   ,corr=0.0;
   double   tau,omega,phi_1,eta;
   double   a,b,d=0.0,m=0.0,u,v,s;
   double   x=0.0,y=0.0,xp,yp;
   
   //__ Initialisation des variables de sortie.
   memset(ext,0,sizeof(loc_circonst_t)*2);
   
   for(i=0;i<2;i++)
     {
      
      //__ omega.
      omega=1/sqrt(1-0.006694385*cos(bess->d[0])*cos(bess->d[0]));
      
      //__ 'u','v','a','b' & 'vr'.
      u=bess->x[0];
      v=omega*bess->y[0];
      a=bess->x[1];
      b=omega*bess->y[1];
      vr=(tmp=a*a+b*b)>0.0 ? sqrt(tmp) : 0.0;
      
      for(j=0;j<2;j++)
        {
         //__ First approximate times.
         s=(a*v-u*b)/vr;
         //__ Calcul du maximum de l'eclipse.
         tau=-(u*a+v*b)/(vr*vr);
         offset=sqrt(1-s*s)/vr;
         corr=tau+(i?offset:-offset);
         t[i]+=corr;
         
         //__ Besselian elements 'x','y','d' & 'm'.
         x=bess->x[0]+t[i]*(bess->x[1]+t[i]*(bess->x[2]+t[i]*bess->x[3]));
         y=bess->y[0]+t[i]*(bess->y[1]+t[i]*(bess->y[2]+t[i]*bess->y[3]));
         d=bess->d[0]+t[i]*(bess->d[1]+t[i]*bess->d[2]);
         m=bess->m[0]+t[i]*(bess->m[1]+t[i]*(bess->m[2]+t[i]*bess->m[3]))+bess->m[4]*dt;
         
         //__ Hourly variations.
         xp=bess->x[1]+t[i]*(2*bess->x[2]+t[i]*3*bess->x[3]);
         yp=bess->y[1]+t[i]*(2*bess->y[2]+t[i]*3*bess->y[3]);
         
#ifdef VERBOSE_EXTREME_POINTS
          {
           char * time[]={"BEGIN","END"};
           printf("___  FOR %s TIME, ACCURATE: %d ___\n",time[i],j);
           printf("OMEGA      : %15.8f\n",omega);
           printf("u          : %15.8f\n",u);
           printf("v          : %15.8f\n",v);
           printf("a          : %15.8f\n",a);
           printf("b          : %15.8f\n",b);
           printf("Vr         : %15.8f\n",vr);
           printf("s          : %15.8f\n",s);
           printf("TAU        : %15.8f\n",tau);
           printf("OFFSET     : %15.8f\n",offset);
           printf("Correction : %15.8f\n",corr);
           printf("T%d         : %15.8f\n",i,t[i]);
           printf("X          : %15.8f\n",x);
           printf("Y          : %15.8f\n",y);
           printf("D          : %15.8f\n",d/RDTODG);
           printf("M          : %15.8f\n",m/RDTODG);
           printf("Xp         : %15.8f\n",xp);
           printf("Yp         : %15.8f\n\n",yp);
          }
#endif
         
         //__ Corrected values for 'u','v','a','b' & 'vr'.
         omega=sqrt(1/(1-0.006694385*cos(d)*cos(d)));
         u=x;
         v=omega*y;
         a=xp;
         b=omega*yp;
         vr=(tmp=a*a+b*b)>0.0 ? sqrt(tmp) : 0.0;
        }
      
      //__ Sun altitude at central eclipse.
      phi_1=0.99664719*omega*omega*y*cos(d);
      eta=atan2(x,-omega*omega*y*sin(d));
      
      ext[i].jd=bess->t0+(t[i]/HOURS_IN_DAY-BESSEL_SHIFT_TIME);
      ext[i].latitude=atan(1.00336409*tan(asin(phi_1)));
      ext[i].longitude=m-eta+1.0027379*(M_PI*2)*dt/86400;
      
#ifdef VERBOSE_EXTREME_POINTS
       {
        char * time[]={"BEGIN","END"};
        printf("___  RESULT FOR %s TIME ___\n",time[i]);
        printf("t          : %15.8f\n",t[i]);
        printf("Sq Omega   : %15.8f\n",omega*omega);
        printf("sin Phi 1  : %15.8f\n",phi_1);
        printf("Eta        : %15.8f\n",eta/RDTODG);
        printf("Latitude   : %15.8f\n",ext[i].latitude/RDTODG);
        printf("Longitude  : %15.8f\n\n",ext[i].longitude/RDTODG);
       }
#endif
      
     }
   return EXIT_SUCCESS;
}


/*
 * Maximum of the central line.
 */

int eclipseMaximumPoint(bessel_elem_t * bess, loc_circonst_t * max, double dt)
{
   
   int      i;
   double   vr;
   double   t=0.0, tau;
   double   omega,sin_phi_1,eta;
   double   a,b,d=0.0,m=0.0,u,v;
   double   x=0.0,y,xp,yp;
   double   kcnt,b1,b2;
   
   //__ Initialisation des variables de sortie.
   memset(max,0,sizeof(loc_circonst_t));
   
   //__ Omega.
   omega=1/sqrt(1-0.006694385*cos(bess->d[0])*cos(bess->d[0]));
   
   //__ 'u','v','a','b' & 'vr'.
   u=bess->x[0];
   v=omega*bess->y[0];
   a=bess->x[1];
   b=omega*bess->y[1];
   vr=sqrt(a*a+b*b);
   
   for(i=0;i<2;i++)
     {
      //__ First approximate times.
      tau=-(u*a+v*b)/(vr*vr);
      t+=tau;
      
      //__ Besselian elements 'x','y','d' & 'm'.
      x=bess->x[0]+t*(bess->x[1]+t*(bess->x[2]+t*bess->x[3]));
      y=bess->y[0]+t*(bess->y[1]+t*(bess->y[2]+t*bess->y[3]));
      d=bess->d[0]+t*(bess->d[1]+t*bess->d[2]);
      m=bess->m[0]+t*(bess->m[1]+t*(bess->m[2]+t*bess->m[3]))+bess->m[4]*dt;
      
      //__ Hourly variations.
      xp=bess->x[1]+t*(2*bess->x[2]+t*3*bess->x[3]);
      yp=bess->y[1]+t*(2*bess->y[2]+t*3*bess->y[3]);
      
#ifdef VERBOSE_MAXIMUM_POINT
      printf("___  ACCURATE: %d ___\n",i);
      printf("OMEGA      : %15.8f\n",omega);
      printf("u          : %15.8f\n",u);
      printf("v          : %15.8f\n",v);
      printf("a          : %15.8f\n",a);
      printf("b          : %15.8f\n",b);
      printf("Vr         : %15.8f\n",vr);
      printf("tau        : %15.8f\n",tau);
      printf("T          : %15.8f\n",t);
      printf("X          : %15.8f\n",x);
      printf("Y          : %15.8f\n",y);
      printf("D          : %15.8f\n",d/RDTODG);
      printf("M          : %15.8f\n",m/RDTODG);
      printf("Xp         : %15.8f\n",xp);
      printf("Yp         : %15.8f\n\n",yp);
#endif
      
      //__ Corrected values for 'u','v','a','b' & 'vr'.
      omega=sqrt(1/(1-0.006694385*cos(d)*cos(d)));
      u=x;
      v=omega*y;
      a=xp;
      b=omega*yp;
      vr=sqrt(a*a+b*b);
      
     }
   
   //__ Sun altitude at central eclipse.
   kcnt=sqrt(1-u*u-v*v);
   b1=omega*sin(d);
   b2=0.99664719*omega*cos(d);
   eta=atan2(x,(kcnt*b2-v*b1));
   sin_phi_1=kcnt*b1+v*b2;
   
   max->jd=bess->t0+(t/HOURS_IN_DAY-BESSEL_SHIFT_TIME);
   max->latitude=atan(1.00336409*tan(asin(sin_phi_1)));
   max->longitude=m-eta+1.0027379*(M_PI*2)*dt/86400;
   
#ifdef VERBOSE_MAXIMUM_POINT
   printf("___  FINAL RESULT ___\n");
   printf("t          : %15.8f\n",t);
   printf("dt         : %15.8f\n",dt);
   printf("Sq Omega   : %15.8f\n",omega*omega);
   printf("sin Phi 1  : %15.8f\n",sin_phi_1);
   printf("Eta        : %15.8f\n",eta/RDTODG);
   printf("Latitude   : %15.8f\n",max->latitude/RDTODG);
   printf("Longitude  : %15.8f\n\n",max->longitude/RDTODG);
#endif
   
   return EXIT_SUCCESS;
}


/*
 * Local circonstances for a eclipse at noon.
 */

int eclipseNoonPoint(bessel_elem_t * bess, loc_circonst_t * ext, double dt)
{
   
   double t=0.0;
   
   t=-(bess->x[0]+t*t*(bess->x[2]+t*bess->x[3]))/bess->x[1];
   t=-(bess->x[0]+t*t*(bess->x[2]+t*bess->x[3]))/bess->x[1];
   
   eclipseCircunsTime(bess,ext,bess->t0+t,dt);
   
   return EXIT_SUCCESS;
}


/*
 * Circunstances of eclipse for given longitudes.
 */

int eclipseCircunsLongitude(bessel_elem_t * bess, loc_circonst_t * loc, double lambda, double tol
                            , double i, double g, double dt)
{
   
   double e,t,x,xp,y,yp,d,m;
   double upsilon,k1,k2;
   double phi,xi,xi_p,eta,eta_p,zeta,tau=0.0;
   double u,v,a,b,n,k;
   double d_phi,w,q;
   double omega;
   double hour_angle;
   double kcnt,tmp;
   double y1,u1,u2,u1p,u2p;
   
   
   /*
    ----------------------------------------------------------
    Type of curve                     i              G
    ----------------------------------------------------------
    * central line                     0.0         irrelevant
    * northern limit of path of       +1.0           +1.0
    total or annular eclipse
    * southern limit of path of       -1.0           +1.0
    total or annular eclipse
    * northern limit of partial       +1.0            0.0
    eclipse
    * southern limit of partial       -1.0            0.0
    eclipse
    * equal magnitude (north)         +1.0         the given G
    * equal magnitude (south)         -1.0         the given G
    ----------------------------------------------------------
    */
   
   
   //__ Initialisation de la structure d'accueil.
   memset(loc,0,sizeof(loc_circonst_t));
   
   //__ Time 't' in hour from reference hour 't0'.
   t=0.0;
   phi=0.0;
   
   do
     {
      
      //__ Besselian elements 'x','y','d','m','u1','u2' - Elements Of Solar Eclipses, pg 11 .
      x=bess->x[0]+t*(bess->x[1]+t*(bess->x[2]+t*bess->x[3]));
      y=bess->y[0]+t*(bess->y[1]+t*(bess->y[2]+t*bess->y[3]));
      d=bess->d[0]+t*(bess->d[1]+t*bess->d[2]);
      m=bess->m[0]+t*(bess->m[1]+t*(bess->m[2]+t*bess->m[3]))+bess->m[4]*dt;
      u1=bess->u1[0]+t*(bess->u1[1]+t*bess->u1[2]);
      u2=bess->u2[0]+t*(bess->u2[1]+t*bess->u2[2]);
      
      //__ Hourly variations.
      xp=bess->x[1]+t*(2*bess->x[2]+t*3*bess->x[3]);
      yp=bess->y[1]+t*(2*bess->y[2]+t*3*bess->y[3]);
      
      //__ Hour angle.
      hour_angle=dgtord(m-lambda);
      
      //__ Rectangular geocentric coordinates of a place - Elements Of Solar Eclipses, pg 10.
      upsilon=atan(0.99664719*tan(phi));
      
      //__ k1 = rho * sin(phi'), k2 = rho * cos(phi') - Elements Of Solar Eclipses, pg 10.
      k1=0.99664719*sin(upsilon);
      k2=cos(upsilon);
      
      //__ Calculate 'xi', 'eta' & 'zeta' - coordonnees du centre de l'ombre.
      xi=k2*sin(hour_angle);
      eta=k1*cos(d)-k2*cos(hour_angle)*sin(d);
      zeta=k1*sin(d)+k2*cos(hour_angle)*cos(d);
      
      //__ Calculate 'xi_p' & 'eta_p' - variation des coordonnees du centre de l'ombre.
      xi_p=0.01745329*rdtodg(bess->m[1])*k2*cos(hour_angle);
      eta_p=0.01745329*(rdtodg(bess->m[1])*xi*sin(d)-zeta*rdtodg(bess->d[1]));
      
      //__.
      omega=1/sqrt(1-0.006694385*cos(d)*cos(d));
      y1=omega*y;
      kcnt=(tmp=1-x*x-y1*y1) > 0.0 ? sqrt(tmp) : 0.0;
      
      //__ Calculate 'u', 'v, 'a', 'b', 'n' - Elements Of Solar Eclipses, pg 18.
      u=x-xi;
      v=y-eta;
      a=xp-xi_p;
      b=yp-eta_p;
      n=sqrt(a*a+b*b);
      
      //__ Calculate width of the path - Elements Of Solar Eclipses, pg 12.
      k=kcnt*kcnt+pow((x*a+y*b),2)/(n*n);
      
      //__ Rayons des sections circulaires des cones de penombre et d'ombre.
      u1p=u1-kcnt*tan(bess->f1);
      u2p=u2-kcnt*tan(bess->f2);
      
      e=u1p-g*(u1p+u2p);
      
      //__ Calculate 'tau', correction in time - Elements Of Solar Eclipses, pg 18 .
      tau=-(u*a+v*b)/(n*n);
      
      //__ Calculate 'delta phi', correction in latitude - Elements Of Solar Eclipses, pg 18 .
      w=(v*a-u*b)/n;
      q=(b*sin(hour_angle)*k1+a*(cos(hour_angle)*sin(d)*k1+cos(d)*k2))/(57.29578*n);
      d_phi=(w+i*fabs(e))/q;
      
      //__ Time 't' and correction 'tau' - Elements Of Solar Eclipses, pg 18.
      t+=tau;
      phi+=dgtord(d_phi);
      phi=asin(sin(phi));
      
     }
   while(fabs(tau)>tol);
   
   
   //__ MaJ de la structure de sortie.
   loc->jd=bess->t0+t-dt/3600;
   loc->kcnt=kcnt;
   loc->h=asin(sin(d)*sin(phi)+cos(d)*cos(phi)*cos(eta));
   loc->longitude=dgtord(lambda);
   loc->latitude=phi;
   loc->duration=7200.0*u2p/n;
   loc->width=12756.0*fabs(u2p)/sqrt(k);
   loc->moon_ratio=(u1p+u2p)/(u1p-u2p);
   
#ifdef VERBOSE_CIRCUNS_LONGITUDE
   printf("\n\n");
   printf("*****************************************************************\n");
   printf("****  CIRCUNSTANCES OF CENTRAL ECLIPSE FOR GIVEN LONGITUDES  ****\n");
   printf("*****************************************************************\n");
   printf("T              : %15.8f\n",t);
   printf("X              : %15.8f\n",x);
   printf("Y              : %15.8f\n",y);
   printf("D              : %15.8f\n",rdtodg(d));
   printf("M              : %15.8f\n",rdtodg(m));
   printf("Xp             : %15.8f\n",xp);
   printf("Yp             : %15.8f\n",yp);
   printf("H              : %15.8f\n",rdtodg(H));
   printf("rho *sin(phi') : %15.8f\n",k1);
   printf("rho *cos(phi') : %15.8f\n",k2);
   printf("xi             : %15.8f\n",xi);
   printf("eta            : %15.8f\n",eta);
   printf("zeta           : %15.8f\n",zeta);
   printf("xi_p           : %15.8f\n",xi_p);
   printf("eta_p          : %15.8f\n",eta_p);
   printf("u              : %15.8f\n",u);
   printf("v              : %15.8f\n",v);
   printf("a              : %15.8f\n",a);
   printf("b              : %15.8f\n",b);
   printf("n              : %15.8f\n",n);
   printf("Tau            : %15.8f\n",tau);
   printf("w              : %15.8f\n",w);
   printf("q              : %15.8f\n",q);
   printf("d_phi          : %15.8f\n",d_phi);
   printf("u1             : %15.8f\n",u1);
   printf("u2             : %15.8f\n",u2);
   printf("OMEGA          : %15.8f\n",omega);
   printf("y1             : %15.8f\n",y1);
   printf("kcnt           : %15.8f\n",kcnt);
   printf("u1p            : %15.8f\n",u1p);
   printf("u2p            : %15.8f\n",u2p);
#endif
   
   return EXIT_SUCCESS;
}

