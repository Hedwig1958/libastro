//
//  elp2000.c
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


#include "elp2000.h"



/*
 * Lunar data type.
 */

lunar_data_t ** _lunar_data(void)
{
   static lunar_data_t * d;
   return &d;
}


/*
 * Full loading of lunar therory ELP2000-82B.
 */

int loadElp2000(elp_perioterms_t * pt, char * work_directory)
{
   int    i, st[NUM_ELP_PERIOTERMS];
   char   file[256];
   
   /*
    __________________________________
    Chargement des fichiers de donnees
    ELP2000, le path est contruit avec
    de la variable d'environnement :
    WORK_DIRECTORY.
    __________________________________
    */
   
   //__ Initialisation of size table 'st[]'.
   for(i=0;i<3;i++)
      st[i] = sizeof(elp_terms_t1);
   for(i=3;i<9;i++)
      st[i] = sizeof(elp_terms_t2);
   for(i=9;i<21;i++)
      st[i] = sizeof(elp_terms_t3);
   for(i=21;i<36;i++)
      st[i] = sizeof(elp_terms_t2);
   
   //__ Loading all terms file of elp2000 theory.
   for(i=0;i<NUM_ELP_PERIOTERMS;i++)
     {
      sprintf(file,"%selp%d.dat",(work_directory?work_directory:""),i+1);
      pt->n[i] = fileElp2000(file,&pt->e[i])/st[i];
     }
   
   return EXIT_SUCCESS;
}





/*
 * Loading lunar solution ELP2000-82B.
 */

int fileElp2000(char * f, void ** t)
{
   FILE         * data;
   int            size;
   struct stat    s;
   
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
   if(!(*t = malloc(s.st_size)))
     {
      fprintf(__ERROR_STREAM__,"Memory size cannot be correctly allocated in function fileElp2000\n");
      exit(EXIT_FAILURE);
     }
   
   //__ Update structure 'perioterms'.
   size = (int)fread(*t,1,s.st_size,data);
   
   //__ Clossing 'data file'.
   fclose(data);
   
   return size;
}




/*
 * Lunar theory ELP2000-82B.
 *
 *  The functions in this file use the Lunar Solution ELP 2000-82B by
 *  Michelle Chapront-Touze and Jean Chapront. *
 */

sph_circ_t elp2000(elp_perioterms_t * pt, double jd)
{
   double w1_v[5],w1;
   double w2_v[5],w2;
   double w3_v[5],w3;
   double eart_v[5],eart;
   double peri_v[5],peri;
   double mer_v[2],mer;
   double ven_v[2],ven;
   double mar_v[2],mar;
   double jup_v[2],jup;
   double sat_v[2],sat;
   double ura_v[2],ura;
   double nep_v[2],nep;
   double tau1,tau2,tau3,tau4,tau5;
   double delnu,dele,delg,delnp,delep;
   double d,lp,l,f,zeta;
   double precess;
   double pa;
   
   double ath = 384747.9806743165;
   double a0 = 384747.9806448954;
   double am = 0.074801329518;
   double alpha = 0.002571881335;
   double dtasm = 2.0*alpha/(3.0*am);
   double r[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
   
   sph_circ_t ret;
   
   //__ Time in Julian centuries from epoch 2000.
   tau1 = (jd-EPOCH_2000)/JULIAN_CENTURIES;
   tau2 = tau1*tau1;
   tau3 = tau2*tau1;
   tau4 = tau3*tau1;
   tau5 = tau4*tau1;
   
   //__ Precession constant.
   precess=sectord(5029.0966);
   
   //__ Lunar arguments, w1.
   w1_v[0] = dgtord(218.0+18.0/C1+59.95571/C2);
   w1_v[1] = sectord(1732559343.73604);
   w1_v[2] = sectord(-5.8883);
   w1_v[3] = sectord(0.6604e-2);
   w1_v[4] = sectord(-0.3169e-4);
   
   w1=w1_v[0]+w1_v[1]*tau1+w1_v[2]*tau2+w1_v[3]*tau3+w1_v[4]*tau4;
   
   //__ Lunar arguments, w2.
   w2_v[0] = dgtord(83.0+21.0/C1+11.67475/C2);
   w2_v[1] = sectord(14643420.2632);
   w2_v[2] = sectord(-38.2776);
   w2_v[3] = sectord(-0.45047e-1);
   w2_v[4] = sectord(0.21301e-3);
   
   w2=w2_v[0]+w2_v[1]*tau1+w2_v[2]*tau2+w2_v[3]*tau3+w2_v[4]*tau4;
   
   //__ Lunar arguments, w3.
   w3_v[0] = dgtord(125.0+2.0/C1+40.39816/C2);
   w3_v[1] = sectord(-6967919.3622);
   w3_v[2] = sectord(6.3622);
   w3_v[3] = sectord(0.7625e-2);
   w3_v[4] = sectord(-0.3586e-4);
   
   w3=w3_v[0]+w3_v[1]*tau1+w3_v[2]*tau2+w3_v[3]*tau3+w3_v[4]*tau4;
   
   //__ Lunar arguments, peri.
   peri_v[0] = dgtord(102.0+56.0/C1+14.42753/C2);
   peri_v[1] = sectord(1161.2283);
   peri_v[2] = sectord(0.5327);
   peri_v[3] = sectord(-0.138e-3);
   peri_v[4] = 0.0;
   
   peri=peri_v[0]+peri_v[1]*tau1+peri_v[2]*tau2+peri_v[3]*tau3;
   
   //__ Planetary arguments for Mercury.
   mer_v[0] = dgtord(252.0+15.0/C1+3.25986/C2);
   mer_v[1] = sectord(538101628.68898);
   
   mer=mer_v[0]+mer_v[1]*tau1;
   
   //__ Planetary arguments for Venus.
   ven_v[0] = dgtord(181.0+58.0/C1+47.28305/C2);
   ven_v[1] = sectord(210664136.43355);
   
   ven=ven_v[0]+ven_v[1]*tau1;
   
   //__ Lunar arguments, eart.
   eart_v[0] = dgtord(100.0+27.0/C1+59.22059/C2);
   eart_v[1] = sectord(129597742.2758);
   eart_v[2] = sectord(-0.0202);
   eart_v[3] = sectord(0.9e-5);
   eart_v[4] = sectord(0.15e-6);
   
   eart=eart_v[0]+eart_v[1]*tau1+eart_v[2]*tau2+eart_v[3]*tau3+eart_v[4]*tau4;
   
   //__ Planetary arguments for Mars.
   mar_v[0] = dgtord(355.0+25.0/C1+59.78866/C2);
   mar_v[1] = sectord(68905077.59284);
   
   mar=mar_v[0]+mar_v[1]*tau1;
   
   //__ Planetary arguments for Jupiter.
   jup_v[0] = dgtord(34.0+21.0/C1+5.34212/C2);
   jup_v[1] = sectord(10925660.42861);
   
   jup=jup_v[0]+jup_v[1]*tau1;
   
   //__ Planetary arguments for Saturn.
   sat_v[0] = dgtord(50.0+4.0/C1+38.89694/C2);
   sat_v[1] = sectord(4399609.65932);
   
   sat=sat_v[0]+sat_v[1]*tau1;
   
   //__ Planetary arguments for Uranus.
   ura_v[0] = dgtord(314.0+3.0/C1+18.01841/C2);
   ura_v[1] = sectord(1542481.19393);
   
   ura=ura_v[0]+ura_v[1]*tau1;
   
   //__ Planetary arguments for Neptune.
   nep_v[0] = dgtord(304.0+20.0/C1+55.19575/C2);
   nep_v[1] = sectord(786550.32074);
   
   nep=nep_v[0]+nep_v[1]*tau1;
   
   //__ Corrections of the constants (fit to DE200/LE200).
   delnu = sectord(0.55604)/w1_v[1];
   dele = sectord(0.01789);
   delg = sectord(-0.08066);
   delnp = sectord(-0.06424)/w1_v[1];
   delep = sectord(-0.12879);
   
   //__ Delaunay's arguments.
   d = w1-eart+M_PI;
   lp = eart-peri;
   l = w1-w2;
   f = w1-w3;
   
   //__ zeta.
   zeta=w1_v[0]+(w1_v[1]+precess)*tau1;
   
   //__ Main problem.
   for(int i=0;i<3;i++)
     {
      for(int j=0;j<pt->n[i];j++)
        {
         elp_terms_t1 *t = &((elp_terms_t1 *)pt->e[i])[j];
         t->a = (i==2) ? t->a-2.0*t->a*delnu/3.0 : t->a;
         double tgv = t->b1+dtasm*t->b5;
         double x = t->a+tgv*(delnp-am*delnu)+t->b2*delg+t->b3*dele+t->b4*delep;
         
         if(i<2)
            r[i%3] += x*sin(t->i1*d+t->i2*lp+t->i3*l+t->i4*f);
         else
            r[i%3] += x*cos(t->i1*d+t->i2*lp+t->i3*l+t->i4*f);
        }
     }
   
   //__ Earth figure perturbations.
   for(int i=3;i<9;i++)
     {
      for(int j=0;j<pt->n[i];j++)
        {
         elp_terms_t2 *t = &((elp_terms_t2 *)pt->e[i])[j];
         double x = (i>=6) ? t->a*tau1 : t->a;
         r[i%3] += x*sin(t->i1*zeta+t->i2*d+t->i3*lp+t->i4*l+t->i5*f+dgtord(t->phi));
        }
     }
   
   //__ Planetary perturbations, table 1.
   for(int i=9;i<15;i++)
     {
      for(int j=0;j<pt->n[i];j++)
        {
         elp_terms_t3 *t = &((elp_terms_t3 *)pt->e[i])[j];
         double x = (i>=12) ? t->a*tau1 : t->a;
         r[i%3] += x*sin(t->i1*mer+t->i2*ven+t->i3*eart+t->i4*mar+t->i5*jup+t->i6*sat+t->i7*ura+t->i8*nep
                       +t->i9*d+t->i10*l+t->i11*f+dgtord(t->phi));
        }
     }
   
   //__ Planetary perturbations, table 2.
   for(int i=15;i<21;i++)
     {
      for(int j=0;j<pt->n[i];j++)
        {
         elp_terms_t3 *t = &((elp_terms_t3 *)pt->e[i])[j];
         double x = (i>=18) ? t->a*tau1 : t->a;
         r[i%3] += x*sin(t->i1*mer+t->i2*ven+t->i3*eart+t->i4*mar+t->i5*jup+t->i6*sat+t->i7*ura
                       +t->i8*d+t->i9*lp+t->i10*l+t->i11*f+dgtord(t->phi));
        }
     }
   
   //__ Tides, Moon figure perturbations, Relativity, Solar eccentricity.
   for(int i=21;i<36;i++)
     {
      for(int j=0;j<pt->n[i];j++)
        {
         elp_terms_t2 *t = &((elp_terms_t2 *)pt->e[i])[j];
         double x = (i>=24 && i<27) ? t->a*tau1 : (i>=33 && i<36) ? t->a*tau2 :t->a;
         r[i%3] += x*sin(t->i1*zeta+t->i2*d+t->i3*lp+t->i4*l+t->i5*f+dgtord(t->phi));
        }
     }
   
   //__ 'pa' Accumulated precession between J2000 and the date.
   pa=5029.0966*tau1+1.112022*tau2+0.00007732*tau3-0.0000235316*tau4;
   
   r[0] = sectord(r[0]+pa)+w1;
   r[1] = sectord(r[1]);
   r[2] = r[2]*a0/ath;
   
   ret.jd = jd;
   ret.type = MOON;
   ret.pos.l = r[0];
   ret.pos.b = r[1];
   ret.pos.radius = r[2];
   
   return ret;
}

