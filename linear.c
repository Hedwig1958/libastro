//
//  linear.c
//  tstastro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "linear.h"

//#define  __SIMPLEX_DEBUG_MODE__
//#define  __LEASTSQR_DEBUG_MODE__
#define  __QUADRATIC_SEARCH_DEBUG_MODE__


static int cmp_vertex(const void * obj1, const void * obj2);


//__ Resolution systeme d'equation lineaire par la methode de Gauss.
int resolve(double * a, double * b, int n)
{
   
   int    i,j,k,max,m=n+1;
   double t;
   
   //__ Phase de triangulation.
   for(i=0;i<n;i++)
     {
      max=i;
      
      //__ Optimisation.
      for(j=i+1;j<n;j++)
         if(fabs(a[j*m+i])>fabs(a[max*m+i]))
            max=j;
      
      for(k=i;k<=n;k++)
        {
         //__ Permutation.
         t=a[i*m+k];
         a[i*m+k]=a[max*m+k];
         a[max*m+k]=t;
        }
      
      //__ Triangulation.
      for(j=i+1;j<n;j++)
         for(k=n;k>=i;k--)
            if(a[i*m+i]!=0.0)
               a[j*m+k]-=a[i*m+k]*a[j*m+i]/a[i*m+i];
            else
               return EXIT_FAILURE;
     }
   
   //__ Phase de substitution.
   memset(b,0,sizeof(double)*n);
   for(i=n-1;i>=0;i--)
     {
      t=0.0;
      for(k=i+1;k<n;k++)
         t+=a[i*m+k]*b[k];
      
      if(a[i*m+i]!=0.0)
         b[i]=(a[i*m+n]-t)/a[i*m+i];
      else
         return EXIT_FAILURE;
     }
   
   return EXIT_SUCCESS;
   
}




/*--- Methode des moindres carres ---*/
int leastsqr(double * coef, int m, double * expe, int n, double * res)
{
   
   /*
    ---------------------------------------------------------
    coef : tableau des coefficients du polynome.
    {c1,c2,c3,c4, etc ...}
    m    : taille du tableau des coefficients du polynome.
    
    expe : tableau des points experimentaux.
    {x1,y1,x2,y2,x3,y3,x4,y4,etc ...}
    n    : taille du tableau des points experimentaux.
    
    res  : tableau des resulats
    {r1,r2,r3,r4,etc ...}
    ---------------------------------------------------------
    */
   
   int     i,j,k,ret;
   double  t, * mtx, * vect;
   
   //__ Tableau d'accueil pour les vecteurs.
   if(!(vect=malloc(sizeof(double)*(m+1)*n)))
      return EXIT_FAILURE;
   
   //__ Tableau d'accueil pour la matrice.
   if(!(mtx=malloc(sizeof(double)*m*(m+1))))
      return EXIT_FAILURE;
   
   //__ Construction des vecteurs de points et d'observations.
   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
        {
         vect[i*n+j]=coef[i]*pow(expe[j*2],i);
#ifdef __LEASTSQR_DEBUG_MODE__
         printf("%9.5f%c",vect[i*n+j],j==n-1?'\n':' ');
#endif
        }
   
   //__ Traitement des 'y'.
   for(j=0;j<n;j++)
     {
      vect[i*n+j]=expe[j*2+1];
#ifdef __LEASTSQR_DEBUG_MODE__
      printf("%9.5f%c",vect[i*n+j],j==n-1?'\n':' ');
#endif
     }
   
   //__Mise en matrice.
   for(i=0;i<m;i++)
      for(j=0;j<=m;j++)
        {
         t=0.0;
         for(k=0;k<n;k++)
            t+=vect[i*n+k]*vect[j*n+k];
         mtx[i*(m+1)+j]=t;
        }
   
#ifdef __LEASTSQR_DEBUG_MODE__
   for(i=0;i<m;i++)
      for(j=0;j<=m;j++)
         printf("%15.5f%c",mtx[i*(m+1)+j],j==m?'\n':' ');
#endif
   
   
   //__ Resolution de la matrice par la methode de Gauss.
   ret=resolve(mtx,res,m);
   
   //__ Produit des resultats avec les valeurs des coefficients.
   for(i=0;i<m;i++)
      res[i]*=coef[i];
   
   free(vect);
   free(mtx);
   
   return ret;
   
}


/*
 * Comparaison de double.
 */

static int cmp_vertex(const void * obj1, const void * obj2)
{
   if((*(double*)obj1)>(*(double*)obj2))
      return 1;
   if((*(double*)obj1)==(*(double*)obj2))
      return 0;
   return -1;
}


/*
 * Controle de depassement des limites.
 */

fbool check_limits(struct SIMPLEX_T * s, int i)
{
   fbool  low, high;
   
   //__ Check low limit.
   if((low=s->mat[LABEL_S+1+i]) < s->lim[i].low && s->lim[i].low_enable)
     {
      printf("WARNING, low limit (%f) reached for parameter %d with value (%f).\n",s->lim[i].low,i,s->mat[LABEL_S+1+i]);
      return low;
     }
   
   //__ Check high limit.
   if((high=s->mat[LABEL_S+1+i]) > s->lim[i].high && s->lim[i].high_enable)
     {
      printf("WARNING, high limit (%f) reached for parameter %d with value (%f).\n",s->lim[i].high,i,s->mat[LABEL_S+1+i]);
      return high;
     }
   
   return false;
}



/*
 * Controle de depassement des limites.
 */

int set_lowlimits(struct SIMPLEX_T * s, int i, double v)
{
   s->lim[i].low=v;
   s->lim[i].low_enable=true;
   return EXIT_SUCCESS;
}


/*
 * Controle de depassement des limites.
 */

int reset_lowlimits(struct SIMPLEX_T * s, int i)
{
   s->lim[i].low_enable=false;
   return EXIT_SUCCESS;
}

/*
 * Controle de depassement des limites.
 */

int set_highlimits(struct SIMPLEX_T * s, int i, double v)
{
   s->lim[i].high=v;
   s->lim[i].high_enable=true;
   return EXIT_SUCCESS;
}


/*
 * Controle de depassement des limites.
 */

int reset_highlimits(struct SIMPLEX_T * s, int i)
{
   s->lim[i].high_enable=false;
   return EXIT_SUCCESS;
}



/*
 * Initialisation d'un simplex.
 */

int init_simplex(simplex_t * s, int n)
{
   
   int i;
   
   //__ Init. constantes of simplex.
   s->n=n;
   s->size=sizeof(double)*BASE*(n+5);
   s->iter=0;
   s->max_iter=100;
   s->tolerance=STD_TOLERANCE;
   s->k_refl=1.0;
   s->k_exp=2.0;
   s->k_cont=0.5;
   
   s->lim=malloc(sizeof(limits_t)*n);
   if(s->lim)
     {
      memset(s->lim,0,sizeof(limits_t)*n);
      for(i=0;i<n;i++)
         s->lim[i].fnc=check_limits;
     }
   else
      return EXIT_FAILURE;
   
   
   //__ Init. labels.
   s->lbl_s=0;
   s->lbl_h=BASE;
   s->lbl_g=BASE*n;
   s->lbl_x=BASE*(n+1);
   s->lbl_r=BASE*(n+2);
   s->lbl_ec=BASE*(n+3);
   s->lbl_dx=BASE*(n+4);
   
   
   //__ Init. Fixed & Expandable Variables.
   s->mat=malloc(s->size);
   if(s->mat)
     {
      memset(s->mat,0,s->size);
      return EXIT_SUCCESS;
     }
   else
      return EXIT_FAILURE;
   
}


/*
 * Positionne 'tolerance' du simplex.
 */

int set_tolerance(simplex_t * s, double tol)
{
   s->tolerance=tol;
   return EXIT_SUCCESS;
}


/*
 * Positionne 'k_refl' du simplex.
 */

int set_reflection(simplex_t * s, double k)
{
   s->k_refl=k;
   return EXIT_SUCCESS;
}


/*
 * Positionne 'k_exp' du simplex.
 */

int set_expansion(simplex_t * s, double k)
{
   s->k_exp=k;
   return EXIT_SUCCESS;
}


/*
 * Positionne 'k_cont' du simplex.
 */

int set_contraction(simplex_t * s, double k)
{
   s->k_cont=k;
   return EXIT_SUCCESS;
}


/*
 * Liberation of zones allocate with 'init'.
 */

int free_simplex(simplex_t * s)
{
   free(s->mat);
   free(s->lim);
   return EXIT_SUCCESS;
}



/*
 * Methode Simplex de Nelder & Mead.
 */

int simplex(simplex_t * s, double * v, double * r, double (*fnc)(double*))
{
   
   int       i,j,k,n=s->n;
   double  * mat=s->mat;
   double    module, max_module;
   
   /*
    --------------------------------------------------------
    Arguments
    
    s      : pointeur de structure de type simplex_t.
    v      : tableau de sommets (vertices).
    r      : tableau de resultat (copie de 'S').
    fnc    : fonction de cout.
    
    
    n=4     | f | x(0)  ->  x(n-1)
    --------+---+------------------
    S       | . | v/r v/r v/r v/r
    H(0)    | . |  v   v   v   v
    .       | . |  v   v   v   v
    H(n-2)  | . |  v   v   v   v
    G       | . |  v   v   v   v
    X       | . |  .   .   .   .
    R       | . |  .   .   .   .
    EC      | . |  .   .   .   .
    dx      | . |  .   .   .   .
    
    --------------------------------------------------------
    */
   
   //__ Chargement des sommets.
   for(i=0;i<HEIGHT;i++)
      memcpy(&mat[1+i*BASE],&v[i*n],sizeof(double)*n);
   
   //__ Initialisation de la zone de resultats.
   memset(r,0,sizeof(double)*n);
   
   //__ Calcul du cout pour chaque sommets.
   for(i=0;i<HEIGHT;i++)
      mat[i*BASE]=fnc(&mat[1+i*BASE]);
   
   do
     {
      s->iter++;
      
      //__ Tri des sommets.
      qsort(mat,HEIGHT,sizeof(double)*BASE,cmp_vertex);
      
      //__ Controles des limites.
      for(i=0;i<n;i++)
         if(s->lim[i].fnc(s,i))
            return EXIT_FAILURE;
      
#ifdef  __SIMPLEX_DEBUG_MODE__
      
      //__ Impression du 'header'.
      printf("               F");
      for(i=0;i<n;i++)
         printf("              X%d",i);
      printf("\n");
      
      //__ Impression de 'S'.
      printf("S ");
      for(i=0;i<BASE;i++)
         printf(" %15.6f",mat[LABEL_S+i]);
      printf("\n");
      
      //__ Impression de 'H'.
      for(j=0;j<n-1;j++)
        {
         printf("H%d",j);
         for(i=0;i<BASE;i++)
            printf(" %15.6f",mat[LABEL_H+(BASE*j)+i]);
         printf("\n");
        }
      
      //__ Impression de 'G'.
      printf("G ");
      for(i=0;i<BASE;i++)
         printf(" %15.6f",mat[LABEL_G+i]);
      printf("\n");
      
      
#endif
      
      
      //__ Calcul du milieu de segment de contenant pas 'G'.
      memset(&mat[LABEL_X],0,sizeof(double)*BASE);
      
      //__ Somme des points sans 'G'.
      for(j=1;j<BASE;j++)
         for(i=0;i<HEIGHT-1;i++)
            mat[LABEL_X+j]+=mat[j+BASE*i];
      
      //__ Moyenne des points sans 'G'.
      for(i=1;i<BASE;i++)
         mat[LABEL_X+i]/=(HEIGHT-1);
      
#ifdef  __SIMPLEX_DEBUG_MODE__
      printf("X ");
      for(i=0;i<BASE;i++)
         printf(" %15.6f",mat[LABEL_X+i]);
      printf("\n");
#endif
      
      //__ Evaluation de 'R'.
      for(i=1;i<BASE;i++)
         mat[LABEL_R+i]=mat[LABEL_X+i]+s->k_refl*(mat[LABEL_X+i]-mat[LABEL_G+i]);
      
      //__ Calcul le cout de 'R'.
      mat[LABEL_R]=fnc(&mat[LABEL_R+1]);
      
      
#ifdef  __SIMPLEX_DEBUG_MODE__
      printf("R ");
      for(i=0;i<BASE;i++)
         printf(" %15.6f",mat[LABEL_R+i]);
      printf("\n");
#endif
      
      //__ Copie 'R' vers 'G'.
      if(mat[LABEL_R]<mat[LABEL_G])
         memcpy(&mat[LABEL_G],&mat[LABEL_R],sizeof(double)*BASE);
      
      
      //__ Expansion de 'R'.
      if(mat[LABEL_R]<mat[LABEL_S])
        {
         
         //__ Calcul de 'E'.
         for(i=1;i<BASE;i++)
            mat[LABEL_EC+i]=mat[LABEL_X+i]+s->k_exp*(mat[LABEL_R+i]-mat[LABEL_X+i]);
         
         //__ Calcul le cout de 'E'.
         mat[LABEL_EC]=fnc(&mat[LABEL_EC+1]);
         
#ifdef  __SIMPLEX_DEBUG_MODE__
         printf("E ");
         for(i=0;i<BASE;i++)
            printf(" %15.6f",mat[LABEL_EC+i]);
         printf("\n");
#endif
         
         //__ Copie 'EC' vers 'G'.
         if(mat[LABEL_EC]<mat[LABEL_R])
            memcpy(&mat[LABEL_G],&mat[LABEL_EC],sizeof(double)*BASE);
         
        }
      
      //__ Compression de 'R'.
      if(mat[LABEL_R]>=mat[LABEL_H])
        {
         
         //__ Calcul de 'C'.
         for(i=1;i<BASE;i++)
            mat[LABEL_EC+i]=mat[LABEL_X+i]+s->k_cont*(mat[LABEL_G+i]-mat[LABEL_X+i]);
         
         //__ Calcul le cout de 'C'.
         mat[LABEL_EC]=fnc(&mat[LABEL_EC+1]);
         
#ifdef  __SIMPLEX_DEBUG_MODE__
         printf("C ");
         for(i=0;i<BASE;i++)
            printf(" %15.6f",mat[LABEL_EC+i]);
         printf("\n");
#endif
         
         //__ Copy de 'EC' vers 'G'.
         if(mat[LABEL_EC]<mat[LABEL_H])
            memcpy(&mat[LABEL_G],&mat[LABEL_EC],sizeof(double)*BASE);
         else
           {
            
            //__ Calcul des nouveaux sommets.
            for(j=1;j<BASE;j++)
               for(i=0;i<HEIGHT;i++)
                  mat[j]=(mat[j+BASE*i]+mat[LABEL_S+j])/2;
            
#ifdef  __SIMPLEX_DEBUG_MODE__
            printf("RECALCUL DE NOUVEAUX SOMMETS\n");
#endif
            
           }
        }
      
      //__ Conditions de sortie.
      max_module=0.0;
      for(k=0;k<HEIGHT;k++)
         for(j=1;j<BASE;j++)
           {
            
            //__ Difference entre deux points contigus.
            mat[LABEL_DX+j]=mat[j+(BASE*k)]-mat[j+(BASE*((k+1)%HEIGHT))];
            
            //__ Calcul le module.
            mat[LABEL_DX]=0.0;
            for(i=1;i<BASE;i++)
               mat[LABEL_DX]+=(mat[LABEL_DX+i]*mat[LABEL_DX+i]);
            
            module=sqrt(mat[LABEL_DX]);
            if(module>max_module)
               max_module=module;
            
#ifdef  __SIMPLEX_DEBUG_MODE__
            printf("MAX MODULE: %15.6f %d %d\n",max_module,j+(BASE*k),j+(BASE*((k+1)%HEIGHT)));
#endif
            
           }
      
     } while(max_module>s->tolerance && k<s->max_iter);
   
   
   //__ Transfert du resultat.
   memcpy(r,&mat[1],sizeof(double)*n);
   
   
#ifdef  __SIMPLEX_DEBUG_MODE__
   printf("ITERATIONS: %d \n",s->iter);
#endif
   
   return EXIT_SUCCESS;
   
}



/*
 * Single variable minimization 'Quadratic Search'.
 */

double quadratic_search(double x0, double x2, double epsilon, double (*fnc)(double), fbool verbose)
{
   
   double x,x1,h,f0,f1,f2;
   
   /*
    ------------------------------------------
    x0, x2  : Bracket for the minimum of f(x).
    epsilon : Required tolerance.
    fnc     : fonction single variable.
    ------------------------------------------
    */
   
   x1=(x0+x2)/2;
   h=x1-x0;
   
   do
     {
      f0=fnc(x0);
      f1=fnc(x1);
      f2=fnc(x2);
      
      if(f0<f1)
        {
         //__ Shift left.
         x2=x1; x1=x0; x0=x1-h;
        }
      
      if(f2<f1)
        {
         //__ Shift right.
         x0=x1; x1=x2; x2=x1+h;
        }
      
      x=x1+(h*(f0-f2)/(2*(f0-2*f1+f2)));
      
      if(fnc(x)<f1)
         x1=x;
      
      h=h/2;
      x0=x1-h;
      x2=x1+h;
      
      if (verbose)
         printf("QUADRATIC SEARCH x1: %f, h: %f, epsilon: %f\n", x1, h, epsilon);
      
     } while(h>epsilon);
   
   return x1;
}


/*
 * Cosinus of complex.
 */

complex_t zcos(complex_t a)
{
   complex_t z={0.0,0.0};
   return z;
}


/*
 * Sinus of complex.
 */

complex_t zsin(complex_t a)
{
   complex_t z;
   z.r=(sin(a.r)*(exp(a.i)+exp(-a.i))/2.0);
   z.i=(cos(a.r)*(exp(-a.i)-exp(a.i))/-2.0);
   return z;
}




