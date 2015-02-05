//
//  linear.h
//  tstastro
//
//  Created by Corto Maltese on 29/11/14.
//  Copyright (c) 2014 Corto Maltese. All rights reserved.
//





#ifndef __linear__
#define __linear__

#include "fbool.h"
#include "common.h"

//__ Definition de type '__COMPLEX__'.

#ifndef   __COMPLEX__
#define   __COMPLEX__

typedef struct COMPLEX_T {double r,i;} complex_t;

#endif /* defined(__COMPLEX__) */




#define STD_TOLERANCE    0.001


#define LABEL_S   (s->lbl_s)
#define LABEL_H   (s->lbl_h)
#define LABEL_G   (s->lbl_g)
#define LABEL_X   (s->lbl_x)
#define LABEL_R   (s->lbl_r)
#define LABEL_EC  (s->lbl_ec)
#define LABEL_DX  (s->lbl_dx)

#define BASE      (s->n+1)
#define HEIGHT    (s->n+1)

struct SIMPLEX_T;

//__ Structure de type 'LIMITS_T'.
typedef struct LIMITS_T
{

  fbool   (*fnc)(struct SIMPLEX_T *,int);
  fbool   low_enable,  high_enable;
  double  low, high;

} limits_t;


//__ Structure de type 'SIMPLEX_T'.
typedef struct SIMPLEX_T
{

  /*
  --------------------------------------------------------
  Parametres

  tolerance :
  k_ref  : coeficient de reflection,    0 < k_ref  < 1
  k_exp  : coefficient d'expansion ,        k_exp  > 1
  k_cont : coefficient de contraction,  0 < k_cont < 1
  --------------------------------------------------------
  */

  //__ Size of the Simplex.
  int          n;
  int          iter;
  int          max_iter;
  size_t       size;
  limits_t   * lim;


  //__ Constantes of the Simplex.
  double       tolerance;
  double       k_refl;
  double       k_exp;
  double       k_cont;

  //__ Labels.
  int          lbl_s;
  int          lbl_h;
  int          lbl_g;
  int          lbl_x;
  int          lbl_r;
  int          lbl_ec;
  int          lbl_dx;

  //__ Fixed & Expandable Variables.
  double     * mat;

} simplex_t;



//__ Least square.
int   resolve(double * a, double * b, int n);
int   leastsqr(double * coef, int m, double * expe, int n, double * res);

//__ Simplex method.
fbool check_limits(simplex_t * s, int n);
int   set_lowlimits(struct SIMPLEX_T * s, int i, double v);
int   reset_lowlimits(struct SIMPLEX_T * s, int i);
int   set_highlimits(struct SIMPLEX_T * s, int i, double v);
int   reset_highlimits(struct SIMPLEX_T * s, int i);
int   init_simplex(simplex_t * s, int n);
int   free_simplex(simplex_t * s);
int   set_tolerance(simplex_t * s, double tol);
int   set_reflection(simplex_t * s, double k);
int   set_expansion(simplex_t * s, double k);
int   set_contraction(simplex_t * s, double k);
int   simplex(simplex_t * s, double * v, double * r, double (*fnc)(double *));

//__ Quadratic Search.
double quadratic_search(double x0, double x2, double epsilon, double (*fnc)(double), fbool verbose);

//__ Complex.
complex_t zcos(complex_t a);
complex_t zsin(complex_t a);


#endif /* defined(__linear__) */











