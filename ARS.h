#ifndef ARS_H
#define ARS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include "bistree.h"
#include "dlist.h"


/* the flag of Interval position  */
typedef enum Position_
{
    LEFT1,
    LEFT2,
    RIGHT1,
    RIGHT2,
    MIDDLE
} Position;

/* a sigment in upper envelope  */
typedef struct Segment_
{
    double x1, x2;
    double alpha, beta;
    double area;
} Segment;

/* an interval of Sn, including a segment in the lower envelope 
   and two segments in the upper envelope  */
typedef struct Interval_
{
    double x1, x2;
    double alpha, beta;
    Position flag;
    Segment leftseg, rightseg;
} Interval;

#define seg_init(seg, xa, xb, k, b, S) \
    (seg).x1 = xa;                     \
    (seg).x2 = xb;                     \
    (seg).alpha = k;                   \
    (seg).beta = b;                    \
    (seg).area = S

#define Int_init(Int, xa, xb, k, b, f) \
    (Int)->x1 = xa;                    \
    (Int)->x2 = xb;                    \
    (Int)->alpha = k;                  \
    (Int)->beta = b;                   \
    (Int)->flag = f

#define INF (1.0 / 0.0)
#define NINF (-1.0 / 0.0)
#define ALPHA(x1, x2) ((log_density(x2) - log_density(x1)) / (x2 - x1))
#define BETA(x, alpha) (log_density(x) - x * alpha)
#define AREA(x1, x2, alpha, beta) ((exp(alpha * x2 + beta) - exp(alpha * x1 + beta)) / alpha)

/*update the area of a segment  */
#define area_update(seg) \
    seg.area = (exp(seg.alpha * seg.x2 + seg.beta) - exp(seg.alpha * seg.x1 + seg.beta)) / seg.alpha

#define uniform() ((rand() + 0.0) / RAND_MAX)

/* evaluate lower_fn(x)  */
#define lower_fn(x, pdlElmt) exp((((Interval *)(pdlElmt)->data)->alpha * x + ((Interval *)(pdlElmt)->data)->beta))

/* evaluate upper_fn(x)  */
double upper_fn(double x, DListElmt *pdlElmt);

/* f(x)  */
#define f(x) (exp(log_density(x)))

/* intersection of two lines  */
#define INTERSECTION(k1, b1, k2, b2) ((b1 - b2) / (k2 - k1))

#define get_AvlNode_TreeNode(node) ((AvlNode *)(node->data))
#define get_DListElmet_AvlNode(avl) ((DListElmt *)((AvlNode *)avl)->data)

#define get_Int_DListElmt(elmt) (((Interval *)((DListElmt *)elmt)->data))
#define get_DListElmt_TreeNode(node) ((DListElmt *)((AvlNode *)node->data)->data)
#define get_Int_TreeNode(node) ((Interval *)((DListElmt *)((AvlNode *)node->data)->data)->data)
#define get_Int_AvlNode(avl) ((Interval *)((DListElmt *)((AvlNode *)avl)->data)->data)

#define get_leftseg_DListElmt(elmt) (((Interval *)((DListElmt *)elmt)->data)->leftseg)
#define get_rightseg_DListElmt(elmt) (((Interval *)((DListElmt *)elmt)->data)->rightseg)

/* global factor varpi  */
double varpi;
double suppl, suppu;
double raw_varpi[7];
unsigned int SEED;

/* max length of Sn[]  */
#define MAX_Sn 2000

/* update varpi  */
#define varpi_endat(i) raw_varpi[i] = INF
#define varpi_update()                    \
    for (pf = raw_varpi; *pf < INF; pf++) \
    varpi += *pf

/* piecewise exonential density r.v. generator  */
double rpe(DList *Int);

/* insert X to Sn[]  */
void Sn_ins(double X, BisTree *Int_Tree, DList *Int_DList, DListElmt *DLhit);

void destroy(void *data);
int compare(const void *key1, const void *key2);
double log_density(double x);

double *rars(size_t n, double (*log_density)(double x), const double *S0, const size_t n0);

/* wrapped for normal r.v.  */
double *rnorm(size_t n, double mu, double sigma);

/* interface to R  */
void rnorm_ars(size_t *n, double *mu, double *sigma, double *out);

#endif
