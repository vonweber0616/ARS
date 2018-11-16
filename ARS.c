
#include "ARS.h"

/* object density  */
double log_density(double x)
{
    return (-x * x / 2);
}

/* compare two disjointed intervals  */
int compare(const void *key1, const void *key2)
{
    /* keys point to AvlNode->data, i.e. DListElmt  */
    Interval *node1 = (Interval *)dlist_data((DListElmt *)key1);
    Interval *node2 = (Interval *)dlist_data((DListElmt *)key2);

    if ((node1->x2) <= (node2->x1))
    {
        return -1;
    }
    else if ((node1->x1) >= (node2->x2))
    {
        return 1;
    }
    return 0;
}

void destroy(void *data)
{
    free(data);
    return;
}

/* ARS r.v. generator  */
double *rars(size_t n,                        /* number  */
             double (*log_density)(double x), /* desity  */
             const double *S0,                /* initial seperators, started with suppl, ended with suppu  */
             const size_t n0                  /* length of Sn[]  */
)
{
    int i;
    //int n0 = 0;
    double *pf, x, alpha, beta, area;
    double alpha1, alpha2, beta1, beta2;
    Interval *pint;
    //DListElmt *pdlElmt;

    extern double suppl, suppu, varpi, raw_varpi[7];
    extern unsigned int SEED;

    /* length of S0[]  */
    suppl = S0[0];
    suppu = S0[n0 - 1];
    //for(pf = S0; *pf != INF; pf++) {}
    //n0 = (pf - S0)/sizeof(double) + 1;

    /* BisTreeElmt->AvlNode->DlistElmt->Interval  */

    /* the set of Intervals as a DList  */
    DList *Int_DList;
    Int_DList = (DList *)malloc(sizeof(DList));
    dlist_init(Int_DList, destroy);

    /* the set of Intervals as a BisTree  */
    BisTree *Int_Tree;
    Int_Tree = (BisTree *)malloc(sizeof(BisTree));
    bistree_init(Int_Tree, compare, NULL);

    /*******************************/
    /* initialize the Intervals */
    /*******************************/

    for (i = 0; i < n0 - 1; i++)
    {
        pint = (Interval *)malloc(sizeof(Interval));

        if (i == 0)
        {
            /* Interval [suppl, x0]  */
            Int_init(pint, S0[0], S0[1], 0, NINF, LEFT1);

            /* left segment  */
            alpha = ALPHA(S0[1], S0[2]);
            beta = BETA(S0[1], alpha);
            area = AREA(S0[0], S0[1], alpha, beta);
            seg_init(pint->leftseg, S0[0], S0[1], alpha, beta, area);

            /* new weights to update varpi  */
            raw_varpi[0] = area;
            varpi_endat(1);
        }
        else if (i == 1)
        {
            /* Interval [x0, x1]  */
            alpha = ALPHA(S0[1], S0[2]);
            beta = BETA(S0[1], alpha);
            Int_init(pint, S0[1], S0[2], alpha, beta, LEFT2);

            /* left segment  */
            alpha = ALPHA(S0[2], S0[3]);
            beta = BETA(S0[2], alpha);
            area = AREA(S0[1], S0[2], alpha, beta);
            seg_init(pint->leftseg, S0[1], S0[2], alpha, beta, area);

            /* new weights to update varpi  */
            raw_varpi[0] = area;
            varpi_endat(1);
        }
        else if (i > 1 && i < n0 - 3)
        {
            /* middle Intervals  */

            /* lower envelope  */
            alpha = ALPHA(S0[i], S0[i + 1]);
            beta = BETA(S0[i], alpha);
            Int_init(pint, S0[i], S0[i + 1], alpha, beta, MIDDLE);

            alpha1 = ALPHA(S0[i - 1], S0[i]);
            beta1 = BETA(S0[i], alpha1);
            alpha2 = ALPHA(S0[i + 1], S0[i + 2]);
            beta2 = BETA(S0[i + 1], alpha2);
            x = (beta2 - beta1) / (alpha1 - alpha2); /* intersection  */

            /* left segment  */
            area = AREA(S0[i], x, alpha1, beta1);
            seg_init(pint->leftseg, S0[i], x, alpha1, beta1, area);

            /* new weights to update varpi  */
            raw_varpi[0] = area;

            /* right segment  */
            area = AREA(x, S0[i + 1], alpha2, beta2);
            seg_init(pint->rightseg, x, S0[i + 1], alpha2, beta2, area);

            /* new weights to update varpi  */
            raw_varpi[1] = area;
            varpi_endat(2);
        }
        else if (i == n0 - 3)
        {
            /* Inteval [x(n-3), x(n-2)]  */
            alpha = ALPHA(S0[n0 - 3], S0[n0 - 2]);
            beta = BETA(S0[n0 - 2], alpha);
            Int_init(pint, S0[n0 - 3], S0[n0 - 2], alpha, beta, RIGHT2);

            /* left segment  */
            alpha = ALPHA(S0[n0 - 4], S0[n0 - 3]);
            beta = BETA(S0[n0 - 3], alpha);
            area = AREA(S0[n0 - 3], S0[n0 - 2], alpha, beta);
            seg_init(pint->leftseg, S0[n0 - 3], S0[n0 - 2], alpha, beta, area);

            /* new weights to update varpi  */
            raw_varpi[0] = area;
            varpi_endat(1);
        }
        else if (i == n0 - 2)
        {
            /* Interval [x(n-2), suppu]  */
            Int_init(pint, S0[n0 - 2], S0[n0 - 1], 0, NINF, RIGHT1);

            /* left segment  */
            alpha = ALPHA(S0[n0 - 3], S0[n0 - 2]);
            beta = BETA(S0[n0 - 2], alpha);
            area = AREA(S0[n0 - 2], S0[n0 - 1], alpha, beta);
            seg_init(pint->leftseg, S0[n0 - 2], S0[n0 - 1], alpha, beta, area);

            /* new weights to update varpi  */
            raw_varpi[0] = area;
            varpi_endat(1);
        }

        /* insert Int_DList  */
        dlist_ins_next(Int_DList, Int_DList->tail, pint);

        /* insert Int_Tree  */
        bistree_insert(Int_Tree, Int_DList->tail);

        /* update varpi  */
        varpi_update();
    }

    /***********************************/
    /* main loop  **********************/
    /***********************************/

    /* count of the object random numbers  */
    size_t rv_count = 0;

    /* output array  */
    double *rv_out = malloc(sizeof(double) * n);

    /* X ~ gn(x)  */
    double X;

    int status;
    BiTreeNode *treeNode;
    pint = (Interval *)malloc(sizeof(Interval));
    DListElmt *pdlElmt = (DListElmt *)malloc(sizeof(DListElmt));
    DListElmt *bak = pdlElmt;
    pdlElmt->data = (void *)pint;

    /* set random seed  */
    SEED = (unsigned int)time(NULL);
    srand(123);

    /* uniform r.v.  */
    double U;

    /* the hit interval  */
    DListElmt **Int_hit;

    /* current value of _fn(x) (upper) and f_n(x) (lower)  */
    double _fn, f_n;

    while (rv_count < n)
    {
        /* generate X ~ gn(x)  */
        X = rpe(Int_DList);
        U = uniform();

        /* lookup X (a fake DListElmt) in the tree  */
        pdlElmt = bak; /*reset pdlElmt  */
        ((Interval *)pdlElmt->data)->x1 = X;
        ((Interval *)pdlElmt->data)->x2 = X;
        Int_hit = &pdlElmt;

        /* lookup changed the value of *data i.e. pdlElmt  */
        if (bistree_lookup(Int_Tree, (void **)Int_hit) == -1)
            continue;

        /* evaluate lower_fn(x) and upper_fn(x)  */
        _fn = upper_fn(X, *Int_hit);
        f_n = lower_fn(X, *Int_hit);

        if (U * _fn <= f_n)
        {
            rv_out[++rv_count] = X;
            continue;
        }
        else if (U * _fn <= f(x))
        {
            rv_out[++rv_count] = X;

            /* condition to insert new X  */
            if (Int_Tree->size < MAX_Sn && get_Int_DListElmt(*Int_hit)->flag == MIDDLE)
                Sn_ins(X, Int_Tree, Int_DList, (DListElmt *)(*Int_hit));
            continue;
        }
    }

    /* temporarl memery  */
    free(pint);
    free(bak);

    /* free the dlist and tree  */
    dlist_destroy(Int_DList);
    bistree_destroy(Int_Tree);

    return rv_out;
}

double rpe(DList *Int)
{
    extern double varpi;
    extern unsigned int SEED;

    double x1, x2, alpha, beta;
    Segment *left, *right;

    double sum = 0;
    double u = uniform();
    for (DListElmt *pdlElmt = Int->head; pdlElmt != NULL; pdlElmt = pdlElmt->next)
    {
        left = &get_leftseg_DListElmt(pdlElmt);

#define choseSeg(seg)             \
    {                             \
        sum += (seg)->area;       \
        if (sum > u * varpi)      \
        {                         \
            x1 = (seg)->x1;       \
            x2 = (seg)->x2;       \
            alpha = (seg)->alpha; \
            break;                \
        }                         \
    }
        choseSeg(left);

        if (get_Int_DListElmt(pdlElmt)->flag == MIDDLE)
        {
            right = &get_rightseg_DListElmt(pdlElmt);
            choseSeg(right);
        }
    }

    /* update u  */
    u = uniform();
    return log((1 - u) * exp(alpha * x1) + u * exp(alpha * x2)) / alpha;
}

double upper_fn(double x, DListElmt *pdlElmt)
{
    Segment *left = &get_leftseg_DListElmt(pdlElmt);

    if (get_Int_DListElmt(pdlElmt)->flag == MIDDLE &&
        x > left->x2) /* right half  */
    {
        Segment *right = &get_rightseg_DListElmt(pdlElmt);
        return exp(right->alpha * x + right->beta);
    }
    return exp(left->alpha * x + left->beta);
}

void Sn_ins(double X, BisTree *Int_Tree, DList *Int_DList, DListElmt *DLhit)
{
    /* change one DListElmt and insert another one      */
    /* X was found in thre Interval (DListElmt *)DLhit  */

    extern double raw_varpi[7];
    extern double varpi;

    Interval *Int = (Interval *)DLhit->data;

    /* copy old Interval  */
    Interval *newInt = (Interval *)malloc(sizeof(Interval));
    memcpy(newInt, Int, sizeof(Interval));

    double alpha, beta, x0, area, *pf;
    double k1, k2, b1, b2;

    /* copy the left new Interval to Int               */
    /* corrected Interval keeps the order in the Tree  */

    raw_varpi[0] = -(Int->leftseg.area + Int->rightseg.area);

    /* Interval  */
    alpha = ALPHA(Int->x1, X);
    beta = BETA(X, alpha);
    Int_init(Int, Int->x1, X, alpha, beta, MIDDLE);

    /* left Segment  */

    k2 = ALPHA(X, newInt->x2);
    b2 = BETA(X, k2);
    x0 = INTERSECTION(Int->leftseg.alpha, Int->leftseg.beta, k2, b2);

    Int->leftseg.x2 = x0;
    area_update(Int->leftseg);

    /* right segment  */
    area = AREA(x0, X, k2, b2);
    seg_init(Int->rightseg, x0, X, k2, b2, area);

    raw_varpi[1] = Int->leftseg.area + Int->rightseg.area;

    /* make another new Interval  */
    /********************************/

    /* Interval  */
    Int_init(newInt, X, newInt->x2, k2, b2, MIDDLE);

    /* left segment  */
    x0 = INTERSECTION(alpha, beta, newInt->rightseg.alpha, newInt->rightseg.beta);
    area = AREA(X, x0, alpha, beta);
    seg_init(newInt->leftseg, X, x0, alpha, beta, area);

    /* right segment  */
    newInt->rightseg.x1 = x0;
    area_update(newInt->rightseg);

    raw_varpi[2] = newInt->leftseg.area + newInt->rightseg.area;
    varpi_endat(3);

    /* insert at the tail of Int_DList  */
    dlist_ins_next(Int_DList, Int_DList->tail, (void *)newInt);

    /* insert to the tree  */
    bistree_insert(Int_Tree, (void *)Int_DList->tail);

    /* update varpi  */
    varpi_update();
}

/* wrapped for noraml r.v.  */
double *rnorm(size_t n, double mu, double sigma)
{
    double S0[] = {NINF, -6, -5, 0, 5, 6, INF};
    double *rv = rars(n, log_density, S0, 7);

    for (int i = 0; i < n; i++)
    {
        rv[i] = rv[i] * sigma + mu;
    }

    return rv;
}

/* output to R  */
void rnorm_ars(size_t *n, double *mu, double *sigma, double *out)
{
    double *rv = rnorm(*n, *mu, *sigma);
    memcpy(out, rv, sizeof(double) * (*n));
    free(rv);
    return;
}

void main(void)
{
    double S0[] = {NINF, -6, -5, 0, 5, 6, INF};
    int n = 1000000;
    struct timeval starttime, endtime;
    long s, ms, us, timeuse;
    gettimeofday(&starttime, 0);

    //double *rv = rars(n, log_density, S0, 7);
    double *rv = rnorm(n, 0, 1);
    /*
    for (int i = 0; i < n; i+=100)
        printf("%6.4f ", rv[i]);  */

    /* elapsed time  */
    gettimeofday(&endtime, 0);
    timeuse = 1000000 * (endtime.tv_sec - starttime.tv_sec) +
              endtime.tv_usec - starttime.tv_usec;
    s = timeuse / 1000000;
    ms = timeuse / 1000 % 1000;
    us = timeuse % 1000;
    printf("%lds %ldms %ldus\n", s, ms, us);

}