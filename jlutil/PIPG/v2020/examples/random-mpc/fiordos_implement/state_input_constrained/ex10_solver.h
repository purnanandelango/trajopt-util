/**
 * This file is generated by FiOrdOs, a program licensed under GPL
 * by copyright holder Automatic Control Laboratory, ETH Zurich.
 * 
 * If you are interested in using this file commercially,
 * please contact the copyright holder.
 */

#ifndef EX10_SOLVER_H
#define EX10_SOLVER_H

#ifdef USE_REALTYPE_SINGLE
  typedef float realtype;
#else
  typedef double realtype;
#endif

/* <<< struct Params >>> */
#if defined(USE_PARAMSMODE_PTR)  && !(defined(MATLAB_MEX_FILE) || defined(RT))
#  define PARAMS_MACRO(name,len) const realtype *name;
#else /* the default */
#  define PARAMS_MACRO(name,len) realtype name[len];
#endif
typedef struct {
    PARAMS_MACRO(g,120)
    PARAMS_MACRO(c,1)
    PARAMS_MACRO(be,80)
} ex10_Params;
#undef PARAMS_MACRO

/* <<< struct Settings >>> */
typedef struct {
    struct {
        int warmstartInner;
    } approach;
    struct {
        realtype init[120];
        int maxit;
        realtype stopgEps;
        int stopgStride;
    } algoInner;
    struct {
        realtype init[80];
        int maxit;
        realtype stopgEps;
        int stopgStride;
    } algoOuter;
} ex10_Settings;

/* <<< struct Result >>> */
typedef struct {
    realtype la[80];
    realtype x[120];
    realtype d;
    int iter;
    int exitflag;
} ex10_Result;

/* <<< struct Work >>> */
typedef struct {
    struct {
        realtype H[120];
        realtype Ae[9600];
        realtype inner_g[120];
        realtype inner_c[1];
        struct {
            realtype l[4];
            realtype u[4];
        } X1;
        struct {
            realtype l[4];
            realtype u[4];
        } X2;
        struct {
            realtype l[4];
            realtype u[4];
        } X3;
        struct {
            realtype l[4];
            realtype u[4];
        } X4;
        struct {
            realtype l[4];
            realtype u[4];
        } X5;
        struct {
            realtype l[4];
            realtype u[4];
        } X6;
        struct {
            realtype l[4];
            realtype u[4];
        } X7;
        struct {
            realtype l[4];
            realtype u[4];
        } X8;
        struct {
            realtype l[4];
            realtype u[4];
        } X9;
        struct {
            realtype l[4];
            realtype u[4];
        } X10;
        struct {
            realtype l[8];
            realtype u[8];
        } X11;
        struct {
            realtype l[8];
            realtype u[8];
        } X12;
        struct {
            realtype l[8];
            realtype u[8];
        } X13;
        struct {
            realtype l[8];
            realtype u[8];
        } X14;
        struct {
            realtype l[8];
            realtype u[8];
        } X15;
        struct {
            realtype l[8];
            realtype u[8];
        } X16;
        struct {
            realtype l[8];
            realtype u[8];
        } X17;
        struct {
            realtype l[8];
            realtype u[8];
        } X18;
        struct {
            realtype l[8];
            realtype u[8];
        } X19;
        struct {
            realtype l[8];
            realtype u[8];
        } X20;
    } Prob;
    struct {
        realtype init[120];
        realtype glob_Linv;
        realtype glob_mu;
        realtype res_z[120];
        int res_stopcode;
        int res_iter;
    } algoInner;
    struct {
        realtype init[80];
        realtype glob_Linv;
        realtype res_z[80];
        int res_stopcode;
        int res_iter;
    } algoOuter;
} ex10_Work;

/* <<< solver INTERFACE >>> */
void ex10_init(ex10_Params *params, ex10_Settings *settings, ex10_Result *result, ex10_Work *work);
void ex10_solve(ex10_Params *params, ex10_Settings *settings, ex10_Result *result, ex10_Work *work);
#endif
