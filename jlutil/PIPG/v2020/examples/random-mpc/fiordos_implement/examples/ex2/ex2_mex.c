/**
 * This file is generated by FiOrdOs, a program licensed under GPL
 * by copyright holder Automatic Control Laboratory, ETH Zurich.
 * 
 * If you are interested in using this file commercially,
 * please contact the copyright holder.
 */

#include <string.h>
#include "ex2_solver.h"
#include "mex.h"

/* input and output arguments */
#define MPARAMS prhs[0]
#define MSETGS prhs[1]
#define MRES plhs[0]

/* copying */
void copyCArrayToM(realtype *src, double *dest, int dim) {
    while (dim--) {
        *dest++ = (double)*src++;
    }
}
void copyMArrayToC(double *src, realtype *dest, int dim) {
    while (dim--) {
        *dest++ = (realtype)*src++;
    }
}


/* THE mex-function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )  {

    /* define variables for structs */
    ex2_Work work ={0};
    ex2_Params params ={0};
    ex2_Settings settings = {0};
    ex2_Result result = {0};

    /* initialize structs */
    ex2_init(&params, &settings, &result, &work);

    /* Check for proper number of arguments */
    if (nrhs != 1 && nrhs!=2) {
        mexErrMsgTxt("One or two input arguments required.");
    }
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }

    /* set params according to mparams from prhs */
    {
    
    mxArray *tmp1;
    int i;
    const char *fname;
    
    if (!mxIsStruct(MPARAMS)) {
        mexErrMsgTxt("MPARAMS must be a structure.");
    }
    for (i=0; i<mxGetNumberOfFields(MPARAMS); ++i) {
        fname= mxGetFieldNameByNumber(MPARAMS, i);
        if (0==strcmp(fname,"g"));
        else
        if (0==strcmp(fname,"c"));
        else
        if (0==strcmp(fname,"be"));
        else
            mexPrintf("Warning: MPARAMS contains the superfluous field '%s'.\n",fname);
    }
    
        tmp1=mxGetField(MPARAMS,0,"g");
        if (NULL==tmp1){
            mexErrMsgTxt("MPARAMS.g is missing.");
        }
        if (!mxIsDouble(tmp1)) {
            mexErrMsgTxt("MPARAMS.g must be a double.");
        }
        if ( mxGetM(tmp1) != 30 || mxGetN(tmp1) != 1 ) {
            mexErrMsgTxt("MPARAMS.g must be matrix of dimension 30x1.");
        }
        copyMArrayToC(mxGetPr(tmp1),params.g,30 );
        
        tmp1=mxGetField(MPARAMS,0,"c");
        if (NULL==tmp1){
            mexErrMsgTxt("MPARAMS.c is missing.");
        }
        if (!mxIsDouble(tmp1)) {
            mexErrMsgTxt("MPARAMS.c must be a double.");
        }
        if ( mxGetM(tmp1) != 1 || mxGetN(tmp1) != 1 ) {
            mexErrMsgTxt("MPARAMS.c must be matrix of dimension 1x1.");
        }
        copyMArrayToC(mxGetPr(tmp1),params.c,1 );
        
        tmp1=mxGetField(MPARAMS,0,"be");
        if (NULL==tmp1){
            mexErrMsgTxt("MPARAMS.be is missing.");
        }
        if (!mxIsDouble(tmp1)) {
            mexErrMsgTxt("MPARAMS.be must be a double.");
        }
        if ( mxGetM(tmp1) != 20 || mxGetN(tmp1) != 1 ) {
            mexErrMsgTxt("MPARAMS.be must be matrix of dimension 20x1.");
        }
        copyMArrayToC(mxGetPr(tmp1),params.be,20 );
        
    
    }

    /* change settings according to msetgs in prhs (if given)*/
    if (nrhs==2) {
    {
    
    mxArray *tmp1;
    int i1;
    const char *fname1;
    mxArray *tmp2;
    int i2;
    const char *fname2;
    
    if (!mxIsStruct(MSETGS)) {
        mexErrMsgTxt("MSETGS must be a structure.");
    }
    for (i1=0; i1<mxGetNumberOfFields(MSETGS); ++i1) {
        fname1 = mxGetFieldNameByNumber(MSETGS, i1);
        tmp1 = mxGetField(MSETGS,0,fname1);
        if ( 0 == strcmp(fname1,"approach")) {
            if (!mxIsStruct(tmp1)) {
                mexErrMsgTxt("MSETGS.approach must be a structure.");
            }
            for (i2=0; i2<mxGetNumberOfFields(tmp1); ++i2) {
                fname2 = mxGetFieldNameByNumber(tmp1, i2);
                tmp2 = mxGetField(tmp1,0,fname2);
                if ( 0 == strcmp(fname2,"warmstartInner")) {
                    if (!mxIsDouble(tmp2)) {
                        mexErrMsgTxt("MSETGS.approach.warmstartInner must be a double.");
                    }
                    if ( mxGetM(tmp2) != 1 || mxGetN(tmp2) != 1 ) {
                        mexErrMsgTxt("MSETGS.approach.warmstartInner must be matrix of dimension 1x1.");
                    }
                    settings.approach.warmstartInner = (int) *mxGetPr(tmp2);
                }
                else
                {
                    mexPrintf("Warning: MSETGS.approach contains the superfluous field '%s'.\n",fname2);
                }
            }
        }
        else
        if ( 0 == strcmp(fname1,"algoInner")) {
            if (!mxIsStruct(tmp1)) {
                mexErrMsgTxt("MSETGS.algoInner must be a structure.");
            }
            for (i2=0; i2<mxGetNumberOfFields(tmp1); ++i2) {
                fname2 = mxGetFieldNameByNumber(tmp1, i2);
                tmp2 = mxGetField(tmp1,0,fname2);
                if ( 0 == strcmp(fname2,"init")) {
                    if (!mxIsDouble(tmp2)) {
                        mexErrMsgTxt("MSETGS.algoInner.init must be a double.");
                    }
                    if ( mxGetM(tmp2) != 30 || mxGetN(tmp2) != 1 ) {
                        mexErrMsgTxt("MSETGS.algoInner.init must be matrix of dimension 30x1.");
                    }
                    copyMArrayToC(mxGetPr(tmp2),settings.algoInner.init,30);
                }
                else
                if ( 0 == strcmp(fname2,"maxit")) {
                    if (!mxIsDouble(tmp2)) {
                        mexErrMsgTxt("MSETGS.algoInner.maxit must be a double.");
                    }
                    if ( mxGetM(tmp2) != 1 || mxGetN(tmp2) != 1 ) {
                        mexErrMsgTxt("MSETGS.algoInner.maxit must be matrix of dimension 1x1.");
                    }
                    settings.algoInner.maxit = (int) *mxGetPr(tmp2);
                }
                else
                if ( 0 == strcmp(fname2,"stopgEps")) {
                    if (!mxIsDouble(tmp2)) {
                        mexErrMsgTxt("MSETGS.algoInner.stopgEps must be a double.");
                    }
                    if ( mxGetM(tmp2) != 1 || mxGetN(tmp2) != 1 ) {
                        mexErrMsgTxt("MSETGS.algoInner.stopgEps must be matrix of dimension 1x1.");
                    }
                    settings.algoInner.stopgEps = (realtype) *mxGetPr(tmp2);
                }
                else
                if ( 0 == strcmp(fname2,"stopgStride")) {
                    if (!mxIsDouble(tmp2)) {
                        mexErrMsgTxt("MSETGS.algoInner.stopgStride must be a double.");
                    }
                    if ( mxGetM(tmp2) != 1 || mxGetN(tmp2) != 1 ) {
                        mexErrMsgTxt("MSETGS.algoInner.stopgStride must be matrix of dimension 1x1.");
                    }
                    settings.algoInner.stopgStride = (int) *mxGetPr(tmp2);
                }
                else
                {
                    mexPrintf("Warning: MSETGS.algoInner contains the superfluous field '%s'.\n",fname2);
                }
            }
        }
        else
        if ( 0 == strcmp(fname1,"algoOuter")) {
            if (!mxIsStruct(tmp1)) {
                mexErrMsgTxt("MSETGS.algoOuter must be a structure.");
            }
            for (i2=0; i2<mxGetNumberOfFields(tmp1); ++i2) {
                fname2 = mxGetFieldNameByNumber(tmp1, i2);
                tmp2 = mxGetField(tmp1,0,fname2);
                if ( 0 == strcmp(fname2,"init")) {
                    if (!mxIsDouble(tmp2)) {
                        mexErrMsgTxt("MSETGS.algoOuter.init must be a double.");
                    }
                    if ( mxGetM(tmp2) != 20 || mxGetN(tmp2) != 1 ) {
                        mexErrMsgTxt("MSETGS.algoOuter.init must be matrix of dimension 20x1.");
                    }
                    copyMArrayToC(mxGetPr(tmp2),settings.algoOuter.init,20);
                }
                else
                if ( 0 == strcmp(fname2,"maxit")) {
                    if (!mxIsDouble(tmp2)) {
                        mexErrMsgTxt("MSETGS.algoOuter.maxit must be a double.");
                    }
                    if ( mxGetM(tmp2) != 1 || mxGetN(tmp2) != 1 ) {
                        mexErrMsgTxt("MSETGS.algoOuter.maxit must be matrix of dimension 1x1.");
                    }
                    settings.algoOuter.maxit = (int) *mxGetPr(tmp2);
                }
                else
                if ( 0 == strcmp(fname2,"stopgEps")) {
                    if (!mxIsDouble(tmp2)) {
                        mexErrMsgTxt("MSETGS.algoOuter.stopgEps must be a double.");
                    }
                    if ( mxGetM(tmp2) != 1 || mxGetN(tmp2) != 1 ) {
                        mexErrMsgTxt("MSETGS.algoOuter.stopgEps must be matrix of dimension 1x1.");
                    }
                    settings.algoOuter.stopgEps = (realtype) *mxGetPr(tmp2);
                }
                else
                if ( 0 == strcmp(fname2,"stopgStride")) {
                    if (!mxIsDouble(tmp2)) {
                        mexErrMsgTxt("MSETGS.algoOuter.stopgStride must be a double.");
                    }
                    if ( mxGetM(tmp2) != 1 || mxGetN(tmp2) != 1 ) {
                        mexErrMsgTxt("MSETGS.algoOuter.stopgStride must be matrix of dimension 1x1.");
                    }
                    settings.algoOuter.stopgStride = (int) *mxGetPr(tmp2);
                }
                else
                {
                    mexPrintf("Warning: MSETGS.algoOuter contains the superfluous field '%s'.\n",fname2);
                }
            }
        }
        else
        {
            mexPrintf("Warning: MSETGS contains the superfluous field '%s'.\n",fname1);
        }
    }
    
    }
    }

    /* call the solver */
    ex2_solve(&params, &settings, &result, &work);

    /*write result to mresult in plhs */
    {
    
    mxArray *tmp1;
    
    {
    const char *fnames[] = {"la", "x", "d", "iter", "exitflag"};
    MRES = mxCreateStructMatrix(1,1, 5, fnames);
        tmp1 = mxCreateDoubleMatrix(20,1,mxREAL);
        copyCArrayToM( result.la, mxGetPr(tmp1), 20 );
        mxSetField(MRES, 0, "la", tmp1);
    
        tmp1 = mxCreateDoubleMatrix(30,1,mxREAL);
        copyCArrayToM( result.x, mxGetPr(tmp1), 30 );
        mxSetField(MRES, 0, "x", tmp1);
    
        tmp1 = mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(tmp1) = (double) result.d;
        mxSetField(MRES, 0, "d", tmp1);
    
        tmp1 = mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(tmp1) = (double) result.iter;
        mxSetField(MRES, 0, "iter", tmp1);
    
        tmp1 = mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(tmp1) = (double) result.exitflag;
        mxSetField(MRES, 0, "exitflag", tmp1);
    
    }
    
    }

}