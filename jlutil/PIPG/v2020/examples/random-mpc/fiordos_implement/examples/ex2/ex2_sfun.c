/**
 * This file is generated by FiOrdOs, a program licensed under GPL
 * by copyright holder Automatic Control Laboratory, ETH Zurich.
 * 
 * If you are interested in using this file commercially,
 * please contact the copyright holder.
 */

#include <stdlib.h>

#include "ex2_solver.c"

#define S_FUNCTION_NAME   ex2_sfun   /**< Name of the S function. */
#define S_FUNCTION_LEVEL  2    /**< S function level. */

#include "simstruc.h"


#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
static void mdlCheckParameters(SimStruct *S)
{
    if ( !mxIsDouble(ssGetSFcnParam(S,0)) || mxGetNumberOfElements(ssGetSFcnParam(S,0)) != 1) {
        ssSetErrorStatus(S,"Parameter 'approach.warmstartInner' to S-function must be a scalar.");
        return;
    }
    if ( !mxIsDouble(ssGetSFcnParam(S,1)) || mxGetNumberOfElements(ssGetSFcnParam(S,1)) != 1) {
        ssSetErrorStatus(S,"Parameter 'algoInner.stopgStride' to S-function must be a scalar.");
        return;
    }
    if ( !mxIsDouble(ssGetSFcnParam(S,2)) || mxGetNumberOfElements(ssGetSFcnParam(S,2)) != 1) {
        ssSetErrorStatus(S,"Parameter 'algoOuter.stopgStride' to S-function must be a scalar.");
        return;
    }
    if ( !mxIsDouble(ssGetSFcnParam(S,3)) || mxGetNumberOfElements(ssGetSFcnParam(S,3)) != 1 || mxGetPr(ssGetSFcnParam(S,3))[0] <= 0 ) {
        ssSetErrorStatus(S,"Parameter 'Ts' to S-function must be a postive scalar.");
        return;
    }
}
#endif /* MDL_CHECK_PARAMETERS */


static void mdlInitializeSizes (SimStruct *S)
{
    int_T i;
    
    /* parameters */
    ssSetNumSFcnParams(S, 4);  /* Number of expected parameters */ 
    #if defined(MATLAB_MEX_FILE)
        if(ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S)) {
            mdlCheckParameters(S);
            if(ssGetErrorStatus(S) != NULL) { return; }
        } else {
            return; /* The Simulink engine reports a mismatch error. */
        }
    #endif
    
    
    /* number of continuous and discrete states */
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);
    
    /* INPUT ports */
    /* number */
    if ( !ssSetNumInputPorts(S, 9) )
    	return;    
    /* dimensions */
    ssSetInputPortMatrixDimensions(S, 0, 30, 1);
    ssSetInputPortMatrixDimensions(S, 1, 1, 1);
    ssSetInputPortMatrixDimensions(S, 2, 20, 1);
    ssSetInputPortMatrixDimensions(S, 3, 30, 1);
    ssSetInputPortMatrixDimensions(S, 4, 1, 1);
    ssSetInputPortMatrixDimensions(S, 5, 1, 1);
    ssSetInputPortMatrixDimensions(S, 6, 20, 1);
    ssSetInputPortMatrixDimensions(S, 7, 1, 1);
    ssSetInputPortMatrixDimensions(S, 8, 1, 1);
    /* direct feedthrough */
    for (i=0; i<9; ++i) {
        ssSetInputPortDirectFeedThrough(S, i, 1);
    }
    
    /* OUTPUT ports */
    /* number */
    if ( !ssSetNumOutputPorts(S, 5) )
    	return;     
    /* dimensions */
    ssSetOutputPortMatrixDimensions(S, 0, 20, 1);
    ssSetOutputPortMatrixDimensions(S, 1, 30, 1);
    ssSetOutputPortMatrixDimensions(S, 2, 1, 1);
    ssSetOutputPortMatrixDimensions(S, 3, 1, 1);
    ssSetOutputPortMatrixDimensions(S, 4, 1, 1);
    

    /* One sample time */
    ssSetNumSampleTimes(S, 1);
    /* size of the block's pointerwork vector */
    ssSetNumPWork(S, 4);
}


static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, *mxGetPr(ssGetSFcnParam(S,3)));
    ssSetOffsetTime(S, 0, 0.0);
}


#define MDL_START /**< Activate call to mdlStart. */
static void mdlStart(SimStruct *S)
{
    /* allocate memory */
    ssGetPWork(S)[0] = (void *) calloc( 1, sizeof(ex2_Params) );
    ssGetPWork(S)[1] = (void *) calloc( 1, sizeof(ex2_Settings) );
    ssGetPWork(S)[2] = (void *) calloc( 1, sizeof(ex2_Result) );
    ssGetPWork(S)[3] = (void *) calloc( 1, sizeof(ex2_Work) );
    
    /* call init */
    ex2_init( (ex2_Params *) ssGetPWork(S)[0], (ex2_Settings *) ssGetPWork(S)[1], (ex2_Result *) ssGetPWork(S)[2], (ex2_Work *) ssGetPWork(S)[3] );
    
    /* adapt settings (from mask) */
    {
        ex2_Settings *settings = (ex2_Settings *) ssGetPWork(S)[1];
        settings->approach.warmstartInner = (int) *mxGetPr(ssGetSFcnParam(S,0));
        settings->algoInner.stopgStride = (int) *mxGetPr(ssGetSFcnParam(S,1));
        settings->algoOuter.stopgStride = (int) *mxGetPr(ssGetSFcnParam(S,2));
    }
}


static void mdlOutputs(SimStruct *S, int_T tid)
{
    int_T i; 
    InputRealPtrsType in_pptr;
    real_T *out_ptr;
    
    ex2_Params *params = (ex2_Params *) ssGetPWork(S)[0];
    ex2_Settings *settings = (ex2_Settings *) ssGetPWork(S)[1];
    ex2_Result *result = (ex2_Result *) ssGetPWork(S)[2];
    ex2_Work *work = (ex2_Work *) ssGetPWork(S)[3];
    
    /* copy params from block input */

    in_pptr = ssGetInputPortRealSignalPtrs(S, 0);
    for ( i=0; i<30; ++i ) {
        params->g[i] = (realtype) *in_pptr[i];
    }
    
    in_pptr = ssGetInputPortRealSignalPtrs(S, 1);
    for ( i=0; i<1; ++i ) {
        params->c[i] = (realtype) *in_pptr[i];
    }
    
    in_pptr = ssGetInputPortRealSignalPtrs(S, 2);
    for ( i=0; i<20; ++i ) {
        params->be[i] = (realtype) *in_pptr[i];
    }
    

    /* adapt settings */

    in_pptr = ssGetInputPortRealSignalPtrs(S, 3);
    for ( i=0; i<30; ++i ) {
        settings->algoInner.init[i] = (realtype) *in_pptr[i];
    }

    in_pptr = ssGetInputPortRealSignalPtrs(S, 4);
    settings->algoInner.maxit = (int) *in_pptr[0];

    in_pptr = ssGetInputPortRealSignalPtrs(S, 5);
    settings->algoInner.stopgEps = (realtype) *in_pptr[0];

    in_pptr = ssGetInputPortRealSignalPtrs(S, 6);
    for ( i=0; i<20; ++i ) {
        settings->algoOuter.init[i] = (realtype) *in_pptr[i];
    }

    in_pptr = ssGetInputPortRealSignalPtrs(S, 7);
    settings->algoOuter.maxit = (int) *in_pptr[0];

    in_pptr = ssGetInputPortRealSignalPtrs(S, 8);
    settings->algoOuter.stopgEps = (realtype) *in_pptr[0];



    /* solve problem */
    ex2_solve(params, settings, result, work);


    /* copy result to block output */

    out_ptr = ssGetOutputPortRealSignal(S, 0);
    for ( i=0; i<20; ++i ) {
        out_ptr[i] = (real_T) result->la[i];
    }
    
    out_ptr = ssGetOutputPortRealSignal(S, 1);
    for ( i=0; i<30; ++i ) {
        out_ptr[i] = (real_T) result->x[i];
    }
    
    out_ptr = ssGetOutputPortRealSignal(S, 2);
    out_ptr[0] = (real_T) result->d;
    
    out_ptr = ssGetOutputPortRealSignal(S, 3);
    out_ptr[0] = (real_T) result->iter;
    
    out_ptr = ssGetOutputPortRealSignal(S, 4);
    out_ptr[0] = (real_T) result->exitflag;
    

}


static void mdlTerminate(SimStruct *S)
{
    int i;
    for ( i=0; i<4; ++i )
    {
        if ( ssGetPWork(S)[i] != NULL )
            free( ssGetPWork(S)[i] );
    }
}


#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif


/*
 *	end of file
 */
