#include "matlab.h"

extern struct mpc_ctl ctl;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs > 0) {
        mexErrMsgTxt("Function doesn't accept any input arguments.");
        return;
    }

    if (nlhs != 1) {
        mexErrMsgTxt("Function accepts only one output variable.");
        return;
    }

    plhs[0] = matlab_create_ctl(&ctl);
    return;
}