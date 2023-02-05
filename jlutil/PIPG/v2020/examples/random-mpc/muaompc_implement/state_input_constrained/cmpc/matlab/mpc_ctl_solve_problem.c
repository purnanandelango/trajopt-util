#include "matlab.h"

extern struct mpc_ctl ctl;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgTxt("Function needs two input variables.");
        return;
    }

    if (!mxIsStruct(prhs[0])) {
        mexErrMsgTxt("Input argument has to be a mpc structure.");
        return;
    }

    if (nlhs != 1) {
        mexErrMsgTxt("Function accepts only one output variable.");
        return;
    }

    matlab_apply_ctl(&ctl, prhs[0]);

    mpc_ctl_solve_problem(&ctl, mxGetPr(prhs[1]));

    plhs[0] = matlab_create_ctl(&ctl);
    return;
}
