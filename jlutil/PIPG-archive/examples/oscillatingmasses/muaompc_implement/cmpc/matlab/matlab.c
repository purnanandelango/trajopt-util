#include "matlab.h"

uint32_t matlab_cell_index(
    uint32_t row,
    uint32_t col,
    uint32_t cols
) {
    return row * cols + col;
}

void matlab_move_matrix(
    real_t* destination,
    const real_t* source,
    uint32_t rows,
    uint32_t cols
) {
    unsigned int row = 0;
    unsigned int col = 0;

    /* loops are neccessary in order to transpose matrices */
    for (row = 0; row < rows; row++) {
        for (col = 0; col < cols; col++) {
            destination[matlab_cell_index(col, row, rows)] =
                (double)source[matlab_cell_index(row, col, cols)];
        }
        col = 0;
    }
}

void matlab_add_matrix(
    mxArray* structure,
    const char* name,
    const real_t* data,
    uint32_t rows,
    uint32_t cols
) {
    mxArray* matrix = mxCreateDoubleMatrix(rows, cols, 0);
    matlab_move_matrix(mxGetPr(matrix), data, rows, cols);
    mxAddField(structure, name);
    mxSetField(structure, 0, name, matrix);
}

void matlab_add_scalar(
    mxArray* structure,
    const char* name,
    real_t data
) {
    mxArray* scalar = mxCreateDoubleScalar(data);
    mxAddField(structure, name);
    mxSetField(structure, 0, name, scalar);
}

void matlab_read_matrix(
    const mxArray* structure,
    const char* name,
    real_t* data,
    uint32_t rows,
    uint32_t cols
) {
    mxArray* matrix = mxGetField(structure, 0, name);
    matlab_move_matrix(data, mxGetPr(matrix), rows, cols);
}

void matlab_read_scalar(
    const mxArray* structure,
    const char* name,
    real_t* data
) {
    mxArray* scalar = mxGetField(structure, 0, name);
    *data = mxGetScalar(scalar);
}

void matlab_read_uint(
    const mxArray* structure,
    const char* name,
    uint32_t* data
) {
    double temporary;
    matlab_read_scalar(structure, name, &temporary);
    *data = (uint32_t)temporary;
}

mxArray* matlab_create_conf(const struct mpc_ctl *data) {
    mxArray* conf = mxCreateStructMatrix(1, 1, 0 , 0);

    matlab_add_scalar(conf, "in_iter", data->conf->in_iter);
    matlab_add_scalar(conf, "ex_iter", data->conf->ex_iter);
    matlab_add_scalar(conf, "warmstart", data->conf->warmstart);

    return conf;
}

mxArray* matlab_create_sys(const struct mpc_ctl *data) {
    mxArray* sys = mxCreateStructMatrix(1, 1, 0, 0);

    matlab_add_matrix(sys, "Ad", data->sys->Ad,
        data->alm->STATES,
        data->alm->fgm->STATES
    );

    matlab_add_matrix(sys, "Bd", data->sys->Bd,
        data->alm->fgm->STATES,
        data->alm->fgm->INPUTS
    );

    matlab_add_scalar(sys, "dt", *(data->sys->dt));

    return sys;
}

mxArray* matlab_create_wmx(const struct mpc_ctl *data) {
    mxArray* wmx = mxCreateStructMatrix(1, 1, 0, 0);

    matlab_add_matrix(wmx, "Q", data->wmx->Q,
        data->alm->fgm->STATES,
        data->alm->fgm->STATES
    );

    matlab_add_matrix(wmx, "R", data->wmx->R,
        data->alm->fgm->INPUTS,
        data->alm->fgm->INPUTS
    );

    matlab_add_matrix(wmx, "P", data->wmx->P,
        data->alm->fgm->STATES,
        data->alm->fgm->STATES
    );

    return wmx;
}

mxArray* matlab_create_qpx(const struct mpc_ctl *data) {
    mxArray* qpx = mxCreateStructMatrix(1, 1, 0, 0);

    matlab_add_matrix(qpx, "HoL", data->qpx->HoL,
        data->qpx->HOR_INPUTS,
        data->qpx->HOR_INPUTS
    );

    matlab_add_matrix(qpx, "gxoL", data->qpx->gxoL,
        data->qpx->HOR_INPUTS, 1
    );

    matlab_add_matrix(qpx, "E", data->qpx->E,
        data->qpx->HOR_MXCONSTRS,
        data->qpx->HOR_INPUTS
    );

    matlab_add_matrix(qpx, "u_lb", data->qpx->u_lb,
        data->qpx->HOR_INPUTS, 1
    );

    matlab_add_matrix(qpx, "u_ub", data->qpx->u_ub,
        data->qpx->HOR_INPUTS, 1
    );

    matlab_add_matrix(qpx, "zx_lb", data->qpx->zx_lb,
        data->qpx->HOR_MXCONSTRS, 1
    );

    matlab_add_matrix(qpx, "zx_ub", data->qpx->zx_ub,
        data->qpx->HOR_MXCONSTRS, 1
    );

    matlab_add_scalar(qpx, "HOR_INPUTS", data->qpx->HOR_INPUTS);
    matlab_add_scalar(qpx, "HOR_MXCONSTRS", data->qpx->HOR_MXCONSTRS);

    return qpx;
}

mxArray* matlab_create_fgm(const struct mpc_ctl *data) {
    mxArray* fgm = mxCreateStructMatrix(1, 1, 0, 0);

    matlab_add_matrix(fgm, "u_0", data->alm->fgm->u_0,
        data->alm->fgm->HOR_INPUTS, 1
    );

    matlab_add_matrix(fgm, "gxoL", data->alm->fgm->gxoL,
        data->alm->fgm->HOR_INPUTS, 1
    );

    matlab_add_matrix(fgm, "groL", data->alm->fgm->groL,
        data->alm->fgm->HOR_INPUTS, 1
    );

    matlab_add_scalar(fgm, "j_in", *(data->alm->fgm->j_in));

    matlab_add_matrix(fgm, "HoL", data->alm->fgm->HoL,
        data->alm->fgm->HOR_INPUTS,
        data->alm->fgm->HOR_INPUTS
    );

    matlab_add_matrix(fgm, "GoL", data->alm->fgm->GoL,
        data->alm->fgm->HOR_INPUTS,
        data->alm->fgm->STATES
    );

    matlab_add_matrix(fgm, "Bh_T", data->alm->fgm->Bh_T,
        data->alm->fgm->HOR_INPUTS,
        data->alm->fgm->HOR_STATES
    );

    matlab_add_matrix(fgm, "u_lb", data->alm->fgm->u_lb,
        data->alm->fgm->HOR_INPUTS, 1
    );

    matlab_add_matrix(fgm, "u_ub", data->alm->fgm->u_ub,
        data->alm->fgm->HOR_INPUTS, 1
    );

    matlab_add_scalar(fgm, "nu", *(data->alm->fgm->nu));

    matlab_add_scalar(fgm, "HOR", data->alm->fgm->HOR);
    matlab_add_scalar(fgm, "STATES", data->alm->fgm->STATES);
    matlab_add_scalar(fgm, "INPUTS", data->alm->fgm->INPUTS);
    matlab_add_scalar(fgm, "HOR_INPUTS", data->alm->fgm->HOR_INPUTS);
    matlab_add_scalar(fgm, "HOR_STATES", data->alm->fgm->HOR_STATES);

    return fgm;
}

mxArray* matlab_create_alm(const struct mpc_ctl *data) {
    mxArray* alm = mxCreateStructMatrix(1, 1, 0, 0);

    mxArray* fgm = matlab_create_fgm(data);
    mxAddField(alm, "fgm");
    mxSetField(alm, 0, "fgm", fgm);

    matlab_add_matrix(alm, "l_0", data->alm->l_0,
        data->alm->HOR_MXCONSTRS, 1
    );

    matlab_add_matrix(alm, "zx_lb", data->alm->zx_lb,
        data->alm->HOR_MXCONSTRS, 1
    );

    matlab_add_matrix(alm, "zx_ub", data->alm->zx_ub,
        data->alm->HOR_MXCONSTRS, 1
    );

    matlab_add_scalar(alm, "i_ex", *(data->alm->i_ex));
    matlab_add_scalar(alm, "mu", *(data->alm->mu));

    matlab_add_matrix(alm, "E", data->alm->E,
        data->alm->HOR_MXCONSTRS,
        data->alm->HOR_INPUTS
    );

    matlab_add_matrix(alm, "Kx_Ai", data->alm->Kx_Ai,
        data->alm->HOR_MXCONSTRS,
        data->alm->STATES
    );

    matlab_add_matrix(alm, "e_lb", data->alm->e_lb,
        data->alm->HOR_MXCONSTRS, 1
    );

    matlab_add_matrix(alm, "e_ub", data->alm->e_ub,
        data->alm->HOR_MXCONSTRS, 1
    );

    matlab_add_scalar(alm, "Linv", *(data->alm->Linv));

    matlab_add_scalar(alm, "STATES", data->alm->STATES);
    matlab_add_scalar(alm, "MXCONSTRS", data->alm->MXCONSTRS);
    matlab_add_scalar(alm, "HOR_INPUTS", data->alm->HOR_INPUTS);
    matlab_add_scalar(alm, "HOR_MXCONSTRS", data->alm->HOR_MXCONSTRS);

    return alm;
}

mxArray* matlab_create_ctl(const struct mpc_ctl *data) {
    mxArray* ctl;
    mxArray* conf;
    mxArray* alm;
    mxArray* sys;
    mxArray* wmx;
    mxArray* qpx;

    ctl = mxCreateStructMatrix(1, 1, 0, 0);

    conf = matlab_create_conf(data);
    mxAddField(ctl, "conf");
    mxSetField(ctl, 0, "conf", conf);

    alm = matlab_create_alm(data);
    mxAddField(ctl, "alm");
    mxSetField(ctl, 0, "alm", alm);

    sys = matlab_create_sys(data);
    mxAddField(ctl, "sys");
    mxSetField(ctl, 0, "sys", sys);

    wmx = matlab_create_wmx(data);
    mxAddField(ctl, "wmx");
    mxSetField(ctl, 0, "wmx", wmx);

    qpx = matlab_create_qpx(data);
    mxAddField(ctl, "qpx");
    mxSetField(ctl, 0, "qpx", qpx);

    matlab_add_matrix(ctl, "u_opt", data->u_opt,
        data->alm->fgm->HOR_INPUTS, 1
    );

    matlab_add_matrix(ctl, "l_opt", data->l_opt,
        data->alm->HOR_MXCONSTRS, 1
    );

    matlab_add_matrix(ctl, "u_ini", data->u_ini,
        data->alm->fgm->HOR_INPUTS, 1
    );

    matlab_add_matrix(ctl, "l_ini", data->l_ini,
        data->alm->HOR_MXCONSTRS, 1
    );


    /* u_ref and x_ref are by default zero pointers. The if statements may
     * never execute (left for completeness/defensive programming).
     */
    if (data->u_ref != 0) {
        matlab_add_matrix(ctl, "u_ref", data->u_ref,
            data->alm->fgm->HOR_INPUTS, 1
        );
    }

    if (data->x_ref != 0) {
        matlab_add_matrix(ctl, "x_ref", data->x_ref,
            data->alm->fgm->HOR_STATES, 1
        );
    }

    return ctl;
}

void matlab_apply_conf(struct mpc_ctl *data, const mxArray* conf) {
    matlab_read_uint(conf, "in_iter", &(data->conf->in_iter));
    matlab_read_uint(conf, "ex_iter", &(data->conf->ex_iter));
    matlab_read_uint(conf, "warmstart", &(data->conf->warmstart));
}

void matlab_apply_fgm(struct mpc_ctl *data, const mxArray* fgm) {
    matlab_read_matrix(fgm, "u_0", data->alm->fgm->u_0, data->alm->fgm->HOR_INPUTS, 1);
    matlab_read_matrix(fgm, "gxoL", data->alm->fgm->gxoL, data->alm->fgm->HOR_INPUTS, 1);
    matlab_read_matrix(fgm, "groL", data->alm->fgm->groL, data->alm->fgm->HOR_INPUTS, 1);
}

void matlab_apply_alm(struct mpc_ctl *data, const mxArray* alm) {
    mxArray* fgm = mxGetField(alm, 0, "fgm");
    matlab_apply_fgm(data, fgm);

    matlab_read_matrix(alm, "l_0", data->alm->l_0, data->alm->HOR_MXCONSTRS, 1);
    matlab_read_matrix(alm, "zx_lb", data->alm->zx_lb, data->alm->HOR_MXCONSTRS, 1);
    matlab_read_matrix(alm, "zx_ub", data->alm->zx_ub, data->alm->HOR_MXCONSTRS, 1);
}

void matlab_apply_ctl(struct mpc_ctl *data, const mxArray* ctl) {
    mxArray* conf;
    mxArray* alm;

    conf = mxGetField(ctl, 0, "conf");
    matlab_apply_conf(data, conf);

    alm = mxGetField(ctl, 0, "alm");
    matlab_apply_alm(data, alm);

    matlab_read_matrix(ctl, "u_opt", data->u_opt, data->alm->fgm->HOR_INPUTS, 1);
    matlab_read_matrix(ctl, "l_opt", data->l_opt, data->alm->HOR_MXCONSTRS, 1);
    matlab_read_matrix(ctl, "u_ini", data->u_ini, data->alm->fgm->HOR_INPUTS, 1);
    matlab_read_matrix(ctl, "l_ini", data->l_ini, data->alm->HOR_MXCONSTRS, 1);

    if (mxGetField(ctl, 0, "u_ref")) {
        data->u_ref = mxCalloc(data->alm->fgm->HOR_INPUTS, sizeof(real_t));
        matlab_read_matrix(ctl, "u_ref", data->u_ref, data->alm->fgm->HOR_INPUTS, 1);
    }

    if (mxGetField(ctl, 0, "x_ref")) {
        data->x_ref = mxCalloc(data->alm->fgm->HOR_STATES, sizeof(real_t));
        matlab_read_matrix(ctl, "x_ref", data->x_ref, data->alm->fgm->HOR_STATES, 1);
    }
}

