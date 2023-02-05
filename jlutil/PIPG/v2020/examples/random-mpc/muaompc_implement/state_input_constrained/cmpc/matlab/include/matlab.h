#ifndef MPC_MATLAB_H
#define MPC_MATLAB_H

#include "mex.h"
#include "mpc.h"

uint32_t matlab_cell_index(
    uint32_t row, 
    uint32_t col, 
    uint32_t cols
);

void matlab_move_matrix(
    real_t* destination, 
    const real_t* source, 
    uint32_t rows,
    uint32_t cols
);

void matlab_add_matrix(
    mxArray* structure,
    const char* name,
    const real_t* data, 
    uint32_t rows,
    uint32_t cols
);

void matlab_add_scalar(
    mxArray* structure, 
    const char* name, 
    real_t data
);

void matlab_read_matrix(
    const mxArray* structure,
    const char* name,
    real_t* data,
    uint32_t rows,
    uint32_t cols
);

void matlab_read_scalar(
    const mxArray* structure,
    const char* name,
    real_t* data
);

void matlab_read_uint(
    const mxArray* structure,
    const char* name,
    uint32_t* data
);

/**
 * These functions create matlab structures
 */
mxArray* matlab_create_conf(const struct mpc_ctl *data);
mxArray* matlab_create_sys(const struct mpc_ctl *data);
mxArray* matlab_create_wmx(const struct mpc_ctl *data);
mxArray* matlab_create_qpx(const struct mpc_ctl *data);
mxArray* matlab_create_fgm(const struct mpc_ctl *data);
mxArray* matlab_create_alm(const struct mpc_ctl *data);
mxArray* matlab_create_ctl(const struct mpc_ctl *data);

/**
 * These functions copy volatile fields from matlab to C structures
 */
void matlab_apply_conf(struct mpc_ctl *data, const mxArray* conf);
void matlab_apply_fgm(struct mpc_ctl *data, const mxArray* fgm);
void matlab_apply_alm(struct mpc_ctl *data, const mxArray* alm);
void matlab_apply_ctl(struct mpc_ctl *data, const mxArray* ctl);

#endif
