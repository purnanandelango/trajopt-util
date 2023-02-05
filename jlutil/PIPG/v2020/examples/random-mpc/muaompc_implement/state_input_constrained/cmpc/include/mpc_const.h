#ifndef MPC_CONST_H
#define MPC_CONST_H

#include "mpc_base.h"

#define STATE_CONSTR  /* State constrained problem. */ 


enum { 
MPC_HOR = 10,  /**< MPC prediction horizon. */
MPC_STATES = 8,  /**< Number of system states. */
MPC_INPUTS = 4,  /**< Number of system inputs. */
MPC_MXCONSTRS = 8, /**< Number of mixed stage constraints. */
MPC_HOR_INPUTS = 40,  /**< Horizon times number of inputs. */
MPC_HOR_STATES = 80,  /**< Horizon times number of states. */
MPC_HOR_MXCONSTRS = 88  /**< Horizon times number of mixed constrained
plus the number of end state constraints. */
}; 

#endif /* MPC_CONST_H */
