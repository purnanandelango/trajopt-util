#ifndef MPC_CONST_H
#define MPC_CONST_H

#include "mpc_base.h"

#define STATE_CONSTR  /* State constrained problem. */ 


enum { 
MPC_HOR = 10,  /**< MPC prediction horizon. */
MPC_STATES = 20,  /**< Number of system states. */
MPC_INPUTS = 10,  /**< Number of system inputs. */
MPC_MXCONSTRS = 20, /**< Number of mixed stage constraints. */
MPC_HOR_INPUTS = 100,  /**< Horizon times number of inputs. */
MPC_HOR_STATES = 200,  /**< Horizon times number of states. */
MPC_HOR_MXCONSTRS = 220  /**< Horizon times number of mixed constrained
plus the number of end state constraints. */
}; 

#endif /* MPC_CONST_H */
