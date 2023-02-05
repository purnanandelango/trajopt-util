"""Random MPC
"""

import muaompc

mpc = muaompc.ltidt.setup_mpc_problem('problem_data.mat')
mpc.generate_c_files(matlab=True)

