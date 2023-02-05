mex -outdir ./@mpc_ctl -I./../include -I./include ./mpc_ctl_get_data.c ./../mpc_const.c ./matlab.c
mex -outdir ./@mpc_ctl -I./../include -I./include ./mpc_ctl_solve_problem.c ./../mpc.c ./../mpc_const.c ./../mpc_inc.c ./../mpc_ref.c ./../mpc_stc.c ./../mtx_ops.c ./matlab.c
mex -outdir ./@mpc_ctl -I./../include -I./include ./mpc_ctl_form_qp.c ./../mpc.c ./../mpc_const.c ./../mpc_inc.c ./../mpc_ref.c ./../mpc_stc.c ./../mtx_ops.c ./matlab.c
