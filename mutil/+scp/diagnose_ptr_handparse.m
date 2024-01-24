function rel_err = diagnose_ptr_handparse(xbar,ubar,prb,sys_cnstr_cost_fun,flag)


    prb.scp_iters = 10;

    switch flag

        case 'affine-var'

            [xbar1,ubar1] = scp.run_ptr_noparam(xbar,ubar,prb,sys_cnstr_cost_fun);
            [xbar2,ubar2] = scp.run_ptr_handparse_noparam(xbar,ubar,prb);
            
            rel_err = norm([xbar1(:);ubar1(:)]-[xbar2(:);ubar2(:)])/norm([xbar1(:);ubar1(:)]);            

        case 'deviation-var'            

            [xbar1,ubar1] = scp.run_ptr_dvar_noparam(xbar,ubar,prb,@sys_cnstr_cost);            
            [xbar2,ubar2] = scp.run_ptr_dvar_handparse_noparam(xbar,ubar,prb);            

            rel_err = norm([xbar1(:);ubar1(:)]-[xbar2(:);ubar2(:)])/norm([xbar1(:);ubar1(:)]);                        

        otherwise

            error("Incorrect choice of PTR.")

    end

end