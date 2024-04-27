function rel_err = diagnose_ctscvx_noparam(xbar,ubar,prb,sys_cnstr_cost_fun,alg1,alg2)
% Compare output of two ct-SCvx implementations
% alg1  = {x,y}
% x     = "" or "dvar_"
% y     = [] or "handparse"

    prb.scp_iters = 10;

    if isempty(alg1{2})
        [xbar1,ubar1] = feval("scp.ctscvx_"+alg1{1}+"noparam",xbar,ubar,prb,sys_cnstr_cost_fun);    
    elseif alg1{2} == "handparse"
        [xbar1,ubar1] = feval("scp.ctscvx_"+alg1{1}+"handparse_noparam",xbar,ubar,prb);
    else
        error("Incorrect choice of first algorithm");
    end

    if isempty(alg2{2})
        [xbar2,ubar2] = feval("scp.ctscvx_"+alg2{1}+"noparam",xbar,ubar,prb,sys_cnstr_cost_fun);    
    elseif alg2{2} == "handparse"
        [xbar2,ubar2] = feval("scp.ctscvx_"+alg2{1}+"handparse_noparam",xbar,ubar,prb);
    else
        error("Incorrect choice of second algorithm");
    end  

    rel_err = norm([xbar1(:);ubar1(:)]-[xbar2(:);ubar2(:)])/norm([xbar1(:);ubar1(:)]);

end