clearvars
close all

load recent_solution

figure
subplot(2,2,1)
plot3(r(1,:),r(2,:),r(3,:),'-r');
hold on 
plot3(xbar(1,:),xbar(2,:),xbar(3,:),'or');
title('Position')

nrm_T(prb.Kfine) = 0;
nrm_v(prb.Kfine) = 0;
for j = 1:prb.Kfine
    nrm_T(j) = norm(u(1:3,j));
    nrm_v(j) = norm(v(1:3,j));
end
nrm_Tbar(prb.K) = 0;
nrm_vbar(prb.K) = 0;
for j = 1:prb.K
    nrm_Tbar(j) = norm(ubar(1:3,j));
    nrm_vbar(j) = norm(xbar(4:6,j));
end

subplot(2,2,2)
plot(tvec,nrm_v,'-m');
hold on 
plot(tvecbar,nrm_vbar,'om');
title('Velocity')

subplot(2,2,3)
plot(tvec,nrm_T,'-b');
hold on 
plot(tvecbar,nrm_Tbar,'ob');
title('Thrust');

subplot(2,2,4)
plot(tau,tvec,'-k');
hold on
plot(tau,u(4,:),'--g');