clearvars
close all

load recent_solution

figure
subplot(2,2,1)
plot3(r1(1,:),r1(2,:),r1(3,:),'-r');
hold on 
plot3(r2(1,:),r2(2,:),r2(3,:),'--b');
plot3(xbar(1,:),xbar(2,:),xbar(3,:),'or');
plot3(xbar(7,:),xbar(8,:),xbar(9,:),'ob');
title('Position')

nrm_T1(prb.Kfine) = 0;
nrm_v1(prb.Kfine) = 0;
nrm_T2(prb.Kfine) = 0;
nrm_v2(prb.Kfine) = 0;
for j = 1:prb.Kfine
    nrm_T1(j) = norm(T1(:,j));
    nrm_v1(j) = norm(v1(:,j));
    nrm_T2(j) = norm(T2(:,j));
    nrm_v2(j) = norm(v2(:,j));    
end
nrm_T1bar(prb.K) = 0;
nrm_v1bar(prb.K) = 0;
nrm_T2bar(prb.K) = 0;
nrm_v2bar(prb.K) = 0;
for j = 1:prb.K
    nrm_T1bar(j) = norm(ubar(1:3,j));
    nrm_v1bar(j) = norm(xbar(4:6,j));
    nrm_T2bar(j) = norm(ubar(5:7,j));
    nrm_v2bar(j) = norm(xbar(10:12,j));
end

subplot(2,2,2)
plot(tvec1,nrm_v1,'-m');
hold on 
plot(tvec2,nrm_v2,'--g');
plot(tvecbar1,nrm_v1bar,'om');
plot(tvecbar2,nrm_v2bar,'og');
title('Velocity')

subplot(2,2,3)
plot(tvec1,nrm_T1,'-b');
hold on 
plot(tvec2,nrm_T2,'--r');
plot(tvecbar1,nrm_T1bar,'ob');
plot(tvecbar2,nrm_T2bar,'or');
title('Thrust');

subplot(2,2,4)
plot(tau,tvec1,'-r');
hold on
plot(tau,tvec2,'--m');
plot(tau,s1,'--g');
plot(tau,s2,'--b');