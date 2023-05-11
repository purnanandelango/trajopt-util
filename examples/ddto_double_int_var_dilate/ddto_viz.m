clearvars
close all
clc

path2results = "results/051123/ex7/";

figure
for stage = 1:2
    dat = load(path2results+"solve_"+num2str(stage));
    plot_ddto_stage_traj(dat.rbar,dat.r,dat.prb)
    if stage == 1
        prb = dat.prb;
    end
end

if prb.n == 2
    for j = 1:prb.ntarg
        plot(prb.rK(1,j),prb.rK(2,j),'ks','MarkerSize',16);
    end
    th = linspace(0,2*pi);
    for j = 1:prb.nobs
        pobs = prb.robs(:,j) + prb.aobs(j)*[cos(th);sin(th)];
        plot(pobs(1,:),pobs(2,:),'-r','LineWidth',2);
    end
end