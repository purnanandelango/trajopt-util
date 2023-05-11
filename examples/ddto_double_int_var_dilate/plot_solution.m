clearvars
close all

load recent_solution_3
% load solve_2
cval = {[0.6,0.4,0],[0.6,0,0.4],[0,0.4,0.6],[0.6,0.4,0.6]};

figure
subplot(2,2,1)
if prb.n == 3    
    for j = 1:prb.ntarg
        plot3(r(1,1:prb.Kstrfine,j),r(2,1:prb.Kstrfine,j),r(3,1:prb.Kstrfine,j),'-k');
        hold on 
        plot3(rbar(1,1:prb.Kstr,j),rbar(2,1:prb.Kstr,j),rbar(3,1:prb.Kstr,j),'+k');        
        plot3(r(1,prb.Kstrfine:end,j),r(2,prb.Kstrfine:end,j),r(3,prb.Kstrfine:end,j),'-','Color',cval{j});
        plot(rbar(1,prb.Kstr+1:end,j),rbar(2,prb.Kstr+1:end,j),rbar(3,prb.Kstr+1:end,j),'o','Color',cval{j});        
    end
elseif prb.n == 2    
    for j = 1:prb.ntarg
        plot(r(1,1:prb.Kstrfine,j),r(2,1:prb.Kstrfine,j),'-k');
        hold on 
        plot(rbar(1,1:prb.Kstr,j),rbar(2,1:prb.Kstr,j),'+k');        
        plot(r(1,prb.Kstrfine:end,j),r(2,prb.Kstrfine:end,j),'-','Color',cval{j});
        plot(rbar(1,prb.Kstr+1:end,j),rbar(2,prb.Kstr+1:end,j),'o','Color',cval{j});
    end
else
    error('Only 2D and 3D cases are possible.')
end
title('Position');

if isfield(prb,'robs')
    if prb.n == 2
        th = linspace(0,2*pi);
        for j = 1:prb.nobs
            pobs = prb.robs(:,j) + prb.aobs(j)*[cos(th);sin(th)];
            plot(pobs(1,:),pobs(2,:),'-k','LineWidth',0.5);
        end
    end
end

ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];


subplot(2,2,2)
for j = 1:prb.ntarg
    plot(tvec(1:prb.Kstrfine,j),nrm_v(1:prb.Kstrfine,j),'-k')
    hold on
    plot(tvecbar(1:prb.Kstr,j),nrm_vbar(1:prb.Kstr,j),'ok');    
    plot(tvec(prb.Kstrfine:end,j),nrm_v(prb.Kstrfine:end,j),'-','Color',cval{j});
    plot(tvecbar(prb.Kstr+1:end,j),nrm_vbar(prb.Kstr+1:end,j),'o','Color',cval{j});
end
title('Velocity')
xlabel('$t$');

subplot(2,2,3)
for j = 1:prb.ntarg
    plot(tvec(1:prb.Kstrfine,j),nrm_T(1:prb.Kstrfine,j),'-k');
    hold on
    plot(tvecbar(1:prb.Kstr,j),nrm_Tbar(1:prb.Kstr,j),'ok');    
    plot(tvec(prb.Kstrfine:end,j),nrm_T(prb.Kstrfine:end,j),'-','Color',cval{j});
    plot(tvecbar(prb.Kstr:end,j),nrm_Tbar(prb.Kstr:end,j),'o','Color',cval{j});
end
title('Thrust');
xlabel('$t$');

subplot(2,2,4)
for j = 1:prb.ntarg
    plt_t = plot(prb.tau,tvecbar(:,j),'d','Color',cval{j});
    hold on
    plot(tau,tvec(:,j),'-','Color',cval{j});    
    plt_dt = plot(prb.tau(1:end-1),diff(tvecbar(:,j)),'-s','Color',cval{j});    
    plt_s = plot(prb.tau,sbar(:,j),'+','Color',cval{j});
    plot(tau,s(:,j),'-','Color',cval{j});    
end
legend([plt_t,plt_dt,plt_s],{'$t$','$\Delta t$','$s$'},'Location','bestoutside')