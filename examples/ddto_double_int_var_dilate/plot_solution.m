clearvars
close all

load recent_solution
cval = {[0.6,0.4,0],[0.6,0,0.4],[0,0.4,0.6],[0.6,0.4,0.6]};


figure
subplot(2,2,1)
if prb.n == 3
    plot3(r(1,1:prb.Kstrfine,1),r(2,1:prb.Kstrfine,1),r(3,1:prb.Kstrfine,1),'-k');
    hold on 
    for j = 1:prb.ntarg
        plot3(r(1,prb.Kstrfine:end,j),r(2,prb.Kstrfine:end,j),r(3,prb.Kstrfine:end,j),'-','Color',cval{j});
        hold on 
    end
else
    plot(r(1,1:prb.Kstrfine,1),r(2,1:prb.Kstrfine,1),'-k');
    hold on 
    for j = 1:prb.ntarg
        plot(r(1,prb.Kstrfine:end,j),r(2,prb.Kstrfine:end,j),'-','Color',cval{j});
        hold on 
    end
end
title('Position');
ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];


subplot(2,2,2)
plot(tvecbar(1:prb.Kstr,1),nrm_vbar(1:prb.Kstr,1),'-ok');
hold on
for j = 1:prb.ntarg
    plot(tvec(prb.Kstrfine:end,j),nrm_v(prb.Kstrfine:end,j),'-','Color',cval{j});
    plot(tvecbar(prb.Kstr:end,j),nrm_vbar(prb.Kstr:end,j),'o','Color',cval{j});
end
title('Velocity')
xlabel('$t$');

subplot(2,2,3)
plot(tvecbar(1:prb.Kstr,1),nrm_Tbar(1:prb.Kstr,1),'-ok');
hold on
for j = 1:prb.ntarg
    plot(tvec(prb.Kstrfine:end,j),nrm_T(prb.Kstrfine:end,j),'-','Color',cval{j});
    plot(tvecbar(prb.Kstr:end,j),nrm_Tbar(prb.Kstr:end,j),'o','Color',cval{j});
end
title('Thrust');
xlabel('$t$');

subplot(2,2,4)
for j = 1:prb.ntarg
    plot(tau,tvec(:,j),'-','Color',cval{j});
    hold on
    plot(tau,s(:,j),'--','Color',cval{j});    
end