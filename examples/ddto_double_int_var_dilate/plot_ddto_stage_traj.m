function [] = plot_ddto_stage_traj(rbar,r,prb)
    cval = {[0.6,0.4,0],[0.6,0,0.4],[0,0.4,0.6],[0.6,0.4,0.6]}; 
    mrksz = 16;
    if prb.n == 2
        for j = 1:prb.ntarg
            plot(rbar(1,1:prb.Kstr,j),rbar(2,1:prb.Kstr,j),'.','MarkerSize',mrksz,'Color',[0,0,0]);
            hold on
            plot(r(1,1:prb.Kstrfine,j),r(2,1:prb.Kstrfine,j),'-','MarkerSize',mrksz,'Color',[0,0,0]);            
            plot(rbar(1,prb.Kstr:end,j),rbar(2,prb.Kstr:end,j),'.','MarkerSize',mrksz,'Color',cval{j});
            plot(r(1,prb.Kstrfine:end,j),r(2,prb.Kstrfine:end,j),'-','MarkerSize',mrksz,'Color',cval{j});            
        end
        plot(rbar(1,1,1),rbar(2,1,1),'d','MarkerSize',mrksz,'Color','black');
        plot(rbar(1,prb.Kstr,1),rbar(2,prb.Kstr,1),'d','MarkerSize',mrksz,'Color','black');
        ax = gca;
        ax.DataAspectRatio = [1,1,1];
        ax.PlotBoxAspectRatio = [1,1,1];   
        ax.XLim = [-10,80];
        ax.YLim = [-10,80];
    end
end