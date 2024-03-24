function plot_vehicle_forces(u,rI,vI,qBI,bdscl,thrust,drag,idx)
    K = size(rI,2);

    xBI = plant.rocket6DoF.compute_bodyaxis(qBI);
    uBI = plant.rocket6DoF.compute_thrustdir(u,qBI);
    if ~isnan(drag.scl)
        AI = plant.rocket6DoF.compute_dragdir(drag,vI,qBI);
    end
    
    j1 = idx{1}(1);
    j2 = idx{1}(2);
    j3 = idx{1}(3);
    
    % plot3(rI(j1,:),rI(j2,:),rI(j3,:),'.b')

    xlabel(horzcat('$',idx{2}{1},'$'))
    ylabel(horzcat('$',idx{2}{2},'$'))
    zlabel(horzcat('$',idx{2}{3},'$'))
    
    % rylim = max(abs(rI(2,:)));
    % rzlim = max(abs(rI(3,:)));
    % ryzlim = 1.1*max([rylim,rzlim]);
    % rxlim = 1.1*max(abs(rI(1,:)));
    
    % axis([-ryzlim,ryzlim,-ryzlim,ryzlim,-0.08*rxlim,rxlim])
    title('Trajectory');

    hold on 

    for k=1:K
        if k<K
        quiver3(rI(j1,k),rI(j2,k),rI(j3,k),...
            bdscl*xBI(j1,k),bdscl*xBI(j2,k),bdscl*xBI(j3,k),...
            'Color','red','LineWidth',4);
        quiver3(rI(j1,k),rI(j2,k),rI(j3,k),...
            -thrust*uBI(j1,k),-thrust*uBI(j2,k),-thrust*uBI(j3,k),...
            'Color',[0,0.7,0.3],'LineWidth',4);
        if ~isnan(drag.scl)
            quiver3(rI(j1,k),rI(j2,k),rI(j3,k),...
                drag.scl*AI(j1,k),drag.scl*AI(j2,k),drag.scl*AI(j3,k),...
                'Color',[0.9,0.7,0],'LineWidth',4);
        end 
        end
        if k==K
        qplot_body = quiver3(rI(j1,k),rI(j2,k),rI(j3,k),...
            bdscl*xBI(j1,k),bdscl*xBI(j2,k),bdscl*xBI(j3,k),...
            'Color','red','LineWidth',4);
        qplot_thrust = quiver3(rI(j1,k),rI(j2,k),rI(j3,k),...
            -thrust*uBI(j1,k),-thrust*uBI(j2,k),-thrust*uBI(j3,k),...
            'Color',[0,0.7,0.3],'LineWidth',4);
        if ~isnan(drag.scl)
            qplot_drag = quiver3(rI(j1,k),rI(j2,k),rI(j3,k),...
                drag.scl*AI(j1,k),drag.scl*AI(j2,k),drag.scl*AI(j3,k),...
                'Color',[0.9,0.7,0],'LineWidth',4);
        end
        end
    end
    if ~isnan(drag.scl)
        legend([qplot_body,qplot_thrust,qplot_drag],{horzcat('Body axis : $C_{\mathcal{I}\leftarrow \mathcal{B}} ',idx{3},'_{\mathcal{B}}$')...
            ,'Thrust : $-C_{\mathcal{I}\leftarrow \mathcal{B}}T_{\mathcal{B}}$','Drag vector : $C_{\mathcal{I}\leftarrow \mathcal{B}}A_{\mathcal{B}}$'},'Location','northwest');
    else
        legend([qplot_body,qplot_thrust],{horzcat('Body axis : $C_{\mathcal{I}\leftarrow \mathcal{B}} ',idx{3},'_{\mathcal{B}}$')...
            ,'Thrust : $-C_{\mathcal{I}\leftarrow \mathcal{B}}T_{\mathcal{B}}$'},'Location','northwest');    
    end
    view(-19,18);

    ax = gca;
    ax.DataAspectRatio = [1,1,1];
    ax.PlotBoxAspectRatio = [1,1,1];
end