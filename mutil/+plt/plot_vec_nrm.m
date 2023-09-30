function plot_vec_nrm(tvec,v,vbnd,norm_type,title_str,varargin)
% Plot the norm of vector-valued signal
    if nargin >= 7 % Caution: plot line can be set only after setting ax.YScale
        line_color = varargin{2};
    else
        line_color = [0,0.2,0.7];
    end
    K = size(v,2);
    nrm_v(K) = 0;
    for k=1:K
        nrm_v(k) = norm(v(:,k),norm_type);
    end
    if nargin == 8
        display_name = varargin{3};
        plot(tvec,nrm_v,'-','Color',line_color,'DisplayName',display_name);
        legend('AutoUpdate','off');
    else
        plot(tvec,nrm_v,'-','Color',line_color)        
    end    
    hold on
    if length(vbnd) == 2
        plot(tvec,vbnd(1)*ones(1,K),'--r')
        plot(tvec,vbnd(2)*ones(1,K),'--r')
        ylim([0.5*min(vbnd),1.1*max(vbnd)]);
    else
        plot(tvec,vbnd*ones(1,K),'--r')
        ylim([0,1.1*vbnd]);
    end
    title(title_str);
    xlim([0,tvec(end)]);
    if nargin >= 6
        ax = gca;
        ax.YScale = varargin{1};
    end
end