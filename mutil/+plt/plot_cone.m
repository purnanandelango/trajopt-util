function plot_cone(cone_angle,cone_axis_length,DCM,color_val,face_alpha)
% Plot 3D cone
% cone_angle is half of the cone angle in deg
% cone_axis_length is the length along the cone axis considered for plotting
% DCM rotates the cone axis with respect to some frame
% color_val specifies the color of the cone with RGB triplet
% face_alpha specifies the transparency of the cone
    r_cone = linspace(0,cone_axis_length) ;
    th_cone = linspace(0,2*pi) ;
    [R_cone,Th_cone] = meshgrid(r_cone,th_cone) ;
    X_cone = R_cone;
    Y_cone = tand(cone_angle)*R_cone .* cos(Th_cone);
    Z_cone = tand(cone_angle)*R_cone .* sin(Th_cone);
    C_cone = zeros(100,100,3);
    for i = 1:100
        for j = 1:100
            x_cone = X_cone(i,j);
            y_cone = Y_cone(i,j);
            z_cone = Z_cone(i,j);
            rot_pos = DCM*[x_cone;y_cone;z_cone];
            X_cone(i,j) = rot_pos(1);
            Y_cone(i,j) = rot_pos(2);
            Z_cone(i,j) = rot_pos(3);
            C_cone(i,j,:) = color_val; 
        end
    end
    hold on % make sure that the existing plot information is not overwritten
    surf(X_cone,Y_cone,Z_cone,C_cone,'EdgeColor','none','FaceColor','flat','FaceAlpha',face_alpha)
end