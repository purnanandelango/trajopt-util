function rot_mat = rot_mat_2D(thet)
% 2D rotation matrix: rotation about z-axis by thet degrees 
    rot_mat = [cosd(thet), -sind(thet); sind(thet), cosd(thet)];
end