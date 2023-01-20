function [R,q] = ned2seu()
% Compute rotation matrix which transforms from North-East-Down (NED)
% coordinates to the South-East-Up (SEU)

    R = [-1   0   0;
          0   1   0;
          0   0  -1];

    % assert(det(R)==1,"Not a valid DCM.")

    q = q_from_dcm(R);

end