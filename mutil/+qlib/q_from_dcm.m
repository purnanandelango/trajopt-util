function q = q_from_dcm(dcm)
% Compute quaternion from direction cosine matrix (DCM) 
% SCALAR-LAST convention 

    q_temp = dcm2quat(dcm); % Employs scalar-first convention
    q = reshape(q_temp([2,3,4,1]),[4,1]);
  
% Note that q and -q are both valid solutions    
end