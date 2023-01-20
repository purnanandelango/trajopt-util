function eul_ang = q_to_eulzyx(q)
% Compute Tait-Bryan angles in deg from quaternion
% Scalar-last convention in quaternion
    q0 = q(4);
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);

    heading_angle = atan2(2*(q0*q3+q1*q2),1-2*(q2^2+q3^2))*180/pi;     % psi (first rotation about z)
    pitch_angle = asin(2*(q0*q2-q3*q1))*180/pi;                        % theta (second rotation about y)
    bank_angle = atan2(2*(q0*q1+q2*q3),1-2*(q1^2+q2^2))*180/pi;        % phi (third rotation about x)
    eul_ang = [heading_angle,pitch_angle,bank_angle]';
end