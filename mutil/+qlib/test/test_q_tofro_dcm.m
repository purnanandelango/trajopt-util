clearvars
close all
clc

M = 10;
for j = 1:M
    q = qlib.q_rand_unit();
    q2 = qlib.q_from_dcm(qlib.q_dcm(q)); % Could be -ve of q
    % norm(q - q2)
    if norm(qlib.q_dcm(q) - qlib.q_dcm(q2))>1e-6
        error("DCM estimates don't match.")
    end
end