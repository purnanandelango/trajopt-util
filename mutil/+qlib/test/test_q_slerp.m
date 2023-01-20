clearvars
clc
close all

q1 = qlib.q_rand_unit();
q3 = qlib.q_rand_unit();
t = 0.5;

q2_v1 = qlib.q_slerp(q1, q3, t);
q2_v2 = quatinterp(q1([4,1,2,3])',q3([4,1,2,3])',t,'slerp');
q2_v2 = q2_v2([2,3,4,1])';

diff(q2_v2-q2_v1)
