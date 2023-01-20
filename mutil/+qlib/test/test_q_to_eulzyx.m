clearvars
close all
clc

q = qlib.q_rand_unit();
s1 = quat2eul(q([4,1,2,3])',"ZYX")*180/pi;
s2 = qlib.q_to_eulzyx(q)';

norm(s1-s2)