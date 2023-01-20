function [ v_skew ] = skew( v )
% Compute the 3x3 skew symmetric cross product operator of a vector

v_skew  = [ 0       -v(3)   v(2);
            v(3)    0       -v(1);
            -v(2)   v(1)    0 ];

end