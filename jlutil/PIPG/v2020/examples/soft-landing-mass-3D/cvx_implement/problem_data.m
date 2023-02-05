% 03/28/21
% simplified 3DoF pdg problem without the pointing constraints


nx = 7;                                                 % state dim
nu = 4;                                                 % input dim    

N = 25;                                                 % horizon length
Del = 80/(N-1);                                         % sampling time (s)
tvec = 0:Del:(N-1)*Del;                                 % time grid (s)

mwet = 2000;                                            % dry weight (kg)
mdry = 150;                                             % wet weight (kg)

gvec = [0;0;-3.71];                                     % accl. due to gravity (m/s^2)

Tmax = 24000;                                           % maximum thrust (N)
rho1 = 0.2*Tmax;                                        % upper bound on thrust (N)
rho2 = 0.8*Tmax;                                        % lower bound on thrust (N)    
Vmax = 300;                                             % upper bound on speen (m/s)
tangamgs = tand( 85 );                                  % tangent of maximum glide-slope angle 

alph = 0.0005;                                          % mas depletion rate (s^-1/m)

x0 = [0;            %4500;
     4000;          %-330;
     1500;          %2400;
      0;            %-40;
      100;          %10;
     -75;           %-10;
      log(mwet)];
rf = zeros(3,1);
vf = zeros(3,1);
xf = [rf;vf;log(mdry)];

A = [zeros(3),eye(3),zeros(3,1);zeros(4,7)];
B = [zeros(3,4);eye(3),zeros(3,1);zeros(1,3),-alph];
g = [zeros(3,1);gvec;0];

% ZOH
Ad = expm(Del*A);                               
Bd = (Del*eye(nx)+0.5*Del*Del*A)*B;             
gd = (Del*eye(nx)+0.5*Del*Del*A)*g;             

nhat = [0;0;1];
E = [1,0,0;
     0,1,0];
 
z1 = log(mwet - alph*rho1*tvec);
z0 = log(mwet - alph*rho2*tvec);
mu1 = rho1*exp(-z0);
mu2 = rho2*exp(-z0);

Q = diag([1e-3,1e-3,1e-3,1e-2,1e-2,1e-2,1e-1]);
R = diag([1e-1,1e-1,1e-1,1e-1]);

yref(nx,N) = 0;
for j = 0:N-1
    yref(:,j+1) = (1-j/(N-1)) * x0 + (j/(N-1)) * xf;
end