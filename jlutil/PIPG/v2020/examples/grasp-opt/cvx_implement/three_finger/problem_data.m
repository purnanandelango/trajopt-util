% 01/26/21
% grasp optimization for a 3D block with a three-finger manipulator

nx = 6;                       % no. of states
nu = 9;                       % no. of inputs
N = 30;                       % no. of discretization points (horizon length considered at each MPC instant)
Nbar = 60;                    % horizon length for full MPC simulation (Nbar-N MPC instants)  

a = 0.01;                     % half the width of block (m)
m = 0.2;                      % mass of the block (kg)
ge = 9.806;                   % accl. due to gravity at sea-level (m s^-2)
% ge = 0.0;
mu1 = 4.5;                      % friction coefficient
mu2 = 4.5;
mu3 = 4.5;
              
Del = 0.4;                    % sampling time (s)
tvec = 0:Del:(Del*(N-1));     % time vector (s)

% system matrices
Ad = [eye(3),   Del*eye(3);
      zeros(3), eye(3)];

Bd = [repmat(Del*Del*0.5*eye(3)/m,[1,3]);
      repmat(Del*eye(3)/m,[1,3])];
  
gd = [0;0;-0.5*Del*Del*ge;0;0;-Del*ge];

% reference trajectory
ra = 0.2;                                     % radius of arc traversed by the block in the reference solution
tau = @(t) ra*(2*t/(N-1)-1);                  % convenient function fo reference trajectory definition (t = 0,...,N-1)
tau1_z = @(t) sqrt(ra*ra-tau(t)^2);
tau2_z = @(t) 3*(ra-tau(t)+sqrt(ra*ra-tau(t)^2));
ybar = zeros(nx*Nbar,1);
y = zeros(nx*N,1);
for j=1:N
%     y((j-1)*nx+1:j*nx) = [tau(j-1);2*tau(j-1);0.5*tau1_z(j-1);0;0;0];
    y((j-1)*nx+1:j*nx) = [tau(j-1);tau(j-1);tau2_z(j-1);0;0;0];
end
yy = reshape(y,[nx,N]);
ybar(1:nx*N) = y(1:nx*N);
ybar(nx*N+1:end) = repmat(y(nx*(N-1)+1:nx*N),[Nbar-N,1]);
yybar = reshape(ybar,[nx,Nbar]);

% initial condition
x0 = y(1:nx);

% position vector of contact points in block frame
s1 = [-a;-tand(10)*a;tand(35)*a];
s2 = [-a;0;-tand(25)*a];
s3 = [a;-tand(10)*a;-tand(30)*a];
% s1 = [-a;0;0];
% s2 = [-a;0;0];
% s3 = [a;0;0];

% matrix in linear equality constraints due to rotational equlibrium
Gam = [skew(s1), skew(s2), skew(s3)];
   
Vmax = 100;                    % max l2-norm of velocity   
F1max = 100;                   % max l2-norm of the grasping force        
F2max = 100;                  % max l2-norm of the grasping force         
F3max = 100;                  % max l2-norm of the grasping force        

% Vmax = 2;                   % max l2-norm of velocity   
% F1max = 2;                  % max l2-norm of the grasping force        
% F2max = 2;                  % max l2-norm of the grasping force         
% F3max = 2;                  % max l2-norm of the grasping force        

% diagnostics
% figure
% plot(yy(1,:),yy(3,:),'-k');
% axis equal

q_wt = 1.0;
qf_wt = 10.0;
r_wt = 1.0;

function A = skew(a)
A = [0,-a(3),a(2);
     a(3),0,-a(1);
    -a(2),a(1),0];
end