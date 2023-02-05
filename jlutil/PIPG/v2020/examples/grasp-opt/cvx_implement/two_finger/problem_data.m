% 02/01/21
% grasp optimization for a 3D block with a two-finger manipulator

nx = 6;                       % no. of states
nu = 6;                       % no. of inputs
N = 20;                       % no. of discretization points

a = 0.1;                      % half the width of block (m)
m = 2;                        % mass of the block (kg)
% ge = 9.806;                   % accl. due to gravity at sea-level (m s^-2)
ge = 0;
mu1 = tand(55);               % friction coefficient
mu2 = tand(55);

Del = 0.9;                    % sampling time (s)
tvec = 0:Del:(Del*(N-1));     % time vector (s)

% system matrices
Ad = [eye(3),   Del*eye(3);
      zeros(3), eye(3)];

Bd = [repmat(Del*Del*0.5*eye(3)/m,[1,2]);
      repmat(Del*eye(3)/m,[1,2])];
  
gd = [0;0;-0.5*Del*Del*ge;0;0;-Del*ge];

% reference trajectory
ra = 2;                       % radius of arc traversed by the block in the reference solution
tau = @(t) ra*(2*t/(N-1)-1);  % convenient function fo reference trajectory definition (t = 0,...,N-1)
y = zeros(nx*N,1);
for j=1:N
%     y((j-1)*nx+1:j*nx) = [tau(j-1);2*tau(j-1);sqrt(ra*ra-tau(j-1)^2);0;0;0];
    y((j-1)*nx+1:j*nx) = zeros(6,1);
end
yy = reshape(y,[nx,N]);

% position vector of contact points in block frame
s1 = [-a;a/2;a/2];
s2 = [a;a/2;a/2];

% matrix in linear equality constraints due to rotational equlibrium
Gam = [skew(s1), skew(s2)];
   
Vmax = 100;                     % max l2-norm of velocity   
F1max = 500;                  % max l2-norm of the grasping force        
F2max = 500;                  % max l2-norm of the grasping force         

% diagnostics
% figure
% plot(yy(1,:),yy(3,:),'-k');
% axis equal

% projection onto ker Gam
proj_ker_Gam = eye(nu) - (Gam')*inv(Gam*(Gam'))*Gam;

function A = skew(a)
A = [0,-a(3),a(2);
     a(3),0,-a(1);
    -a(2),a(1),0];
end