% Rocket soft landing
% Thrust lower bound is approximated using halfspace
% Mass dynamics is excluded
% Thrust pointing constraint is convex i.e. \theta \in (0,\pi/2)
% Glide slope constraint is included

% problem parameters
nx = 6;
nu = 3;

m = 10.0;                      % mass
tfin = 15.0;                   % duration of flight
N = 20;                        % no. of discretization points (horizon length considered at each MPC instant)
Nbar = 40;                     % horizon length for full MPC simulation (Nbar-N MPC instants)

Del = tfin/(N-1);              % sampling time
tvec = 0:Del:tfin;             % time grid
ge = 9.806;                    % magnitude of acceleration due to gravity

% LQR weighting matrices
Q = eye(6);
Qf = eye(6);
R = 0.1*eye(3);

% continuous time
Ac = [zeros(3,3) eye(3);
      zeros(3,3) zeros(3,3)];
Bc = [zeros(3,3);
      eye(3)/m];
gc = [zeros(5,1);-ge];         % gravity vector expressed in East(x)-North(y)-Up(z) frame

% ZOH discretization
Ad = expm(Del*Ac);
Bd = Del*(eye(6)+0.5*Del*Ac)*Bc;
gd = Del*(eye(6)+0.5*Del*Ac)*gc;

% boundary conditions
x0 = [4.4;13.5;23;-3;2;-5];
xf = zeros(6,1);

% constraints
gamgs = 50*pi/180;            % maximum glide-slope cone angle
thetmax = 7*pi/180;           % maximum gimbal angle
Vmax = 8;                     % maximum velocity
umax = 1.2*m*ge;              % maximum thrust magnitude
umin = 0.6*m*ge;              % thrust lower bound approximation with halfspace
cvec = [0,0,1]*tan(gamgs);
dvec = [0,0,1]*tan(thetmax);

% straight-line interpolation between x0 and xf
ybar = zeros(6,Nbar);

% for i = 1:Nbar
%     ybar(:,i) = (1-(i-1)/(Nbar-1))*x0 + (i-1)*xf/(Nbar-1);
% end
% y = ybar(:,1:N);

y = zeros(nx,N);
for i = 1:N
    y(:,i) = (1-(i-1)/(N-1))*x0 + (i-1)*xf/(N-1);
end
ybar(:,1:N) = y;
ybar(:,N+1:end) = repmat(y(:,N),[1,Nbar-N]);