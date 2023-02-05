% solve 3D grasp optimization problem

clear all
clc
close all

problem_data;

cvx_solver mosek

cvx_begin
    variables x(nx*N) u(nu*(N-1))
    minimize  q_wt*sum_square(x(1:(N-1)*nx)-y(1:(N-1)*nx)) + qf_wt*sum_square(x((N-1)*nx+1:N*nx)-y((N-1)*nx+1:N*nx)) + r_wt*sum_square(u)
    subject to
        for j=1:N-1
            x(j*nx+1:(j+1)*nx) == Ad*x((j-1)*nx+1:j*nx) + Bd*u((j-1)*nu+1:j*nu) + gd; % dynamics constraint
            Gam*u((j-1)*nu+1:j*nu) == 0;                                              % rotational equilibrium
            norm(u((j-1)*nu+2:(j-1)*nu+3)) <= mu1*u((j-1)*nu+1);                      % friction cone
            norm(u((j-1)*nu+5:(j-1)*nu+6)) <= mu2*u((j-1)*nu+4);                      % friction cone            
            norm(u((j-1)*nu+8:(j-1)*nu+9)) <= -mu3*u((j-1)*nu+7);                     % friction cone            
            norm(u((j-1)*nu+1:(j-1)*nu+3)) <= F1max;                                  % grasp force magnitude
            norm(u((j-1)*nu+4:(j-1)*nu+6)) <= F2max;                                  % grasp force magnitude
            norm(u((j-1)*nu+7:(j-1)*nu+9)) <= F3max;                                  % grasp force magnitude
            norm(x(j*nx+4:j*nx+6)) <= Vmax;                                           % velocity magnitude                  
        end
        x(1:nx) == x0; 
cvx_end

% constraint activity
Dyn(1,N-1) = 0;
RotEq(1,N-1) = 0;
Fcone(3,N-1) = 0;
Fmag(3,N-1) = 0;
Vmag(1,N-1) = 0;
for j=1:N-1
    Dyn(j) = norm(x(j*nx+1:(j+1)*nx) - Ad*x((j-1)*nx+1:j*nx) - Bd*u((j-1)*nu+1:j*nu) - gd,Inf); % dynamics constraint
    RotEq(j) = norm(Gam*u((j-1)*nu+1:j*nu),Inf);                                                % rotational equilibrium
    Fcone(1,j) = norm(u((j-1)*nu+2:(j-1)*nu+3)) - mu1*u((j-1)*nu+1);                             % friction cone
    Fcone(2,j) = norm(u((j-1)*nu+5:(j-1)*nu+6)) - mu2*u((j-1)*nu+4);                             % friction cone            
    Fcone(3,j) = norm(u((j-1)*nu+8:(j-1)*nu+9)) + mu3*u((j-1)*nu+7);                             % friction cone            
    Fmag(1,j) = norm(u((j-1)*nu+1:(j-1)*nu+3)) - F1max;                                         % grasp force magnitude
    Fmag(2,j) = norm(u((j-1)*nu+4:(j-1)*nu+6)) - F2max;                                         % grasp force magnitude
    Fmag(3,j) = norm(u((j-1)*nu+7:(j-1)*nu+9)) - F3max;                                         % grasp force magnitude
    Vmag(j) = norm(x(j*nx+4:j*nx+6)) - Vmax;                                                    % velocity magnitude                  
end
fprintf("\nConstraint Activity:\n");
fprintf("dynamics constr: %.3f\n",max(abs(Dyn)));
fprintf("rotational eqb.: %.3f\n",max(abs(RotEq)));
fprintf("friction cone 1: %.3f\n",max(Fcone(1,:)));
fprintf("friction cone 2: %.3f\n",max(Fcone(2,:)));
fprintf("friction cone 3: %.3f\n",max(Fcone(3,:)));
fprintf("friction mag. 1: %.3f\n",max(Fmag(1,:)));
fprintf("friction mag. 2: %.3f\n",max(Fmag(2,:)));
fprintf("friction mag. 3: %.3f\n",max(Fmag(3,:)));
fprintf("velocity mag.  : %.3f\n",max(Vmag));

x = reshape(x,[nx,N]);
u = reshape(u,[nu,N-1]);

plot_soln(x,u,yy,N,s1,s2,s3,a)