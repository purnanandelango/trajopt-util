clearvars
close all
clc

N = 100;
x1 = 1;
xN = 10;
tau = linspace(0,1,N);

y_lin   = grid.ends2interp(x1,xN,tau,'poly',1);
y_quad  = grid.ends2interp(x1,xN,tau,'poly',2);
y_cub   = grid.ends2interp(x1,xN,tau,'poly',3);
y_quar  = grid.ends2interp(x1,xN,tau,'poly',4);
y_exp   = grid.ends2interp(x1,xN,tau,'exp',10);
y_log   = grid.ends2interp(x1,xN,tau,'log',1);

plot(tau,y_lin,'-','DisplayName','$m=1$');
hold on 
plot(tau,y_quad,'-','DisplayName','$m=2$');
plot(tau,y_cub,'-','DisplayName','$m=3$');
plot(tau,y_quar,'-','DisplayName','$m=4$');
plot(tau,y_exp,'-','DisplayName','exp');
plot(tau,y_log,'-','DisplayName','log');
legend('Location','best');