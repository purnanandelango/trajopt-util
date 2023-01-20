# Nonlinear system based on van der Pol oscillator

$$
\begin{aligned}
\dot{x} ={} & f = \left[ \begin{array}{c} p_1x_2u_1 \\ p_1u_2x_2 - p_1u_2x_1^2x_2 - p_1p_2x_1 \end{array}\right]\\
\frac{\partial f}{\partial x} ={} & \left[\begin{array}{cc} 0 & p_1u_1 \\ -p_1p_2 - 2p_1u_2x_1x_2 & p_1u_2(1-x_1^2) \end{array}\right]\\
\frac{\partial f}{\partial u} ={} & \left[\begin{array}{cc} p_1x_2 & 0 \\ 0 & p_1x_2 - p_1 x_1^2x_2 \end{array}\right]\\
\frac{\partial f}{\partial p} ={} & \left[\begin{array}{cc} x_2u_1 & 0 \\ u_2x_2 - u_2x_1^2x_2 - p_2x_1 & -p_1x_2 \end{array}\right]
\end{aligned}
$$