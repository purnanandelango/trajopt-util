# PIPG

Discrete-time optimal control problem

$$
\begin{aligned}\underset{u_t}{\operatorname{minimize}}~~~&~\sum_{t=1}^{N} x_t^\top Q_t x_t + q_t^\top x_t + u_t^\top R_t u_t + r_t^\top & & \\
\operatorname{subject~to}~~&~x_{t+1} = A_t x_t + B^-_t u_t + B_{t+1}^+u_{t+1}+g_t, & & t = 1,\ldots,N-1,\\
  &~x_t\in\mathbb{D}^x_t, & & t = 1,\ldots,N,\\
  &~u_t\in\mathbb{D}^u_t, & & t = 1,\ldots,N.\\
\end{aligned}
$$ (1)

To track known state reference $x_t^{\text{ref}}$ and/or a control reference $u^{\text{ref}}_t$, choose $q_t = -2x_t^{\text{ref}}$ and $r_t = -2u^{\text{ref}}_t$. The boundary conditions on states and control   