# DDTO with varying dilation

## Observations

- Large $\Delta t$ causes large single shooting propagation error.
- With variable dilation, $x^i_k = x^j_k$ is not enough.
- Shrinking the constraint boundary for suboptimality is beneficial.
- Choice of $K^\star$ is crucial.
- ECOS and Gurobi are able to drive slack variables down to $10^{-8}$ while Mosek only manages $10^{-5}$,

## Scenario Ideas

 - Necessity of the constraint $u^j_t = u^i_t$.
 - Different types of suboptimality constraints and the consequence to SOCP: several small cones vs single large cones.
 - Vehicle slows down initially to gain deferrability.
 - Increase $K^\star$ to increase deferrability and reduce quality of trajectories?

