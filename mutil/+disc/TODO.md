# TODO

 - Support hybrid discretization:
    - `Impulse` or `FBP` for actual control input and `ZOH` for dilation factors.
 - Add support for `FBP` and `Impulse` discretizations when parameters are decision variables.
 - Parallelize `compute_foh_v...` and `compute_zoh_v...` via `parfor`.
 - Remove `time_grid_via_b` and `time_of_maneuver_via_b`. 