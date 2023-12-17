function x = callpiqp(P, c, A, b, G, h, x_lb, x_ub)
    solver = piqp('dense');
    solver.update_settings('verbose', false, 'compute_timings', false);
    solver.setup(P, c, A, b, G, h, x_lb, x_ub);
    result = solver.solve();
    x = result.x;
end