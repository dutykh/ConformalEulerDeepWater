% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Stand-in for scripts/Plot.m, used only by ../conformal_readme_data.m. The script run_Euler.m
% calls Plot(eta, t) at every output time; this version draws nothing: it records the highest
% elevation of the free surface and the current time step dt of run_Euler.m, keeps the initial
% profile and the profile with the highest crest met so far (in physical coordinates), and stops
% the run with the error 'readme:stop' at t = 26 s. The records live in global variables, created
% here, since run_Euler.m clears the workspace when it starts.

function Plot (eta, t)

    global readme_t readme_max readme_dt readme_first readme_peak
    t_stop = 26.0;

    eta = eta(:);
    x = Conf2Real(eta);
    readme_t(end+1, 1) = t;
    readme_max(end+1, 1) = max(eta);
    readme_dt(end+1, 1) = evalin('base', 'dt');
    if isempty(readme_first)
        readme_first = [x eta];
    end
    if isempty(readme_peak) || max(eta) > max(readme_peak(:, 2))
        readme_peak = [x eta repmat(t, numel(eta), 1)];
    end
    if t >= t_stop
        error('readme:stop', 'The run for the figure of the README stops at t = %g s.', t);
    end

end % Plot ()
