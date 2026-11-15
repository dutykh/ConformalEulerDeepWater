% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Stand-in for scripts/Plot.m, used only by ../conformal_readme_data.m. The script run_Euler.m
% calls Plot(eta, t) at every output time (about sixteen times per carrier period); this version
% draws nothing and records, in global variables (run_Euler.m clears the workspace when it starts):
%   readme_rec    one row per output time: t, the highest elevation and its abscissa x, and the
%                 elevation at the wave gauges at 2.25 m, 6.25 m and 9 m from the wavemaker;
%   readme_snap   the free surface (x, eta) at t = 0 and at the output times closest to the
%                 passage of the centre of the group at these gauges;
%   readme_max    for every x of a uniform grid, the highest elevation met there so far, which
%                 a wave gauge at x would record as its highest crest.
% The elevation, known on the conformal grid xi, is interpolated onto the uniform grid of the
% physical abscissa x = xi - H[eta]. A crest passes a given x between two output times: there
% the elevation is the cubic Hermite interpolant of its values and of its time derivatives at
% the two output times, the derivative at fixed x being d eta/dt = -psi_xi/x_xi (the kinematic
% condition in conformal variables), and its maximum is taken on 16 sub-steps.

function Plot (eta, t)

    global readme_rec readme_snap readme_max readme_last ka xi N
    dtw = evalin('base', 'dtw'); cg = evalin('base', 'cg'); v = evalin('base', 'v');
    gauges = [2.25 6.25 9.0];                 % gauges of Figure 1 of the article (m)

    eta = eta(:);
    x = Conf2Real(eta);
    L = N*(xi(2) - xi(1));                    % length of the periodic domain
    m = 64;                                   % periodic extension for the interpolation
    xe = [x(end-m+1:end) - L; x; x(1:m) + L];
    ee = [eta(end-m+1:end); eta; eta(1:m)];
    eu = interp1(xe, ee, xi, 'spline');       % elevation on the uniform grid of x
    rate = -real(ifft(-ka.*v(:,2)))./(1 + real(ifft(ka.*v(:,1))));
    du = interp1(xe, [rate(end-m+1:end); rate; rate(1:m)], xi, 'spline');   % d eta/dt at fixed x
    [top, jc] = max(eta);
    readme_rec(end+1, :) = [t top x(jc) interp1(xe, ee, gauges, 'spline')];

    % highest elevation at every x of the uniform grid
    if isempty(readme_max)
        readme_max = eu;
    else
        t1 = readme_last(1, 3); h = t - t1;
        e1 = readme_last(:, 1); d1 = h*readme_last(:, 2); d2 = h*du;
        a = 2*e1 - 2*eu + d1 + d2;
        b = -3*e1 + 3*eu - 2*d1 - d2;
        for s = (1:16)/16
            readme_max = max(readme_max, ((a*s + b)*s + d1)*s + e1);
        end
    end
    readme_last = [eu du t*ones(N, 1)];

    % snapshots: t = 0, and the centre of the group at each gauge (t = X/cg)
    instants = [0 gauges/cg];
    if any(abs(t - instants) < dtw/2)
        if isempty(readme_snap), readme_snap = {}; end
        readme_snap{end+1} = [t*ones(N, 1) x eta];
    end

end % Plot ()
