% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Data of the figure of the README: the run of scripts/run_Euler.m as it stands, the Peregrine
% breather of the experiment of Shemer & Alperovich (2013) from the wavemaker to the distance
% X0 = 9 m where the Peregrine solution has its maximum. The script is run from the folder
% recording/, whose Plot.m takes the place of scripts/Plot.m and records, at every output time,
% the highest elevation and the elevation at the gauges at 2.25 m, 6.25 m and 9 m from the
% wavemaker, the highest elevation met at every x, and the free surface when the centre of the
% group leaves the wavemaker and passes these gauges. This script writes, all in SI units,
%   record.csv       t, max eta, its x, eta at the three gauges;
%   maxima.csv       x, highest elevation met at x during the run;
%   snapshots.csv    t, x, eta, for the four instants one after the other;
%   invariants.csv   t, mass, momentum, energy (the lists Time, Mass, Moment, Energy of run_Euler.m).
%
% Run it from this folder with GNU Octave or MATLAB (about 80 s):
%   octave --no-gui conformal_readme_data.m        or, in MATLAB,   run conformal_readme_data
% plot_readme_figure.py then draws conformal-breather.png from the four files.

set(0, 'defaultfigurevisible', 'off');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'scripts'));
cd(fullfile(fileparts(mfilename('fullpath')), 'recording'));
run_Euler

% run_Euler.m clears the workspace when it starts: the records come back through the globals,
% and its own variables (Time, Mass, Moment, Energy, ...) are the ones left by the run.
global readme_rec readme_snap readme_max
here = fileparts(mfilename('fullpath'));
cd(here);
dlmwrite(fullfile(here, 'record.csv'), readme_rec, 'precision', '%.16e');
dlmwrite(fullfile(here, 'maxima.csv'), [xi readme_max], 'precision', '%.16e');
dlmwrite(fullfile(here, 'snapshots.csv'), vertcat(readme_snap{:}), 'precision', '%.16e');
dlmwrite(fullfile(here, 'invariants.csv'), [Time Mass Moment Energy], 'precision', '%.16e');
[top, n] = max(readme_rec(:, 2));
fprintf('%d output times up to t = %.3f s; highest crest %.6f m at t = %.3f s, x = %.3f m\n', ...
        size(readme_rec, 1), readme_rec(end, 1), top, readme_rec(n, 1), readme_rec(n, 3));
fprintf('largest change over the run: mass %.2e m^2, momentum %.2e (relative), energy %.2e (relative)\n', ...
        max(abs(Mass - Mass(1))), max(abs(Moment - Moment(1)))/abs(Moment(1)), ...
        max(abs(Energy - Energy(1)))/Energy(1));
