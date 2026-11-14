% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Data of the figure of the README: the focusing of the Peregrine breather computed by
% scripts/run_Euler.m, run as it stands over its first 26 seconds. The script is run from the
% folder recording/, whose Plot.m takes the place of scripts/Plot.m: at every output time
% (0.1 s) it records the highest elevation of the free surface and the time step, it keeps the
% initial profile and the profile with the highest crest, and it stops the run at t = 26 s. This
% script writes the records to crest.csv (t, max eta, dt), initial.csv (x, eta) and peak.csv
% (x, eta, t), all in SI units.
%
% Run it from this folder with GNU Octave or MATLAB (about two minutes):
%   octave --no-gui conformal_readme_data.m        or, in MATLAB,   run conformal_readme_data
% plot_readme_figure.py then draws conformal-breather.png from the three files.

set(0, 'defaultfigurevisible', 'off');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'scripts'));
cd(fullfile(fileparts(mfilename('fullpath')), 'recording'));
try
    run_Euler
catch err
    if ~strcmp(err.identifier, 'readme:stop'), rethrow(err); end
end

% run_Euler.m clears the workspace when it starts: the records come back through the globals.
global readme_t readme_max readme_dt readme_first readme_peak
here = fileparts(mfilename('fullpath'));
cd(here);
dlmwrite(fullfile(here, 'crest.csv'), [readme_t readme_max readme_dt], 'precision', '%.16e');
dlmwrite(fullfile(here, 'initial.csv'), readme_first, 'precision', '%.16e');
dlmwrite(fullfile(here, 'peak.csv'), readme_peak, 'precision', '%.16e');
[top, n] = max(readme_max);
fprintf('%d output times up to t = %.2f s; highest crest %.6f m at t = %.2f s (initially %.6f m)\n', ...
        numel(readme_t), readme_t(end), top, readme_t(n), readme_max(1));
fprintf('time step: %.2e s at the start, %.2e s at the end\n', readme_dt(1), readme_dt(end));
