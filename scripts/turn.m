%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, Khalifa University           %%%
%%% of Science and Technology, Abu Dhabi, UAE          %%%
%%% Date:   2025-05-23                                 %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

%%% Description: Advance solution through one linear rotation step with spectral filter.
%%% Inputs:
%%%    w  - Nx2 numeric array of Fourier coefficients [eta_hat, phi_hat].
%%%    dt - time increment for rotation (scalar).
%%% Globals:
%%%    kf, omg, gom, omega, N
%%% Outputs:
%%%    v  - Nx2 rotated and filtered Fourier solution.

function v = turn (w, dt)

    global kf omg gom omega N
    
    % Validate inputs
    if ~isnumeric(w) || size(w,2)~=2, error('turn:w','w must be an Nx2 numeric array'); end
    if ~isnumeric(dt), error('turn:dt','dt must be numeric'); end
    
    v = zeros(N, 2);
    
    cosw = cos(omega*dt);
    sinw = sin(omega*dt);
    
    v(:,1) = kf.*(cosw.*w(:,1) - omg.*sinw.*w(:,2));
    v(:,2) = kf.*(gom.*sinw.*w(:,1) + cosw.*w(:,2));
    
end % turn ()