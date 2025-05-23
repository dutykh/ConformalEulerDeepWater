%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, Khalifa University           %%%
%%% of Science and Technology, Abu Dhabi, UAE          %%%
%%% Date:   2025-05-23                                 %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

%> Description: Compute nonlinear right-hand side for conformal Euler evolution.
%> Inputs:
%>    w  - Nx2 matrix of [eta_hat, phi_hat] in Fourier space.
%>    dt - time increment for linear rotation step.
%> Globals:
%>    N, sk, k, ka (grid size, spectral operators).
%> Outputs:
%>    rhs - Nx2 matrix of right-hand sides in Fourier space.

function rhs = RHS (w, dt)

    global N sk k ka

    % Validate inputs
    if ~isnumeric(w) || size(w,2)~=2, error('RHS:w','w must be an Nx2 numeric array'); end
    if ~isnumeric(dt), error('RHS:dt','dt must be numeric'); end

    % declaration of the result and memory preallocation
    rh = zeros(N,2);    % nonlinear part of the RHS
    
    % v = turn(w, -dt);   % we integrate back the linear terms
    v = w;
    
    gam_hat = v(:,1);   % Fourier transform of the free surface elevation
    phi_hat = v(:,2);   % Fourier transform of the velocity potential
    
    gam_xi_hat = 1i*k.*gam_hat;         % compute Fourier transform of the derivative of \eta == \gamma
    gam_xi = real(ifft(gam_xi_hat));    % and in the real space
    
    gamH_hat = -ka.*gam_hat;            % Fourier transform of the Hilbert transform of the derivative of \gamma
    gamH = real(ifft(gamH_hat));        % the same in the real space
    
    chi_xi = 1 - gamH;                  % dX/d\xi expressed in terms of H[\gam_xi]
    J = chi_xi.^2 + gam_xi.^2;          % Jacobian of the transform
    J(J<eps) = eps;                     % avoid division by zero
    
    phi_hat_xi = 1i*k.*phi_hat;         % differentiate the potential in the Fourier space
    phi_xi = real(ifft(phi_hat_xi));    % and come back to the real space
    
    psi_xi_hat = -ka.*phi_hat;          % reconstruct the stream function derivative \psi_\xi with the same operator H[d_\xi]
    psi_xi = real(ifft(psi_xi_hat));    % and in the real space
    
    % Here we compute some nonlinear terms appearing at the RHS:
    psi_xiJ = psi_xi./J;
    psi_xiJ_hat = fft(psi_xiJ);
    psi_xiJ = real(ifft(psi_xiJ_hat));
    HpsiJ_hat = 1i*sk.*psi_xiJ_hat;
    HpsiJ = real(ifft(HpsiJ_hat));
    
    % Finally we can assemble the right-hand sides of evolution equations:
    rh(:,1) = fft(gam_xi.*HpsiJ + chi_xi.*psi_xiJ.*(J-1) + psi_xi.*gamH);
    rh(:,2) = fft(0.5*(psi_xi.^2-phi_xi.^2)./J + phi_xi.*HpsiJ);
    
    rhs = turn(rh, dt);  % we rotate the RHS and return the result

end % RHS ()