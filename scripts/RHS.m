%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

% Right-hand side function
function rhs = RHS (w, dt)

    global N sk k ka

    % declaration of the result and memory preallocation
    rh = zeros(N,2);    % nonlinear part of the RHS
    
    v = turn(w, -dt);   % we integrate back the linear terms
    
    gam_hat = v(:,1);   % Fourier transform of the free surface elevation
    phi_hat = v(:,2);   % Fourier transform of the velocity potential
    
    gam_xi_hat = 1i*k.*gam_hat;         % compute Fourier transform of the derivative of \eta == \gamma
    gam_xi = real(ifft(gam_xi_hat));    % and in the real space
    
    gamH_hat = -ka.*gam_hat;            % Fourier transform of the Hilbert transform of the derivative of \gamma
    gamH = real(ifft(gamH_hat));        % the same in the real space
    
    chi_xi = 1 - gamH;                  % dX/d\xi expressed in terms of H[\gam_xi]
    J = chi_xi.^2 + gam_xi.^2;          % Jacobian of the transform
    
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