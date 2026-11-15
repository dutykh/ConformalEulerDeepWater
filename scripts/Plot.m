%%% -------------------------------------------------- %%%
%%% This function plots the solution in physical space %%%
%%% along with  its Fourier spectrum                   %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, Khalifa University           %%%
%%% of Science and Technology, Abu Dhabi, UAE          %%%
%%% Date:   2025-05-23                                 %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Input parameters:                                  %%%
%%%    eta : free surface elevation                    %%%
%%%      t : current instant of time (for the title)   %%%
%%% -------------------------------------------------- %%%

%%% Description: Plot the free-surface elevation in real space and its Fourier spectrum.
%%% Inputs:
%%%    eta - column vector of free-surface elevation.
%%%    t   - current time (scalar) for title annotation.
%%% Globals:
%%%    et0, k, k0, l, N (initial amplitude, wavenumbers, domain half-length, grid size).

function Plot (eta, t)

    % Clear figure and validate inputs
    clf;
    if ~isvector(eta), error('Plot:eta','eta must be a vector.'); end
    if ~isnumeric(t), error('Plot:t','t must be numeric.'); end
    eta = eta(:);
    global et0 k k0 l N

    %%% We plot the initial condition:
    h = subplot(2,1,1);
    p = get(h, 'Position'); p(4) = p(4) - 0.03; p(2) = p(2) + 0.02;
    set(h, 'Position', p);

    x = Conf2Real(eta);

    plot(x, eta, 'b-');
    xlim([-l l]); ylim([-3.5*et0 3.5*et0])
    xlabel('$x$', 'interpreter', 'latex');
    ylabel('$\eta(x,t)$', 'interpreter', 'latex');
    title(['Peregrine breather at t = ', num2str(t, '%3.2f')], 'FontAngle', 'italic');

    h = subplot(2,1,2);
    p = get(h, 'Position'); p(4) = p(4) - 0.03; p(2) = p(2) + 0.02;
    set(h, 'Position', p);

    Power = abs(fft(eta)/N).^2;
    kk = k(2:N/2)/k0; Power = Power(2:N/2);  % positive wavenumbers, in units of k0
    keep = Power > 0;                        % the filter sets the last modes exactly to zero
    loglog(kk(keep), Power(keep), 'k.');
    xlabel('$k/k_0$', 'interpreter', 'latex');
    ylabel('$|\hat\eta|^2(k)$', 'interpreter', 'latex');
    xlim([kk(1) kk(end)]); ylim([1e-30 1]);
    set(gca, 'YTick', [1e-30 1e-20 1e-10 1])
    title('Fourier spectrum', 'FontAngle', 'italic');

    set(gcf, 'Color', [0.95, 0.95, 0.95]);
    drawnow
    
end % Plot ()