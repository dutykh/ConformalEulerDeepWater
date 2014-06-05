%%% -------------------------------------------------- %%%
%%% This function plots the solution in physical space %%%
%%% along with  its Fourier spectrum                   %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Input parameters:                                  %%%
%%%    eta : free surface elevation                    %%%
%%%      t : current instant of time (for the title)   %%%
%%% -------------------------------------------------- %%%

function Plot (eta, t)

	global a0 l N

	%%% We plot the initial condition:
	h = subplot(2,1,1);
	p = get(h, 'pos'); p(4) = p(4) - 0.03; p(2) = p(2) + 0.02;
	set(h, 'pos', p);

	x = Conf2Real(eta);

	plot(x, eta, 'b-');
	xlim([-l l]); ylim([-1.5*a0 1.5*a0])
	xlabel('$x$', 'interpreter', 'latex');
	ylabel('$\eta(x,t)$', 'interpreter', 'latex');
	title(['Peregrine breather at t = ', num2str(t, '%3.2f')], 'FontAngle', 'oblique');

	h = subplot(2,1,2);
	p = get(h, 'pos'); p(4) = p(4) - 0.03; p(2) = p(2) + 0.02;
	set(h, 'pos', p);

	Power = abs(fft(eta)/N).^2;
	loglog(1:N/2, Power(1:N/2), 'k-');
	xlabel('$\log\, (k)$', 'interpreter', 'latex');
	ylabel('$|\hat\eta|^2(k)$', 'interpreter', 'latex');
	xlim([0 N/2+1]); ylim([1e-30 1e2]);
	set(gca, 'YTick', [1e-30 1e-20 1e-10 1])
	title('Fourier spectrum', 'FontAngle', 'oblique');

	set(gcf, 'Color', [0.95, 0.95, 0.95]);
	drawnow
	
end % Plot ()