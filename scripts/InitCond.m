%%% -------------------------------------------------- %%%
%%% We set-up the Peregrine breather as the initial    %%%
%%% condition. Please, note that the breather is       %%%
%%% windowed to enforce the periodic boundary          %%%
%%% conditions. The windowing is done far from the     %%%
%%% focusing point in order to annihilate its effect.  %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Input parameters:                                  %%%
%%%		t0 : initial moment of time                    %%%
%%%		 x :                                           %%%
%%% Output parameters:                                 %%%
%%%	  eta0 : initial free surface elevation at t = t0  %%%
%%%	  phi0 : velocity potential at the free surface    %%%
%%% -------------------------------------------------- %%%

function [eta0, phi0] = InitCond (t0, x0)

	global cg0 e0 g k0 k ka om0 sk tol win

	% a simple iterative procedure to obtain the initial condition in
	% the conformal space:
	err = inf;
	eta0 = FreeSurf(x0);
	while (err > tol)
		x = Conf2Real(eta0);
		eta = FreeSurf(x);
		err = norm(eta - eta0, inf);
		eta0 = eta;
	end % while ()

	phi0_hat = fft(eta0).*sk.*sqrt(g./ka); phi0_hat(1) = 0.0;
	phi0 = real(ifft(phi0_hat));

	% Function which returns the shape of the free surface as a function of x:
	function eta0 = FreeSurf (x)
		T = 0.25*e0*e0*om0*t0;
		X = sqrt(2)*e0*k0*(x - cg0*t0);
	
		Amp = ((4 + 16i*T)./(1 + 4*X.^2 + 16*T^2) - 1)*exp(2i*T);
		eta0 = 0.5*real(e0/k0*Amp.*exp(-1i*(k0*x - om0*t0))).*win;
		eta0 = eta0 - mean(eta0);	% we remove the mean to be consistent with the definition of \eta
	end % FreeSurf ()

end % InitCond ()