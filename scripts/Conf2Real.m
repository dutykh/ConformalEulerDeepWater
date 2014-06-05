%%% -------------------------------------------------- %%%
%%% Transformation to the physical horiz. coordinate   %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function x = Conf2Real (eta)
	global sk xi
	
	x = xi + real(ifft(1i*sk.*fft(eta)));
end % Conf2Real ()