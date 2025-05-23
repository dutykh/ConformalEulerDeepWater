%%% -------------------------------------------------- %%%
%%% Transformation to the physical horiz. coordinate   %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, Khalifa University           %%%
%%% of Science and Technology, Abu Dhabi, UAE          %%%
%%% Date:   2025-05-23                                 %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

%%% Description: Convert surface elevation in conformal space to real horizontal coordinate.
%%% Inputs:
%%%    eta - column vector of surface elevation in conformal coordinates.
%%% Outputs:
%%%    x   - column vector of physical horizontal positions.
%%% Globals:
%%%    sk, xi (spectral variables and grid coordinates).

function x = Conf2Real (eta)
	
	% Input validation
	if ~isvector(eta)
		error('Conf2Real:eta','Input eta must be a vector.');
	end
	eta = eta(:);

	global sk xi
	
	x = xi - real(ifft(1i*sk.*fft(eta)));
	
end % Conf2Real ()