%%% -------------------------------------------------- %%%
%%% Numerical solution of full Euler equations deep w  %%%
%%%     using the conformal mapping technique          %%%
%%% Simulation of the Peregrine breather (NLS)         %%%
%%% -------------------------------------------------- %%%
%%% For the initialisation of the breather, please,    %%%
%%% refer to L. Shemer & L. Alperovich. 25, 051701,    %%%
%%% "Peregrine breather revisited", Phys. Fluids, 2013 %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Acknowledgements:                                  %%%
%%% * Prof. Didier Clamond, Univ.Nice Sophia Antipolis %%%
%%% * Dr. Bernard Ee, Tel Aviv University, Israel      %%%
%%% * Prof. Lev Shemer, Tel Aviv University, Israel    %%%
%%% -------------------------------------------------- %%%

% We clear the workspace
close all, clear all
format longE

% Declaration of global variables
global et0 gom k ka kf k0 l omg omega sk xi N

%%% Physical parameters of the breather:
g    = 9.81;					% gravity acceleration
T0   = 0.75;					% carrier wave period
et0  = 0.01;					% characteristic wave amplitude
x0   = -17.74;					% initial position (center) of the breather
om0  = 2*pi/T0;					% wave frequency
k0   = om0^2/g;					% wavenumber found from the dispersion relation
l0   = 2*pi/k0;					% wavelength
a0   = et0/sqrt(2);		    	% carrier wave amplitude
epsi = a0*k0;					% characteristic wave steepness
cg   = 0.5*sqrt(g/k0);			% wave group speed in deep water
Tmax = -30.30;					% time where max amplitude is reached 
								% (computed from the analytical solution)

%%% Numerical parameters:
N = 16384;               		% number of Fourier modes in discrete solution

%%% Definition of the computational domain:
l = 35.0*l0;            		% half-length of the domain in \xi space
dxi = 2*l/N;            		% distance between two \xi-points
xi = (1-N/2:N/2)'*dxi;  		% \xi-space discretization

%%% Discretization in the Fourier space:
% vector of wavenumbers:
k = [0:N/2 1-N/2:-1]'*pi/l;

% Window function for the free surface elevation:
win = zeros(N,1);
for j=1:length(xi), % now create window function which will be multiplied with the amplitude function
    if (abs(xi(j)) <= l) & (abs(xi(j)) > 0.5*l)
        win(j) = 0.0;
    elseif (abs(xi(j)) <= 0.5*l) & (abs(xi(j)) > 0.5*l - 2.0*l0)
        if xi(j) < 0,
            coef = 1/(2*l0);
        else 
        	coef = -1/(2*l0);
        end;
        win(j) = coef*(xi(j)) + 0.25*l/l0;
    else 
    	win(j) = 1.0;
    end
end

% Exponential de-aliasing treatment proposed by T. Hou et al.:
Lambd = 15;
kam = 0.60*max(abs(k));
kf = exp(-5.0*(abs(k)/kam).^Lambd); % dealiasing function

%%% Pseudo-differential operators:
sk = sign(k);   % needed for the Hilbert transform = -i*sign(k)
ka = abs(k);    % Hilbert transform of the derivative

% Needed for the integration of linear terms:
omega = sqrt(g*ka);
omg = omega/g;
gom = g./omega; gom(1) = 0.0;

%%% Initial condition specification in the conformal space:
xi1 = x0 + xi;					% shifted space variable
Xit = epsi*om0*(xi1/cg - Tmax); % scaled time variable
Xis  = epsi^2*k0*xi1;			% scaled space variable
% Peregrine's analytical solution:
A   = sqrt(2)*(1 - 4*(1-4*1i*Xis)./(1 + 4*Xit.^2 + 16*Xis.^2)).*exp(-2.0*1i*Xis);
Env = a0*A.*exp(1i*k0*xi1);		% wave envelope
eta0 = real(Env).*win;			% windowed free surface elevation

err = inf;
tol = 1e-14;
while (err > tol)
	chi0 = Conf2Real(eta0);

	xi1 = x0 + chi0;				% shifted space variable
	Xit = epsi*om0*(xi1/cg - Tmax); % scaled time variable
	Xis  = epsi^2*k0*xi1;			% scaled space variable
	% Peregrine's analytical solution:
	A   = sqrt(2)*(1 - 4*(1-4*1i*Xis)./(1 + 4*Xit.^2 + 16*Xis.^2)).*exp(-2.0*1i*Xis);
	Env = a0*A.*exp(1i*k0*xi1);		% wave envelope
	eta1 = real(Env).*win;			% windowed free surface elevation
	
	err = norm(eta1 - eta0, inf);
	eta0 = eta1;
end % while ()
eta0 = eta0 - mean(eta0);		% we remove the mean

phi0_hat = -1i*fft(eta0).*sk.*sqrt(g./ka);
phi0_hat(1) = 0.0;
phi0 = real(ifft(phi0_hat));

v = zeros(N,2);
v(:,1) = fft(eta0);
v(:,2) = fft(phi0);

% FFTW plan:
fftw('planner', 'hybrid');

%%% Time-stepping parameters:
t    = 0.0;      			% the discrete time variable
Tf   = 95.0*T0*sqrt(g);		% final simulation time
dtw  = 0.1;   				% write time of the results
tol  = 5e-10;  				% tolerance parameter for time stepping
dt   = 0.01;    			% initial guess of the time step
err1 = 1.0;   				% local error guess
tloc = 0.0;   				% local time variable
pf   = false;   			% plot flag

%%% Define the graphic window:
FigHandle = figure(1);
set(FigHandle, 'Position', [100, 100, 849, 495]);
% and plot the initial condition:
Plot(eta0, t);

% Parameters for the adaptive time-stepping (PI control):
% taken according to the book of Hairer (Part II)
al = 0.7/5; be = 0.4/5;

%%% In order to check the accuracy of computations we can follow
%%% the following quantities:
Mass = [];    % List of the masses at different moments of time
Moment = [];  % List of the momenta at different moments of time
Energy = [];  % List of the (total) energy values at different moments of time
Time = [];    % finally the list of times when we compute the invariants

while (t < Tf) % main loop in time
    t = t + dt; tloc = tloc + dt;
    
    % We implement the embedded Cash-Karp method of the order 5(4):
    k1 = dt*RHS(v, 0);
    k2 = dt*RHS(v + 0.2*k1, 0.2*dt);
    k3 = dt*RHS(v + 3/40*k1 + 9/40*k2, 0.3*dt);
    k4 = dt*RHS(v + 0.3*k1 - 0.9*k2 + 1.2*k3, 0.6*dt);
    k5 = dt*RHS(v - 11/54*k1 + 2.5*k2 - 70/27*k3 + 35/27*k4, dt);
    k6 = dt*RHS(v + 1631/55296*k1 + 175/512*k2 + 575/13824*k3 + 44275/110592*k4 + 253/4096*k5, 7/8*dt);
    
    w = 37/378*k1 + 250/621*k3 + 125/594*k4 + 512/1771*k6;
    w4 = 2825/27648*k1 + 18575/48384*k3 + 13525/55296*k4 + 277/14336*k5 + 0.25*k6;
    
    % we choose the new time step according to PI control:
    err2 = norm(real(ifft(w - w4)), inf);
    dt = abs(dt)*((tol/err2)^al)*((err1/tol)^be);
    err1 = err2;
    
    v = turn(v + w, -dt); % 5th order solution after changing the variables

    % adjust the time step 
    if ((tloc + dt > dtw) && (~pf))
        dt = dtw - tloc;
        pf = true;
    end % if (tloc)
    
    % if the time to plot has come:
    if ((tloc >= dtw - tol) && (pf))
        tloc = 0.0; pf = false;   % we reset the counters
        eta = real(ifft(v(:,1))); % find the free surface elevation in real space
        Plot(eta, t);             % plot the result
        
        phi = real(ifft(v(:,2)));             % velocity potential
        u = real(ifft(1i*k.*v(:,2)));         % its horizontal derivative
        
        chi_xi = 1 + real(ifft(ka.*v(:,1)));
        eta_xi = real(ifft(1i*k.*v(:,1)));
        psi_xi = real(ifft(-ka.*v(:,2)));

        mass = dxi*sum(eta.*chi_xi);          % integral representing the total mass
        moment = dxi*sum(eta.*u);             % total momentum
        kin = 0.5*dxi*sum(-psi_xi.*phi);      % kinetic energy
        pot = 0.5*dxi*g*sum((eta.^2).*chi_xi);% potential energy
        energy = kin + pot;                   % total energy
        
        % we store all the values in the lists:
        Time = [Time; t];
        Mass = [Mass; mass];
        Moment = [Moment; moment];
        Energy = [Energy; energy];
    end % if ()

end % while (t)