%%% -------------------------------------------------- %%%
%%% Numerical solution of full Euler equations deep w  %%%
%%%     using the conformal mapping technique          %%%
%%% Simulation of the Peregrine breather (NLS)         %%%
%%% -------------------------------------------------- %%%
%%% For the initialisation of the breather, please,    %%%
%%% refer to L. Shemer & L. Alperovich. 25, 051701,    %%%
%%% "Peregrine breather revisited", Phys. Fluids, 2013 %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, Khalifa University           %%%
%%% of Science and Technology, Abu Dhabi, UAE          %%%
%%% Date:   2025-05-23                                 %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Acknowledgements:                                  %%%
%%% * Prof. Didier Clamond, Univ.Nice Sophia Antipolis %%%
%%% * Dr. Bernard Ee, Tel Aviv University, Israel      %%%
%%% * Prof. Lev Shemer, Tel Aviv University, Israel    %%%
%%% * Prof. Francesco Fedele, Professor in the School  %%%
%%% of Civil and Environmental Engineering and in      %%%
%%% the School of Electrical and Computer Engineering  %%%
%%% at Georgia Institute of Technology, GA, USA        %%%
%%% -------------------------------------------------- %%%

% We clear the workspace
close all;
clearvars;
format longE

% Declaration of global variables
global et0 gom k ka kf k0 l omg omega sk xi N

%%% Physical parameters of the breather: the experiment of Shemer & Alperovich (2013) in a
%%% wave tank 18 m long, with water 0.6 m deep (k0*h = 7: deep water for the carrier), the
%%% wavemaker at x = 0 and wave gauges from 1.85 m to 13.3 m. All quantities in SI units.
g    = 9.81;					% gravity acceleration (m/s^2)
T0   = 0.587;					% carrier wave period (s)
et0  = 0.01;					% amplitude zeta0 of the carrier, far from the hump (m)
X0   = 9.0;					% distance from the wavemaker to the maximum of the breather (m)
Xend = X0;					% distance travelled by the centre of the group at the end (m)
om0  = 2*pi/T0;					% wave frequency (rad/s)
k0   = om0^2/g;					% wavenumber found from the dispersion relation (1/m)
l0   = 2*pi/k0;					% wavelength (m)
a0   = et0/sqrt(2);		    	% characteristic amplitude of the envelope (m)
epsi = a0*k0;					% characteristic wave steepness
cg   = 0.5*sqrt(g/k0);			% wave group speed in deep water (m/s)

% With these values, k0 = 11.68 1/m, l0 = 0.538 m, cg = 0.458 m/s and epsi = 0.0826. Up to X0
% the slope of the surface stays below 24 degrees. Followed beyond X0, towards the last gauge
% at 13.3 m (Xend = 13.3), the crests passing through the peak of the group steepen beyond
% 45 degrees near 11 m and, with N = 32768, overturn near 12 m: the computation enters the onset
% of breaking, where the exponential filter removes energy (3 to 4 per cent by 13.3 m) and the
% crests depend on the number of Fourier modes.

%%% Numerical parameters:
N = 16384;               		% number of Fourier modes in discrete solution

%%% Definition of the computational domain:
l = 35.0*l0;            		% half-length of the domain in \xi space (m)
dxi = 2*l/N;            		% distance between two \xi-points
xi = (1-N/2:N/2)'*dxi;  		% \xi-space discretization

%%% Discretization in the Fourier space:
% vector of wavenumbers:
k = [0:N/2 1-N/2:-1]'*pi/l;

% Window of the wave group. The wavemaker of the experiment runs for 70 periods, tapered over
% the first and the last two; the group moves at cg = l0/(2*T0), so in space it occupies
% 35 wavelengths, the half of the domain, with tapers one wavelength long at both ends:
window = @(x) min(1, max(0, (0.5*l - abs(x))/l0));

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

%%% Initial condition: the wave train at the time t = 0 when the centre of the group leaves the
%%% wavemaker. At the distance x from the wavemaker and the time t, the Peregrine breather of the
%%% spatial NLS equation is, after Eqs. (1), (2) and (4) of Shemer & Alperovich (2013),
%%%   zeta = Re[a0*A*exp(i*(k0*x - om0*t))],
%%%   A = -sqrt(2)*(1 - 4*(1 - 4i*Xis)/(1 + 4*Xit^2 + 16*Xis^2))*exp(-2i*Xis),
%%% with the scaled time Xit = epsi*om0*(x/cg - t) and the scaled distance Xis = epsi^2*k0*(x - X0)
%%% (xi and eta in the article): the envelope grows from 1.37*et0 at the wavemaker (Xis = -0.717)
%%% to 3*et0 at x = X0, where the Peregrine solution has its maximum at the time X0/cg.
Xis = @(x) epsi^2*k0*(x - X0);	% scaled space variable
Xit = @(x) epsi*om0*x/cg;		% scaled time variable at t = 0
A   = @(x) -sqrt(2)*(1 - 4*(1 - 4i*Xis(x))./(1 + 4*Xit(x).^2 + 16*Xis(x).^2)).*exp(-2i*Xis(x));
surface = @(x) real(a0*A(x).*exp(1i*k0*x)).*window(x);

% The elevation is a function of the physical abscissa x = \xi - H[eta]: the fixed point
% eta(\xi) = surface(x(\xi)) is found by iteration. The calm water outside the group stays at
% eta = 0, the still water level of the tank.
eta0 = surface(xi);
err = inf;
tol = 1e-14;
while (err > tol)
	eta1 = surface(Conf2Real(eta0));
	err = norm(eta1 - eta0, inf);
	eta0 = eta1;
end % while ()

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
Tf   = Xend/cg;				% final time (s): the group travels Xend at the speed cg
dtw  = Tf/round(16*Tf/T0);	% write time of the results: about sixteen times per period (s)
tol  = 5e-10;  				% tolerance parameter for time stepping
dt   = 0.01;    			% initial guess of the time step (s)
err1 = 1.0;   				% local error guess
tloc = 0.0;   				% local time variable
pf   = false;   			% plot flag

%%% Post-processing control:
compute_velocities = false;	% toggle for computing velocities in post-processing

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

while (t < Tf - 1e-9) % main loop in time, up to the last write time Tf

    h = dt;                    % the time step taken now
    t = t + h; tloc = tloc + h;
    
    % We implement the embedded Cash-Karp method of the order 5(4) for the variable
    % w(s) = exp(-s*L)*u(t + s), where u = [eta_hat, phi_hat] and L is the linear part of the
    % equations: RHS(w, s) evaluates the nonlinear terms at the state exp(s*L)*w of the stage
    % and rotates them back by exp(-s*L) (method of Lawson):
    k1 = h*RHS(v, 0);
    k2 = h*RHS(v + 0.2*k1, 0.2*h);
    k3 = h*RHS(v + 3/40*k1 + 9/40*k2, 0.3*h);
    k4 = h*RHS(v + 0.3*k1 - 0.9*k2 + 1.2*k3, 0.6*h);
    k5 = h*RHS(v - 11/54*k1 + 2.5*k2 - 70/27*k3 + 35/27*k4, h);
    k6 = h*RHS(v + 1631/55296*k1 + 175/512*k2 + 575/13824*k3 + 44275/110592*k4 + 253/4096*k5, 7/8*h);
    
    w = 37/378*k1 + 250/621*k3 + 125/594*k4 + 512/1771*k6;
    w4 = 2825/27648*k1 + 18575/48384*k3 + 13525/55296*k4 + 277/14336*k5 + 0.25*k6;
    
    % we choose the next time step according to PI control:
    err2 = norm(real(ifft(w - w4)), inf);
    dt = h*((tol/err2)^al)*((err1/tol)^be);
    err1 = err2;
    
    v = turn(v + w, -h); % 5th order solution, rotated by the linear flow over the step taken

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
        
        % Compute velocities if toggle is enabled:
        if (compute_velocities)

            %%% -------------------------------------------------- %%%
            %%% Surface velocity computation                       %%%
            %%% Contributed by: Prof. Francesco Fedele             %%%
            %%% Associate Professor, School of Civil &             %%%
            %%% Environmental Engineering, Georgia Institute       %%%
            %%% of Technology, Atlanta, Georgia, USA               %%%
            %%% -------------------------------------------------- %%%
            
            % Compute required derivatives in conformal space:
            phi_xi = real(ifft(1i*k.*v(:,2)));      % xi derivative of velocity potential
            psi_xi = -real(ifft(ka.*v(:,2)));       % xi derivative of stream function
            chi_xi = 1 + real(ifft(ka.*v(:,1)));    % xi derivative of conformal mapping
            chi_xixi = real(ifft(1i*k.*ka.*v(:,1))); % second xi derivative
            eta_xi = real(ifft(1i*k.*v(:,1)));      % xi derivative of free surface
            eta_xixi = -real(ifft(k.^2.*v(:,1)));   % second xi derivative
            
            % Compute physical space derivatives:
            % First x-derivative of free surface elevation
            eta_x = eta_xi ./ chi_xi;
            
            % Second x-derivative of free surface elevation
            eta_xx = (eta_xixi .* chi_xi - eta_xi .* chi_xixi) ./ (chi_xi.^3);
            
            % Jacobian of the conformal transformation:
            J = chi_xi.^2 + eta_xi.^2;
            
            % Compute surface velocities (U horizontal, V vertical):
            % Using the complex velocity formulation
            U = (phi_xi .* chi_xi + psi_xi .* eta_xi) ./ J;  % horizontal velocity
            V = -(psi_xi .* chi_xi - phi_xi .* eta_xi) ./ J; % vertical velocity
            
            % Compute xi derivatives of velocities:
            U_xi = real(ifft(1i*k.*fft(U)));
            V_xi = real(ifft(1i*k.*fft(V)));
            
            % Transform to physical x-derivatives:
            U_x = (U_xi .* chi_xi - V_xi .* eta_xi) ./ J; % horizontal derivative of U
            V_x = (U_xi .* eta_xi + V_xi .* chi_xi) ./ J; % horizontal derivative of V
            
            % Note: The velocities U, V and their derivatives U_x, V_x
            % are now available for further analysis or storage
        end % if (compute_velocities)

    end % if ()

end % while (t)

% Conservation of the mass and of the energy (the Hamiltonian) over the run:
fprintf('t = %.2f s: mass %.3e m^2 (largest change %.1e m^2), energy %.6e (largest relative change %.1e)\n', ...
        t, Mass(1), max(abs(Mass - Mass(1))), Energy(1), max(abs(Energy - Energy(1)))/Energy(1));