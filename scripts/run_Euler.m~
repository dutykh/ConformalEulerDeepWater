%%% -------------------------------------------------- %%%
%%% Numerical solution of full Euler equations deep w  %%%
%%%     using the conformal mapping technique          %%%
%%% Simulation of the Peregrine breather               %%%
%%% -------------------------------------------------- %%%
%%% Computation of the potential:                      %%%
%%%   \phi_t + g*\eta = O(\eps^2), z = 0               %%%
%%%   -i*\omega*\phi + g\eta = 0,  z = 0               %%%
%%%   \phi = -i*g / \omega * \eta, z = 0               %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

% We clear the workspace
close all, clear all
format longE

% Parameters needed in other functions:
global a0 cg0 e0 g gom k ka kf k0 l l0 om0 omg omega sk win xi N

%%% Physical parameters:
g = 9.8;                % gravity acceleration

%%% Physical parameters related to the test-case of Peregrine's breather:
l0 = 0.538;             % characteristic wavelength
k0 = 2*pi/l0;           % characteristic wavenumber
zeta0 = 0.01;           % carrier amplitude
a0 = zeta0/sqrt(2);     % characteristic wave amplitude
e0 = a0*k0;             % wave steepness
om0 = sqrt(g*k0);       % angular frequency
cg0 = 0.5*sqrt(g/k0);   % group velocity
T0 = 0.587;             % wave period
L = 19.0;               % domain is specified to be [-L*l0, L*l0]
t0 = -2.0*l0*L;         % initial moment of time

%%% Numerical parameters:
N = 4096;               % number of Fourier modes in discrete solution

%%% Definition of the computational domain:
l = 19.0*l0;            % half-length of the domain in \xi space
dxi = 2*l/N;            % distance between two \xi-points
xi = (1-N/2:N/2)'*dxi;  % \xi-space discretization

% Change the 2nd argument to tune the window size:
win = tukeywin(N, 0.25);% Tukey windowing function

%%% Discretization in the Fourier space:
% vector of wavenumbers:
k = [0:N/2 1-N/2:-1]'*pi/l;

% Exponential de-aliasing treatment proposed by T. Hou et al.:
Lambd = 15;
kam = 0.70*max(abs(k));
kf = exp(-5.0*(abs(k)/kam).^Lambd); % dealiasing function

%%% Pseudo-differential operators:
sk = sign(k);   % needed for the Hilbert transform = -i*sign(k)
ka = abs(k);    % Hilbert transform of the derivative

% Needed for the integration of linear terms:
omega = sqrt(g*ka);
omg = omega/g;
gom = g./omega; gom(1) = 0.0;

%%% Initial condition computation:
[eta0, phi0] = InitCond(t0, xi);

v = zeros(N,2);
v(:,1) = fft(eta0);
v(:,2) = fft(phi0);

%%% Define the graphic window:
FigHandle = figure(1);
set(FigHandle, 'Position', [100, 100, 849, 495]);
% and plot the initial condition:
Plot(eta0, t0);
pause;

% FFTW plan:
fftw('planner', 'hybrid');

%%% Time-stepping parameters:
t = t0;       % the discrete time variable
Tf = 6.0;     % final simulation time
dtw = 0.05;   % write time of the results
tol = 5e-14;  % tolerance parameter
dt = 0.01;    % initial guess of the time step
err1 = 1.0;   % local error guess
tloc = 0.0;   % local time variable
pf = false;   % plot flag

% Parameters for the adaptive time-stepping (PI control):
% taken according to the book of Hairer (Part II)
al = 0.7/5; be = 0.4/5;

% %%% In order to check the accuracy of computations we can follow
% %%% the following quantities:
% Mass = [];    % List of the masses at different moments of time
% Moment = [];  % List of the momenta at different moments of time
% Energy = [];  % List of the (total) energy values at different moments of time
% Time = [];    % finally the list of times when we compute the invariants

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
        
        % phi = real(ifft(v(:,2)));             % velocity potential
        % u = real(ifft(1i*k.*v(:,2)));         % its horizontal derivative
        
        % chi_xi = 1 + real(ifft(ka.*v(:,1)));
        % eta_xi = real(ifft(1i*k.*v(:,1)));
        % psi_xi = real(ifft(-ka.*v(:,2)));

        % mass = dxi*sum(eta.*chi_xi);          % integral representing the total mass
        % moment = dxi*sum(eta.*u);             % total momentum
        % kin = 0.5*dxi*sum(-psi_xi.*phi);      % kinetic energy
        % pot = 0.5*dxi*g*sum((eta.^2).*chi_xi);% potential energy
        % energy = kin + pot;                   % total energy
        
        % % we store all the values in the lists:
        % Time = [Time; t];
        % Mass = [Mass; mass];
        % Moment = [Moment; moment];
        % Energy = [Energy; energy];
    end % if()

end % while (t)