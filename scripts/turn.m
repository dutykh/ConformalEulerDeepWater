%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function v = turn (w, dt)

    global kf omg gom omega N
    
    v = zeros(N, 2);
    
    cosw = cos(omega*dt);
    sinw = sin(omega*dt);
    
    v(:,1) = kf.*(cosw.*w(:,1) - omg.*sinw.*w(:,2));
    v(:,2) = kf.*(gom.*sinw.*w(:,1) + cosw.*w(:,2));
end