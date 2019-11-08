function Aftershock = f_EQ_Conc(Aftershock,tmp)
%==========================================================================
% Aftershock magnitude
Aftershock.Mw = [Aftershock.Mw; tmp.Mw];

% Aftershock spatial pryoperties
Aftershock.x  = [Aftershock.x; tmp.x];
Aftershock.y  = [Aftershock.y; tmp.y];
Aftershock.z  = [Aftershock.z; tmp.z];

% Aftershock time properties
Aftershock.g  = [Aftershock.g; tmp.g];
Aftershock.t  = [Aftershock.t; tmp.t];
Aftershock.dt = [Aftershock.dt; tmp.dt];

% Aftershock properties
Aftershock.L  = [Aftershock.L; tmp.L];
Aftershock.W  = [Aftershock.W; tmp.W];
Aftershock.D  = [Aftershock.D; tmp.D];


