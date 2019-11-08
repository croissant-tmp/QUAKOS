function Mainshock = f_EQ_Mainshock()
%==========================================================================
% Generate mainshock sequence
%==========================================================================
% Load relevant parameters
global parEQ parFAULT;

if (strcmp(parEQ.type,'Periodic'))
 
    % Mainshock magnitude
    Mainshock.Mw = ones(1,parEQ.Nmain).*parEQ.Mwevent;
    %Mainshock.Mw = (7.9-4.4)*rand(1,500)+4.4;
    
    % Mainshock time
    Mainshock.t  = 0:parEQ.period:parEQ.T;
    %%Mainshock.t  = 1:500;
    % Mainshock.t = cumsum(round(((263+68)-(263-68)).*rand(1,parEQ.Nmain) + (263-68)) * 365);
    % Mainshock.t = cumsum(round((330-30).*rand(1,parEQ.Nmain)+30) * 365);
    % Mainshock.t = Mainshock.t - min(Mainshock.t);
    
elseif (strcmp(parEQ.type,'Gutenberg-Richter')) 
    % Mainshock magnitude
    Mainshock.Mw  = [];
    for i = 1:length(parEQ.mw)
        N            = round(10^(parEQ.a-parEQ.b*parEQ.mw(i)))-round(10^(parEQ.a-parEQ.b*(parEQ.mw(i)+parEQ.mw_dm)));
        Mainshock.Mw = [Mainshock.Mw parEQ.mw(i) + parEQ.mw_dm.*(rand(1,N)-0.5)];
    end
    % Shuffle vector
    Mainshock.Mw  = Mainshock.Mw(randperm(length(Mainshock.Mw)));
    % Mainshock time
    Mainshock.t   = parEQ.T.*rand(1,length(Mainshock.Mw));
end    

% Compute the rupture geometry of the Earthquake (L, W, D) - Leonard 2010
[Mainshock.L,Mainshock.W,Mainshock.D] = fRuptureGeometry(Mainshock.Mw);

% Mainshock Coordinates epicenter
dx            = 0;
Mainshock.x   = round(dx + ((parFAULT.W-dx) - dx).*rand(1,length(Mainshock.Mw)));
Mainshock.y   = round(dx + ((parFAULT.L-dx) - dx).*rand(1,length(Mainshock.Mw)));
Mainshock.z   = zeros(1,length(Mainshock.Mw));

% Initialize - Mainshock aftershock generation number (by definition 0)
Mainshock.g   = zeros(1,length(Mainshock.Mw));

% Initialize -  Mainshock dt
Mainshock.dt  = zeros(1,length(Mainshock.Mw));
