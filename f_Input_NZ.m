function f_Input_NZ()
%==========================================================================
%                           QUAKOS Parameters
%==========================================================================
% Define all the parameter values necessary for a QUAKOS simulation.
% These parameters are global
%
% Dev: P. Steer and T. Croissant
% Last update: 10/2018
%==========================================================================
global parFAULT parRUP parEQ parPGA parLSD parSURF;
%==========================================================================
%% Parameters - grid DEM
parSURF.res      = 100;                                                    % DEM resolution high 
parSURF.reslarge = 1000;                                                   % DEM resolution low 
%==========================================================================
%% Parameters - Fault dimension
% Load Fault planar geometry and DEM of the area of interest
% Define left and right (in longitude sense) tips of the fault (lat/lon)
parFAULT.Latl   = -44.51;                                                  % left x-coordinates <decim degree>
parFAULT.Lonl   = 167.83;                                                  % left y-coordinates <decim degree>
parFAULT.Latr   = -42.50;                                                  % right x-coordinates <decim degree>
parFAULT.Lonr   = 171.85;                                                  % right y-coordinates <decim degree>
% Find UTM zone '##A'
[~,~,utm_temp]  = deg2utm((parFAULT.Latl+parFAULT.Latr)/2,(parFAULT.Lonl+parFAULT.Lonr)/2);
parFAULT.utm    = [utm_temp(1:2) utm_temp(4)];
% Fault geometry
parFAULT.dip    = 60;                                                      % dipping angle of the fault plane (0<dip<=90) on the right of the fault (when going from start to end)
parFAULT.rake   = 172;                                                     % orientation of the slip on the fault plane (-180<=rake<=180) reverse=90, normal=-90, left-lat=0, right-lat=180
parFAULT.ztop   = 0;                                                       % depth of the fault trace (m)
parFAULT.ztip   = 19000;                                                   % depth of the fault tip (m)
parFAULT.res    = 100;                                                     % fault discretization step (in x and y)
% x-y coordinates of the start point (defined on the left side)
[parFAULT.xll,parFAULT.yll] = deg2utmsimple(parFAULT.Latl,parFAULT.Lonl,str2double(parFAULT.utm(1:2)),parFAULT.utm(3));
% x-y coordinates of the second point
[parFAULT.xrr,parFAULT.yrr] = deg2utmsimple(parFAULT.Latr,parFAULT.Lonr,str2double(parFAULT.utm(1:2)),parFAULT.utm(3));
% fault trace length
parFAULT.L      = sqrt((parFAULT.xrr-parFAULT.xll).^2+(parFAULT.yrr-parFAULT.yll).^2);
% strike angle between the fault trace and the north defined clockwise (0<=strike<360)
parFAULT.strike = atan2d((parFAULT.xrr-parFAULT.xll),(parFAULT.yrr-parFAULT.yll));
if parFAULT.strike<0
    parFAULT.strike=parFAULT.strike+360;
end
% Fault width (m)
parFAULT.W = (parFAULT.ztip-parFAULT.ztop)/sind(parFAULT.dip);

%==========================================================================
% Parameters - Fault dimension
%==========================================================================
%% Rupture parameters [Leonard, 2010]
% Define fault style, 2 choices: 'Strike-Slip' or 'Reverse'
parRUP.type = 'Strike-Slip';

if (strcmp(parRUP.type,'Strike-Slip'))
    % Strike-slip fault - For L > 45 km
    parRUP.mu       = 33e9;                                                % Shear modulus (Pa)
    parRUP.Beta     = 2/3;                                                 % parameter Leonard or Beta=1
    parRUP.C1       = 15;                                                  % parameter Leonard
    parRUP.C10      = parFAULT.W;                                          % Crust mean thickness - SHOULD WE SET IT EQUAL TO FAAULt WIDTH ????
    parRUP.C11      = 1;                                                   % Crust mean thickness - SHOULD WE SET IT EQUAL TO FAAULt WIDTH ????
    parRUP.C2       = 3.6e-5;                                              % parameter Leonard
    % Define maximum and minimum earthquake magnitude
    parRUP.Lmax     = parFAULT.L;                                          % Max. length
    parRUP.Wmax     = parFAULT.W;                                          % Max. width
    if parRUP.Lmax< 5000
        parRUP.MoLmax   = parRUP.mu*parRUP.C11^(3/2)*parRUP.C2.*parFAULT.L.^(3);
        % parRUP.MoWmax   = parRUP.mu*parRUP.C11^(3/2)*parRUP.C2.*parFAULT.L.^(3/2)                                          % Max. moment (based on width)
    elseif parRUP.Lmax>=5000 && parRUP.Lmax<45000
        parRUP.MoLmax   = parRUP.mu*parRUP.C1^(3/2)*parRUP.C2.*parFAULT.L.^(5/2);
    else
        parRUP.MoLmax   = parRUP.mu*parRUP.C10^(3/2)*parRUP.C2.*parFAULT.L.^(3/2);
    end
    %   parRUP.MoWmax   = parRUP.C10;                                      % Max. moment (based on width)
    %   parRUP.Momax    = min(parRUP.MoLmax,parRUP.MoWmax);                % Max. moment
    parRUP.Momax    = parRUP.MoLmax;                                       % Max. moment
    parRUP.Mwmax    = 2/3*log10(parRUP.Momax)-6.07;                        % Max. magnitude
    parRUP.Lmin     = parFAULT.res;                                        % Min. length
    parRUP.Wmin     = parFAULT.res;                                        % Min. width
    parRUP.MoLmin   = parRUP.mu*parRUP.C1^(3/2)*parRUP.C2*parRUP.Lmin.^(3/2*(1+parRUP.Beta));
    parRUP.MoWmin   = parRUP.mu*parRUP.C1^(-3/(2*parRUP.Beta))*parRUP.C2*parRUP.Wmin.^(3/2*(1+1/parRUP.Beta));
    parRUP.Momin    = max(parRUP.MoLmin,parRUP.MoWmin);
    parRUP.Mwmin    = 2/3*log10(parRUP.Momin)-6.07;
    
elseif (strcmp(parRUP.type,'Reverse'))
    % Reverse fault
    parRUP.mu       = 30e9;                                                % Shear modulus (Pa)
    parRUP.Beta     = 2/3;                                                 % parameter Leonard or Beta=1
    parRUP.C1       = 17.5;                                                % parameter Leonard
    parRUP.C2       = 3.8e-5;                                              % parameter Leonard
    % Define maximum and minimum earthquake magnitude
    parRUP.Lmax     = parFAULT.L;
    parRUP.Wmax     = parFAULT.W;
    parRUP.MoLmax   = parRUP.mu*parRUP.C1^(3/2)*parRUP.C2*parFAULT.L.^(3/2*(1+parRUP.Beta));
    parRUP.MoWmax   = parRUP.mu*parRUP.C1^(-3/(2*parRUP.Beta))*parRUP.C2*parFAULT.W.^(3/2*(1+1/parRUP.Beta));
    parRUP.Momax    = min(parRUP.MoLmax,parRUP.MoWmax);
    parRUP.Mwmax    = 2/3*log10(parRUP.Momax)-6.07;
    parRUP.Lmin     = parFAULT.res;
    parRUP.Wmin     = parFAULT.res;
    parRUP.MoLmin   = parRUP.mu*parRUP.C1^(3/2)*parRUP.C2*parRUP.Lmin.^(3/2*(1+parRUP.Beta));
    parRUP.MoWmin   = parRUP.mu*parRUP.C1^(-3/(2*parRUP.Beta))*parRUP.C2*parRUP.Wmin.^(3/2*(1+1/parRUP.Beta));
    parRUP.Momin    = max(parRUP.MoLmin,parRUP.MoWmin);
    parRUP.Mwmin    = 2/3*log10(parRUP.Momin)-6.07;
end

%==========================================================================
%% Tectonic input scenario
% Type of mainshocks - 'Periodic' or 'Gutenberg-Richter'
parEQ.type = 'Periodic';

parEQ.mw_min    = parRUP.Mwmin;                                            % Minimum Magnitude (deterministic)
parEQ.mw_max    = parRUP.Mwmax;                                            % Maximum Magnitude (deterministic)
parEQ.mw_dm     = 0.1;                                                     % Magnitude steps
parEQ.mw        = parEQ.mw_min:parEQ.mw_dm:parEQ.mw_max;
parEQ.T         = 1*270*365-1;                                             % Total time duration (days)
parEQ.b         = 1;                                                       % GR param for Mainshocks and Aftershocks
if (strcmp(parEQ.type,'Periodic'))
    parEQ.period    = 270*365;                                             % Periodicity of the Main shocks
    parEQ.Mwevent   = parEQ.mw_max;                                        % Magnitude of the reccurent event
    parEQ.Nmain     = floor(parEQ.T./parEQ.period)+1;                      % Number of mainshocks
elseif (strcmp(parEQ.type,'Gutenberg-Richter'))
    parEQ.srate     = 0.1;                                                 % Number of mainshock per day
    parEQ.Nmain     = parEQ.T*parEQ.srate;                                 % Number of mainshocks
    parEQ.a         = log10(parEQ.T.*parEQ.srate) + parEQ.b*parEQ.mw_min;  % GR param for Mainshocks - for Mw>=Mwmin
end
parEQ.Dm        = 1.25;                                                    % GR-Bath param
parEQ.c         = 0.1;                                                     % temporal Omori param
parEQ.p         = 1.25;                                                    % temporal Omori param
parEQ.d         = 4.0;                                                     % spatial Omori param
parEQ.q         = 1.35;                                                    % spatial Omori param
parEQ.mw_c      = 4.4;                                                     % Minimum critical Magnitude (by default parEQ.mw_min)

%==========================================================================
%% Parameters - PGA model
% PGA - [O'Brien et al. 2016]
parPGA.b1       = 2.175;                                                   % PGA empirical law
parPGA.b2       = 0.765;                                                   % PGA empirical law
parPGA.b3       = 1.398;                                                   % PGA empirical law
% PGA - [Campbell & Borzognia 2008] (Page 10), PGA unit <g>
parPGA.c        =  1.880;                                                  % PGA empirical law
parPGA.n        =  1.180;                                                  % PGA empirical law
parPGA.c0       = -1.715;                                                  % PGA empirical law
parPGA.c1       =  0.500;                                                  % PGA empirical law
parPGA.c2       = -0.530;                                                  % PGA empirical law
parPGA.c3       = -0.262;                                                  % PGA empirical law
parPGA.c4       = -2.118;                                                  % PGA empirical law
parPGA.c5       =  0.170;                                                  % PGA empirical law
parPGA.c6       =  5.600;                                                  % PGA empirical law
parPGA.c7       =  0.280;                                                  % PGA empirical law
parPGA.c8       = -0.120;                                                  % PGA empirical law
parPGA.c9       =  0.490;                                                  % PGA empirical law
parPGA.c10      =  1.058;                                                  % PGA empirical law
parPGA.c11      =  0.040;                                                  % PGA empirical law
parPGA.c12      =  0.610;                                                  % PGA empirical law
parPGA.k1       =  865.0;                                                  % PGA empirical law
parPGA.k2       = -1.186;                                                  % PGA empirical law
parPGA.k3       =  1.839;                                                  % PGA empirical law
parPGA.VS30     =  180;                                                    % S wave velocity in the top 30 m - VS30 (m/s)
parPGA.A1100    =  0.10;                                                   % reference PGA (in g) at VS=1100 m/s
parPGA.Z2_5     =  0;                                                      % Sediment depth - Z2.5 (m)

%==========================================================================
%% Parameters - Landslide model
parLSD.type     = 'bedrock';
% PGA - Landslide density relationship [Meunier et al. 2007]
parLSD.alpha    = 4;                                                       % Dependency of Pls on PGA, taiwan = 2.7
parLSD.beta     = 0.5;                                                     % Background Pls
parLSD.PGAcr    = 1.7;                                                     % Threshold of PGA for LSD <m.s-2>
% PDF of landslide area [Malamud et al, 2004]
parLSD.a        = 2000;                                                    % position of max <m2>
parLSD.s        = -200;                                                    % roll over <m2>
parLSD.rho      = 1.4;                                                     % tail slope -(rho + 1)
parLSD.Al_min   = 50;                                                      % min landslide area <m2>
parLSD.Al_max   = 1.5E6;                                                     % max landslide area <m2>

% Landslide volume-area relationship [Larsen et al, 2010]
if (strcmp(parLSD.type,'bedrock'))
    % Parameterization bedrock landslides
    parLSD.alpha_l  = 0.05;
    parLSD.gamma    = 1.5;
else
    % Parameterization mixed bedrock-soil landslides
    parLSD.alpha_l  = 0.15;
    parLSD.gamma    = 1.33;
end

%% Out
disp('Parameters loaded ...')
