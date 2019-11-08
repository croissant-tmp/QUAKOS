function Earthquake = f_Earthquakes_Generation()
global parFAULT parEQ;
%==========================================================================
% Generation of earthquakes series (Mainshocks + Aftershocks)
% 
% Input parameters are loaded with the function f_input.m 
%       - parFAULT, parEQ
% 
% Output
%       - Earthquake : structure containing all the properties associated with each
% earthquakes
%
% Associated functions:
%      - f_EQ_Mainshock:  Series of mainshocks
%      - f_EQ_BASS_local: Series of aftershocks associated to each
%                         mainshocks
%
% Dev: T. Croissant and P. Steer
% Last update: 10/2018
%==========================================================================
% Generation of mainshock series
Mainshock  = f_EQ_Mainshock();
%==========================================================================
% Generation of aftershock based on the mainshocks
Aftershock = f_EQ_BASS(Mainshock);
%==========================================================================
% Earthquakes series properties - concatenation of Main and Aftershocks
Earthquake.Mw = horzcat(Mainshock.Mw,vertcat(Aftershock.Mw)');             % Magnitude
Earthquake.x  = horzcat(Mainshock.x,vertcat(Aftershock.x)');               % x-coord
Earthquake.y  = horzcat(Mainshock.y,vertcat(Aftershock.y)');               % y-coord
Earthquake.z  = horzcat(Mainshock.z,vertcat(Aftershock.z)');               % depth
Earthquake.g  = horzcat(Mainshock.g,vertcat(Aftershock.g)');               % generation
Earthquake.t  = horzcat(Mainshock.t,vertcat(Aftershock.t)');               % time
Earthquake.dt = horzcat(Mainshock.dt,vertcat(Aftershock.dt)');             % time 
Earthquake.L  = horzcat(Mainshock.L',vertcat(Aftershock.L)');              % rupture length
Earthquake.W  = horzcat(Mainshock.W',vertcat(Aftershock.W)');              % rupture width
Earthquake.D  = horzcat(Mainshock.D',vertcat(Aftershock.D)');              % surface displacement
%--------------------------------------------------------------------------
% Sort by earthquake time
[~,inew] = sort(Earthquake.t,'ascend');
Earthquake.Mw = Earthquake.Mw(inew);
Earthquake.x  = Earthquake.x(inew);
Earthquake.y  = Earthquake.y(inew);
Earthquake.z  = Earthquake.z(inew);
Earthquake.g  = Earthquake.g(inew);
Earthquake.t  = Earthquake.t(inew);
Earthquake.dt = Earthquake.dt(inew);
Earthquake.L  = Earthquake.L(inew);
Earthquake.W  = Earthquake.W(inew);
Earthquake.D  = Earthquake.D(inew);
%--------------------------------------------------------------------------
% Remove Earthquakes outside the fault bounding box
inew = find(Earthquake.x>0 & Earthquake.x<parFAULT.W & Earthquake.y>0 & Earthquake.y<parFAULT.L & Earthquake.t<=parEQ.T);
Earthquake.Mw = Earthquake.Mw(inew);
Earthquake.x  = Earthquake.x(inew);
Earthquake.y  = Earthquake.y(inew);
Earthquake.z  = Earthquake.z(inew);
Earthquake.g  = Earthquake.g(inew);
Earthquake.t  = Earthquake.t(inew);
Earthquake.dt = Earthquake.dt(inew);
Earthquake.L  = Earthquake.L(inew);
Earthquake.W  = Earthquake.W(inew);
Earthquake.D  = Earthquake.D(inew);
%--------------------------------------------------------------------------
% Convert distances in km
Earthquake.L_km  = Earthquake.L ./ 1000;
Earthquake.W_km  = Earthquake.W ./ 1000;
%--------------------------------------------------------------------------
% Save the initial coordinates of the earthquakes
Earthquake.xini = Earthquake.x;
Earthquake.yini = Earthquake.y;
Earthquake.zini = Earthquake.z;
%--------------------------------------------------------------------------
% Translate and rotate earthquakes to match the orientation of the fault plane geometry (position, strike and dip)
[Earthquake] = frelocateonfault(Earthquake);
%--------------------------------------------------------------------------
Earthquake.Ztor   = zeros(1,length(Earthquake.Mw));
Earthquake.ind_mw = find(Earthquake.Mw >= parEQ.mw_c);
%--------------------------------------------------------------------------
disp('Earthquake frequency-magnitude computed ...')









