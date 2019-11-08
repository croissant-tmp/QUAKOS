function [V_time,V_time_a] = Vs_Reconstruction(V_Qt,Vini,tcon,tG,time)
%==========================================================================
% V_Qt = 50;
% Vini = 500000;
% time = run.time;
% tcon = 10;
% tG   = 10;
%==========================================================================
% Eros to Quakos post-seismic volume evolution for each % interval
V_perc = [0 0.2:0.1:1];
V_t    = Vini - Vini .* V_perc;

t80  = 23.0.*(0.033*V_Qt).*(1+0.033*V_Qt.^0.79).^((0.1-1)/0.79);

t = [0 ...
     4*(0.02*V_Qt).*(1+0.02*V_Qt.^1.05).^((0.01-1)/1.05)  ...              % 20%
     7.*(0.02.*V_Qt).*(1+0.02.*V_Qt.^1.0).^((0.01-1)/1.0) ...              % 30%
     07.0 .* (0.033.*V_Qt) .* (1+0.033.*V_Qt.^1.01).^((0.1-1)/1.01) ...    % 40%
	 10.8 .* (0.033.*V_Qt) .* (1+0.033.*V_Qt.^1).^((0.1-1)/1)       ...    % 50%
	 13.0 .* (0.033.*V_Qt) .* (1+0.033.*V_Qt.^0.95).^((0.1-1)/0.95) ...    % 60%
	 17.0 .* (0.033.*V_Qt) .* (1+0.033.*V_Qt.^0.93).^((0.1-1)/0.93) ...    % 70%
	 23.0.*(0.033*V_Qt).*(1+0.033*V_Qt.^0.79).^((0.1-1)/0.79)       ...    % 80%
     V_Qt ...                                                              % 90%
     -(Vini.*0.2 - ((Vini.*0.1-Vini.*0.2)./(V_Qt-t80)).*t80) ./ ((Vini.*0.1-Vini.*0.2)./(V_Qt-t80))]; % 100%
%==========================================================================
V_tmp = interp1(t,V_t,time,'linear');
V_tmp(isnan(V_tmp)) = 0;

tmp = V_tmp;
tmp(tmp==0) = [];
    
V_time(1:tG + tcon + length(tmp)) = tmp(1);
V_time(tG + tcon + 1:tG + tcon + length(tmp)) = tmp;
V_time(1:tG) = 0;

if length(V_time) < length(time)
    V_time(tG + tcon + length(tmp)+1:length(time)) = 0;
else
    V_time(length(time)+1:end) = [];
end

V_time_a = V_time;
V_time_a(1:tG + tcon) = 0;
if length(V_time_a) > length(time)
    V_time_a(length(time)+1:end) = [];
end



% run.V_t(1:LandslideDB.t(p)+LandslideDB.tcon(p)+length(tmp),p) = tmp(1);
% run.V_t(LandslideDB.t(p)+LandslideDB.tcon(p)+1:LandslideDB.t(p)+LandslideDB.tcon(p)+length(tmp),p) = tmp;
% run.V_t(1:LandslideDB.t(p),p) = 0;
