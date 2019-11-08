clearvars -global; clear all; clc; close all;
addpath('utilities','-end')
%==========================================================================
%                 ____                       ____   ____
%                |    | |    |    /\    | / |    | |
%                |    | |    |   /__\   |/  |    | |____
%                |    | |    |  /    \  |\  |    |      |
%                |____\ |____| /      \ | \ |____|  ____| 
%                      \
%==========================================================================
% Developped by Thomas Croissant and Philippe Steer
% Géosciences Rennes, CNRS/Université Rennes 1
% mail: thomas.croissant@durham.ac.uk or philippe.steer@univ-rennes.fr
% Last update: 01/2019
%==========================================================================
%                           INPUT PARAMETERS
%==========================================================================
f_Input_NZ();
%==========================================================================
%                              TOPOGRAPHY
%==========================================================================
% topogrid  = f_LoadDEM();
load('data\DEM\topogrid.mat');
%==========================================================================
%                                FAULT
%==========================================================================
[faultgrid,par] = f_FaultGrid(topogrid);
%==========================================================================
%                             EARTHQUAKES
%==========================================================================
Earthquake = f_Earthquakes_Generation();
%==========================================================================
%                             LANDSLIDES
%==========================================================================
% Generate Landslide population
[Earthquake,Landslide,LandslideDB] = f_Landslides_Generation(topogrid,faultgrid,Earthquake,par);

% Connection to river network
LandslideDB = f_LS_Connect2River(LandslideDB,topogrid);

%==========================================================================
%                        POST-SEISMIC EVACUATION
%==========================================================================
run.time = 0:1:2000;                                        % time vector
run.vcon = [100000 10 1 0.1];                                            % connection velocity

% Compute transport capacity associated to each landslide
LandslideDB.Qt     = f_QT(LandslideDB.dA,'eff',1,0,7.5);
LandslideDB.Vls_Qt = LandslideDB.V./LandslideDB.Qt;
LandslideDB.Vls_Qt(LandslideDB.Vls_Qt<1) = 1;
LandslideDB.t      = round(LandslideDB.t);

% Post-seismic evacuation of landslide clusters
for v = 1:length(run.vcon)
    % Processing
    LandslideDB.tcon_V = round(LandslideDB.d2s./run.vcon(v));
    LandslideDB.tcon_V(LandslideDB.conINI==1) = 0;
    for p = 1:length(LandslideDB.V)
        [run.Vt{v}(:,p),run.Vta{v}(:,p)] = Vs_Reconstruction(LandslideDB.Vls_Qt(p),LandslideDB.V(p),LandslideDB.tcon_V(p),LandslideDB.t(p),run.time);
    end
    % Compute output
    out.V(:,v)     = sum(run.Vt{v},2);
    out.Vnorm(:,v) = sum(run.Vt{v},2)./sum(run.Vt{v}(1,:));
    out.Qs(:,v)    = -gradient(out.V(:,v));
    for i = 1:length(run.time)
        out.Nls(i,v)  = numel(find(run.Vt{v}(i,:)>0));
        out.Nlsa(i,v) = numel(find(run.Vta{v}(i,:)>0));
    end
end

% See PLOTS.m for outputs!








