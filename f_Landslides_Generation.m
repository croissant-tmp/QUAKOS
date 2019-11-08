function [Earthquake,Landslide,LandslideDB] = f_Landslides_Generation(surfgrid,faultgrid,Earthquake,par)
%==========================================================================
%                   Generate landslide population
%==========================================================================
% Assign landslide population to each seismic event
%
% Input 
%       - surfgrid : strucure containing the geographic and topographic 
%       information of the study area
%       - faultgrid : structure containing all the properties of the fault
%       grid
%       - Earthquake : structure containing all the properties associated 
%       with each earthquakes
%       - par : structure containing information relevant to the PGA
%       computation
% 
% Output
%       - Earthquake : structure containing all the properties associated 
%       with each earthquakes
%       - Landslide : structure containing information relevant to the PGA
%       computation
%       - LandslideDB : same than Landslide but filtered for particular
%       catchments (DB: drainage basin)
%
% Called functions: f_PGA; 
%                   f_EQ_rupture_and_distance; 
%                   f_LS_Sample; 
%                   f_LS_Init; 
%                   f_LS_Conc;
%
% Dev.: T Croissant and P. Steer
% Last update: 10/2018
%==========================================================================
% Initialize Landslide structures
Landslide   = f_LS_Init();
LandslideDB = f_LS_Init();

for i = 1:length(Earthquake.ind_mw)
    % Rupture patch
    [par,Earthquake] = f_EQ_rupture_and_distance(Earthquake,faultgrid,surfgrid,par,Earthquake.ind_mw(i));
    
    % Sophisticated NGA model for PGA (Campbell & Bozorgnia 2008)
    Earthquake.PGA{i} = f_PGA(Earthquake.Mw(Earthquake.ind_mw(i)),Earthquake.Ztor(Earthquake.ind_mw(i)),par);
    
    % Landslides generation
    [ls_tmp,lsDB_tmp] = f_LS_Sample(Earthquake.PGA{i},Earthquake.t(Earthquake.ind_mw(i)),Earthquake.dt(Earthquake.ind_mw(i)),Earthquake.g(Earthquake.ind_mw(i)),surfgrid);
    % Concatenation of landslide population
    Landslide   = f_LS_Conc(Landslide,ls_tmp);
    LandslideDB = f_LS_Conc(LandslideDB,lsDB_tmp);
end

Landslide.Atot(1) = []; Landslide.Alstot(1) = []; Landslide.Vlstot(1) = []; Landslide.N(1) = []; 
LandslideDB.Atot(1) = []; LandslideDB.Alstot(1) = []; LandslideDB.Vlstot(1) = []; LandslideDB.N(1) = []; 

disp('Landslides generated ...')

%==========================================================================


