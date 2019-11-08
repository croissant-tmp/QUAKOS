function [faultgrid,par] = f_FaultGrid(surfgrid)
global parFAULT parPGA;
%==========================================================================
% Generation of the fault plane grid
% 
% Input 
%       - parFAULT, parPGA [global parameters]
%       - surfgrid : strucure containing the geographic and topographic 
%       information of the study area
% 
% Output
%       - faultgrid : structure containing all the properties of the fault
%       grid
%       - par : structure containing information relevant to the PGA
%       computation
%
% Associated functions:
%      - frelocateonfault:  Series of mainshocks
%
% Dev: T. Croissant and P. Steer
% Last update: 10/2018
%==========================================================================
% Define the fault plane grid
[faultgrid.x,faultgrid.y] = meshgrid(0:parFAULT.res:parFAULT.W,0:parFAULT.res:parFAULT.L);
faultgrid.z = zeros(size(faultgrid.x));
%==========================================================================
% Save the initial coordinates of the earthquakes and fault grid
faultgrid.xini = faultgrid.x;  faultgrid.yini  = faultgrid.y;  faultgrid.zini  = faultgrid.z;
%==========================================================================
% Translate and rotate the fault grid to match the orientation of the fault plane geometry (position, strike and dip)
faultgrid = frelocateonfault(faultgrid);
%==========================================================================
% VS30 map initialization
par.VS30  = ones(size(surfgrid.xl)).*parPGA.VS30;
% A1100 map initialization
par.A1100 = ones(size(surfgrid.xl)).*parPGA.A1100;
% Z2_5 map initialization
par.Z2_5  = ones(size(surfgrid.xl)).*parPGA.Z2_5;
%==========================================================================
% Generate matrixes of cumulated displacement
par.Dcum  = zeros(size(faultgrid.x));
%==========================================================================
% Identify surface point above the rupture
par.indover = inpolygon(surfgrid.xl,surfgrid.yl,[faultgrid.x(1,1) faultgrid.x(1,end) faultgrid.x(end,end) faultgrid.x(end,1)],[faultgrid.y(1,1) faultgrid.y(1,end) faultgrid.y(end,end) faultgrid.y(end,1)]);
par.indover = find(par.indover==1);

%==========================================================================
disp('Fault grid created ...')