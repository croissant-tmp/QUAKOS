function [var] = frelocateonfault(var)
%==========================================================================
% Parameters
global parFAULT;

% Stock coordinate
x     = var.x;
y     = var.y;
z     = var.z; 
% rotate the coordinates to respect the dip
var.x = (x.*cosd(parFAULT.dip))-((z-parFAULT.ztop).*sind(parFAULT.dip)); 
var.y = y; 
var.z = (x.*sind(parFAULT.dip))+((z-parFAULT.ztop).*cosd(parFAULT.dip))+parFAULT.ztop; 

% rotate the coordinates to respect the strike
% Stock coordinate
x     = var.x;
y     = var.y;
z     = var.z; 
var.x = (x.*cosd(-parFAULT.strike))-(y.*sind(-parFAULT.strike)); 
var.y = (x.*sind(-parFAULT.strike))+(y.*cosd(-parFAULT.strike)); 
var.z = var.z;

% Translate coordinates
% Stock coordinate
x     = var.x;
y     = var.y;
z     = var.z; 
var.x = x + parFAULT.xll;
var.y = y + parFAULT.yll;
var.z = z;