function [par,Earthquake] = f_EQ_rupture_and_distance(Earthquake,faultgrid,surfgrid,par,ii)
global parFAULT parRUP parSURF;
%==========================================================================
%
%
%
% Dev.: T Croissant P. Steer
% Last update: 10/2018
%==========================================================================
% Determine pixel that include the earthquake
Dist    = sqrt((Earthquake.x(ii)-faultgrid.x).^2+(Earthquake.y(ii)-faultgrid.y).^2);
[~,ind] = min(Dist(:));
[ni,nj] = size(Dist);
[i,j]   = ind2sub(size(Dist),ind);

% Compute the Earthquake width and length in pixels
nL = min(floor(Earthquake.L(ii)/parFAULT.res),floor(parRUP.Lmax/parFAULT.res));
nW = min(floor(Earthquake.W(ii)/parFAULT.res),floor(parRUP.Wmax/parFAULT.res));

% Distribute this width and length around the central pixel
% Start with a random repartition
nLdown = round(rand(1).*nW);
nLup   = nL-nLdown;
nWdown = round(rand(1).*nW);
nWup   = nW-nWdown;
% Readjust if it violates the grid borders - L
DnL = 0;
if i+nLup > ni
    DnL = ni-(i+nLup);
elseif i-nLdown < 1
    DnL = 1-(i-nLdown);
end
nLdown = nLdown-DnL;
nLup   = nLup+DnL;
% Readjust if it violates the grid borders - W
DnW=0;
if j+nWup>nj
    DnW=nj-(j+nWup);
elseif j-nWdown<1
    DnW=1-(j-nWdown);
end
nWdown = nWdown-DnW;
nWup   = nWup+DnW;

% Vertical distance between the Earth'surface and the top of the rupture patch
% Ztor=(nj-(j+nWup)).*parFAULT.res.*tand(parFAULT.dip)+parFAULT.ztop;
% Zdor=(nj-(j-nWdown)).*parFAULT.res.*tand(parFAULT.dip)+parFAULT.ztop;
Earthquake.Ztor(ii) = faultgrid.z(i,j-nWdown);                                            % Depth to the top of the rupture
Earthquake.Zdor(ii) = faultgrid.z(i,j+nWup);                                              % Depth to the tip of the rupture

% Co-seismic displacement on the fault plane (raster)
par.Disp = zeros(size(faultgrid.x));
par.Disp(i-nLdown:i+nLup,j-nWdown:j+nWup) = Earthquake.D(ii);

% Compute distance to the horizontal projection of the rupture plane
iUL = i-nLdown ; jUL = j+nWup;   % Corner up-left
iUR = i+nLup   ; jUR = j+nWup;   % Corner up-right
iDL = i-nLdown ; jDL = j-nWdown; % Corner down-left
iDR = i+nLup   ; jDR = j-nWdown; % Corner down-right
xv  = [faultgrid.x(iUL,jUL),faultgrid.x(iUR,jUR),faultgrid.x(iDR,jDR),faultgrid.x(iDL,jDL)];  % rupture edges polygon x
yv  = [faultgrid.y(iUL,jUL),faultgrid.y(iUR,jUR),faultgrid.y(iDR,jDR),faultgrid.y(iDL,jDL)];  % rupture edges polygon y
in  = inpolygon(surfgrid.xl(par.indover),surfgrid.yl(par.indover),xv,yv);                     % Find surface point above the rupture
in  = par.indover(in);

if isempty(in)==0
    I = zeros(size(surfgrid.xl));                                          % Create a mask
    I(in) = 1;
    % Distance to the mask
    par.Rjb = bwdist(I).*parSURF.reslarge;
else
    par.Rjb = sqrt((surfgrid.xl-Earthquake.x(ii)).^2+(surfgrid.yl-Earthquake.y(ii)).^2);
end
% Take an approximate value of the distance to the rupture plane
par.Rrup = sqrt(par.Rjb.^2+((Earthquake.Ztor(ii)+Earthquake.Zdor(ii))./2).^2);

% Cumulated slip on the fault plan
par.Dcum = par.Dcum + par.Disp;

