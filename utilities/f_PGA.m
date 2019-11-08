function PGA = f_PGA(MW,ZTOR,par)
%==========================================================================
% --- OUTPUT
% PGA:   Peak ground acceleration raster (g) - geometrical mean of PGAh
%
% --- INPUTS
% --- scalars
% MW:    Earthquake magnitude
% RAKE:  orientation of the slip on the fault plane (-180<=rake<=180) reverse=90, normal=-90, left-lat=0, right-lat=180 (degrees)
% DIP:   Fault plane dip angle (0<dip<=90) on the right of the fault (when going from start to end) (degrees)
% ZTOR:  Vertical distance between the Earth'surface and the top of the rupture patch (m)
%
% --- rasters
% RRUP:  Closest distance to the rupture patch (m)
% RJB:   Closest distance to the surface projection of the rupture patch (m)
%
% Based on Campbell and Bozorgnia, 2008, (Journal: Earthquake Spectra), hereafter referred as CB2008 
% Written by P. Steer
% Last update: 08/2017.
%==========================================================================
% Parameters
global parPGA parFAULT;
%--------------------------------------------------------------------------
% Convert all distances in km
ZTOR     = ZTOR     ./ 1000;
par.Rrup = par.Rrup ./ 1000;
par.Rjb  = par.Rjb  ./ 1000;
par.Z2_5 = par.Z2_5 ./ 1000;
%--------------------------------------------------------------------------
% Magnitude dependency function (scalar)
% CB2008, eq. [2]
if MW<=5.5
    fmag = parPGA.c0+parPGA.c1*MW;
elseif MW>5.5 && MW<=6.5
    fmag = parPGA.c0+parPGA.c1*MW+parPGA.c2*(MW-5.5);
else
    fmag = parPGA.c0+parPGA.c1*MW+parPGA.c2*(MW-5.5)+parPGA.c3*(MW-6.5);
end
%--------------------------------------------------------------------------
% Distance dependency function (raster)
% CB2008, eq. [3]
fdis = (parPGA.c4+parPGA.c5*MW).*log(sqrt((par.Rrup).^2+parPGA.c6^2));
%--------------------------------------------------------------------------
% Fault mechanism dependency function (scalar)
% CB2008, eq. [4] & [5]
if ZTOR < 1
    ffltZ = ZTOR;
else
    ffltZ = 1;
end

if parFAULT.rake > 30 && parFAULT.rake < 150
    FRV = 1;
else
    FRV = 0;
end

if parFAULT.rake > -150 && parFAULT.rake < -30
    FNM = 1;
else
    FNM = 0;
end

fflt = parPGA.c7.*FRV.*ffltZ + parPGA.c8.*FNM;
%--------------------------------------------------------------------------
% Hanging-wall term dependency (raster - beacuse of FHNGR)
% CB2008, eq. [6], [7], [8], [9] & [10]
FHNGR = zeros(size(par.Rrup)); 
FHNGR(par.Rjb == 0) = 1;
% eq. [7]
ind = find(par.Rjb>0);
if ZTOR < 1
    ind11   = find(par.Rrup >= sqrt(par.Rjb.^2+1)-par.Rjb);
    ind12   = find(par.Rrup <  sqrt(par.Rjb.^2+1)-par.Rjb);
    ind21   = find(par.Rrup >= sqrt(par.Rjb.^2+1));
    ind22   = find(par.Rrup <  sqrt(par.Rjb.^2+1));
    ind1121 = intersect(ind11,ind21);
    ind1122 = intersect(ind11,ind22);
    ind1221 = intersect(ind21,ind21);
    ind1222 = intersect(ind22,ind22);
    FHNGR(intersect(ind,ind1121)) = par.Rrup(intersect(ind,ind1121))./par.Rrup(intersect(ind,ind1121));
    FHNGR(intersect(ind,ind1122)) = par.Rrup(intersect(ind,ind1122))./sqrt(par.Rjb(intersect(ind,ind1122)).^2+1);
    FHNGR(intersect(ind,ind1221)) = (sqrt(par.Rjb(intersect(ind,ind1221)).^2+1)-par.Rjb(intersect(ind,ind1221)))./par.Rrup(intersect(ind,ind1221));
    FHNGR(intersect(ind,ind1222)) = (sqrt(par.Rjb(intersect(ind,ind1222)).^2+1)-par.Rjb(intersect(ind,ind1222)))./sqrt(par.Rjb(intersect(ind,ind1222)).^2+1);
else
    FHNGR(ind) = (par.Rrup(ind)-par.Rjb(ind))./par.Rrup(ind);
end
% eq. [8]
if MW<=6.0
    FHNGM = 0;
elseif MW>6.0 && MW<=6.5
    FHNGM = 2*(MW-6.0);
else
    FHNGM = 1;
end
% eq. [9]
if ZTOR>=20
    FHNGZ=0;
else
    FHNGZ=(20-ZTOR)/20;
end
% eq. [10]
if parFAULT.dip <= 70
    FHNGD = 1;
else
    FHNGD = (90-parFAULT.dip)/90;
end
% eq. [6]
fhng = parPGA.c9.*FHNGR.*FHNGM.*FHNGZ.*FHNGD;
%--------------------------------------------------------------------------
% Shallow site response dependency (Raster)
% CB2008, eq. [11]
fsite = zeros(size(par.Rrup));
fsite(par.VS30 < parPGA.k1) = parPGA.c10.*log(par.VS30(par.VS30<parPGA.k1)./parPGA.k1) + ...
                              parPGA.k2.*( log(par.A1100(par.VS30<parPGA.k1) + ...
                              parPGA.c.*(par.VS30(par.VS30<parPGA.k1)./parPGA.k1).^parPGA.n) - ...
                              log(par.A1100(par.VS30<parPGA.k1)+parPGA.c));
                        
fsite(par.VS30>=parPGA.k1 & par.VS30<1100) = (parPGA.c10+parPGA.k2.*parPGA.n) .* ...
                                             log(par.VS30(par.VS30>=parPGA.k1 & par.VS30<1100)./parPGA.k1);
                                         
fsite(par.VS30 >= 1100)     = (parPGA.c10+parPGA.k2.*parPGA.n).*log(1100./parPGA.k1);
%--------------------------------------------------------------------------
% Basin response dependency (raster)
% CB2008, eq. [12]
fsed = zeros(size(par.Rrup));
fsed(par.Z2_5<1) = parPGA.c11.*(par.Z2_5(par.Z2_5<1)-1);
fsed(par.Z2_5>=1 & par.Z2_5<=3) = 0;
fsed(par.Z2_5>3) = parPGA.c12.*parPGA.k3.*exp(-0.75).*(1-exp(-0.25.*(par.Z2_5(par.Z2_5>3)-3)));
%--------------------------------------------------------------------------
% COMPUTE PGA BY COMBINING ALL PREVIOUS TERMS
% CB2008, eq. [1]
PGA = exp(fmag + fdis + fflt + fhng + fsite + fsed);                       % PGA (g)

%--------------------------------------------------------------------------
% subplot(3,2,1);imagesc(fmag); axis equal tight;colorbar;title('fmag');
% subplot(3,2,2);imagesc(fdis); axis equal tight;colorbar;title('fdis');
% subplot(3,2,3);imagesc(fflt); axis equal tight;colorbar;title('fflt');
% subplot(3,2,4);imagesc(fhng); axis equal tight;colorbar;title('fhng');
% subplot(3,2,5);imagesc(fsite);axis equal tight;colorbar;title('fsite');
% subplot(3,2,6);imagesc(fsed); axis equal tight;colorbar;title('fsed');




