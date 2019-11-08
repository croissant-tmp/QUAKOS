function [Landslide,LandslideDB] = f_LS_Sample(PGA,t,dt,g,surfgrid)
global parLSD parSURF;
%==========================================================================
%
% Output:
% Landslide: landslide properties: area, volume, coordinates - format: structure
% nLSD:      total landslide number                          - format: single value
% Atot:      total area affected by landslides               - format: single value
% Alstot:    total landslide area                            - format: single value
% Vlstot:    total landslide volume                          - format: single value
%
% Input:
% PGA:   Peak ground acceleration (g) (geometrical mean of PGAh) - format RASTER
% t:
% dt:
% g:
% surfgrid
%
% references: Malamud et al, 2004; Meunier et al, 2007; Larsen et al, 2010.
%
% Dev.: T. CRoissant & P. Steer
% Last update: 10/2018.
%==========================================================================
% Compute landslide spatial density
%==========================================================================
PGA_ms2 = PGA .* 9.81;
% Landslide spatial density map 
Pls = zeros(size(PGA_ms2));
Pls(PGA_ms2 > parLSD.PGAcr) = parLSD.alpha.*PGA(PGA_ms2>parLSD.PGAcr) - parLSD.beta;
% Apply Pls map corrections
Pls(Pls<0) = 0;
Pls(surfgrid.S_lr20.Z == 0) = 0;
%Pls(surfgrid.DB_lr.Z  == 0) = 0;

%==========================================================================
Landslide   = f_LS_Init();
LandslideDB = f_LS_Init();

if (max(Pls(:)) > 0)
    % Total area affected by landsliding (m2) 
    Landslide.Atot   = sum(sum(Pls>0)).*parSURF.reslarge^2;                % <m2>
    % Total landslide area (m2)
    Alstot = sum(sum(Pls./100)).*parSURF.reslarge^2;             % <m2>
    %======================================================================
    % Compute PDF Area of landslides
    Al.Al     = logspace(log10(parLSD.Al_min),log10(parLSD.Al_max),2E3);   % Landslide area <m2>
    Al.dAl    = diff(Al.Al); Al.dAl = [Al.dAl Al.dAl(end)];
    Al.pdf_Al = 1./(parLSD.a.*gamma(parLSD.rho)).*(parLSD.a./(Al.Al-parLSD.s)).^(parLSD.rho+1).*exp(-parLSD.a./(Al.Al-parLSD.s));
    if (sum(Al.pdf_Al.*Al.dAl)<0.99 || sum(Al.pdf_Al.*Al.dAl)>1.01); disp('pdf ls. Area error!'); end
    Al.cdf_Al = cumtrapz(Al.Al,Al.pdf_Al);
    
    %======================================================================
    % Sample landslide from the PDF
    mc.nl     = 100000;
    mc.cdfAl  = rand(1,mc.nl);
    mc.Al_ini = interp1(Al.cdf_Al, Al.Al, mc.cdfAl); mc.Al_ini(isnan(mc.Al_ini)) = [];
    mc.cumAl  = cumsum(mc.Al_ini);
    mc.indAl  = find(mc.cumAl > Alstot,1 ,'first');
    mc.Al_vec = mc.Al_ini(1:mc.indAl);
    ALS       = mc.Al_vec;
    
    %======================================================================
    
    nLSD = numel(ALS); % Total number of expected landslides
    
    % Map of number of landslides per pixel    
    nLSDMAP = round(Pls./sum(sum(Pls)).*nLSD);
    nLSD2   = sum(sum(nLSDMAP));
    
    if (nLSD2 > nLSD)
        ALS(nLSD:nLSD2) = ones(1,nLSD2-nLSD+1) .* parLSD.Al_min;
    end
    
    nLSD = nLSD2;
    VLS  = parLSD.alpha_l .* ALS .^ parLSD.gamma;
    ind  = find(nLSDMAP > 0);
    
    k = 0;
    for i = 1:numel(ind)
        for j = 1:nLSDMAP(ind(i))
            k = k+1;
            Landslide.x(k)  = surfgrid.xl(ind(i))+(rand-0.5).*parSURF.reslarge;
            Landslide.y(k)  = surfgrid.yl(ind(i))+(rand-0.5).*parSURF.reslarge;
            Landslide.z(k)  = surfgrid.DEM_lr.Z(ind(i));
            Landslide.t(k)  = t;
            Landslide.dt(k) = dt;
            Landslide.g(k)  = g;
            Landslide.A(k)  = ALS(k);
            Landslide.V(k)  = VLS(k);
        end
    end
        
    Landslide = f_LS_correction(Landslide,surfgrid);
    
    Landslide.IX = coord2ind(surfgrid.x,surfgrid.y,Landslide.x,Landslide.y)';
    % Assign landslide to their catchments
    Landslide.DB = surfgrid.DB.Z(Landslide.IX);
    % Assign initial connectivity to the landslide cluster
    Landslide = f_LS_Connectivity(Landslide);
    % Cluster properties
    Landslide.Alstot = sum(Landslide.A);
    Landslide.Vlstot = sum(Landslide.V);   
    Landslide.N      = numel(Landslide.V);
    %==========================================================================
    % Filter out the landslides out of the West Coast catchments
    %==========================================================================
    LandslideDB = f_LS_indDB(Landslide,surfgrid);
end





