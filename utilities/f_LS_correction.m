function Landslide = f_LS_correction(Landslide,surfgrid)
%==========================================================================
% Correction of Landslide
% Delete landslides that are located on pixels with z = 0 or S < 11.3°
indLScor = [];
for i = 1:length(Landslide.V)
    
    indx = find(surfgrid.xl(1,:)>Landslide.x(i),1,'first');
    indy = find(surfgrid.yl(:,1)<Landslide.y(i),1,'first');
    
    indS = surfgrid.S_lr.Z(indy,indx);
    indZ = surfgrid.DEM_lr.Z(indy,indx);
    
    if(indS < 11.3 || indZ == 0)
        indLScor = [indLScor i];
    end
end

% Apply corrections
Landslide.x(indLScor)  = []; Landslide.y(indLScor)  = []; Landslide.z(indLScor)  = [];
Landslide.t(indLScor)  = []; Landslide.dt(indLScor) = []; Landslide.g(indLScor)  = [];
Landslide.A(indLScor)  = []; Landslide.V(indLScor)  = [];

