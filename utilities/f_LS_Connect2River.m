function LandslideDB = f_LS_Connect2River(LandslideDB,surfgrid)
%==========================================================================
% Apply corrections to Landslides
%
% dev: T. Croissant
% Last update: 10/2018
%==========================================================================
nlsDB = unique(LandslideDB.DB);
inc = 1;
for n = 1:length(nlsDB)
    % Worflow to extract only the stream of the considered catchment
    load(strcat('data\basin\basin_',num2str(nlsDB(n)),'.mat'));
    
    % landslide population
    ixLS_db = LandslideDB.IX(LandslideDB.DB==nlsDB(n));
    
    for i = 1:length(ixLS_db)
        [IX,distance] = flowpathextract(basin.FD,ixLS_db(i));
        tmp = basin.A.Z(IX);
        if sum(find(tmp>5E5,1)) > 0
            LandslideDB.ixS(inc) = IX(find(tmp>5E5,1)); 
            LandslideDB.d2s(inc) = distance(find(tmp>5E5,1));
            LandslideDB.dA(inc)  = basin.A.Z(LandslideDB.ixS(inc));
        else
            [xtmp,ytmp] = ind2coord(surfgrid.DEM_hr,ixLS_db(i));
            dist_tmp    = sqrt(sum(bsxfun(@minus, [basin.S.x basin.S.y],[xtmp,ytmp]).^2,2));
            LandslideDB.ixS(inc) = basin.S.IXgrid(find(dist_tmp == min(dist_tmp),1),:);
            LandslideDB.d2s(inc) = dist_tmp(find(dist_tmp == min(dist_tmp),1));
            LandslideDB.dA(inc)  = basin.A.Z(LandslideDB.ixS(inc));
        end
        inc = inc + 1;
    end
    disp(n./length(nlsDB)*100)
end

rndD = 100.*rand(1,length(LandslideDB.IX)) - 50;
rndD(LandslideDB.d2s==0) = 0;

LandslideDB.d2s = double(LandslideDB.d2s+rndD);
LandslideDB.d2s = LandslideDB.d2s - 50;
LandslideDB.d2s = max(LandslideDB.d2s,0);

atmp = surfgrid.A;
atmp.Z(surfgrid.DB.Z == 0) = NaN;
St   = STREAMobj(surfgrid.FD,atmp>5E5);

LandslideDB.Sx  = St.x;
LandslideDB.Sy  = St.y;





% plot(LandslideDB.d2s,d2s,'.')
%==========================================================================
% Assign initial connectivity 
%LandslideDB.conINI_bu = LandslideDB.conINI;
%LandslideDB = f_LS_reassignConnectivity(LandslideDB);


% D    = flowdistance(surfgrid.FD,St.IXgrid);
% D.Z = double(D.Z);
% 
% d2s2 = D.Z(LandslideDB.IX);
% 
% rndD = 100.*rand(1,length(LandslideDB.IX)) - 50;
% rndD(d2s==0) = 0;
% 
% LandslideDB.d2s = double(d2s+rndD);
% LandslideDB.d2s = LandslideDB.d2s - 50;
% LandslideDB.d2s = max(LandslideDB.d2s,0);