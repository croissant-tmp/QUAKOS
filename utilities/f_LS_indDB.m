function LandslideDB = f_LS_indDB(Landslide,surfgrid)
%==========================================================================
% Correction of Landslide
nlsDB = unique(Landslide.DB);

ind_badDB = [];
for i = 1:length(nlsDB)
    aDB = sum(sum(surfgrid.DB.Z==nlsDB(i)))*surfgrid.DB.cellsize.^2;
    if aDB<5E5
        ind_badDB = [ind_badDB find(Landslide.DB == nlsDB(i))];
    end   
end

% index the landslide within the west coast catchments
indDB  = find(Landslide.DB>0);
indDB  = setdiff(indDB,ind_badDB);
% Sorting the indexes by catchments number
[~,is] = sort(Landslide.DB(indDB));
indDB  = indDB(is);

%==========================================================================
% Define the properties of Landslide belonging to the West Coast
%==========================================================================
% Reduce the properties of the whole population to the one only present in the catchments
LandslideDB.x  = Landslide.x(indDB);                                       % x-axis coordinate
LandslideDB.y  = Landslide.y(indDB);                                       % y-axis coordinate
LandslideDB.IX = coord2ind(surfgrid.DEM_hr,LandslideDB.x,LandslideDB.y)';  % index coordinate
LandslideDB.z  = Landslide.z(indDB);
LandslideDB.t  = round(Landslide.t(indDB)./365);                           % Landslide trigger time
LandslideDB.dt = zeros(1,length(LandslideDB.t));                           % Landslide trigger time
LandslideDB.g  = Landslide.g(indDB);                                       % Landslide generation
LandslideDB.A  = Landslide.A(indDB);                                       % Landslide area [m2]
LandslideDB.V  = Landslide.V(indDB);                                       % Landslide volume [m3]
LandslideDB.conINI = Landslide.conINI(indDB);                              % Landslide initial connectivity status [0 1]
LandslideDB.DB = Landslide.DB(indDB);                                      % Landslide catchment

LandslideDB.Atot   = Landslide.Atot ;
LandslideDB.Alstot = sum(LandslideDB.A);
LandslideDB.Vlstot = sum(LandslideDB.V);
LandslideDB.N      = numel(LandslideDB.A);
