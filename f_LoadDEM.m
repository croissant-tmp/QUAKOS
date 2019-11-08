function topogrid  = f_LoadDEM()
%==========================================================================
% Load DEM and compute relevant geomophic properties of the landscape
%
% Output: topogrid <structure>
%
% supplementary functions used:
%
% topotoolbox_v2.2: - reproject2utm.m
%                   - getcoordinates.m
%                   - gradient8.m
%                   - fillsinks.m
%                   - GRIDobj.m
%                   - FLOWobj.m
%                   - flowacc.m
%
% sources:          - https://topotoolbox.wordpress.com/
%          
%
% dev: T. Croissant & P. Steer
% Last update: 04/2017
%==========================================================================
global parSURF
%==========================================================================
% Load DEM 
load('DEM.mat');    % High resolution DEM
load('DEMf.mat');   % Low resolution DEM
%==========================================================================
% Process high-res DEM
topogrid.DEM_hr               = DEM;
[topogrid.xvec,topogrid.yvec] = getcoordinates(DEM); topogrid.yvec = topogrid.yvec';
[topogrid.x,topogrid.y]       = meshgrid(topogrid.xvec,topogrid.yvec);

% Process low-res DEM
topogrid.DEM_lr               = DEMf;
[parSURF.x,parSURF.y]         = getcoordinates(DEMf);
[topogrid.xl,topogrid.yl]     = meshgrid(parSURF.x,parSURF.y);

%==========================================================================
% Compute grid slope (D8 algorithm)
topogrid.S_hr   = gradient8(DEM,'deg');                                
topogrid.S_hr.Z = double(topogrid.S_hr.Z);
topogrid.S_hr20 = topogrid.S_hr;
topogrid.S_hr20.Z(topogrid.S_hr.Z < atand(0.2)) = 0;

topogrid.S_lr   = gradient8(DEMf,'deg');                                   
topogrid.S_lr.Z = double(topogrid.S_lr.Z);
topogrid.S_lr20 = topogrid.S_lr;
topogrid.S_lr20.Z(topogrid.S_lr.Z < atand(0.2)) = 0;

%==========================================================================
% Compute drainage area
topogrid.FD   = FLOWobj(DEM,'preprocess','fill');
topogrid.A    = flowacc(topogrid.FD).*(topogrid.FD.cellsize^2);
topogrid.A.Z(topogrid.A.Z<5E5) = 0;
topogrid.A.Z(DEM.Z==0) = 0;

%surfgrid.SO   = streamorder(surfgrid.FD,surfgrid.A.Z>5E5);

% Stream extraction
% surfgrid.Str   = surfgrid.A;
% surfgrid.Str.Z = false(surfgrid.A.size);
% surfgrid.Str.Z(surfgrid.A.Z>5E5) = true;

%==========================================================================
% Extract catchments 
[xDB,yDB]     = f_findOulet(topogrid,5E6);
topogrid.DB   = drainagebasins(topogrid.FD,topogrid.xvec(xDB),topogrid.yvec(yDB));
topogrid.DB.Z = double(topogrid.DB.Z);
topogrid.DB   = shufflelabel(topogrid.DB);
% nrDB          = numel(unique(surfgrid.DB.Z(:)))-1;                         % nr of drainage basins 
STATS         = regionprops(topogrid.DB.Z,'PixelIdxList','Area','Centroid');

for i = 1:length(STATS)
    STATS(i).Area = str2double(num2str(round(STATS(i).Area * topogrid.DB.cellsize^2/1e6)));
    DB_area_tmp(i) = STATS(i).Area;
    DB_cent_tmp(i,:) = STATS(i).Centroid;
end

indDB   = find(DB_area_tmp>=10);
DB_area = DB_area_tmp(indDB);
DB_cent = DB_cent_tmp(indDB,:);

DB_cent(:,3) = sqrt(DB_cent(:,1).^2+DB_cent(:,2).^2)*topogrid.DB.cellsize;

pxl_c = []; pxl_x = []; pxl_y = []; pxl_a = []; pxl_n = []; pxl_o = [];

for i = 1:length(indDB)
    tmp_px = STATS(indDB(i)).PixelIdxList;
    [x,y]  = ind2coord(topogrid.DEM_hr,tmp_px);
    pxl_c = [pxl_c tmp_px'];
    pxl_x = [pxl_x x'];
    pxl_y = [pxl_y y'];
    pxl_a = [pxl_a DB_area(i).*ones(1,length(tmp_px))];
    pxl_n = [pxl_n i.*ones(1,length(tmp_px))];
    pxl_o = [pxl_o DB_cent(i,3).*ones(1,length(tmp_px))];
end

topogrid.mCatch(1,:) = pxl_c;
topogrid.mCatch(2,:) = pxl_x;
topogrid.mCatch(3,:) = pxl_y;
topogrid.mCatch(4,:) = pxl_a;
topogrid.mCatch(5,:) = pxl_n;
topogrid.mCatch(6,:) = pxl_o;

%==========================================================================
res =[];
ris =[];
for i = 1:max(max(topogrid.DB.Z))
    atmp = topogrid.A;
    atmp.Z = zeros(DEM.size);
    atmp.Z(topogrid.DB.Z==i) = topogrid.A.Z(topogrid.DB.Z==i);
    St = STREAMobj(topogrid.FD,atmp>5E5);
    res = [res; St.IXgrid];
    ris = [ris; i.*ones(length(St.IXgrid),1)];
end

topogrid.StreamDB(1,:) = res;
topogrid.StreamDB(2,:) = ris;


disp('DEM loaded ...')