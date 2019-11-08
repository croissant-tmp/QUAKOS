function [FD,A,Gt,S] = f_loadDB(surfgrid,N)

DEM   = surfgrid.DEM_hr;
DEM.Z = nan.*ones(DEM.size);
DEM.Z(surfgrid.DB.Z==N) = surfgrid.DEM_hr.Z(surfgrid.DB.Z==N);

bw   = bwboundaries(surfgrid.DB.Z==N);

img.xdb_bound = DEM.refmat(3,2)-bw{1}(:,1)*100;
img.ydb_bound = DEM.refmat(3,1)+bw{1}(:,2)*100;

FD    = FLOWobj(DEM,'preprocess','fill');
A     = flowacc(FD).*(FD.cellsize^2);
A.Z(isnan(DEM.Z)) = NaN;
Gt    = gradient8(DEM,'tan');
S     = STREAMobj(FD,A>5E5);