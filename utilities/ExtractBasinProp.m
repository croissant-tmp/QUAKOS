nlsDB = unique(LandslideDB.DB);
for n = 343%1:length(nlsDB)
    % Worflow to extract only the stream of the considered catchment
    [basin.FD,basin.A,basin.Gt,basin.S] = f_loadDB(surfgrid,n);
    save(strcat('data\basin\basin_',num2str(n),'.mat'),'basin');
end
