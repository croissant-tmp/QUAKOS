%==========================================================================
%                                 PLOTS
%==========================================================================
%              Earthquake series - Mainshocks + Aftershocks
%--------------------------------------------------------------------------
img.V = log10(Landslide.Vlstot);
figure
scatter(Earthquake.t/365,Earthquake.Mw,5,[0.6 0.6 0.6],'filled')
hold on
scatter(Earthquake.t(Earthquake.ind_mw)/365,Earthquake.Mw(Earthquake.ind_mw),10,img.V,'filled')
scatter(Earthquake.t(Earthquake.Mw==max(Earthquake.Mw))/365,Earthquake.Mw(Earthquake.Mw==max(Earthquake.Mw)),100,'rp','filled')
ylim([2.5 8])
xlabel('Time (yr)'); ylabel('Earthquake magnitude');
legend('Mainshock','Aftershocks - trigerring landslides','Aftershocks')
%==========================================================================
%                          Co-seismic Landslide Map
%--------------------------------------------------------------------------
topogrid.DEM_hr.Z(topogrid.DEM_hr.Z==0)=NaN;
img.DB = topogrid.DB;
img.DB.Z(img.DB.Z>1) = 1;
img.DB.Z(img.DB.Z==0) = NaN;
atmp = topogrid.A;
atmp.Z(topogrid.DB.Z == 0) = NaN;
St   = STREAMobj(topogrid.FD,atmp>5E5);
db   = unique(LandslideDB.DB);

img.pls  = 1*(10.^(0.301+0.28316*log10(LandslideDB.V)));
img.plsc = 1*(10.^(0.301+0.28316*log10(Landslide.V)));

figure
imageschs(topogrid.DEM_hr,img.DB,'colormap',[238/255 223/255 204/255],'colorbar',false,'ticklabels','none','usepermanent',true)
hold on
scatter(LandslideDB.Sx,LandslideDB.Sy,10,[30/255 149/255 237/255],'filled')
for i = 1:length(db)
    bw   = bwboundaries(topogrid.DB.Z==db(i));
    img.xdb_bound = topogrid.DEM_hr.refmat(3,2)-bw{1}(:,1)*100;
    img.ydb_bound = topogrid.DEM_hr.refmat(3,1)+bw{1}(:,2)*100;
    plot(img.ydb_bound,img.xdb_bound,'r','LineWidth',2)
end
plotstreamorder(St,'colormap',[30/255 149/255 237/255],'LineWidth',max([1 2 3 4 5]/2,1),'legend',false);
scatter(Landslide.x,Landslide.y,img.plsc,'w','filled','MarkerEdgeColor','k');
scatter(LandslideDB.x,LandslideDB.y,img.pls,'k','filled','MarkerEdgeColor','w');
scatter(LandslideDB.x(LandslideDB.conINI==1),LandslideDB.y(LandslideDB.conINI==1),img.pls(LandslideDB.conINI==1),'r','filled','MarkerEdgeColor','k','MarkerFaceColor',[255/255 236/255 102/255]);
xlim([4.2E5 4.9E5]); ylim([5.17E6 5.23E6]);

%==========================================================================
%                 Post-seismic evolution of the cluster
%--------------------------------------------------------------------------
figure
subplot(2,2,1); plot(run.time,out.V);  xlabel('Time (yr)'); ylabel('Norm. V_l_s_,_t_o_t'); %ylim([0 1]);
subplot(2,2,2); semilogy(run.time,out.Qs); ylim([1E4 1E8]); xlabel('Time (yr)'); ylabel('Mobilisation flux (m^3.yr^-^1)');
subplot(2,2,3); semilogy(run.time,out.Nls); ylim([1 100000]); xlabel('Time (yr)'); ylabel('Total landslide number');
subplot(2,2,4); semilogy(run.time,out.Nlsa); ylim([1 10000]); xlabel('Time (yr)'); ylabel('Total active landslide number');

%==========================================================================
%                              Animations
%--------------------------------------------------------------------------
% f_PGA_gif_lowRes(topogrid,Earthquake,'test.gif',1)
