function [lsActive,Landslide] = f_LS_ActiveInit(Landslide)
%==========================================================================
%
% 
%
%==========================================================================
[~,I] = sort(Landslide.A);

lsActive.V1(1,:)  = Landslide.V;
lsActive.V2(1,:)  = Landslide.V;
lsActive.Vtot1(1) = sum(Landslide.V);
lsActive.Vtot2(1) = sum(Landslide.V);
lsActive.Ntot1(1) = numel(Landslide.V);
lsActive.Ntot2(1) = numel(Landslide.V);

Landslide.Qt   = Landslide.Qt;
Landslide.t2dn = Landslide.t2dn;
Landslide.Qt2  = repmat(Landslide.Qt,1000,1);

for i = 1:length(Landslide.Qt)
    if Landslide.t2dn(i) >= 1
        Landslide.Qt2(1:Landslide.t2dn(i),i) = 0;
    end
end

edge = 0.5:6.5;
lsActive.SO_1(:,1) = histcounts(double(Landslide.SO),edge);
lsActive.SO_2(:,1) = histcounts(double(Landslide.SO),edge);

lsActive.SOd_1 = double(Landslide.SO);
lsActive.SOd_2 = double(Landslide.SO);


% lsActive.Vini = [];
% lsActive.V    = [];
% lsActive.Qt = [];
% lsActive.VsQt = [];
% lsActive.x    = [];
% lsActive.y    = [];




