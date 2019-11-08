function Landslide = f_LS_Connectivity(Landslide)
%==========================================================================
% Assign an initial connectivity to landlslides triggered by an event
%
% Dev: T. Croissant
% Last update: 06/2017
%==========================================================================
% Bins of landslide area necessary to compute the connectivity function
edgesA = logspace(2,7,50);
expA   = linspace(2,7,50);
bincA  = 10.^((expA(2:end) + expA(1:end-1))/2);
% Connectivity function here derived from Li et al, [2016]
conIni   = 10.^(-0.06+0.34*log10(bincA)); conIni(bincA > 1E6) = 100;

%==========================================================================
% Assign full connectivity to the first x% of landslide in respect with the
% percentage for every landslide generation.

% Find the succession of landslide triggering events
gen   = unique(Landslide.t);
conLS = [];
for g = 1:length(gen)
    
    Atmp     = Landslide.A(Landslide.t == gen(g));
    conLStmp = zeros(1,length(Atmp));
    % Count the numbers of landslides inside each bins of connectivity
    nA       = histcounts(Atmp,edgesA);
    nAc      = round(nA .* (conIni/100));
    % sum(nAc)./sum(nA)
    
    for i = 1:length(nAc)
        if nAc(i) > 0
            ind = find(Atmp > edgesA(i) & Atmp <= edgesA(i+1));
            ind = ind(randperm(length(ind)));
            conLStmp(ind(1:nAc(i))) = 1;
        else
            ind = [];
        end
    end
    conLS = [conLS conLStmp];
    clear conLStmp Atmp;
end

Landslide.conINI = conLS;

% ckc = numel(find(Landslide.conINI==1))./length(Landslide.conINI);

%disp('Connectivity computed ...')

%==========================================================================
% figure
% semilogx(bincA,conIni,'k.')
% 
% figure
% semilogx(bincA,nA,bincA,nAc)





















