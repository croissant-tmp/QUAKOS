function [L,W,D] = fRuptureGeometry(Mw)

global parRUP;


%==========================================================================
% Determine the geometrical properties of earthquakes (Leonard 2010)
Mo = 10.^(3/2.*Mw+9.105);

if (strcmp(parRUP.type,'Strike-Slip'))
    
    % Leonard (2010)- Strike-slip fault - General case (trilinear relationship)
    % if L < 5 km
    ind1=find(Mo<parRUP.mu*parRUP.C11^(3/2)*parRUP.C2*5000.^(3));
        L(ind1)  = (Mo(ind1)./(parRUP.mu*parRUP.C11^(3/2)*parRUP.C2)).^(1/3);
        W(ind1)  = parRUP.C11.*L(ind1);
        D(ind1)  = parRUP.C2*parRUP.C11^(1/2)*L(ind1); 
    % if 5 km < L < 45km
    ind2=find(Mo>=parRUP.mu*parRUP.C11^(3/2)*parRUP.C2*5000.^3 & Mo <parRUP.mu*parRUP.C10^(3/2)*parRUP.C2*45000.^(3/2));
        L(ind2)  = (Mo(ind2)./(parRUP.mu*parRUP.C1^(3/2)*parRUP.C2)).^(2/5);
        W(ind2)  = parRUP.C1.*L(ind2).^(2/3);
        D(ind2)  = parRUP.C2*parRUP.C1^(1/2)*L(ind2).^(5/6);
    % if 45 km < L 
    ind3=find(Mo>=parRUP.mu*parRUP.C10^(3/2)*parRUP.C2*45000.^(3/2));
        L(ind3)  = (Mo(ind3)./(parRUP.mu*parRUP.C10^(3/2)*parRUP.C2)).^(2/3);
        W(ind3)  = parRUP.C10.*ones(1,length(L(ind3)));
        D(ind3)  = parRUP.C2*parRUP.C10^(1/2)*L(ind3).^(1/2);    
    
elseif (strcmp(parRUP.type,'Reverse'))
    
    % Leonard (2010)- Reverse fault - General case
    L  = (Mo./(parRUP.mu*parRUP.C1^(3/2)*parRUP.C2)).^(2/(3*(1+parRUP.Beta)));
    W  = parRUP.C1.*L.^parRUP.Beta;
    D  = parRUP.C2*parRUP.C1^(1/2)*L.^(0.5*(1+parRUP.Beta));
end

L=L';
W=W';
D=D';