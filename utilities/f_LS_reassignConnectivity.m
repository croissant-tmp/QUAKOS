function LandslideDB = f_LS_reassignConnectivity(LandslideDB)
%==========================================================================
% Assign the connectivity status to landslide presenting the lowest
% connection time
%==========================================================================
nlsDB  = unique(LandslideDB.DB);

for n = 1:length(nlsDB)
    
    ixLS_db = find(LandslideDB.DB==nlsDB(n));
    
    if sum(LandslideDB.conINI(ixLS_db)) >= 1
        
        ixLS_con      = find(LandslideDB.conINI(ixLS_db)==1) + min(ixLS_db)-1;
        [~,ixLS_Tcon] = sort(LandslideDB.tcon(ixLS_db));
        
        ixLS_Tcon = ixLS_Tcon + min(ixLS_db) - 1;
        ixLS_Tcon = ixLS_Tcon(1:length(ixLS_con));
        
        ixLS_Tcon1 = setdiff(ixLS_Tcon,ixLS_con);
        ixLS_con1 = setdiff(ixLS_con,ixLS_Tcon);
        
        Atmp1 = LandslideDB.A(ixLS_con1);
        Atmp2 = LandslideDB.A(ixLS_Tcon1);
        Vtmp1 = LandslideDB.V(ixLS_con1);
        Vtmp2 = LandslideDB.V(ixLS_Tcon1);
        
        LandslideDB.conINI(ixLS_Tcon1) = 1;
        LandslideDB.conINI(ixLS_con1) = 0;
        LandslideDB.A(ixLS_Tcon1) = Atmp1;
        LandslideDB.A(ixLS_con1) = Atmp2;
        LandslideDB.V(ixLS_Tcon1) = Vtmp1;
        LandslideDB.V(ixLS_con1) = Vtmp2;
        
    end
end





