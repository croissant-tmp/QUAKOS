function Qt = f_QT(dA,Qt_tp,Qt_ts,runT,p_effv)
%==========================================================================
% Compute alongstream transport capacity
%
% dev: T. Croissant
% Last update: 04/2017
%==========================================================================

% dA     = LandslideDB.dA;
% Qt_tp  = 'stoch';
% Qt_ts  = 1;

% Parameters
kwn    = 0.005;                                                            % n. width index
ksn    = 180;                                                              % n. steepness index
D50    = 0.3;                                                              % Median grain size
tau_c  = 0.035*9.81*1700*D50;                                              % critical shear stress
p_eff  = p_effv/(365*86400);                                                    % mean annual runoff
k      = 1.0;                                                              % Q variability
W      = kwn * dA.^0.5;                                                    % alongS width
S      = ksn * dA.^-0.45;                                                  % alongS slope
Qm     = p_eff .* dA;                                                      % alongS mean Q

%==========================================================================
% Along stream geometry and discharge
Qn     = logspace(-2,log10(200),2E4);
dQn    = diff(Qn); dQn = [dQn dQn(end)];
pdf_Q  = (k^(k+1))./gamma(k+1) .* exp(-k./Qn) .* Qn.^-(2+k);
cdf_Q  = cumsum(pdf_Q.*dQn);
recT   = (1-cumsum(pdf_Q.*dQn)).^-1;
Qn_1yr = Qn(find(recT>365,1));

%==========================================================================
% Compute daily bedload transport capacity

switch Qt_tp
    
    case 'eff'   %=========================================================
        if Qt_ts == 1
            Q     = Qm .* Qn_1yr;
            tau   = 9810 * (0.035.*Q./W).^0.6 .* S.^0.7;
            tau(tau<tau_c) = tau_c;
            Qt   = 86400 .* W .* 1.5E-5 .* (tau - tau_c).^1.5;
        else
            Q     = Qm .* Qn_1yr;
            for i = 1:length(Qm)
                W{i} = f(t);
                tau   = 9810 * (0.035.*Q(i)./W{i}).^0.6 .* S(i).^0.7;
                tau(tau<tau_c) = tau_c;
                Qt{i}   = 86400 .* W{i} .* 1.5E-5 .* (tau - tau_c).^1.5;
            end
        end
        
    case 'stoch' %=========================================================
        
        if Qt_ts == 1 % Stochatic-based annual transport capacity
            nQ = max(cdf_Q)*rand(1,runT*365);
            qq = interp1(cdf_Q, Qn, nQ);
            for i = 1:length(Qm)
                Q   = Qm(i) .* qq;
                tau   = 9810 * (0.035.*Q./W(i)).^0.6 .* S(i).^0.7;
                tau(tau<tau_c) = tau_c;
                Qt_d{i} = 86400 .* W(i) .* 1.5E-5 .* (tau - tau_c).^1.5;
                inc = 0;
                for j = 1:runT
                    Qt{i}(j) = sum(Qt_d{i}(inc+1:365*j));
                    inc = inc + 365;
                end
            end
  
        else % Long-term bedload transport capacity (Lague et al, [2005])
            for i = 1:length(Qm)
                Q   = Qm(i) .* Qn;
                tau   = 9810 * (0.035.*Q./W(i)).^0.6 .* S(i).^0.7;
                tau(tau<tau_c) = tau_c;
                Qt_eb{i} = 86400.* W(i) .* 1.5E-5 .* (tau - tau_c).^1.5;
                Qt(i)    = 365*sum(Qt_eb{i}.*pdf_Q.*dQn);  
            end
        end
        
    otherwise
        
        disp('Error check Qt type')
end








