function Aftershock = f_EQ_BASS(Mainshock)
global parEQ;
%==========================================================================
%                               BASS model
%
% Generate a aftershocks series associated to every mainshock
%
% Output : Aftershock - format: structure
% Input  : Mainshock  - format: structure
%
% Dev: T. Croissant & P. Steer
% Last update: 10/2018
%
% Based on Turcotte et al, 2007, GRL
%==========================================================================
Aftershock.Mw   = [];
Aftershock.x    = [];
Aftershock.y    = [];
Aftershock.z    = [];
Aftershock.g    = [];
Aftershock.t    = [];
Aftershock.dt   = [];
Aftershock.L    = [];
Aftershock.W    = [];
Aftershock.D    = [];

for i = 1:length(Mainshock.Mw)
    % First generation of daughter EQ.
    Nd{1}    = round(10.^(parEQ.b*(Mainshock.Mw(i) - parEQ.Dm - parEQ.mw_min)));
    % Random vectors generation
    rnd.Pm   = rand(Nd{1},1);
    rnd.Pt   = rand(Nd{1},1);
    rnd.Pr   = rand(Nd{1},1);
    % Angular sampling - Ellipsoidal
    alpha    = linspace(0,2*pi,360);
    r        = sqrt(1./((cos(alpha)./Mainshock.W(i)).^2+(sin(alpha)./Mainshock.L(i)).^2));
    r        = r./sum(r);
    rnd.Pphi = randsample(alpha,Nd{1},'true',r)';
    
    % Magnitude, Occurrence Time, Coord and Generation number of first generation
    mw_d{1}  = (log10(rnd.Pm) - parEQ.b * parEQ.mw_min)/-parEQ.b;
    r_d{1}   = ((1 ./ rnd.Pr).^(1/(parEQ.q-1)) - 1).*parEQ.d.*10.^(0.5*Mainshock.Mw(i));
    x_d{1}   = Mainshock.x(i) + r_d{1} .* cos(rnd.Pphi);
    y_d{1}   = Mainshock.y(i) + r_d{1} .* sin(rnd.Pphi);
    z_d{1}   = Mainshock.z(i) + zeros(Nd{1},1);
    g_d{1}   = Mainshock.g(i) + ones(Nd{1},1);
    t_inc    = parEQ.c.*((1 ./ rnd.Pt).^(1/(parEQ.p-1)) - 1);
    t_d{1}   = Mainshock.t(i) + t_inc;
    dt_d{1}  = Mainshock.dt(i) + t_inc;
    
    % Compute the rupture geometry of the Earthquake (L, W, D) - Leonard 2010
    [L_d{1},W_d{1},D_d{1}] = fRuptureGeometry(mw_d{1});
    
    % Variable initialization
    Aftershock_tmp.Mw   = [mw_d{1}];
    Aftershock_tmp.x    = [x_d{1}];
    Aftershock_tmp.y    = [y_d{1}];
    Aftershock_tmp.z    = [z_d{1}];
    Aftershock_tmp.g    = [g_d{1}];
    Aftershock_tmp.t    = [t_d{1}];
    Aftershock_tmp.dt   = [dt_d{1}];
    Aftershock_tmp.L    = [L_d{1}];
    Aftershock_tmp.W    = [W_d{1}];
    Aftershock_tmp.D    = [D_d{1}];
    
    inc = 1;
    while Nd{inc} > 0
        
        % Number of aftershocks of the inc-generation
        Nd_tmp{inc+1} = round(10.^(parEQ.b*(mw_d{inc} - parEQ.Dm - parEQ.mw_min)));
        Nd{inc+1}     = sum(Nd_tmp{inc+1});
        indNd         = propIndex(Nd_tmp{inc+1});
        
        if isempty(indNd)==0
            
            % Random vectors generation
            rnd.Pm      = rand(Nd{inc+1},1);
            rnd.Pt      = rand(Nd{inc+1},1);
            rnd.Pr      = rand(Nd{inc+1},1);
            % Angular sampling - ellipsoidal
            alpha = repmat(linspace(0,2*pi,360),Nd{inc},1);
            Wmat  = repmat(W_d{inc},1,360);Lmat=repmat(L_d{inc},1,360);
            % matrix of probablity of ellispoidal distance
            rmat  = sqrt(1./((cos(alpha)./Wmat).^2+(sin(alpha)./Lmat).^2));
            rmat  = rmat./repmat(sum(rmat,2),1,360);
            Ntosample = zeros(Nd{inc},1);
            Ntosample(1:max(indNd)) = accumarray(indNd,1);
            indtemp   = find(Ntosample>0);
            l1 = 1;
            rnd.Pphi  = zeros(Nd{inc+1},1);
            for j = 1:numel(indtemp)
                k  = indtemp(j);
                l2 = l1+Ntosample(k)-1;
                rnd.Pphi(l1:l2) = randsample(alpha(k,:),Ntosample(k),'true',rmat(k,:))';
                l1 = l2+1;
            end
            
            % Magnitude, Occurrence Time and Coord first generation
            mw_d{inc+1} = (log10(rnd.Pm) - parEQ.b * parEQ.mw_min)/-parEQ.b;
            r_d{inc+1}  = ((1 ./ rnd.Pr).^(1/(parEQ.q-1)) - 1).* parEQ.d.*10.^(0.5*mw_d{inc+1});
            x_d{inc+1}  = x_d{inc}(indNd) + r_d{inc+1} .* cos(rnd.Pphi);
            y_d{inc+1}  = y_d{inc}(indNd) + r_d{inc+1} .* sin(rnd.Pphi);
            z_d{inc+1}  = z_d{inc}(indNd) + zeros(Nd{inc+1},1);
            g_d{inc+1}  = g_d{inc}(indNd) + ones(Nd{inc+1},1);
            t_inc       = parEQ.c.*((1 ./ rnd.Pt).^(1/(parEQ.p-1)) - 1);
            t_d{inc+1}  = t_d{inc}(indNd) + t_inc;
            dt_d{inc+1} = dt_d{inc}(indNd)+ t_inc;
            
            % Compute the rupture geometry of the Earthquake (L, W, D) - Leonard 2010
            [L_d{inc+1},W_d{inc+1},D_d{inc+1}] = fRuptureGeometry(mw_d{inc+1});
            
            % Variable update [push-back]
            Aftershock_tmp.Mw = [Aftershock_tmp.Mw; mw_d{inc+1}];
            Aftershock_tmp.x  = [Aftershock_tmp.x;  x_d{inc+1}];
            Aftershock_tmp.y  = [Aftershock_tmp.y;  y_d{inc+1}];
            Aftershock_tmp.z  = [Aftershock_tmp.z;  z_d{inc+1}];
            Aftershock_tmp.g  = [Aftershock_tmp.g;  g_d{inc+1}];
            Aftershock_tmp.t  = [Aftershock_tmp.t;  t_d{inc+1}];
            Aftershock_tmp.dt = [Aftershock_tmp.dt; dt_d{inc+1}];
            Aftershock_tmp.L  = [Aftershock_tmp.L;  L_d{inc+1}];
            Aftershock_tmp.W  = [Aftershock_tmp.W;  W_d{inc+1}];
            Aftershock_tmp.D  = [Aftershock_tmp.D;  D_d{inc+1}];
            
        end
        inc = inc + 1;
    end
    Aftershock = f_EQ_Conc(Aftershock,Aftershock_tmp);
end
