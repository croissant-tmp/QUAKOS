function Landslide = f_LS_Conc(Landslide,tmp)
%==========================================================================
% Landslides coordinates and elevation
Landslide.x  = [Landslide.x tmp.x];
Landslide.y  = [Landslide.y tmp.y];
Landslide.IX = [Landslide.IX tmp.IX];
Landslide.z  = [Landslide.z tmp.z];

% Landslides time properties
Landslide.t  = [Landslide.t tmp.t];
Landslide.dt = [Landslide.dt tmp.dt];
Landslide.g  = [Landslide.g tmp.g];

% Landslides properties
Landslide.A      = [Landslide.A tmp.A];
Landslide.V      = [Landslide.V tmp.V];
Landslide.DB     = [Landslide.DB tmp.DB];
Landslide.conINI = [Landslide.conINI tmp.conINI];

% Landslides cluster properties
Landslide.Atot   = [Landslide.Atot tmp.Atot];
Landslide.Alstot = [Landslide.Alstot tmp.Alstot];
Landslide.Vlstot = [Landslide.Vlstot tmp.Vlstot];
Landslide.N      = [Landslide.N tmp.N];

    