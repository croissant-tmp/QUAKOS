function Landslide = f_LS_Init()

% Landslides coordinates and elevation
Landslide.x  = [];                                                         % x-coord
Landslide.y  = [];                                                         % y-coord
Landslide.IX = [];                                                         % index-coord
Landslide.z  = [];                                                         % Elevation

% Landslides time properties
Landslide.t  = [];                                                         % Landslide occurence time
Landslide.dt = [];                                                         % Landslide time compared to mainshock
Landslide.g  = [];                                                         % Landslide generation compared to mainshock

% Landslides properties
Landslide.A  = [];                                                         % Landslide area
Landslide.V  = [];                                                         % Landslide volume
Landslide.conINI = [];                                                     % Landslide catchment
Landslide.DB  = [];                                                        % Landslide catchment

% Landslides cluster properties
Landslide.Atot   = 0;
Landslide.Alstot = 0;                                                      % Total area of the landslides
Landslide.Vlstot = 0;                                                      % Total Volume of the landslides
Landslide.N      = 0;                                                      % Total number of  landslides