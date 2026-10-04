function P = default_params()
% DEFAULT_PARAMS  Baseline parameter set for the kirigami wingbox study.
%
%   P = default_params() returns a struct with every geometric, material and
%   loading parameter used by the simulation. ALL UNITS ARE SI.
%   Override any field to run your own design, e.g.
%       P = default_params(); P.chord = 0.45; P.n_cells = 6;
%
% ---------------------------------------------------------------------
% GEOMETRIC PARAMETERS  (metres)
% ---------------------------------------------------------------------
P.chord    = 0.30;   % wingbox chord c                    [m]
P.height   = 0.08;   % wingbox height h (web height)      [m]
P.n_cells  = 4;      % number of cells N along the chord  [-]
P.cell_len = 0.075;  % nominal spanwise cell length       [m]  (info only:
                     %  chord/N is what the 2-D section actually uses)
P.seam_pos = 0.40;   % seam location gamma in each cell: seam at
                     %  x_left_web + gamma*cell_wide, 0<gamma<1  [-]
                     %  gamma ~= 0.5 makes the section asymmetric -> Ixy ~= 0
P.flap_ratio = 0.35; % seam-lift amplitude factor beta [-]
                     %  (kept for future variants; lift uses d*sin(theta))

% ---------------------------------------------------------------------
% WALL THICKNESSES  (metres)  - thin-walled idealisation
% ---------------------------------------------------------------------
P.t_up  = 1.0e-3;    % upper skin thickness               [m]
P.t_lo  = 1.0e-3;    % lower skin thickness               [m]
P.t_web = 1.2e-3;    % side-wall (web/rib) thickness      [m]

% ---------------------------------------------------------------------
% MORPHING SWEEP
% ---------------------------------------------------------------------
P.theta_min  = 0;    % closed configuration               [deg]
P.theta_max  = 90;   % fully deployed configuration       [deg]
P.theta_step = 15;   % sweep resolution                   [deg]

% ---------------------------------------------------------------------
% MATERIAL (aluminium 2024-T3 typical wingbox values)
% ---------------------------------------------------------------------
P.E  = 72e9;         % Young's modulus                    [Pa]
P.G  = 27e9;         % shear modulus                      [Pa]
P.nu = 0.33;         % Poisson ratio (informational, E & G are used)

% ---------------------------------------------------------------------
% LOADING
% ---------------------------------------------------------------------
P.Vy = 1000;         % transverse shear force (upward)    [N]
P.Vx = 0;            % chordwise shear force              [N]
P.load_pt = [0.30, 0.04];  % line of action of the applied force [m,m]
                          %  used for the eccentricity / torque study
P.L_ref   = 0.50;    % reference beam length for twist theta_t = T*L/Kt [m]

% ---------------------------------------------------------------------
% VISUALISATION / EXPORT
% ---------------------------------------------------------------------
P.save_figs = false; % save PNGs of every figure
P.fig_dir   = 'figures';
P.verbose   = true;  % print intermediate calculations to the console
end
