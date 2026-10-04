function [R, T] = morphing_simulation(P)
% MORPHING_SIMULATION  Parametric sweep of the kirigami wingbox over theta.
%
%   [R, T] = morphing_simulation(P) loops over theta_min:theta_step:theta_max
%   and for each morphing configuration computes:
%       - cross-section properties   (section_properties.m)
%       - shear-flow distribution    (shear_flow.m, applied at the CENTROID)
%       - shear-center location      (shear_center.m)
%       - torsional properties       (bredt_torsion.m)
%   and finally the bending-torsion coupling for a force applied along the
%   VERTICAL LINE x = P.load_pt(1):
%       T(theta) = V * e(theta),      e = x_loadline - xSC
%       theta_t  = T / Kt             (uniform-torsion estimate, L_ref)
%
%   OUTPUTS:
%       R  struct of arrays, one entry per theta (centroids, Ix..., SC...)
%       T  MATLAB table with one row per morphing angle - the "results table"
%          required by the project statement.
%
%   UNITS: SI throughout.  With verbose=true the three validation
%   configurations (closed / intermediate / deployed) print their full
%   intermediate calculations for the presentation.

if nargin < 1, P = default_params(); end
verbose = isfield(P, 'verbose') && P.verbose;

thetas = P.theta_min : P.theta_step : P.theta_max;
nT = numel(thetas);

% preallocate result arrays
R = struct();
R.theta   = thetas(:);
R.xbar    = zeros(nT,1);  R.ybar   = zeros(nT,1);
R.A       = zeros(nT,1);
R.Ix      = zeros(nT,1);  R.Iy     = zeros(nT,1);  R.Ixy = zeros(nT,1);
R.xSC     = zeros(nT,1);  R.ySC    = zeros(nT,1);
R.ex      = zeros(nT,1);  R.ey     = zeros(nT,1);
R.J       = zeros(nT,1);  R.Kt     = zeros(nT,1);
R.Torque  = zeros(nT,1);  R.theta_t = zeros(nT,1);
R.G_geo   = cell(nT,1);            % keep geometry of each state
R.shear   = cell(nT,1);            % keep shear-flow solution of each state

fprintf('\n===== MORPHING SWEEP: %d configurations, theta = %g:%g:%g deg =====\n', ...
        nT, P.theta_min, P.theta_step, P.theta_max);
fprintf('%6s %10s %10s %11s %11s %11s %10s %10s %11s %11s\n', ...
        'theta', 'xbar[m]', 'ybar[m]', 'Ix[m^4]', 'Iy[m^4]', 'Ixy[m^4]', ...
        'xSC[m]', 'ySC[m]', 'e[m]', 'theta_t[deg]');

for it = 1:nT
    th = thetas(it);

    % --- geometry + section properties -----------------------------------
    G  = kirigami_geometry(P, th);
    SP = section_properties(G.E);

    % --- shear flow for the reference load (applied at the CENTROID) -----
    sh = shear_flow(G, SP, P.Vy, P.Vx);

    % --- shear center -----------------------------------------------------
    [ex, ey, SC] = shear_center(G, SP, abs(P.Vy) + abs(P.Vx) + eps);
    % scale-normalised: V only appears as a divisor, any V > 0 works

    % --- torsion ----------------------------------------------------------
    [J, Kt, ~] = bredt_torsion(G, SP, P);

    % --- bending-torsion coupling for the eccentric load ------------------
    e_load = P.load_pt(1) - SC.x;          % signed eccentricity [m]
    Torque = P.Vy * e_load;                % T = V*e  [N*m]
    theta_t = Torque / Kt;                 % twist [rad]  (uniform torsion)

    % --- store ------------------------------------------------------------
    R.xbar(it) = SP.xbar;   R.ybar(it) = SP.ybar;   R.A(it) = SP.A;
    R.Ix(it) = SP.Ix;  R.Iy(it) = SP.Iy;  R.Ixy(it) = SP.Ixy;
    R.xSC(it) = SC.x;  R.ySC(it) = SC.y;
    R.ex(it) = ex;     R.ey(it) = ey;
    R.J(it) = J;       R.Kt(it) = Kt;
    R.Torque(it) = Torque;  R.theta_t(it) = theta_t;
    R.G_geo{it} = G;       R.shear{it} = sh;

    fprintf('%6.1f %10.4f %10.4f %11.4e %11.4e %11.4e %10.4f %10.4f %11.4f %11.4f\n', ...
            th, SP.xbar, SP.ybar, SP.Ix, SP.Iy, SP.Ixy, SC.x, SC.y, e_load, theta_t*180/pi);

    % --- verbose dumps for the three validation configurations ------------
    if verbose
        if abs(th - P.theta_min) < 1e-9
            fprintf('\n########## CONFIG 1: CLOSED (theta = %g deg) ##########\n', th);
            dump_config(P, th);
        elseif abs(th - (P.theta_min+P.theta_max)/2) < P.theta_step/2
            fprintf('\n########## CONFIG 2: INTERMEDIATE (theta = %g deg) ##########\n', th);
            dump_config(P, th);
        elseif abs(th - P.theta_max) < 1e-9
            fprintf('\n########## CONFIG 3: DEPLOYED (theta = %g deg) ##########\n', th);
            dump_config(P, th);
        end
    end
end

% ---- results table -------------------------------------------------------
T = table(R.theta, R.xbar, R.ybar, R.Ix, R.Iy, R.Ixy, ...
          R.xSC, R.ySC, R.ex, R.ey, R.J, R.Kt, R.Torque, R.theta_t, ...
          'VariableNames', {'theta_deg', 'xbar_m', 'ybar_m', 'Ix_m4', 'Iy_m4', ...
          'Ixy_m4', 'xSC_m', 'ySC_m', 'ex_m', 'ey_m', 'J_m4', 'Kt_Nm_per_rad', ...
          'Torque_Nm', 'twist_deg'});
end

% ============================================================================
function dump_config(P, th)
% Full intermediate-calculations dump for one configuration (presentation aid).
G  = kirigami_geometry(P, th);
SP = section_properties(G.E, true);
sh = shear_flow(G, SP, P.Vy, P.Vx, true);
shear_center(G, SP, 1000, true);
bredt_torsion(G, SP, P, true);

% per-segment shear-flow table (q at both ends of each wall segment)
fprintf('  --- q(s) per wall segment [N/m] (stored n1->n2 direction) ---\n');
fprintf('  elem      x1        y1        x2        y2        q1      q2\n');
for e = 1:size(G.E,1)
    fprintf('  %4d  %8.4f  %8.4f  %8.4f  %8.4f  %8.1f  %8.1f\n', ...
            e, G.E(e,1), G.E(e,2), G.E(e,3), G.E(e,4), sh.q1(e), sh.q2(e));
end
fprintf('  resultant: Fx = %.3f N, Fy = %.3f N\n\n', sh.Fx, sh.Fy);
end
