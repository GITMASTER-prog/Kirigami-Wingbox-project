function figs = visualization(R, P)
% VISUALIZATION  All report figures for the kirigami wingbox study.
%
%   figs = visualization(R, P) creates the 10 required figures from sweep
%   results R (morphing_simulation) and parameters P.  Returns the figure
%   handles.  Each figure is report-ready (labelled axes, legend, grid).

if nargin < 2, P = default_params(); end
figs = struct();

th = R.theta;  n = numel(th);
cmap = turbo(n);   % colour per morphing state (Fig 9 trajectory)

% ---- FIG 1: geometry at several morphing angles --------------------------
figs.f1 = figure('Name','Fig1 Geometry','Color','w','Position',[60 60 1200 420]);
N = P.n_cells;
pick = max(1, round(linspace(1, n, min(n,5))));
for jj = 1:numel(pick)
    subplot(1, numel(pick), jj); hold on; axis equal;
    G = R.G_geo{pick(jj)};
    draw_geometry(G);
    title(sprintf('\\theta = %g\\circ', R.theta(pick(jj))));
    xlabel('x [m]'); ylabel('y [m]'); grid on;
end
sgtitle('Kirigami wingbox cross-section vs morphing angle');

% ---- FIG 2: centroid migration -------------------------------------------
figs.f2 = figure('Name','Fig2 Centroid','Color','w');
plot(th, R.xbar*1e3, 'o-', 'LineWidth',1.4); hold on;
plot(th, R.ybar*1e3, 's-', 'LineWidth',1.4);
grid on; xlabel('\theta [deg]');
legend('x_{bar} [mm]', 'y_{bar} [mm]', 'Location','best');
title('Centroid location vs morphing angle');

% ---- FIG 3: shear-center components vs theta ------------------------------
figs.f3 = figure('Name','Fig3 SC location','Color','w');
plot(th, R.xSC*1e3, 'o-', 'LineWidth',1.4); hold on;
plot(th, R.ySC*1e3, 's-', 'LineWidth',1.4);
plot(th, R.xbar*1e3, '--', 'LineWidth',1.0);
plot(th, R.ybar*1e3, ':', 'LineWidth',1.4);
grid on; xlabel('\theta [deg]'); ylabel('[mm]');
legend('x_{SC}', 'y_{SC}', 'x_{bar}', 'y_{bar}', 'Location','best');
title('Shear-center and centroid vs morphing angle');

% ---- FIG 4: eccentricity vs theta -----------------------------------------
figs.f4 = figure('Name','Fig4 Eccentricity','Color','w');
yyaxis left;  plot(th, R.ex*1e3, 'o-', 'LineWidth',1.4); ylabel('e_x [mm]');
yyaxis right; plot(th, R.ey*1e3, 's-', 'LineWidth',1.4); ylabel('e_y [mm]');
grid on; xlabel('\theta [deg]');
title('Shear-center eccentricity from centroid vs morphing angle');
legend('e_x', 'e_y', 'Location','best');

% ---- FIG 5-7: section properties vs theta ---------------------------------
figs.f5 = figure('Name','Fig5 Ix','Color','w');
plot(th, R.Ix*1e8, 'o-', 'LineWidth',1.4); grid on;
xlabel('\theta [deg]'); ylabel('I_x [10^{-8} m^4]');
title('Second moment of area I_x vs morphing angle');

figs.f6 = figure('Name','Fig6 Iy','Color','w');
plot(th, R.Iy*1e8, 's-', 'LineWidth',1.4); grid on;
xlabel('\theta [deg]'); ylabel('I_y [10^{-8} m^4]');
title('Second moment of area I_y vs morphing angle');

figs.f7 = figure('Name','Fig7 Ixy','Color','w');
plot(th, R.Ixy*1e8, '^-', 'LineWidth',1.4); grid on;
xlabel('\theta [deg]'); ylabel('I_{xy} [10^{-8} m^4]');
title('Product of inertia I_{xy} vs morphing angle');

% ---- FIG 8: shear-flow distribution at selected configurations -------------
figs.f8 = figure('Name','Fig8 Shear flow','Color','w','Position',[80 80 1200 420]);
sel = unique([1, round(n/2), n]);
for jj = 1:numel(sel)
    subplot(1, numel(sel), jj); hold on; axis equal;
    plot_shear(R.G_geo{sel(jj)}, R.shear{sel(jj)}, ...
               R.xbar(sel(jj)), R.ybar(sel(jj)), R.xSC(sel(jj)), R.ySC(sel(jj)));
    title(sprintf('q(s), \\theta = %g\\circ', R.theta(sel(jj))));
end
sgtitle(sprintf('Shear flow for V_y = %g N (applied at centroid)', P.Vy));

% ---- FIG 9: shear-center trajectory in the cross-section plane -------------
figs.f9 = figure('Name','Fig9 SC trajectory','Color','w'); hold on; axis equal;
% faint geometry envelope for context
for it = 1:n
    G = R.G_geo{it};
    for e = 1:size(G.E,1)
        plot([G.E(e,1) G.E(e,3)], [G.E(e,2) G.E(e,4)], '-', ...
             'Color', [0.85 0.85 0.85], 'HandleVisibility','off');
    end
end
scatter(R.xSC*1e3, R.ySC*1e3, 42, th, 'filled');
plot(R.xSC*1e3, R.ySC*1e3, '-', 'Color', [0.3 0.3 0.3]);
plot(R.xbar*1e3, R.ybar*1e3, 'k--');      % centroid path for contrast
colormap(gca, cmap); cb = colorbar; cb.Label.String = '\theta [deg]';
plot(R.xSC(1)*1e3, R.ySC(1)*1e3, 'p', 'MarkerSize',16, ...
     'MarkerFaceColor','y', 'MarkerEdgeColor','k');
text(R.xSC(1)*1e3, R.ySC(1)*1e3, '  closed', 'FontWeight','bold');
plot(R.xSC(end)*1e3, R.ySC(end)*1e3, 'p', 'MarkerSize',16, ...
     'MarkerFaceColor','m', 'MarkerEdgeColor','k');
text(R.xSC(end)*1e3, R.ySC(end)*1e3, '  deployed', 'FontWeight','bold');
grid on; xlabel('x [mm]'); ylabel('y [mm]');
title('Shear-center migration (closed \rightarrow deployed)');
legend({'SC trajectory','centroid path'}, 'Location','best');

% ---- FIG 10: 3D morphing + SC migration ------------------------------------
figs.f10 = figure('Name','Fig10 3D morphing','Color','w');
draw_morph_3d(R, P);

% ---- coupling figures (Part 7) ---------------------------------------------
figs.f11 = figure('Name','Fig11 Torque','Color','w');
yyaxis left;  plot(th, R.Torque, 'o-', 'LineWidth',1.4); ylabel('T = V\cdote [N*m]');
yyaxis right; plot(th, R.theta_t*180/pi, 's-', 'LineWidth',1.4);
ylabel('\theta_t [deg]');
grid on; xlabel('\theta [deg]');
title('Bending-torsion coupling: torque and twist vs morphing angle');
legend('torque T [N*m]', 'twist \theta_t [deg]', 'Location','best');

end

% ============================================================================
function draw_geometry(G)
% Draw all wall segments + nodes.
for e = 1:size(G.E,1)
    plot([G.E(e,1) G.E(e,3)], [G.E(e,2) G.E(e,4)], 'k-', 'LineWidth', 1.6);
end
plot(G.nodes(:,1), G.nodes(:,2), 'r.', 'MarkerSize', 8);
end

% ============================================================================
function plot_shear(G, sh, xbar, ybar, xSC, ySC)
% Shear-flow plot: colour-coded segments + direction arrows + key markers.
qmid = (sh.q1 + sh.q2)/2;
qmax = max(abs(qmid));
for e = 1:size(G.E,1)
    x1 = G.E(e,1); y1 = G.E(e,2); x2 = G.E(e,3); y2 = G.E(e,4);
    cval = abs(qmid(e)) / max(qmax, eps);
    col = [1-cval, 1-0.5*cval, 1-0.15*cval];     % red-hot for high |q|
    if cval < 0.15, col = [0.2 0.6 1]; end       % blue for low |q|
    plot([x1 x2], [y1 y2], '-', 'Color', col, 'LineWidth', 1.2 + 5*cval);
    % direction arrow at the midpoint (physical flow direction)
    xmv = (x1+x2)/2; ymv = (y1+y2)/2;
    dx = x2-x1; dy = y2-y1; Ln = hypot(dx,dy);
    u = dx/Ln; v = dy/Ln;
    if qmid(e) < 0, u = -u; v = -v; end
    quiver(xmv, ymv, 0.006*u, 0.006*v, 0, 'Color', col, 'LineWidth', 1.2, ...
           'MaxHeadSize', 2);
end
plot(xbar, ybar, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
text(xbar, ybar, '  C', 'FontWeight','bold');
plot(xSC, ySC, 'bp', 'MarkerSize', 13, 'MarkerFaceColor','g', 'MarkerEdgeColor','k');
text(xSC, ySC, '  SC', 'FontWeight','bold');
grid on;
end

% ============================================================================
function draw_morph_3d(R, P)
% 3D view: each theta gets a spanwise slice; SC path drawn through 3D space.
n = numel(R.theta);
for it = 1:n
    z0 = (it-1) * P.cell_len;
    G = R.G_geo{it};
    for e = 1:size(G.E,1)
        plot3([G.E(e,1) G.E(e,3)], [G.E(e,2) G.E(e,4)], [z0 z0], ...
              'k-', 'LineWidth', 0.8); %#ok<AGROW>
    end
    plot3(R.xSC(it)*ones(1,2), R.ySC(it)*ones(1,2), [z0 z0]+0.004, ...
          'r.', 'MarkerSize', 14); %#ok<AGROW>
end
hold on;
plot3(R.xSC, R.ySC, (0:n-1)*P.cell_len + 0.004, 'r-o', 'LineWidth', 1.4, ...
      'MarkerFaceColor','r');
grid on; view(3);
xlabel('x [m]'); ylabel('y [m]'); zlabel('spanwise z [m]');
title('Wingbox morphing (spanwise stack) and shear-center migration');
axis equal;
end
