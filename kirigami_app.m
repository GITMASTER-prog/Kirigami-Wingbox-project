function app = kirigami_app()
% KIRIGAMI_APP  Simple interactive UI for the kirigami wingbox study.
%
%   app = kirigami_app() opens a uifigure with all geometric / material /
%   loading inputs, the six project buttons, and live outputs.  It reuses
%   exactly the same solver functions as the batch run (main.m), so UI and
%   script results are identical.
%
%   Buttons:
%     GENERATE GEOMETRY        - draw the section at the current theta
%     CALCULATE SHEAR CENTER   - full section + SC analysis at current theta
%     RUN MORPHING SIMULATION  - sweep theta_min..theta_max, plot migration
%     SHOW SHEAR FLOW          - colour-coded q(s) at the current theta
%     START ANIMATION          - morphing animation
%     EXPORT RESULTS           - write CSV + MAT of the last sweep
%
%   State is kept in the returned struct so  app.parameters.chord = 0.5
%   style tweaks (or direct edits in the UI) work between runs.

app = struct();
app.fig = uifigure('Name', 'Kirigami Wingbox - Shear Center Simulator', ...
                   'Position', [80 80 1180 660]);

% ---------------- left panel: inputs -------------------------------------
pL = uipanel(app.fig, 'Title', 'Inputs (SI units)', 'FontWeight','bold', ...
             'Position', [14 14 300 630]);
labels = {'Wing chord c [m]', 'Wingbox height h [m]', 'Number of cells N [-]', ...
          'Cell length [m]', 'Seam position gamma [-]', ...
          'Upper skin t_up [mm]', 'Lower skin t_lo [mm]', 'Web t_web [mm]', ...
          'Morphing angle theta [deg]', 'Shear force Vy [N]', ...
          'theta_min [deg]', 'theta_max [deg]', 'theta_step [deg]'};
defs   = {0.30, 0.08, 4, 0.075, 0.40, 1.0, 1.0, 1.2, 45, 1000, 0, 90, 15};
app.edits = cell(numel(labels),1);
for k = 1:numel(labels)
    uilabel(pL, 'Text', labels{k}, 'Position', [12, 596-44*k+30, 175, 22]);
    app.edits{k} = uieditfield(pL, 'numeric', ...
        'Position', [195, 596-44*k+30, 88, 24], 'Value', defs{k});
end

% ---------------- middle: buttons + outputs -------------------------------
pB = uipanel(app.fig, 'Title', 'Actions', 'FontWeight', 'bold', ...
             'Position', [322 430 300 214]);
btn = {'GENERATE GEOMETRY', 'CALCULATE SHEAR CENTER', 'SHOW SHEAR FLOW', ...
       'RUN MORPHING SIMULATION', 'START ANIMATION', 'EXPORT RESULTS'};
app.buttons = gobjects(numel(btn),1);
for k = 1:numel(btn)
    app.buttons(k) = uibutton(pB, 'push', 'Text', btn{k}, ...
        'Position', [16, 214-34*k-4, 268, 27], ...
        'FontWeight', 'bold', 'BackgroundColor', [0.94 0.94 0.98]);
end

pO = uipanel(app.fig, 'Title', 'Outputs (current theta)', ...
             'FontWeight', 'bold', 'Position', [322 14 300 408]);
out_lab = {'A [cm^2]','xbar [mm]','ybar [mm]','Ix [cm^4]','Iy [cm^4]', ...
           'Ixy [cm^4]','x_SC [mm]','y_SC [mm]','e_x [mm]','e_y [mm]', ...
           'T = Vy*e [N*m]','Kt [N*m/rad]','twist [deg]','q_max [N/m]'};
app.out = uilabel(pO, 'Text', '  (press CALCULATE SHEAR CENTER)', ...
    'Position', [12 20 270 370], 'FontSize', 12, ...
    'HorizontalAlignment','left', 'VerticalAlignment','top');

% ---------------- right: main axes ----------------------------------------
app.ax = uiaxes(app.fig, 'Position', [636 180 530 460]);
title(app.ax, 'Cross-section'); axis(app.ax, 'equal'); grid(app.ax, 'on');
app.ax2 = uiaxes(app.fig, 'Position', [636 14 530 158]);
title(app.ax2, 'SC migration'); grid(app.ax2, 'on');

% ---------------- state ---------------------------------------------------
app.parameters = default_params();
app.results = [];
app.sweep_table = [];

% ---------------- callbacks -----------------------------------------------
app.buttons(1).ButtonPushedFcn = @(src,ev) cb_generate(app);
app.buttons(2).ButtonPushedFcn = @(src,ev) cb_sc(app);
app.buttons(3).ButtonPushedFcn = @(src,ev) cb_flow(app);
app.buttons(4).ButtonPushedFcn = @(src,ev) cb_sweep(app);
app.buttons(5).ButtonPushedFcn = @(src,ev) cb_anim(app);
app.buttons(6).ButtonPushedFcn = @(src,ev) cb_export(app);
end

% ============================================================================
function P = read_params(app)
P = default_params();
P.chord    = app.edits{1}.Value;
P.height   = app.edits{2}.Value;
P.n_cells  = round(app.edits{3}.Value);
P.cell_len = app.edits{4}.Value;
P.seam_pos = app.edits{5}.Value;
P.t_up  = app.edits{6}.Value * 1e-3;     % mm -> m
P.t_lo  = app.edits{7}.Value * 1e-3;
P.t_web = app.edits{8}.Value * 1e-3;
P.theta_deg = app.edits{9}.Value;
P.Vy   = app.edits{10}.Value;
P.theta_min  = app.edits{11}.Value;
P.theta_max  = app.edits{12}.Value;
P.theta_step = app.edits{13}.Value;
P.load_pt = [P.chord, P.height/2];   % load line at the trailing edge
P.verbose = false;
end

% ============================================================================
function G = current_G(app)
P = read_params(app);
G = kirigami_geometry(P, P.theta_deg);
end

% ============================================================================
function cb_generate(app)
cla(app.ax); hold(app.ax, 'on');
G = current_G(app);
for e = 1:size(G.E,1)
    plot(app.ax, [G.E(e,1) G.E(e,3)], [G.E(e,2) G.E(e,4)], 'k-', 'LineWidth', 1.8);
end
plot(app.ax, G.nodes(:,1), G.nodes(:,2), 'r.', 'MarkerSize', 7);
title(app.ax, sprintf('Kirigami geometry, \\theta = %g\\circ', ...
      read_params(app).theta_deg));
axis(app.ax, 'equal'); grid(app.ax, 'on');
end

% ============================================================================
function cb_sc(app)
P = read_params(app);
G  = kirigami_geometry(P, P.theta_deg);
SP = section_properties(G.E);
[ex, ey, SC] = shear_center(G, SP, abs(P.Vy) + 1);
[~, Kt] = bredt_torsion(G, SP, P);
sh = shear_flow(G, SP, P.Vy, 0);
qmax = max(abs((sh.q1 + sh.q2)/2));
e_load = P.load_pt(1) - SC.x;
set(app.out, 'Text', sprintf([ ...
    '  A     = %8.3f cm^2\n  xbar  = %8.2f mm\n  ybar  = %8.2f mm\n' ...
    '  Ix    = %8.3f cm^4\n  Iy    = %8.3f cm^4\n  Ixy   = %8.3f cm^4\n' ...
    '  x_SC  = %8.2f mm\n  y_SC  = %8.2f mm\n  e_x   = %8.2f mm\n' ...
    '  e_y   = %8.2f mm\n  T     = %8.2f N*m\n  Kt    = %8.1f N*m/rad\n' ...
    '  twist = %8.3f deg\n  q_max = %8.0f N/m'], ...
    SP.A*1e4, SP.xbar*1e3, SP.ybar*1e3, SP.Ix*1e8, SP.Iy*1e8, SP.Ixy*1e8, ...
    SC.x*1e3, SC.y*1e3, ex*1e3, ey*1e3, ...
    P.Vy*e_load, Kt, (P.Vy*e_load/Kt)*180/pi, qmax));
% overlay markers on the geometry
cla(app.ax); hold(app.ax, 'on');
for e = 1:size(G.E,1)
    plot(app.ax, [G.E(e,1) G.E(e,3)], [G.E(e,2) G.E(e,4)], 'k-', 'LineWidth', 1.8);
end
plot(app.ax, SP.xbar, SP.ybar, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
plot(app.ax, SC.x, SC.y, 'p', 'MarkerSize', 15, 'MarkerFaceColor', [0.2 0.8 0.3]);
legend(app.ax, {'walls', 'centroid', 'shear center'}, 'Location', 'southwest');
title(app.ax, sprintf('\\theta = %g\\circ: C and SC', P.theta_deg));
axis(app.ax, 'equal'); grid(app.ax, 'on');
end

% ============================================================================
function cb_flow(app)
P = read_params(app);
G  = kirigami_geometry(P, P.theta_deg);
SP = section_properties(G.E);
sh = shear_flow(G, SP, P.Vy, 0);
[~, ~, SC] = shear_center(G, SP, abs(P.Vy)+1);
cla(app.ax); hold(app.ax, 'on');
qmid = (sh.q1 + sh.q2)/2;  qmax = max(abs(qmid)); %#ok<NASGU>
for e = 1:size(G.E,1)
    x1 = G.E(e,1); y1 = G.E(e,2); x2 = G.E(e,3); y2 = G.E(e,4);
    cval = abs(qmid(e))/max(qmax, eps);
    col = [1-cval, 1-0.5*cval, 1-0.15*cval];
    if cval < 0.15, col = [0.2 0.6 1]; end
    plot(app.ax, [x1 x2], [y1 y2], '-', 'Color', col, 'LineWidth', 1.2+5*cval);
    dx = x2-x1; dy = y2-y1; Ln = hypot(dx,dy); u = dx/Ln; v = dy/Ln;
    if qmid(e) < 0, u = -u; v = -v; end
    quiver(app.ax, (x1+x2)/2, (y1+y2)/2, 0.008*u, 0.008*v, 0, ...
           'Color', col, 'LineWidth', 1.1, 'MaxHeadSize', 2);
end
plot(app.ax, SP.xbar, SP.ybar, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
plot(app.ax, SC.x, SC.y, 'p', 'MarkerSize', 15, 'MarkerFaceColor', [0.2 0.8 0.3]);
quiver(app.ax, P.load_pt(1), P.height/2, 0, 0.25*P.height, 0, ...
       'Color', [0 0.3 1], 'LineWidth', 2.4, 'MaxHeadSize', 1.2);
title(app.ax, sprintf('q(s) for V_y = %g N, \\theta = %g\\circ', P.Vy, P.theta_deg));
axis(app.ax, 'equal'); grid(app.ax, 'on');
end

% ============================================================================
function cb_sweep(app)
P = read_params(app);
[R, T] = morphing_simulation(P);
app.results = R;  app.sweep_table = T;
cla(app.ax2); hold(app.ax2, 'on');
plot(app.ax2, R.xSC*1e3, R.ySC*1e3, 'o-', 'LineWidth', 1.4);
plot(app.ax2, R.xbar*1e3, R.ybar*1e3, 'k--');
for k = 1:5:numel(R.theta)
    text(app.ax2, R.xSC(k)*1e3, R.ySC(k)*1e3, ...
         sprintf(' %g\\circ', R.theta(k)), 'FontSize', 8);
end
title(app.ax2, 'SC migration (o) vs centroid path (- -)');
xlabel(app.ax2, 'x [mm]'); ylabel(app.ax2, 'y [mm]');
uialert(app.fig, sprintf('Sweep done: %d configurations.  See command window for the results table.', ...
        numel(R.theta)), 'Morphing simulation', 'Icon', 'success');
end

% ============================================================================
function cb_anim(app)
P = read_params(app);
if ~isempty(app.results)
    kirigami_animation(app.results, P);
else
    [R, ~] = morphing_simulation(P);
    app.results = R;
    kirigami_animation(R, P);
end
end

% ============================================================================
function cb_export(app)
if isempty(app.sweep_table)
    uialert(app.fig, 'Run the morphing simulation first.', 'Export', 'Icon', 'warning');
    return;
end
writetable(app.sweep_table, 'kirigami_app_results.csv');
uialert(app.fig, 'Wrote kirigami_app_results.csv', 'Export', 'Icon', 'success');
end
