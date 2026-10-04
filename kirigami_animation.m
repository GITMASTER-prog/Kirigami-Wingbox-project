function kirigami_animation(R, P)
% KIRIGAMI_ANIMATION  Animate the kirigami wingbox morphing closed -> deployed.
%
%   kirigami_animation(R, P) plays a frame-by-frame animation of the section
%   morphing while showing, at every frame:
%       - the current morphing angle (title + on-figure text)
%       - the centroid (black +)
%       - the shear center (green pentagram, with a trail)
%       - the applied transverse shear force (blue arrow at the load point)
%       - the shear-flow direction on every wall (quiver arrows)
%
%   The sweep in R is interpolated to P.anim_frames smooth frames, so the
%   shear center visibly MIGRATES rather than jumping.  All quantities SI.

if nargin < 2, P = default_params(); end
n_frames = 90;
if isfield(P, 'anim_frames'), n_frames = P.anim_frames; end

thetas  = R.theta;
frames  = min(n_frames, 120);
f_theta = linspace(thetas(1), thetas(end), frames);

% interpolate the SC / centroid paths for smooth motion
xSC  = interp1(thetas, R.xSC,  f_theta, 'pchip');
ySC  = interp1(thetas, R.ySC,  f_theta, 'pchip');
xbar = interp1(thetas, R.xbar, f_theta, 'pchip');
ybar = interp1(thetas, R.ybar, f_theta, 'pchip');

fig = figure('Name', 'Kirigami morphing animation', 'Color', 'w', ...
             'Position', [100 100 900 640]);
ax  = axes(fig); hold(ax, 'on');
txt_box = annotation(fig, 'textbox', [0.72 0.72 0.25 0.22], 'String', '', ...
                     'EdgeColor', 'k', 'BackgroundColor', 'w', ...
                     'Interpreter', 'tex', 'FitBox', 'on');

xL = P.load_pt(1);

for k = 1:frames
    th = f_theta(k);

    % --- exact geometry and shear flow for this angle ---------------------
    G  = kirigami_geometry(P, th);
    SP = section_properties(G.E);
    sh = shear_flow(G, SP, P.Vy, P.Vx);

    cla(ax); hold(ax, 'on');

    % walls (grey)
    for e = 1:size(G.E, 1)
        hl = plot(ax, [G.E(e,1) G.E(e,3)], [G.E(e,2) G.E(e,4)], '-', ...
                  'Color', [0.55 0.55 0.55], 'LineWidth', 2.2);
    end
    % shear-flow arrows (colour = magnitude, direction = physical flow)
    qmid = (sh.q1 + sh.q2)/2;
    qmax = max(abs(qmid));  qmax = max(qmax, eps);
    for e = 1:size(G.E, 1)
        dx = G.E(e,3)-G.E(e,1);  dy = G.E(e,4)-G.E(e,2);
        Ln = hypot(dx, dy);  u = dx/Ln;  v = dy/Ln;
        if qmid(e) < 0, u = -u; v = -v; end
        cval = abs(qmid(e))/qmax;
        hf = quiver(ax, (G.E(e,1)+G.E(e,3))/2, (G.E(e,2)+G.E(e,4))/2, ...
                    0.011*u, 0.011*v, 0, ...
                    'Color', [1-0.7*cval, 1-0.2*cval, 0.15], ...
                    'LineWidth', 1.3, 'MaxHeadSize', 1.5);
    end
    % SC trail up to this frame
    ht = plot(ax, xSC(1:k), ySC(1:k), '-', 'Color', [0.4 0.4 0.4], ...
              'LineWidth', 1.4);
    % centroid and shear center
    hc = plot(ax, xbar(k), ybar(k), '+', 'MarkerSize', 12, ...
              'LineWidth', 2, 'Color', 'k');
    hs = plot(ax, xSC(k), ySC(k), 'p', 'MarkerSize', 16, ...
              'MarkerFaceColor', [0.2 0.8 0.3], 'MarkerEdgeColor', 'k');
    % applied force arrow at the load point
    L0 = 0.18 * P.height;
    hv = quiver(ax, xL, P.height/2, 0, L0, 0, 'Color', [0 0.3 1], ...
                'LineWidth', 2.6, 'MaxHeadSize', 1.2);

    legend(ax, [hl hf ht hc hs hv], {'walls', 'shear-flow q(s)', ...
           'SC trail', 'centroid', 'shear center', 'applied V_y'}, ...
           'Location', 'southwest');

    title(ax, sprintf('\\theta = %5.1f\\circ   (closed \\rightarrow deployed)', th));
    set(txt_box, 'String', sprintf( ...
        '\\theta = %5.1f\\circ\nSC = (%.1f, %.1f) mm\nC  = (%.1f, %.1f) mm\ne = load_x - SC_x = %.2f mm', ...
        th, xSC(k)*1e3, ySC(k)*1e3, xbar(k)*1e3, ybar(k)*1e3, ...
        (xL - xSC(k))*1e3));

    xlabel(ax, 'x [m]'); ylabel(ax, 'y [m]');
    axis(ax, 'equal');
    xlim(ax, [-0.25*P.chord, 1.45*P.chord]);
    ylim(ax, [-0.35*P.height, 1.75*P.height]);
    grid(ax, 'on');

    drawnow;
    pause(0.04);
end

fprintf(['Animation finished (%d frames).  The shear center migrated from ' ...
         '(%.1f, %.1f) mm to (%.1f, %.1f) mm.\n'], frames, ...
        xSC(1)*1e3, ySC(1)*1e3, xSC(end)*1e3, ySC(end)*1e3);
end
