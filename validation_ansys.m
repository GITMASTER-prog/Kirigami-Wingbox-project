function V = validation_ansys(P, R, varargin)
% VALIDATION_ANSYS  Part 10 helper: MATLAB <-> ANSYS Mechanical comparison.
%
%   V = validation_ansys(P, R) does two things:
%
%   (1) EXPORT for ANSYS: for the three validation configurations
%       (closed / intermediate / fully deployed) it writes
%           ansys_validation/cfg_<i>_theta<..>_nodes.csv   node table
%           ansys_validation/cfg_<i>_theta<..>_walls.csv   wall segments
%       so the identical thin-walled model can be rebuilt in ANSYS
%       (DesignModeler > Concept > Lines, or SpaceClaim; then SHELL meshes).
%
%   (2) COMPARE: if you have run the ANSYS models and measured the shear
%       center of each configuration, fill in the SC_ANSYS column in
%       ansys_validation/sc_ansys_template.csv (x_SC [m], y_SC [m]) and call
%           V = validation_ansys(P, R);        % export only
%           V = validation_ansys(P, R, true);  % export + compare
%       The error metric is the project definition:
%           Error(%) = |SC_MATLAB - SC_ANSYS| / |SC_ANSYS| * 100
%
%   OUTPUT: struct V with .config (table), .files (cell of exported paths),
%   and, in compare mode, .comparison (table with the error percentages).

if nargin < 1, P = default_params(); end
if nargin < 2 || isempty(R)
    P.verbose = false;
    [R, ~] = morphing_simulation(P);
end
do_compare = (nargin >= 3) && varargin{1};

outdir = 'ansys_validation';
if ~exist(outdir, 'dir'), mkdir(outdir); end

% ---- select the three validation configurations --------------------------
thetas = R.theta;
i_min  = 1;                          % closed
i_mid  = find(thetas >= (thetas(1)+thetas(end))/2, 1);   % intermediate
i_max  = numel(thetas);              % deployed
idx    = unique([i_min, i_mid, i_max]);

V = struct();
V.config = table(thetas(idx).', R.xSC(idx), R.ySC(idx), ...
                 'VariableNames', {'theta_deg', 'xSC_MATLAB_m', 'ySC_MATLAB_m'});
V.files = cell(numel(idx), 1);

fprintf('\n===== ANSYS VALIDATION EXPORT =====\n');
for k = 1:numel(idx)
    it = idx(k);
    G  = R.G_geo{it};
    tag = sprintf('cfg_%d_theta_%gdeg', k, R.theta(it));

    f_nodes = fullfile(outdir, [tag '_nodes.csv']);
    f_walls = fullfile(outdir, [tag '_walls.csv']);

    % nodes: id, x, y  (ANSYS-friendly 0-based ids)
    n_nodes = size(G.nodes, 1);
    nid = (0:n_nodes-1).';
    T_nodes = table(nid, G.nodes(:,1), G.nodes(:,2), ...
        'VariableNames', {'node_id', 'x_m', 'y_m'});
    writetable(T_nodes, f_nodes);

    % walls: elem id, node1, node2, thickness, also explicit coordinates
    m = size(G.E, 1);
    eid = (0:m-1).';
    T_walls = table(eid, G.C(:,1)-1, G.C(:,2)-1, G.E(:,5), ...
                    G.E(:,1), G.E(:,2), G.E(:,3), G.E(:,4), ...
        'VariableNames', {'elem_id', 'node1', 'node2', 'thickness_m', ...
                          'x1_m', 'y1_m', 'x2_m', 'y2_m'});
    writetable(T_walls, f_walls);

    V.files{k} = f_nodes;
    V.files{end+1} = f_walls; %#ok<AGROW>
    fprintf('  wrote %s  (theta = %g deg, SC_MATLAB = (%.2f, %.2f) mm)\n', ...
            tag, R.theta(it), R.xSC(it)*1e3, R.ySC(it)*1e3);
end

% template for the user's ANSYS results
f_tpl = fullfile(outdir, 'sc_ansys_template.csv');
T_tpl = table(V.config.theta_deg, nan(size(idx)), nan(size(idx)), ...
    'VariableNames', {'theta_deg', 'xSC_ANSYS_m', 'ySC_ANSYS_m'});
writetable(T_tpl, f_tpl);
fprintf('  template for your ANSYS results: %s\n', f_tpl);

% ---- comparison mode -----------------------------------------------------
if do_compare
    f_res = fullfile(outdir, 'sc_ansys_results.csv');
    if ~exist(f_res, 'file')
        error(['comparison requested but %s not found - fill in the ' ...
               'template and save it as sc_ansys_results.csv first'], f_res);
    end
    A = readtable(f_res);
    % match rows by theta
    xM = zeros(size(A.theta_deg));  yM = zeros(size(A.theta_deg));
    for k = 1:height(A)
        j = find(abs(thetas - A.theta_deg(k)) < 1e-6, 1);
        xM(k) = R.xSC(j);  yM(k) = R.ySC(j);
    end
    errX = abs(xM - A.xSC_ANSYS_m) ./ abs(A.xSC_ANSYS_m) * 100;
    errY = abs(yM - A.ySC_ANSYS_m) ./ abs(A.ySC_ANSYS_m) * 100;
    V.comparison = table(A.theta_deg, xM, A.xSC_ANSYS_m, errX, ...
                         yM, A.ySC_ANSYS_m, errY, ...
        'VariableNames', {'theta_deg', 'xSC_MATLAB_m', 'xSC_ANSYS_m', ...
                          'error_x_pct', 'ySC_MATLAB_m', 'ySC_ANSYS_m', ...
                          'error_y_pct'});
    fprintf('\n===== MATLAB vs ANSYS COMPARISON =====\n');
    disp(V.comparison);
end
end
