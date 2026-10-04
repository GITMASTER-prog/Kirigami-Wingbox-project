%MAIN  Orchestrates the whole kirigami-wingbox study.
%
%   Run:   >> main
%
%   Steps:
%     0. self-test of the toolchain (analytic validations)
%     1. morphing sweep over theta      -> results table T, results struct R
%     2. report figures                 -> Figures 1-11
%     3. morphing animation             -> Figure "Kirigami morphing animation"
%     4. export results                 -> CSV + MAT
%
%   Everything is driven by default_params.m - edit the parameters there or
%   override them here, e.g.  P = default_params(); P.chord = 0.45;

clear; clc; close all;

%% 0. self-test --------------------------------------------------------------
self_test();

%% 1. parameters + morphing sweep --------------------------------------------
P = default_params();
P.verbose = true;                 % print intermediate calculations

[R, T] = morphing_simulation(P);
disp(T);                          % the required results table

%% 2. report figures ----------------------------------------------------------
figs = visualization(R, P);

%% 3. animation ----------------------------------------------------------------
P.anim_frames = 90;               % smoothness of the animation
kirigami_animation(R, P);

%% 4. export -------------------------------------------------------------------
out_csv = 'kirigami_results.csv';
out_mat = 'kirigami_results.mat';
writetable(T, out_csv);
save(out_mat, 'R', 'P', 'T', 'figs');
fprintf(['\nExported:\n  %s  (results table)\n' ...
         '  %s  (full workspace: R, P, T, figs)\n'], out_csv, out_mat);

if isfield(P, 'save_figs') && P.save_figs
    if ~exist(P.fig_dir, 'dir'), mkdir(P.fig_dir); end
    fn = fieldnames(figs);
    for k = 1:numel(fn)
        figure(figs.(fn{k}));
        exportgraphics(gcf, fullfile(P.fig_dir, [fn{k} '.png']), 'Resolution', 200);
    end
    fprintf('  figures/*.png  (report figures)\n');
end

fprintf('\nDone. Suggested ANSYS validation configs: theta = %g, %g, %g deg.\n', ...
        P.theta_min, (P.theta_min+P.theta_max)/2, P.theta_max);
