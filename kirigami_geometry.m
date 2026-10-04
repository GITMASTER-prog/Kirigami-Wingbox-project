function G = kirigami_geometry(P, theta_deg, slit_side)
% KIRIGAMI_GEOMETRY  Parametric cross-section of the cellular kirigami wingbox.
%
%   G = kirigami_geometry(P, theta_deg) builds the thin-walled cross-section
%   for morphing angle theta_deg [deg], using parameters in struct P
%   (see default_params.m).  slit_side ('right' default, or 'left') chooses
%   where the analysis slit is placed - the final closed-section solution is
%   independent of this choice (checked by self_test.m).
%
%   MODEL (simplified from the 2023 cellular-kirigami morphable-wingbox papers):
%     - flat BOTTOM skin (t_lo) along the whole chord at y = 0
%     - vertical WEBS   (t_web) at each cell boundary x_i, height h
%     - TOP skin (t_up) per cell: TWO flat panels hinged at the web tops and
%       joined at a raised "seam" (peak).  The seam sits at
%           x = x_left_web + gamma*dc        (gamma = P.seam_pos)
%           y = h + Delta,   Delta = dc*sin(theta)
%       so theta = 0 reproduces the CLOSED rectangular box and larger theta
%       lifts the seams: closed -> deployed.
%
%   OUTPUT struct G:
%       G.E      [n x 5]  x1 y1 x2 y2 t        wall midline segments
%       G.C      [n x 2]  node indices of (n1, n2); each row is stored so that
%                         the analysis graph is a TREE after the slit
%       G.loops  struct-array, one per cell:
%                 .elems  element indices along the cell loop
%                 .signs  +1/-1 loop (CCW) direction vs stored n1->n2
%                 .poly   CCW corner points of the loop (for area)
%       G.nodes           list of node coordinates [x y]
%       G.theta           morphing angle [deg]
%       G.seam_lift       seam lift Delta [m]
%
%   NOTE: wall midlines only - the thin-walled idealisation ignores
%   thickness in the plane of the section.

if nargin < 3, slit_side = 'right'; end

c  = P.chord;  h = P.height;  N = P.n_cells;
gam = P.seam_pos;  t_up = P.t_up;  t_lo = P.t_lo;  t_web = P.t_web;
dc = c / N;                      % cell width along the chord
Delta = dc * sind(theta_deg);    % seam lift for this morphing state
th    = sind(theta_deg);         % (kept for possible slope uses)

xs = linspace(0, c, N+1);        % web x positions

% ----- node table -------------------------------------------------------
% ids (1-based in MATLAB), row k of `nodes` is node id k:
%   bottom i : id i          (i = 1..N+1)   at (x_i, 0)
%   top    i : id N+1+i      (i = 1..N+1)   at (x_i, h)
%   seam   i : id 2N+2+i     (i = 1..N)     at (x_i + gamma*dc, h+Delta)
%   slit   i : id 3N+2+i     (i = 1..N)     duplicate web-top node so the
%                                            cut skin forms an open tree
if strcmp(slit_side, 'right')
    slit_xy = [xs(2:N+1).', h*ones(N,1)];   % dup of web-top i+1
else
    slit_xy = [xs(1:N).',   h*ones(N,1)];   % dup of web-top i
end

nodes = [xs.', zeros(N+1,1); ...          % bottom skin nodes  (ids 1..N+1)
         xs.', h*ones(N+1,1); ...         % web-top nodes      (ids N+2..2N+2)
         xs(1:N).' + gam*dc, (h+Delta)*ones(N,1); ...  % seams
         slit_xy];                        % slit duplicates

bot  = @(i) i;            % bottom node i, ids 1 .. N+1
top  = @(i) N+1+i;        % web-top i,  ids N+2 .. 2N+2
seam = @(i) 2*N+2+i;      % seam i,     ids 2N+3 .. 3N+2
slit = @(i) 3*N+2+i;      % slit dup i, ids 3N+3 .. 4N+2

% ----- elements ---------------------------------------------------------
E = zeros(4*N+1, 5);  C = zeros(4*N+1, 2);
k = 0;

for i = 1:N                          % --- bottom skin panels
    k = k+1;  e_bot(i) = k;
    E(k,:) = [xs(i), 0, xs(i+1), 0, t_lo];
    C(k,:) = [bot(i), bot(i+1)];
end
for i = 1:N+1                        % --- webs
    k = k+1;  e_web(i) = k;
    E(k,:) = [xs(i), 0, xs(i), h, t_web];
    C(k,:) = [bot(i), top(i)];
end
for i = 1:N                          % --- top skin, left panel (web_i -> seam)
    k = k+1;  e_left(i) = k;
    E(k,:) = [xs(i), h, xs(i)+gam*dc, h+Delta, t_up];
    if strcmp(slit_side, 'right')
        C(k,:) = [top(i), seam(i)];
    else
        C(k,:) = [slit(i), seam(i)];
    end
end
for i = 1:N                          % --- top skin, right panel (seam -> web_{i+1})
    k = k+1;  e_right(i) = k;
    E(k,:) = [xs(i)+gam*dc, h+Delta, xs(i+1), h, t_up];
    if strcmp(slit_side, 'right')
        C(k,:) = [seam(i), slit(i)];
    else
        C(k,:) = [seam(i), top(i+1)];
    end
end
assert(k == 4*N+1, 'element counter mismatch');

% ----- cell loops (CCW when viewed with x right / y up) ------------------
% loop i = bottom_i -> web_{i+1} -> right_i (reversed) -> left_i (reversed)
%        -> web_i (reversed)
for i = 1:N
    loops(i).elems = [e_bot(i), e_web(i+1), e_right(i), e_left(i), e_web(i)];
    loops(i).signs = [+1, +1, -1, -1, -1];
    loops(i).poly  = [ xs(i),        0; ...
                       xs(i+1),      0; ...
                       xs(i+1),      h; ...
                       xs(i)+gam*dc, h+Delta; ...
                       xs(i),        h ];
end

G = struct('E', E, 'C', C, 'loops', loops, ...
           'nodes', nodes, 'theta', theta_deg, 'seam_lift', Delta);
end
