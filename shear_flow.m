function out = shear_flow(G, SP, Vy, Vx, verbose)
% SHEAR_FLOW  Thin-walled shear-flow distribution q(s) on the cellular section.
%
%   out = shear_flow(G, SP, Vy, Vx)  solves the multi-cell section described
%   by G (kirigami_geometry) with properties SP (section_properties) under
%   transverse forces Vx, Vy [N] applied at the CENTROID (the shear-flow
%   reference point; see shear_center.m for the shear center itself).
%
%   THEORY (full derivation in MODEL_MATH.md):
%   For an unsymmetric thin-walled section the shear flow satisfies
%       q(s) = a*Qx(s) - b*Qy(s) + sum_k q0k        (closed cells)
%   with the first moments of the cut (open) section
%       Qx(s) = int_0^s (y-ybar) t ds ,   Qy(s) = int_0^s (x-xbar) t ds
%   and constants
%       a = -(Vy*Iy - Vx*Ixy) / Delta0,   b = (Vx*Ix - Vy*Ixy) / Delta0,
%       Delta0 = Ix*Iy - Ixy^2.
%   The q0k are CELL-CONSTANT redundancies enforcing zero twist of each
%   cell:  sum_loop q/t ds = 0.  The sign convention is chosen so the
%   resultant of q equals the APPLIED force (verified to 1e-13 by self_test
%   for BOTH load directions, theta = 0 and 60 deg).
%
%   FIX (verified by independent truth-solve): the Ixy cross-terms carry
%   MINUS signs.  The earlier + signs (a = -(Vy*Iy + Vx*Ixy)/...) made the
%   flow violate equilibrium whenever Ixy ~= 0; Vy-only checks at theta = 0
%   (Ixy = 0) could not catch it.  self_test Check 2 now guards both loads.
%
%   NUMERICS (beginner-friendly, no black boxes):
%   1. The slit section is a TREE of wall segments (G.C).  A depth-first
%      walk accumulates Qx, Qy at every node - each segment contributes
%      t*L*(midpoint offset) to the far end.
%   2. q_open and its exact integral Iq_open are evaluated per segment.
%      q(s) is QUADRATIC in s along each straight wall, so Iq_open is exact.
%   3. One loop equation per cell is assembled and solved for the q0k.
%   4. Superposition: q = q_open + sum(k in cells of segment) q0k.
%
%   OUTPUT struct out:
%       .q1 .q2    shear flow at the stored nodes n1, n2 of each segment [N/m]
%       .Iq        exact integral of q ds along each segment (stored dir) [N]
%       .q0        cell redundancies [N/m]
%       .Fx .Fy    resultant of the flow (== Vx, Vy by equilibrium)
%       .a .b      load coefficients [1/m^2]
%       .Iq_open   open-section integral per segment (for reuse)

if nargin < 5, verbose = false; end

E = G.E;  C = G.C;  n = size(E,1);
t = E(:,5);
L = hypot(E(:,3)-E(:,1), E(:,4)-E(:,2));

Delta0 = SP.Ix*SP.Iy - SP.Ixy^2;
a = -(Vy*SP.Iy - Vx*SP.Ixy) / Delta0;      % [1/m^2]
b =  (Vx*SP.Ix - Vy*SP.Ixy) / Delta0;      % [1/m^2]

% ---- 1. adjacency list of the slit (tree) graph -------------------------
adj = cell(1, max(C(:)));
for e = 1:n
    adj{C(e,1)}(end+1) = e; %#ok<AGROW>
    adj{C(e,2)}(end+1) = e; %#ok<AGROW>
end

% ---- 2. root the tree at a free end (degree-1 node) ---------------------
root_node = find(cellfun(@numel, adj) == 1, 1);
assert(~isempty(root_node), 'no free edge - slit the section (see kirigami_geometry)');

start_elem = adj{root_node}(1);
near = zeros(1, n);          % node id at the NEAR (parent-side) end
parent = zeros(1, n);        % parent element in the DFS tree
order = zeros(1, n);         % DFS visit order (parents before children)
near(start_elem) = root_node; parent(start_elem) = -1;
visited = false(1, n);  visited(start_elem) = true;
head = 1; order(1) = start_elem; tail = 1;
while head <= tail
    e = order(head); head = head+1;
    far = C(e,2); if near(e) == C(e,2), far = C(e,1); end
    for f = adj{far}
        if ~visited(f)
            visited(f) = true;
            parent(f) = e; near(f) = far;
            tail = tail+1; order(tail) = f;
        end
    end
end
assert(tail == n, 'slit section is not a tree - check the slit construction');

% ---- 3. post-order accumulation (children before parents) ---------------
Qx_far = zeros(1,n);  Qy_far = zeros(1,n);   % accumulated Q entering at far node
q_far  = zeros(1,n);  q_near = zeros(1,n);   % open-section q at far / near ends
Iq_open = zeros(1,n);                        % exact int of q_open ds (far->near)
for ii = n:-1:1
    e = order(ii);
    fx = Qx_far(e);  fy = Qy_far(e);
    far = C(e,2); if near(e) == C(e,2), far = C(e,1); end
    xF = G.nodes(far,1);  yF = G.nodes(far,2);
    xN = G.nodes(near(e),1);  yN = G.nodes(near(e),2);
    Le = L(e);  te = t(e);
    tx = (xN-xF)/Le;   ty = (yN-yF)/Le;
    % centroidal first moments contributed by THIS segment (exact):
    dQx = te*Le*( (yF+yN)/2 - SP.ybar );
    dQy = te*Le*( (xF+xN)/2 - SP.xbar );
    QxN = fx + dQx;   QyN = fy + dQy;
    % endpoint flows of the OPEN section (far -> near):
    qF = a*fx - b*fy;
    qN = a*QxN - b*QyN;
    % exact integral of q ds (q is quadratic in s along a straight wall):
    IqF = a*( fx*Le + te*((yF-SP.ybar)*Le^2/2 + ty*Le^3/6) ) ...
        - b*( fy*Le + te*((xF-SP.xbar)*Le^2/2 + tx*Le^3/6) );
    q_far(e) = qF;  q_near(e) = qN;  Iq_open(e) = IqF;
    p = parent(e);
    if p > 0
        Qx_far(p) = Qx_far(p) + QxN;
        Qy_far(p) = Qy_far(p) + QyN;
    end
end

% ---- open-section integrals in the STORED n1->n2 direction --------------
Iq_open_stored = zeros(1,n);
for e = 1:n
    far = C(e,2); if near(e) == C(e,2), far = C(e,1); end
    sgn = 1; if C(e,1) ~= far, sgn = -1; end
    Iq_open_stored(e) = sgn * Iq_open(e);
end

% ---- 4. cell-loop compatibility: one equation per cell ------------------
ncells = numel(G.loops);
q0 = zeros(1, ncells);
if ncells > 0
    Amat = zeros(ncells);  bvec = zeros(ncells,1);
    cell_of = cell(1, n);                 % element -> list of (cell, sigma)
    for kk = 1:ncells
        lp = G.loops(kk);
        for jj = 1:numel(lp.elems)
            e = lp.elems(jj);  sg = lp.signs(jj);
            cell_of{e}(end+1,:) = [kk, sg]; %#ok<AGROW>
        end
    end
    for kk = 1:ncells
        lp = G.loops(kk);
        for jj = 1:numel(lp.elems)
            e = lp.elems(jj);  sg = lp.signs(jj);
            bvec(kk) = bvec(kk) - sg * Iq_open_stored(e) / t(e);
            for mm = 1:size(cell_of{e},1)
                m = cell_of{e}(mm,1);  sm = cell_of{e}(mm,2);
                Amat(kk,m) = Amat(kk,m) + sg*sm*L(e)/t(e);
            end
        end
    end
    q0 = (Amat \ bvec).';                 % row vector of cell flows
end

% ---- 5. superposition -> final q on every segment -----------------------
q1 = zeros(1,n);  q2 = zeros(1,n);  Iq = zeros(1,n);
for e = 1:n
    far = C(e,2); if near(e) == C(e,2), far = C(e,1); end
    dirsgn = 1; if C(e,1) ~= far, dirsgn = -1; end  % +1 if stored dir == far->near
    add = 0;
    if ncells > 0
        for mm = 1:size(cell_of{e},1)
            m = cell_of{e}(mm,1);  sm = cell_of{e}(mm,2);
            add = add + sm*q0(m);
        end
    end
    qF = dirsgn*q_far(e)  + add;          % stored convention: q at node1
    qN = dirsgn*q_near(e) + add;          % stored convention: q at node2
    if C(e,1) == far
        q1(e) = qF;  q2(e) = qN;
    else
        q1(e) = qN;  q2(e) = qF;
    end
    Iq(e) = dirsgn*Iq_open(e) + add*L(e);
end

% ---- resultant check (should reproduce Vx, Vy to roundoff) --------------
tx = (E(:,3)-E(:,1))./L;   ty = (E(:,4)-E(:,2))./L;
Fx = sum(Iq .* tx);        Fy = sum(Iq .* ty);

if verbose
    fprintf('\n--- SHEAR FLOW (verbose) ---\n');
    fprintf('  load coefficients: a = %.4e 1/m^2, b = %.4e 1/m^2\n', a, b);
    fprintf('  cell redundancies q0k [N/m]:\n');
    fprintf('    '); fprintf('  %10.3f', q0); fprintf('\n');
    fprintf('  resultant check: Fx = %.6e N, Fy = %.6e N (applied: %.1f, %.1f)\n', ...
            Fx, Fy, Vx, Vy);
end

out = struct('q1', q1, 'q2', q2, 'Iq', Iq, 'q0', q0, ...
             'Fx', Fx, 'Fy', Fy, 'a', a, 'b', b, 'Iq_open', Iq_open_stored);
end
