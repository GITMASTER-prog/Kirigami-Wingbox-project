function [J, Kt, out] = bredt_torsion(G, SP, P, verbose)
% BREDT_TORSION  Torsional constant J and stiffness Kt of the multi-cell section.
%
%   [J, Kt, out] = bredt_torsion(G, SP, P) solves the multi-cell Bredt-Batho
%   system for a UNIT twist rate and returns the torsion constant.
%
%   THEORY (Bredt-Batho, multicell):
%   Equilibrium:      T = 2 * sum_k q0_k * A0_k            (k = cell)
%   Twist rate:       dphi/dz = q_loop / (2 A0_k G) loop-averaged, i.e.
%                     sum_loop (q ds / t) = 2 A0_k G (dphi/dz)
%   Solving both for a unit torque T = 1 N*m gives the twist rate thp
%   [1/m] (with G from P.G); the torsion constant is
%                     J = T / (G * dphi/dz)                [m^4]
%   and the torsional stiffness of a reference-length shaft is
%                     Kt = G*J / L_ref                     [N*m/rad]
%   (G and J for the current morphing state; L_ref = P.L_ref).
%
%   OUTPUT:
%       J    torsion constant [m^4]
%       Kt   torsional stiffness for length L_ref [N*m/rad]
%       out  struct: .q0 cell flows for T = 1 [N/m], .thp twist rate [1/m],
%                    .A0 cell areas [m^2], .T_check equilibrium check [N*m]
%
%   UNITS: SI throughout.

if nargin < 4, verbose = false; end

E = G.E;  L = hypot(E(:,3)-E(:,1), E(:,4)-E(:,2));  t = E(:,5);
ncells = numel(G.loops);
A0 = zeros(1, ncells);
for kk = 1:ncells
    A0(kk) = shoelace_area(G.loops(kk).poly);
end

% element -> (cell, sigma) incidence
cell_of = cell(1, size(E,1));
for kk = 1:ncells
    lp = G.loops(kk);
    for jj = 1:numel(lp.elems)
        e = lp.elems(jj);  sg = lp.signs(jj);
        cell_of{e}(end+1,:) = [kk, sg]; %#ok<AGROW>
    end
end

% unknowns: [q0_1 .. q0_ncells, thp]      (T = 1 N*m, G from P)
Amat = zeros(ncells+1);  bvec = zeros(ncells+1,1);
for kk = 1:ncells
    lp = G.loops(kk);
    for jj = 1:numel(lp.elems)
        e = lp.elems(jj);  sg = lp.signs(jj);
        for mm = 1:size(cell_of{e},1)
            m = cell_of{e}(mm,1);  sm = cell_of{e}(mm,2);
            Amat(kk,m) = Amat(kk,m) + sg*sm*L(e)/t(e);
        end
    end
    Amat(kk, ncells+1) = -2*A0(kk)*P.G;   % twist-rate coupling
end
for kk = 1:ncells
    Amat(ncells+1, kk) = 2*A0(kk);        % T = 2*sum(q0*A0)
end
bvec(ncells+1) = 1;                       % unit torque

sol  = Amat \ bvec;
q0   = sol(1:ncells).';                   % cell flows for T = 1 [N/m]
thp  = sol(ncells+1);                     % twist rate for T = 1 [1/m]
J    = 1 / (P.G * thp);                   % [m^4]
Kt   = P.G * J / P.L_ref;                 % [N*m/rad]

T_check = 2 * sum(q0 .* A0);              % should be 1.0

out = struct('q0', q0, 'thp', thp, 'A0', A0, 'T_check', T_check);

if verbose
    fprintf('\n--- BREDT TORSION (verbose) ---\n');
    fprintf('  cell areas A0 [m^2]: '); fprintf('  %.6e', A0); fprintf('\n');
    fprintf('  unit-torque cell flows q0 [N/m]: '); fprintf('  %10.3f', q0); fprintf('\n');
    fprintf('  twist rate at T=1: %.6e 1/m  ->  J = %.6e m^4\n', thp, J);
    fprintf('  equilibrium check T = 2*sum(q0*A0) = %.9f N*m (should be 1)\n', T_check);
    fprintf('  torsional stiffness Kt = G*J/L_ref = %.6e N*m/rad\n\n', Kt);
end
end

% ------------------------------------------------------------------------
function A = shoelace_area(poly)
% Shoelace formula for the polygon area of a cell loop (CCW positive).
s = 0;
n = size(poly,1);
for ii = 1:n
    x1 = poly(ii,1); y1 = poly(ii,2);
    x2 = poly(mod(ii, n)+1,1); y2 = poly(mod(ii, n)+1,2);
    s = s + x1*y2 - x2*y1;
end
A = 0.5 * abs(s);
end
