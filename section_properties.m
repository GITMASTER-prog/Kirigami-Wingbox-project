function SP = section_properties(E, verbose)
% SECTION_PROPERTIES  Thin-walled cross-section properties from wall segments.
%
%   SP = section_properties(E)          % E = [x1 y1 x2 y2 t] rows
%   SP = section_properties(E, true)    % also prints the per-element table
%
%   THEORY (thin-walled idealisation, t << all other dimensions):
%   Each wall segment is a thin rectangle of length L and thickness t,
%   so its area is A_e = t*L and all section integrals are evaluated EXACTLY
%   along the straight midline (the t^2 terms neglected by the idealisation
%   are of order (t/L)^2 and are dropped):
%
%       A  = sum_e  t_e * L_e
%       xbar = (1/A) sum_e t_e*L_e*xm_e          xm = midpoint x
%       ybar = (1/A) sum_e t_e*L_e*ym_e          ym = midpoint y
%
%       Ix  = sum_e t_e*L_e*( (ym-ybar)^2 + dy^2/12 )   dy = y2-y1
%       Iy  = sum_e t_e*L_e*( (xm-xbar)^2 + dx^2/12 )   dx = x2-x1
%       Ixy = sum_e t_e*L_e*( (xm-xbar)*(ym-ybar) + dx*dy/12 )
%
%   The dx^2/12 etc. terms are the strip's own second moment about its
%   midline axis (b*h^3/12 for a thin rectangle), which keeps inclined
%   walls exact - not just the parallel-axis approximation.
%
%   OUTPUT struct SP:
%       .A .xbar .ybar .Ix .Iy .Ixy      the usual section properties
%       .L .xm .ym                       per-element length / midpoints
%
%   UNITS: E rows in metres -> A [m^2], I [m^4].

if nargin < 2, verbose = false; end

x1 = E(:,1); y1 = E(:,2); x2 = E(:,3); y2 = E(:,4); t = E(:,5);
dx = x2 - x1;            dy = y2 - y1;
L  = hypot(dx, dy);                      % wall midline length
assert(all(L > 0), 'zero-length wall segment in E');

A  = sum(t .* L);
xm = (x1 + x2)/2;                        % midpoints
ym = (y1 + y2)/2;

% ---- first moments -> centroid ------------------------------------------
Qx = sum(t .* L .* ym);                  % first moment about x axis (raw)
Qy = sum(t .* L .* xm);                  % first moment about y axis (raw)
xbar = Qy / A;
ybar = Qx / A;

% ---- second moments / product of inertia (about the CENTROID) -----------
Ix  = sum(t .* L .* ( (ym - ybar).^2 + dy.^2/12 ));
Iy  = sum(t .* L .* ( (xm - xbar).^2 + dx.^2/12 ));
Ixy = sum(t .* L .* ( (xm - xbar).*(ym - ybar) + dx.*dy/12 ));

SP = struct('A', A, 'xbar', xbar, 'ybar', ybar, ...
            'Ix', Ix, 'Iy', Iy, 'Ixy', Ixy, 'L', L, 'xm', xm, 'ym', ym);

if verbose
    fprintf('\n--- SECTION PROPERTIES (verbose) ---\n');
    fprintf('  elem      L[m]         t[m]       t*L[m^2]     xm[m]      ym[m]\n');
    for e = 1:size(E,1)
        fprintf('  %4d  %10.4f  %10.4f  %11.4e  %9.4f  %9.4f\n', ...
                e, L(e), t(e), t(e)*L(e), xm(e), ym(e));
    end
    fprintf('  A    = %.6e m^2\n', A);
    fprintf('  xbar = %.6f  m,   ybar = %.6f m\n', xbar, ybar);
    fprintf('  Ix   = %.6e m^4\n', Ix);
    fprintf('  Iy   = %.6e m^4\n', Iy);
    fprintf('  Ixy  = %.6e m^4\n', Ixy);
    fprintf('  Delta0 = Ix*Iy - Ixy^2 = %.6e m^8\n\n', Ix*Iy - Ixy^2);
end
end
