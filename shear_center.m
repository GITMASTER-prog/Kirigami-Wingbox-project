function [ex, ey, SC] = shear_center(G, SP, V, verbose)
% SHEAR_CENTER  Shear-center location from the zero-twist condition.
%
%   [ex, ey, SC] = shear_center(G, SP, V) solves TWO shear-flow load cases
%   (pure Vy = V and pure Vx = V) and evaluates the moment of the internal
%   shear flows about the CENTROID.  The shear center S = (xSC, ySC) is the
%   point where applying the resultant produces ZERO net twist:
%
%       0 = M_c + ((c - S) x R)_z = M_c - ex*Vy + ey*Vx
%
%   (c = centroid, ex = xSC - xbar, ey = ySC - ybar, R = applied force).
%   Pure-Vy case  ->  ex = +M_y0 / V
%   Pure-Vx case  ->  ey = -M_x0 / V
%   (derivation and sign convention in MODEL_MATH.md; both were confirmed
%    against the textbook C-channel result e = b^2 h^2 t / (4 Ix) by self_test)
%
%   MOMENT OF THE FLOW: the integrand (r x t_hat)_z is CONSTANT along each
%   straight wall segment, so with Iq_e the exact integral of q ds,
%       M_c = sum_e Iq_e * [ (x1_e - xbar)*ty_e - (y1_e - ybar)*tx_e ]
%
%   OUTPUT:
%       ex, ey   shear-center offsets from the centroid [m]
%       SC       struct: .x .y (absolute), .ex .ey, .M_y0 .M_x0 [N*m]
%
%   UNITS: V [N], positions [m], moment [N*m].

if nargin < 4, verbose = false; end

outY = shear_flow(G, SP, V, 0);
outX = shear_flow(G, SP, 0, V);

L  = hypot(G.E(:,3)-G.E(:,1), G.E(:,4)-G.E(:,2));
tx = (G.E(:,3)-G.E(:,1))./L;   ty = (G.E(:,4)-G.E(:,2))./L;
r1 = (G.E(:,1)-SP.xbar)*ty - (G.E(:,2)-SP.ybar)*tx;   % (r x t)_z, constant per segment

M_y0 = sum(outY.Iq .* r1);      % moment of internal flows, pure-Vy case [N*m]
M_x0 = sum(outX.Iq .* r1);      % moment of internal flows, pure-Vx case [N*m]

ex =  M_y0 / V;
ey = -M_x0 / V;

SC = struct('ex', ex, 'ey', ey, ...
            'x', SP.xbar + ex, 'y', SP.ybar + ey, ...
            'M_y0', M_y0, 'M_x0', M_x0);

if verbose
    fprintf('\n--- SHEAR CENTER (verbose) ---\n');
    fprintf('  moment about centroid, pure-Vy case: M_y0 = %+.6e N*m\n', M_y0);
    fprintf('  moment about centroid, pure-Vx case: M_x0 = %+.6e N*m\n', M_x0);
    fprintf('  ex = M_y0/V = %+.6e m,  ey = -M_x0/V = %+.6e m\n', ex, ey);
    fprintf('  SC (absolute) = (%.6f, %.6f) m\n\n', SC.x, SC.y);
end
end
