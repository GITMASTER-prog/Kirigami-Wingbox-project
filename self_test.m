function self_test()
% SELF_TEST  Analytic validation of the whole toolchain, inside MATLAB.
%
%   self_test() reruns the same checks that were used to verify the solver
%   design (see MODEL_MATH.md section 8).  It prints PASS/FAIL per check and
%   returns nothing.  Run it FIRST after copying the code to a new machine.
%
%   Checks:
%     1. theta = 0 closed box: A, ybar, Ix, Ixy vs closed forms
%     2. shear-flow resultant == applied force (equilibrium), two angles
%     3. symmetric closed box: SC == centroid
%     4. C-channel: SC vs textbook  e = b^2 h^2 t / (4 Ix)
%     5. single-cell Bredt: J = 4 A0^2 / point(d s/t), q0 = T/(2 A0)
%     6. cut-independence: same q from 'right' and 'left' slits

fprintf('==== KIRIGAMI WINGBOX TOOLCHAIN SELF-TEST ====\n');
n_fail = 0;

% ---------------- 1. closed rectangular box -------------------------------
P = default_params();
G = kirigami_geometry(P, 0);
SP = section_properties(G.E);
c = P.chord; h = P.height; N = P.n_cells;
t_up = P.t_up; t_lo = P.t_lo; t_web = P.t_web;

A_cf  = c*t_lo + c*t_up + (N+1)*h*t_web;
Ix_cf = 2*(c*t_up)*(h/2)^2 + (N+1)*t_web*h^3/12;
n_fail = n_fail + chk('A  vs closed form',        abs(SP.A - A_cf)/A_cf,  1e-12);
n_fail = n_fail + chk('ybar vs h/2',              abs(SP.ybar - h/2),     1e-12);
n_fail = n_fail + chk('Ix vs closed form',        abs(SP.Ix - Ix_cf)/Ix_cf, 1e-12);
n_fail = n_fail + chk('Ixy vs 0 (symmetric)',     abs(SP.Ixy),            1e-15);
n_fail = n_fail + chk('xbar vs c/2 (symmetric)',  abs(SP.xbar - c/2),     1e-12);

% ---------------- 2. equilibrium -------------------------------------------
for th = [0, 60]
    G = kirigami_geometry(P, th);
    SP = section_properties(G.E);
    sh = shear_flow(G, SP, 1000, 0);
    n_fail = n_fail + chk(sprintf('Fy = Vy at theta=%g', th), ...
                          abs(sh.Fy - 1000)/1000, 1e-9);
    n_fail = n_fail + chk(sprintf('Fx = 0  at theta=%g', th), ...
                          abs(sh.Fx)/1000, 1e-9);
end

% ---------------- 3. symmetric box: SC == centroid -------------------------
% (rebuild at theta = 0 - G, SP currently hold the theta = 60 case from above)
G = kirigami_geometry(P, 0);
SP = section_properties(G.E);
[ex, ey] = shear_center(G, SP, 1000);
n_fail = n_fail + chk('ex (sym box)', abs(ex), 1e-9);
n_fail = n_fail + chk('ey (sym box)', abs(ey), 1e-9);

% ---------------- 4. C-channel textbook SC ---------------------------------
% web 0,+-h/2 ; flanges to +b ; uniform t.  Analytic:
%   xbar = ..., Ix = t h^3/12 + 2 t b (h/2)^2,  e = b^2 h^2 t / (4 Ix)
hb = 0.10; bb = 0.06; tc = 1.2e-3;
E2 = [0, -hb/2, 0, hb/2, tc;
      0,  hb/2, bb, hb/2, tc;
      0, -hb/2, bb, -hb/2, tc];
% store the C-channel as a minimal G-like struct (one open tree, no loops)
G2 = struct();
G2.E = E2;
G2.C = [1 2; 2 3; 1 4];
G2.nodes = [0, -hb/2; 0, hb/2; bb, hb/2; bb, -hb/2];
G2.loops = struct('elems', {}, 'signs', {}, 'poly', {});
SP2 = section_properties(E2);
Ix2 = tc*hb^3/12 + 2*tc*bb*(hb/2)^2;
e_cf = bb^2 * hb^2 * tc / (4 * Ix2);
% solve the open section by hand here (no cells): reuse shear_flow
outY = shear_flow(G2, SP2, 1000, 0);
outX = shear_flow(G2, SP2, 0, 1000);
L2  = hypot(E2(:,3)-E2(:,1), E2(:,4)-E2(:,2));
tx2 = (E2(:,3)-E2(:,1))./L2;  ty2 = (E2(:,4)-E2(:,2))./L2;
r12 = (E2(:,1)-SP2.xbar)*ty2 - (E2(:,2)-SP2.ybar)*tx2;
M_y0 = sum(outY.Iq .* r12);
M_x0 = sum(outX.Iq .* r12);
ex2 = M_y0/1000;  ey2 = -M_x0/1000;
x_SC2 = SP2.xbar + ex2;
n_fail = n_fail + chk('C-channel Ix', abs(SP2.Ix - Ix2)/Ix2, 1e-12);
n_fail = n_fail + chk('C-channel SC x = -e', abs(x_SC2 - (-e_cf))/e_cf, 1e-10);
n_fail = n_fail + chk('C-channel SC y = 0', abs(ey2), 1e-10);

% ---------------- 5. single-cell Bredt -------------------------------------
P1 = P; P1.n_cells = 1;
G1 = kirigami_geometry(P1, 0);
SP1 = section_properties(G1.E);
[J1, ~, out1] = bredt_torsion(G1, SP1, P1);
A0 = c*h;
dst = 2*c/P1.t_up + 2*h/P1.t_web;
J_cf = 4*A0^2/dst;
n_fail = n_fail + chk('Bredt J = 4A0^2/(ds/t)', abs(J1 - J_cf)/J_cf, 1e-10);
n_fail = n_fail + chk('Bredt q0 = T/(2A0)', abs(out1.q0 - 1/(2*A0)), 1e-9);

% ---------------- 6. cut-independence ---------------------------------------
Ga = kirigami_geometry(P, 60, 'right');
Gb = kirigami_geometry(P, 60, 'left');
Sa = shear_flow(Ga, section_properties(Ga.E), 1000, 0);
Sb = shear_flow(Gb, section_properties(Gb.E), 1000, 0);
qmax = max(abs(Sa.q1));
n_fail = n_fail + chk('q1 slit right == slit left', max(abs(Sa.q1 - Sb.q1))/qmax, 1e-9);
n_fail = n_fail + chk('q2 slit right == slit left', max(abs(Sa.q2 - Sb.q2))/qmax, 1e-9);

% ----------------------------------------------------------------------------
fprintf('================================================\n');
if n_fail == 0
    fprintf('ALL SELF-TESTS PASSED.\n');
else
    fprintf('%d CHECK(S) FAILED - do not trust the results!\n', n_fail);
end
end

% ----------------------------------------------------------------------------
function ok = chk(name, err, tol)
ok = double(err < tol);
if ok
    fprintf('  PASS  %-42s err = %.3e\n', name, err);
else
    fprintf('  FAIL  %-42s err = %.3e (tol %.1e)\n', name, err, tol);
end
end
