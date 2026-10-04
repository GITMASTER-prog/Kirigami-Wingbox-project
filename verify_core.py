# verify_core.py
# ---------------------------------------------------------------------------
# NUMERIC VERIFICATION of the thin-walled solver that will be ported to MATLAB.
# Mirrors: kirigami_geometry.m / section_properties.m / shear_flow.m /
#          shear_center.m / bredt_torsion.m
#
# MODEL — "seam-lift cellular kirigami wingbox" (simplified from the 2023
#   cellular kirigami morphable-wingbox papers):
#   - flat bottom skin (t_lo) along the full chord at y = 0
#   - vertical webs (t_w) at each cell boundary x_i, height h
#   - top skin (t_up) per cell: TWO panels hinged at the web tops and joined
#     at a "seam" (peak) that LIFTS by Delta = d*sin(theta) as the structure
#     deploys. Seam sits at gamma*d from the cell's left web (gamma != 0.5
#     makes the section chordwise-asymmetric -> Ixy != 0 -> SC migrates in x).
#   - theta = 0  -> closed N-cell rectangular box (classic wingbox)
#   - theta > 0  -> cells stay CLOSED but grow taller & asymmetric
#
# THEORY IMPLEMENTED (all documented in MODEL_MATH.md):
#   1. Thin-walled strip section properties (exact for thin rectangles)
#   2. Generalized open-section shear flow  q(s) = a*Qx(s) - b*Qy(s),
#      a = (Vy*Iy + Vx*Ixy)/Delta0, b = (Vx*Ix + Vy*Ixy)/Delta0,
#      Delta0 = Ix*Iy - Ixy^2        [reduces to -Vy*Qx/Ix when Ixy=0]
#   3. Multi-cell closure: one constant redundancy q0_k per cell from
#      compatibility  sum(loop) q/t ds = 0
#   4. Shear center from zero-twist condition using TWO load cases (Vx, Vy)
#   5. Multi-cell Bredt-Batho torsion: J = T/(G*dphi/dz)
# ---------------------------------------------------------------------------
import numpy as np
from collections import deque

# ---------------- 1. geometry ----------------
def build_geometry(c=0.30, h=0.08, N=4, gamma=0.40,
                   t_up=1.0e-3, t_lo=1.0e-3, t_w=1.2e-3, theta_deg=45.0,
                   slit_side='right'):
    """Returns:
       E  (n x 5): x1 y1 x2 y2 t          (element local dir = node1 -> node2)
       C  (n x 2): node indices (0-based) with SLIT duplicates at the
                   right-panel/web-top joints so the open section is a TREE.
       loops: list of dicts per cell: 'elems' (element idx in CCW order),
              'signs' (+1 if CCW traversal matches stored n1->n2),
              'poly'  (CCW vertex list for shoelace area)
       nodes: list of (x, y)
    """
    th = np.radians(theta_deg)
    d = c / N
    Delta = d * np.sin(th)                      # seam lift
    xs = np.linspace(0.0, c, N + 1)

    # node id scheme (0-based):
    #   bottom node i (i=1..N+1): id = i-1
    #   top node i    (i=1..N+1): id = N+i          (at (x_i, h))
    #   peak i        (i=1..N)  : id = 2N+1+i
    #   slit dup of top node i+1  : id = 3N+1+i     (i=1..N)
    nodes = []
    for i in range(1, N + 2):
        nodes.append((xs[i-1], 0.0))            # ids 0..N
    for i in range(1, N + 2):
        nodes.append((xs[i-1], h))              # ids N+1..2N+1
    for i in range(1, N + 1):
        nodes.append((xs[i-1] + gamma*d, h + Delta))   # ids 2N+2..3N+1
    for i in range(1, N + 1):
        if slit_side == 'right':
            nodes.append((xs[i], h))            # ids 3N+2..4N+1 (dup of top i+1)
        else:
            nodes.append((xs[i-1], h))          # ids 3N+2..4N+1 (dup of top i)

    def bot(i): return i - 1          # bottom node i
    def top(i): return N + i          # top node i
    def peak(i): return 2*N + 1 + i   # peak i
    def slit(i): return 3*N + 1 + i   # slit dup (right: of top i+1; left: of top i)

    E, C = [], []
    for i in range(1, N + 1):         # bottom skin
        E.append([xs[i-1], 0.0, xs[i], 0.0, t_lo]); C.append([bot(i), bot(i+1)])
    for i in range(1, N + 2):         # webs
        E.append([xs[i-1], 0.0, xs[i-1], h, t_w]); C.append([bot(i), top(i)])
    for i in range(1, N + 1):         # left panel
        if slit_side == 'right':
            E.append([xs[i-1], h, xs[i-1] + gamma*d, h + Delta, t_up])
            C.append([top(i), peak(i)])       # stored: web-top -> seam
        else:
            E.append([xs[i-1], h, xs[i-1] + gamma*d, h + Delta, t_up])
            C.append([slit(i), peak(i)])      # stored: slit dup (web-top) -> seam
    for i in range(1, N + 1):         # right panel
        if slit_side == 'right':
            E.append([xs[i-1] + gamma*d, h + Delta, xs[i], h, t_up])
            C.append([peak(i), slit(i)])      # stored: seam -> slit dup (web-top i+1)
        else:
            E.append([xs[i-1] + gamma*d, h + Delta, xs[i], h, t_up])
            C.append([peak(i), top(i+1)])     # stored: seam -> real web-top i+1
    E = np.array(E); C = np.array(C)

    # cell loops (CCW pentagons): bottom_i, web_{i+1}, right_i(rev), left_i(rev), web_i(rev)
    loops = []
    for i in range(1, N + 1):
        e_bot = i - 1                 # bottom i        (0-based idx)
        e_webR = N + (i + 1) - 1      # web i+1
        e_right = 3*N + 1 + (i - 1)   # right panel i
        e_left = 2*N + 1 + (i - 1)    # left panel i
        e_webL = N + i - 1            # web i
        elems = [e_bot, e_webR, e_right, e_left, e_webL]
        signs = [+1, +1, -1, -1, -1]  # CCW traversal vs stored n1->n2
        poly = [(xs[i-1], 0.0), (xs[i], 0.0), (xs[i], h),
                (xs[i-1] + gamma*d, h + Delta), (xs[i-1], h)]
        loops.append(dict(elems=elems, signs=signs, poly=poly))
    return E, C, loops, nodes

def shoelace(poly):
    n = len(poly); s = 0.0
    for i in range(n):
        x1, y1 = poly[i]; x2, y2 = poly[(i+1) % n]
        s += x1*y2 - x2*y1
    return 0.5*s

# ---------------- 2. section properties (thin-walled strips) ----------------
def section_properties(E):
    x1, y1, x2, y2, t = E[:,0], E[:,1], E[:,2], E[:,3], E[:,4]
    L = np.hypot(x2-x1, y2-y1)
    A = float(np.sum(t*L))
    xm, ym = (x1+x2)/2, (y1+y2)/2
    dx, dy = x2-x1, y2-y1
    xbar = float(np.sum(t*L*xm)/A)
    ybar = float(np.sum(t*L*ym)/A)
    Ix  = float(np.sum(t*L*((ym-ybar)**2 + dy**2/12.0)))
    Iy  = float(np.sum(t*L*((xm-xbar)**2 + dx**2/12.0)))
    Ixy = float(np.sum(t*L*((xm-xbar)*(ym-ybar) + dx*dy/12.0)))
    return dict(A=A, xbar=xbar, ybar=ybar, Ix=Ix, Iy=Iy, Ixy=Ixy, L=L)

# ---------------- 3. shear flow ----------------
def shear_flow(E, C, SP, Vy, Vx=0.0, loops=None):
    """Open-section baseline on the slit tree + one constant q0_k per cell.
       Returns per element: q1,q2 (endpoint q, stored n1->n2 positive),
       Iq (exact integral of q ds, n1->n2 positive), plus resultant check."""
    n = len(E)
    t = E[:,4]
    L = np.hypot(E[:,2]-E[:,0], E[:,3]-E[:,1])
    Delta0 = SP['Ix']*SP['Iy'] - SP['Ixy']**2
    # SIGN CONVENTION: q = a*Qx - b*Qy with the a,b below gives a wall
    # resultant equal to the APPLIED force (+Vx, +Vy) - verified by the
    # equilibrium checks for BOTH load directions (Check 2).
    # NOTE the Ixy cross-terms carry MINUS signs; an earlier version had
    # + signs (and/or b negated), which violates equilibrium whenever
    # Ixy != 0 - Vy-only checks at theta=0 (Ixy=0) cannot catch it.
    a = -(Vy*SP['Iy'] - Vx*SP['Ixy'])/Delta0
    b =  (Vx*SP['Ix'] - Vy*SP['Ixy'])/Delta0

    # --- adjacency ---
    adj = {}
    for e in range(n):
        adj.setdefault(C[e,0], []).append(e)
        adj.setdefault(C[e,1], []).append(e)

    # --- root the element tree at a free end (degree-1 node) ---
    root_node = None
    for nd, el in adj.items():
        if len(el) == 1:
            root_node = nd; break
    assert root_node is not None, "no free edge - section is closed everywhere"
    start_elem = adj[root_node][0]
    near = {start_elem: root_node}
    parent = {start_elem: None}
    order = []
    dq = deque([start_elem])
    while dq:
        e = dq.popleft(); order.append(e)
        n1, n2 = C[e]
        far = n2 if near[e] == n1 else n1
        for f in adj[far]:
            if f not in parent:
                parent[f] = e; near[f] = far; dq.append(f)
    assert len(order) == n, "open section is not a tree (slit missing?)"

    # --- post-order accumulation (children before parents) ---
    children = {e: [] for e in range(n)}
    for e, p in parent.items():
        if p is not None:
            children[p].append(e)
    Qx_far = np.zeros(n); Qy_far = np.zeros(n)   # accumulated Q entering elem at far node
    q_far  = np.zeros(n); q_near = np.zeros(n)   # open-section q at far/near ends
    Iq_open = np.zeros(n)                        # exact integral of q_open ds (far->near)
    for e in reversed(order):
        fx, fy = Qx_far[e], Qy_far[e]
        xfa, yfa = C[e,0], C[e,1]
        # coordinates of far/near nodes
        fxn = C[e,1] if near[e] == C[e,0] else C[e,0]  # far node id
        nxn = near[e]
        xF, yF = NODE_XY[fxn]; xN, yN = NODE_XY[nxn]
        Le = L[e]; te = t[e]
        tx, ty = (xN-xF)/Le, (yN-yF)/Le
        # own contribution along e (far -> near): CENTROIDAL first moments
        # Qx(s) = int (y - ybar) t ds ; Qy(s) = int (x - xbar) t ds
        dxQ = te*Le*((yF + yN)/2.0 - SP['ybar'])
        dyQ = te*Le*((xF + xN)/2.0 - SP['xbar'])
        QxN = fx + dxQ; QyN = fy + dyQ
        # endpoint q (open section)
        qF = a*fx - b*fy
        qN = a*QxN - b*QyN
        # exact integral of q ds (q is QUADRATIC in s along a straight element):
        # int Qx ds = Qx_far*L + t*((y_far-ybar)*L^2/2 + ty*L^3/6)
        IqF = a*(fx*Le + te*((yF - SP['ybar'])*Le**2/2 + ty*Le**3/6)) \
            - b*(fy*Le + te*((xF - SP['xbar'])*Le**2/2 + tx*Le**3/6))
        q_far[e], q_near[e], Iq_open[e] = qF, qN, IqF
        # deliver to parent (accumulate at near node)
        p = parent[e]
        if p is not None:
            Qx_far[p] += QxN; Qy_far[p] += QyN

    # open-section integral in STORED n1->n2 convention (needed by compatibility)
    Iq_open_stored = np.zeros(n)
    for e in range(n):
        fxn = C[e,1] if near[e] == C[e,0] else C[e,0]
        Iq_open_stored[e] = (1.0 if C[e,0] == fxn else -1.0)*Iq_open[e]

    # --- multi-cell redundancies ---
    ncells = len(loops) if loops else 0
    q0 = np.zeros(0)
    if ncells:
        Amat = np.zeros((ncells, ncells)); bvec = np.zeros(ncells)
        cell_of = {}   # element -> list of (cell, sigma) ; sigma: cell-CCW vs stored dir
        for k, lp in enumerate(loops):
            for e, sg in zip(lp['elems'], lp['signs']):
                cell_of.setdefault(e, []).append((k, sg))
        for k, lp in enumerate(loops):
            for e, sg in zip(lp['elems'], lp['signs']):
                bvec[k] -= sg*Iq_open_stored[e]/t[e]
                for (m, sm) in cell_of[e]:
                    Amat[k, m] += sg*sm*L[e]/t[e]
        q0 = np.linalg.solve(Amat, bvec)

    # --- final q (stored n1->n2 convention) ---
    # q_open was accumulated in the tree's far->near direction; the constant
    # cell redundancies live in the loop (stored) direction and must NOT be
    # flipped by dirsgn.
    q1 = np.zeros(n); q2 = np.zeros(n); Iq = np.zeros(n)
    cell_of2 = {}
    if ncells:
        for k, lp in enumerate(loops):
            for e, sg in zip(lp['elems'], lp['signs']):
                cell_of2.setdefault(e, []).append((k, sg))
    for e in range(n):
        n1, n2 = C[e]
        fxn = C[e,1] if near[e] == C[e,0] else C[e,0]
        dirsgn = 1.0 if n1 == fxn else -1.0   # +1 if stored dir == far->near
        add = 0.0
        for (m, sm) in cell_of2.get(e, []):
            add += sm*q0[m]
        qF = dirsgn*q_far[e] + add            # q at far node, stored convention
        qN = dirsgn*q_near[e] + add           # q at near node, stored convention
        if n1 == fxn:
            q1[e], q2[e] = qF, qN
        else:
            q1[e], q2[e] = qN, qF
        Iq[e] = dirsgn*Iq_open[e] + add*L[e]
    # resultant check
    tx = (E[:,2]-E[:,0])/L; ty = (E[:,3]-E[:,1])/L
    Fx = float(np.sum(Iq*tx)); Fy = float(np.sum(Iq*ty))
    return q1, q2, Iq, q0, dict(Fx=Fx, Fy=Fy, a=a, b=b, Iq_open_stored=Iq_open_stored)

# ---------------- 4. shear center ----------------
def moment_about_centroid(E, Iq, SP):
    """Mz of internal shear flows about the centroid.
       (r x t)_z is CONSTANT along a straight element, so
       M = sum_e Iq_e * [ (x1-xbar)*ty - (y1-ybar)*tx ]."""
    L = np.hypot(E[:,2]-E[:,0], E[:,3]-E[:,1])
    tx = (E[:,2]-E[:,0])/L; ty = (E[:,3]-E[:,1])/L
    c = (E[:,0]-SP['xbar'])*ty - (E[:,1]-SP['ybar'])*tx
    return float(np.sum(Iq*c))

def shear_center(E, C, SP, loops, V=1000.0):
    q1y, q2y, Iqy, _, _ = shear_flow(E, C, SP, V, 0.0, loops)
    q1x, q2x, Iqx, _, _ = shear_flow(E, C, SP, 0.0, V, loops)
    M_y0 = moment_about_centroid(E, Iqy, SP)   # moment from pure Vy
    M_x0 = moment_about_centroid(E, Iqx, SP)   # moment from pure Vx
    # Zero-twist condition about ANY reference P:  0 = M_P + ((P - S) x R)_z
    # About the centroid c, R = (Vx, Vy), S = c + (ex, ey):
    #   0 = M_c + (xbar-xS)*Vy - (ybar-yS)*Vx = M_c - ex*Vy + ey*Vx
    # Vy-only case:  ex = +M_y0/Vy ;  Vx-only case:  ey = -M_x0/Vx
    ex = +M_y0/V     # SC x-offset from centroid
    ey = -M_x0/V     # SC y-offset from centroid
    return ex, ey, dict(M_y0=M_y0, M_x0=M_x0)

# ---------------- 5. multi-cell Bredt torsion ----------------
def bredt_J(E, C, SP, loops, G=1.0, T=1.0):
    """Solve multicell Bredt-Batho for pure torque T (CCW+).
       Unknowns: q0_1..q0_n, twist rate th'.  Returns J = T/(G*th')."""
    ncells = len(loops)
    A0 = np.array([shoelace(lp['poly']) for lp in loops])
    L = np.hypot(E[:,2]-E[:,0], E[:,3]-E[:,1]); t = E[:,4]
    cell_of = {}
    for k, lp in enumerate(loops):
        for e, sg in zip(lp['elems'], lp['signs']):
            cell_of.setdefault(e, []).append((k, sg))
    # unknown vector: [q0_1..q0_ncells, thp]
    Amat = np.zeros((ncells+1, ncells+1)); bvec = np.zeros(ncells+1)
    for k, lp in enumerate(loops):
        for e, sg in zip(lp['elems'], lp['signs']):
            for (m, sm) in cell_of[e]:
                Amat[k, m] += sg*sm*L[e]/t[e]
        Amat[k, ncells] = -2.0*A0[k]*G          # = 2*A0*G*th'
    for m in range(ncells):                     # equilibrium: T = 2*sum(q0*A0)
        Amat[ncells, m] = 2.0*A0[m]
    bvec[ncells] = T
    sol = np.linalg.solve(Amat, bvec)
    q0 = sol[:ncells]; thp = sol[ncells]
    J = T/(G*thp)
    return J, q0, thp

# ======================= VERIFICATION =======================
if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    OK = lambda name, err, tol: print(("  PASS  " if err < tol else "  FAIL  ") +
                                      "%-52s  err = %.3e" % (name, err))

    c, h, N = 0.30, 0.08, 4
    t_up = t_lo = 1.0e-3; t_w = 1.2e-3

    print("="*72)
    print("CHECK 1: theta=0 rectangular multi-cell box vs closed forms")
    print("="*72)
    E, C, loops, nodes = build_geometry(theta_deg=0.0, N=N)
    NODE_XY = nodes
    SP = section_properties(E)
    A_cf = 2*c*t_up + 2*h*t_w + (N-1)*0            # skins + outer webs (inner webs counted in 2h)
    # exact: bottom c*t_lo + top c*t_up + (N+1) webs h*t_w
    A_cf = c*t_lo + c*t_up + (N+1)*h*t_w
    OK("A vs closed form", abs(SP['A']-A_cf)/A_cf, 1e-12)
    OK("ybar vs h/2", abs(SP['ybar']-h/2), 1e-15)
    Ix_cf = 2*(c*t_up)*(h/2)**2 + (N+1)*t_w*h**3/12.0
    OK("Ix vs closed form", abs(SP['Ix']-Ix_cf)/Ix_cf, 1e-12)
    OK("Ixy vs 0 (theta=0 symmetric)", abs(SP['Ixy']), 1e-20)

    print("="*72)
    print("CHECK 2: shear-flow equilibrium (theta=0 and theta=60)")
    print("="*72)
    for thd in (0.0, 60.0):
        E, C, loops, nodes = build_geometry(theta_deg=thd, N=N)
        NODE_XY = nodes
        SP = section_properties(E)
        V = 1000.0
        q1, q2, Iq, q0, info = shear_flow(E, C, SP, V, 0.0, loops)
        OK("theta=%4.1f: Fy - Vy" % thd, abs(info['Fy']-V), 1e-6*V)
        OK("theta=%4.1f: Fx - 0" % thd, abs(info['Fx']), 1e-6*V)
        _, _, _, _, infox = shear_flow(E, C, SP, 0.0, V, loops)
        OK("theta=%4.1f: Fx - Vx" % thd, abs(infox['Fx']-V), 1e-6*V)
        OK("theta=%4.1f: Fy - 0" % thd, abs(infox['Fy']), 1e-6*V)
        print("        ncells=%d, q0=%s" % (len(loops), np.array2string(q0)))

    print("="*72)
    print("CHECK 2b: cut-independence (slit right vs slit left, theta=60)")
    print("="*72)
    E1, C1, loops1, nodes1 = build_geometry(theta_deg=60.0, N=N, slit_side='right')
    E2l, C2l, loops2, nodes2 = build_geometry(theta_deg=60.0, N=N, slit_side='left')
    NODE_XY = nodes1
    SP1 = section_properties(E1)
    qa = shear_flow(E1, C1, SP1, 1000.0, 0.0, loops1)[0:2]
    NODE_XY = nodes2   # left build has its own node-id table - MUST use its own coords
    qb = shear_flow(E2l, C2l, SP1, 1000.0, 0.0, loops2)[0:2]
    scale = max(abs(np.array(qa)).max(), 1e-30)
    OK("q1 slit-right vs slit-left", np.max(np.abs(qa[0]-qb[0]))/scale, 1e-9)
    OK("q2 slit-right vs slit-left", np.max(np.abs(qa[1]-qb[1]))/scale, 1e-9)

    print("="*72)
    print("CHECK 3: shear center of symmetric closed box = centroid")
    print("="*72)
    E, C, loops, nodes = build_geometry(theta_deg=0.0, N=N)
    NODE_XY = nodes
    SP = section_properties(E)
    ex, ey, MM = shear_center(E, C, SP, loops)
    OK("ex (should be 0)", abs(ex), 1e-9)
    OK("ey (should be 0)", abs(ey), 1e-9)

    print("="*72)
    print("CHECK 4: C-channel open section vs analytical SC")
    print("="*72)
    # web (0,-h/2)->(0,h/2); flanges (0,+-h/2)->(b,+-h/2), uniform t
    hb, bb, tc = 0.10, 0.06, 1.2e-3
    E2 = np.array([[0,-hb/2,0,hb/2,tc],
                   [0, hb/2,bb,hb/2,tc],
                   [0,-hb/2,bb,-hb/2,tc]])
    C2 = np.array([[0,1],[1,2],[3,4]])
    # node ids: 0=(0,-h/2) 1=(0,h/2) 2=(b,h/2) 3=(b,-h/2)... need consistency:
    # web: n1=(0,-h/2) n2=(0,h/2) -> ids 0,1
    # top flange: (0,h/2)->(b,h/2) -> ids 1,2
    # bottom flange: (0,-h/2)->(b,-h/2) -> ids 0,3
    C2 = np.array([[0,1],[1,2],[0,3]])
    NODE_XY2 = [(0,-hb/2),(0,hb/2),(bb,hb/2),(bb,-hb/2)]
    # patch shear_flow to use this node table (module-level name, __main__ scope):
    NODE_XY = NODE_XY2
    SP2 = section_properties(E2)
    Ix_cf = tc*hb**3/12 + 2*(tc*bb)*(hb/2)**2
    e_cf = tc*bb**2*hb**2/(4*Ix_cf)
    ex2, ey2, _ = shear_center(E2, C2, SP2, loops=[])
    x_SC_abs = SP2['xbar'] + ex2
    OK("C-channel Ix vs closed form", abs(SP2['Ix']-Ix_cf)/Ix_cf, 1e-12)
    OK("C-channel SC x vs analytic (=-e left of web)", abs(x_SC_abs-(-e_cf))/e_cf, 1e-10)
    OK("C-channel SC y vs analytic (=0)", abs(ey2), 1e-10)
    print("        analytic x_SC = %.6e m, numeric x_SC = %.6e m" % (-e_cf, x_SC_abs))

    print("="*72)
    print("CHECK 5: Bredt torsion, single cell vs J = 4A0^2 / (ds/t)")
    print("="*72)
    E, C, loops, nodes = build_geometry(theta_deg=0.0, N=1, c=0.30, h=0.08)
    NODE_XY = nodes
    SP = section_properties(E)
    J, q0, thp = bredt_J(E, C, SP, loops)
    A0 = 0.30*0.08
    dst = (0.30/1e-3)*2 + (0.08/1.2e-3)*2
    J_cf = 4*A0**2/dst
    OK("J vs closed form", abs(J-J_cf)/J_cf, 1e-10)
    OK("q0 vs T/(2A0)", abs(q0[0]-1.0/(2*A0)), 1e-9)

    print("="*72)
    print("CHECK 6: SC migration with morphing (informational)")
    print("="*72)
    print("  theta[deg]   xbar[mm]   ybar[mm]    Ix[cm4]    Iy[cm4]     Ixy[cm4]   ex[mm]   ey[mm]")
    for thd in (0, 15, 30, 45, 60, 75, 90):
        E, C, loops, nodes = build_geometry(theta_deg=thd, N=N)
        NODE_XY = nodes
        SP = section_properties(E)
        ex, ey, _ = shear_center(E, C, SP, loops)
        print("  %8.1f  %9.3f  %9.3f  %10.4f  %10.4f  %10.4f  %7.3f  %7.3f" %
              (thd, SP['xbar']*1e3, SP['ybar']*1e3, SP['Ix']*1e8, SP['Iy']*1e8,
               SP['Ixy']*1e8, ex*1e3, ey*1e3))
    print("\nAll numeric checks done.")
