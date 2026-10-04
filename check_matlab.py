# check_matlab.py - static checks on the .m files + numeric mirror of the ported code
import re, glob, sys

def strip_comments_strings(src):
    out = []
    for line in src.splitlines():
        # remove strings first (naive but fine for these files)
        line2 = re.sub(r"'[^']*'", "''", line)
        # strip comment after % (not inside a string, which we removed)
        idx = line2.find('%')
        if idx >= 0:
            line2 = line2[:idx]
        out.append(line2)
    return '\n'.join(out)

fails = 0
for f in sorted(glob.glob('*.m')):
    src = open(f, encoding='utf-8').read()
    code = strip_comments_strings(src)
    # remove transposes: ' after ] ) . word digits is transpose in MATLAB, not a string.
    # Our regex above already removed strings; now remove remaining quotes (transposes).
    code_nt = re.sub(r"(?<!\w)'", '', code)
    for op, cl, name in [('(', ')', 'paren'), ('[', ']', 'bracket'), ('{', '}', 'brace')]:
        d = code_nt.count(op) - code_nt.count(cl)
        if d != 0:
            print('FAIL %s: unbalanced %s in %s (%d)' % (f, name, op+cl, d))
            fails += 1
    # keyword balance
    for kw in ('if', 'for', 'while', 'switch'):
        n_open = len(re.findall(r'(?<!\w)%s(?!\w)' % kw, code_nt))
        n_close = len(re.findall(r'(?<!\w)end(?!\w)', code_nt))
        # 'end' used as index (e.g. x(2:end)) - remove those
        n_end_index = len(re.findall(r'\(\s*[^()]*\bend\b[^()]*\)', code_nt))
    # function/end rough check per file: count 'function ' defs vs trailing end
    nfun = len(re.findall(r'^\s*function\b', src, re.M))
    print('checked %s: functions=%d' % (f, nfun))

print('\nparen/bracket/brace balance:', 'FAIL' if fails else 'OK')

# ---------------- numeric mirror of the MATLAB geometry/properties ----------
import numpy as np

def matlab_geometry(P, theta_deg, slit_side='right'):
    """1:1 mirror of kirigami_geometry.m with 1-based ids converted to 0-based."""
    c, h, N = P['chord'], P['height'], P['n_cells']
    gam = P['seam_pos']
    dc = c/N
    Delta = dc*np.sin(np.radians(theta_deg))
    xs = np.linspace(0, c, N+1)
    # 1-based ids: bot(i)=i ; top(i)=N+1+i ; seam(i)=2N+2+i ; slit(i)=3N+2+i
    bot  = lambda i: i-1
    top  = lambda i: N+i
    seam = lambda i: 2*N+1+i
    slit = lambda i: 3*N+1+i
    if slit_side == 'right':
        slit_xy = np.column_stack([xs[1:], h*np.ones(N)])   # dup of web-top i+1
    else:
        slit_xy = np.column_stack([xs[:N], h*np.ones(N)])   # dup of web-top i
    nodes = np.vstack([np.column_stack([xs, np.zeros(N+1)]),
                       np.column_stack([xs, h*np.ones(N+1)]),
                       np.column_stack([xs[:N]+gam*dc, (h+Delta)*np.ones(N)]),
                       slit_xy])
    E, C = [], []
    e_bot = {}; e_web = {}; e_left = {}; e_right = {}
    k = 0
    for i in range(1, N+1):
        k += 1; e_bot[i] = k-1
        E.append([xs[i-1], 0, xs[i], 0, P['t_lo']]); C.append([bot(i), bot(i+1)])
    for i in range(1, N+2):
        k += 1; e_web[i] = k-1
        E.append([xs[i-1], 0, xs[i-1], h, P['t_web']]); C.append([bot(i), top(i)])
    for i in range(1, N+1):
        k += 1; e_left[i] = k-1
        E.append([xs[i-1], h, xs[i-1]+gam*dc, h+Delta, P['t_up']])
        if slit_side == 'right':
            C.append([top(i), seam(i)])
        else:
            C.append([slit(i), seam(i)])
    for i in range(1, N+1):
        k += 1; e_right[i] = k-1
        E.append([xs[i-1]+gam*dc, h+Delta, xs[i], h, P['t_up']])
        if slit_side == 'right':
            C.append([seam(i), slit(i)])
        else:
            C.append([seam(i), top(i+1)])
    assert k == 4*N+1, k
    loops = []
    for i in range(1, N+1):
        loops.append(dict(elems=[e_bot[i], e_web[i+1], e_right[i], e_left[i], e_web[i]],
                          signs=[+1, +1, -1, -1, -1]))
    return np.array(E), np.array(C), loops, nodes

P = dict(chord=0.30, height=0.08, n_cells=4, seam_pos=0.40,
         t_up=1e-3, t_lo=1e-3, t_web=1.2e-3)
for side in ('right', 'left'):
    E, C, loops, nodes = matlab_geometry(P, 45.0, side)
    n = len(E)
    # tree check: BFS from a degree-1 node must reach all n elements
    adj = {}
    for e in range(n):
        adj.setdefault(C[e,0], []).append(e)
        adj.setdefault(C[e,1], []).append(e)
    roots = [nd for nd, el in adj.items() if len(el) == 1]
    root = roots[0]
    seen = {root}; dq = [root]; used = set()
    while dq:
        u = dq.pop()
        for e in adj[u]:
            if e in used: continue
            a, b = C[e]
            v = b if a == u else a
            if v not in seen:
                seen.add(v); used.add(e); dq.append(v)
    print('tree check (%s cut, all %d elements, %d free ends):' % (side, n, len(roots)),
          'OK' if len(used) == n else 'FAIL %d/%d' % (len(used), n))
E, C, loops, nodes = matlab_geometry(P, 45.0)
n = len(E)
print('element count: %d (expected %d)' % (n, 4*P['n_cells']+1))
print('node count: %d (expected %d)' % (len(nodes), 4*P['n_cells']+2))

# theta=0 closed-form section properties with the ported formulas
E0, C0, loops0, nodes0 = matlab_geometry(P, 0.0)
x1,y1,x2,y2,t = E0[:,0],E0[:,1],E0[:,2],E0[:,3],E0[:,4]
L = np.hypot(x2-x1, y2-y1)
A = np.sum(t*L)
xbar = np.sum(t*L*(x1+x2)/2)/A
ybar = np.sum(t*L*(y1+y2)/2)/A
Ix = np.sum(t*L*((y1+y2)/2-ybar)**2 + t*L*(y2-y1)**2/12)
c, h, N = P['chord'], P['height'], P['n_cells']
A_cf = c*P['t_lo'] + c*P['t_up'] + (N+1)*h*P['t_web']
Ix_cf = 2*(c*P['t_up'])*(h/2)**2 + (N+1)*P['t_web']*h**3/12
print('A err: %.2e' % (abs(A-A_cf)/A_cf))
print('ybar err: %.2e' % (abs(ybar-h/2)))
print('Ix err: %.2e' % (abs(Ix-Ix_cf)/Ix_cf))
print('loop element sets: %s' % ([l['elems'] for l in loops0],))
sys.exit(1 if fails else 0)
