"""Three hand-written rules for Q = 3, evaluated by direct simulation on all 24 non-constant inputs (sympy, exact).
A: hear nothing -> own value; otherwise the smallest value heard.
B: as A, except that a lone message carrying 2 is ignored (decide own value).
m: the smallest value in the view (own value included).
Prints the factored agreement polynomial per input class and the worst inputs on each piece of [0,1],
and the crossings A=B, B=m on the all-distinct inputs."""
import itertools
import sympy as sp
from rules3 import P, Qs, PROCS, PATTERNS, view, poly


def A(v):
    h = [a for a in v[1:] if a is not None]
    return v[0] if not h else min(h)


def B(v):
    h = [a for a in v[1:] if a is not None]
    if h == [2]:
        return v[0]
    return v[0] if not h else min(h)


def m(v):
    return min(a for a in v if a is not None)


def analyse(name, f):
    polys = {}
    for x in itertools.product(range(3), repeat=3):
        if len(set(x)) < 2:
            continue
        cnt = [0] * 7
        for pt in PATTERNS:
            ds = {f(view(i, x, pt)) for i in PROCS}
            for i in PROCS:  # strong validity
                assert f(view(i, x, pt)) in [a for a in view(i, x, pt) if a is not None]
            if len(ds) > 1:
                cnt[len(pt)] += 1
        polys.setdefault(tuple(cnt), []).append("".join(map(str, x)))
    print(f"== rule {name}")
    for c, xs in sorted(polys.items()):
        print(f"  counts {list(c)} inputs {' '.join(xs)}  agreement = {sp.factor(1 - poly(c))}")
    cs = list(polys)
    brk = {sp.Integer(0), sp.Integer(1)}
    for a, b in itertools.combinations(cs, 2):
        d = sp.Poly(poly(a) - poly(b), P)
        if not d.is_zero:
            brk |= {r for r in d.real_roots() if 0 < r < 1}
    brk = sorted(brk, key=float)
    last = None
    for lo, hi in zip(brk, brk[1:]):
        best = max(cs, key=lambda c: poly(c).subs(P, (lo + hi) / 2))
        if best != last:
            print(f"  from {float(lo):.6f} ({lo}): worst {' '.join(polys[best])} error = {sp.factor(poly(best))}")
        last = best
    return polys


for name, f in (("A", A), ("B", B), ("m", m)):
    analyse(name, f)
eA = 1 - P**2 * Qs**2 * (1 + Qs + Qs**2)
eB = 1 - P**2 * Qs * (1 + Qs)
print("A-B on all-distinct:", sp.factor(sp.expand(eA - eB)), " roots in (0,1):", [r.evalf(8) for r in sp.Poly(sp.expand(eA - eB), P).real_roots() if 0 < r < 1])
print("B-m on all-distinct:", sp.factor(sp.expand(eB - (1 - P**2))))
