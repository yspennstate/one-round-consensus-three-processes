"""All-p lower bound for deterministic one-round binary consensus, n = 3, independent links.

Claim: for every deterministic algorithm A and every p in [0,1],
    max_x E_x(A)(p) >= E_c(p),
x over the six non-constant inputs, E_c the courteous error (the same on every input).

D_x(A) = c(A,x) - c_court (7 integers, Bernstein form in p^k q^(6-k)).
Filter: some x has all seven entries >= 0, so D_x >= 0 on [0,1].
Residue: the distinct 6-tuples (D_x)_x are checked exactly: the roots in [0,1] of all
six polynomials are isolated together (squarefree product), and max_x D_x is tested
at every isolating-interval endpoint that is not a root and at the midpoint of every
gap; between consecutive roots the sign pattern of the six is constant, and at a root
of D_x the value D_x = 0 already satisfies the inequality.

usage: python cert4.py
"""
import sys
import numpy as np
import cons3 as C
from cert3 import build, COURTEOUS


def residue_ok(tup):
    import sympy as sp
    p = sp.symbols("p")
    L = len(tup[0]) - 1
    polys = [sp.Poly(sum(int(d) * p ** k * (1 - p) ** (L - k) for k, d in enumerate(v)), p) for v in tup]
    F = sp.Poly(1, p)
    for P in polys:
        if not P.is_zero:
            F = F * P
    F = sp.Poly(sp.sqf_part(F.as_expr()), p)
    pts = set()
    ivs = []
    if F.degree() > 0:
        for (a, b), m in F.intervals(inf=0, sup=1, eps=sp.Rational(1, 10 ** 6)):
            ivs.append((sp.Rational(a), sp.Rational(b)))
    ivs.sort()
    edges = [sp.Rational(0)] + [e for iv in ivs for e in iv] + [sp.Rational(1)]
    for e in edges:
        if F.eval(e) != 0:
            pts.add(e)
    for u, v in zip(edges, edges[1:]):
        if v > u:
            pts.add((u + v) / 2)
    worst = None
    for t in sorted(pts):
        if F.eval(t) == 0:
            continue
        m = max(P.eval(t) for P in polys)
        if m < 0:
            return False, t, m
    return True, None, None


def main():
    tabs = build()
    court = None
    for x, T, rid in tabs:
        c = T[rid[0][COURTEOUS], rid[1][COURTEOUS], rid[2][COURTEOUS]].astype(np.int32)
        assert court is None or (c == court).all(), "courteous is not an equalizer"
        court = c
    print("courteous counts (every input):", list(map(int, court)))
    nb = 1 << C.NB
    a0s = range(nb)
    if "neg" in sys.argv[1:]:
        # negative control: raise the target by one pattern at k = 3; courteous itself must now
        # land in the residue and fail the exact check (D_x = -p^3 q^3 < 0 on (0,1) for every x)
        court = court.copy()
        court[3] += 1
        a0s = [COURTEOUS]
        print("NEGATIVE CONTROL: target raised at k=3, scanning a0 = courteous only")
    residue = {}
    nres = 0
    for a0 in a0s:
        good = None
        Ds = []
        for x, T, rid in tabs:
            V = T[rid[0][a0]].astype(np.int32)
            D = V[rid[1][:, None], rid[2][None, :], :] - court[None, None, :]
            Ds.append(D)
            g = (D >= 0).all(axis=-1)
            good = g if good is None else (good | g)
        if not good.all():
            idx = np.argwhere(~good)
            nres += idx.shape[0]
            stack = np.stack([D[~good] for D in Ds], axis=1)  # (m, 6, 7)
            u = np.unique(stack.reshape(stack.shape[0], -1), axis=0)
            for r in map(tuple, u):
                residue[r] = residue.get(r, 0) + 1
        if a0 % 128 == 127:
            print(f"  a0 {a0 + 1}/{nb}: residue algorithms {nres}, distinct residue tuples {len(residue)}")
            sys.stdout.flush()
    print(f"residue algorithms (no single input dominates in Bernstein form): {nres}; distinct tuples {len(residue)}")
    sys.stdout.flush()
    bad = 0
    for k, r in enumerate(sorted(residue)):
        tup = [r[7 * j:7 * j + 7] for j in range(6)]
        ok, t, m = residue_ok(tup)
        if not ok:
            bad += 1
            if bad <= 20:
                print("  FAIL at p =", t, "max D =", m, "tuple", tup)
        if k % 500 == 499:
            print(f"  residue checked {k + 1}/{len(residue)}, failures {bad}")
            sys.stdout.flush()
    print("CERT4 ALL OK" if bad == 0 else f"CERT4 FAILED on {bad} distinct tuples")


if __name__ == "__main__":
    main()
