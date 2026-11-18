"""Describe the q = 3 rules CP-SAT found and the error polynomials they give on every input.

For each rule file: print each process's decisions on its multi-valued views; recompute per-input counts by
direct simulation (no milp3 code); list the distinct error polynomials E_x(p) = sum c_k p^k q^(6-k), factored;
and on [0,1] find which input attains the maximum, with the breakpoints (sympy real roots of the differences).
Also checks 1 - p^2 q (1+q) against the counts [1,6,13,15,11,5,1] and factors (1-p^2) minus it.

usage: python rules3.py out/rule_n3_q3_p1_10.json [...]
"""
import sys
import json
import itertools
import sympy as sp

P = sp.symbols("p")
Qs = 1 - P
PROCS = (0, 1, 2)
LINKS = [(s, r) for r in PROCS for s in PROCS if s != r]
PATTERNS = [frozenset(c) for k in range(7) for c in itertools.combinations(LINKS, k)]


def view(i, x, pat):
    return (x[i],) + tuple((x[j] if (j, i) in pat else None) for j in PROCS if j != i)


def poly(c):
    return sp.expand(sum(ck * P ** k * Qs ** (6 - k) for k, ck in enumerate(c)))


def analyse(path):
    rows = json.load(open(path, encoding="utf-8"))
    rule = {(i, tuple(vw)): v for i, vw, v in rows}
    print(f"== {path}: {len(rule)} multi-valued views")
    for i in PROCS:
        print(f"  process {i}: " + "; ".join(f"{vw}->{v}" for (j, vw), v in sorted(rule.items(), key=str) if j == i))
    polys = {}
    for x in itertools.product(range(3), repeat=3):
        if len(set(x)) < 2:
            continue
        cnt = [0] * 7
        for pt in PATTERNS:
            ds = set()
            for i in PROCS:
                v = view(i, x, pt)
                vs = sorted({a for a in v if a is not None})
                ds.add(vs[0] if len(vs) == 1 else rule[(i, v)])
            if len(ds) > 1:
                cnt[len(pt)] += 1
        polys.setdefault(tuple(cnt), []).append("".join(map(str, x)))
    for c, xs in sorted(polys.items()):
        print(f"  counts {list(c)} inputs {' '.join(xs)}  E = {sp.factor(1 - poly(c))} subtracted from 1")
    # which count vector is the max on which part of [0,1]
    cs = list(polys)
    brk = {sp.Integer(0), sp.Integer(1)}
    for a, b in itertools.combinations(cs, 2):
        d = sp.Poly(poly(a) - poly(b), P)
        if d.is_zero:
            continue
        for r in d.real_roots():
            if 0 < r < 1:
                brk.add(r)
    brk = sorted(brk, key=lambda r: float(r))
    for lo, hi in zip(brk, brk[1:]):
        mid = (lo + hi) / 2
        best = max(cs, key=lambda c: poly(c).subs(P, mid))
        print(f"  on ({float(lo):.6f}, {float(hi):.6f}) worst inputs {' '.join(polys[best])} counts {list(best)}")


def main():
    c = [1, 6, 13, 15, 11, 5, 1]
    e = poly(c)
    print("E[1,6,13,15,11,5,1] - (1 - p^2 q (1+q)) =", sp.simplify(e - (1 - P ** 2 * Qs * (1 + Qs))))
    print("(1-p^2) - E =", sp.factor(1 - P ** 2 - e))
    for path in sys.argv[1:]:
        analyse(path)


if __name__ == "__main__":
    main()
