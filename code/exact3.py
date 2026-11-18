"""Exact optimum of deterministic one-round consensus (n processes, q values,
independent links) at a rational p = a/b, by CP-SAT (integer arithmetic throughout).

At p = a/b a delivery pattern with k of the L links delivered has integer weight
a^k (b-a)^(L-k); the total weight is b^L. For each non-constant input x,
    t >= b^L - sum_pat w(pat) * [all processes decide the same value].
CP-SAT minimises the integer t; OPTIMAL status is a proof. Controls both ways:
the bound t <= OPT-1 must come back INFEASIBLE, t <= OPT FEASIBLE, and the rule
returned is re-evaluated by exact_eval from milp3.py. The rule is written as JSON.

usage: python exact3.py n q a/b [a/b ...]
"""
import sys
import json
from fractions import Fraction
from ortools.sat.python import cp_model
import milp3 as M


def model(n, q, p, cap=None):
    others, links, lidx, inputs = M.setup(n, q)
    L = len(links)
    a, b = p.numerator, p.denominator
    pats = range(1 << L)
    w = [a ** bin(pt).count("1") * (b - a) ** (L - bin(pt).count("1")) for pt in pats]
    m = cp_model.CpModel()
    d = {}
    views = {}
    for x in inputs:
        for pt in pats:
            for i in range(n):
                vw = M.view(i, x, pt, others, lidx)
                sv = M.seen(vw)
                if len(sv) >= 2 and (i, vw) not in views:
                    views[(i, vw)] = sv
                    for v in sv:
                        d[(i, vw, v)] = m.NewBoolVar("")
                    m.AddExactlyOne([d[(i, vw, v)] for v in sv])
    t = m.NewIntVar(0, b ** L, "t")
    for x in inputs:
        terms = []
        for pt in pats:
            vws = [M.view(i, x, pt, others, lidx) for i in range(n)]
            common = set(range(q))
            for vw in vws:
                common &= set(M.seen(vw))
            for v in sorted(common):
                y = m.NewBoolVar("")
                for i, vw in enumerate(vws):
                    if (i, vw, v) in d:
                        m.AddImplication(y, d[(i, vw, v)])
                terms.append((w[pt], y))
        m.Add(t + sum(c * y for c, y in terms) >= b ** L)
    if cap is None:
        m.Minimize(t)
    else:
        m.Add(t <= cap)
    return m, d, views, t, inputs, (others, links, lidx), b ** L


def run(n, q, p, cap=None, workers=4, tl=36000):
    m, d, views, t, inputs, geo, W = model(n, q, p, cap)
    s = cp_model.CpSolver()
    s.parameters.num_search_workers = workers
    s.parameters.max_time_in_seconds = tl
    st = s.Solve(m)
    rule = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        rule = {}
        for (i, vw), sv in views.items():
            for v in sv:
                if s.Value(d[(i, vw, v)]):
                    rule[(i, vw)] = v
    return s.StatusName(st), (s.Value(t) if rule else None), rule, inputs, geo, W


def main():
    n, q = int(sys.argv[1]), int(sys.argv[2])
    for sp in sys.argv[3:]:
        p = Fraction(sp)
        st, tv, rule, inputs, geo, W = run(n, q, p)
        ev = M.exact_eval(n, q, rule, inputs, geo)
        worst = max(M.poly_at(c, p) for _, c in ev)
        print(f"n={n} q={q} p={p}: {st} OPT={Fraction(tv, W)} exact_eval={worst} match={Fraction(tv, W) == worst}")
        st1, _, _, _, _, _ = run(n, q, p, cap=tv - 1)
        st0, _, _, _, _, _ = run(n, q, p, cap=tv)
        print(f"   control t<=OPT-1: {st1} (must be INFEASIBLE); t<=OPT: {st0} (must be FEASIBLE or OPTIMAL)")
        for x, c in ev:
            print(f"   input {''.join(map(str, x))} counts {c}")
        with open(f"out/rule_n{n}_q{q}_p{p.numerator}_{p.denominator}.json", "w", encoding="utf-8") as f:
            json.dump([[i, [None if s is None else s for s in vw], v] for (i, vw), v in sorted(rule.items(), key=str)], f)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
