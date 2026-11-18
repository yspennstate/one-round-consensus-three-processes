"""Optimal one-round rules under symmetry, n = 3 (CP-SAT, exact integers, built on exact3.model).

mode anon : a decision depends only on (own value, multiset of values heard), the same for every process.
mode cyc  : anon, and also commutes with the value shift v -> v+1 mod Q.
For each p prints OPT under the mode, with the controls (OPT-1 infeasible), and the decision table.
Compare OPT with the unrestricted optimum from exact3 (out/exact_q3*.txt).

usage: python sym3.py Q mode a/b [a/b ...]
"""
import sys
from fractions import Fraction
from ortools.sat.python import cp_model
import exact3 as E


def key(vw, v, Q, mode):
    own = vw[0]
    heard = tuple(sorted(h for h in vw[1:] if h is not None))
    if mode == "anon":
        return (own, heard, v)
    return (0, tuple(sorted((h - own) % Q for h in heard)), (v - own) % Q)


def build(Q, p, mode, cap=None):
    m, d, views, t, inputs, geo, W = E.model(3, Q, p, cap)
    first = {}
    for (i, vw, v), var in d.items():
        k = key(vw, v, Q, mode)
        if k in first:
            m.Add(var == first[k])
        else:
            first[k] = var
    return m, d, t, W, first


def solve(Q, p, mode, cap=None):
    m, d, t, W, first = build(Q, p, mode, cap)
    s = cp_model.CpSolver()
    s.parameters.num_search_workers = 4
    s.parameters.max_time_in_seconds = 36000
    st = s.Solve(m)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return s.StatusName(st), None, None, W
    table = {}
    for k, var in first.items():
        if s.Value(var):
            table[k[:2]] = k[2]
    return s.StatusName(st), s.Value(t), table, W


def main():
    Q, mode = int(sys.argv[1]), sys.argv[2]
    for sp in sys.argv[3:]:
        p = Fraction(sp)
        st, tv, table, W = solve(Q, p, mode)
        if tv is None:
            print(f"Q={Q} {mode} p={p}: {st}")
            continue
        st1 = solve(Q, p, mode, cap=tv - 1)[0]
        print(f"Q={Q} {mode} p={p}: {st} OPT={Fraction(tv, W)}  control OPT-1: {st1}")
        print("   " + "; ".join(f"own {k[0]} heard {list(k[1])} -> {v}" for k, v in sorted(table.items())))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
