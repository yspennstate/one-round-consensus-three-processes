"""Decode the Q=3, K=47 SAT model from satcheck.part_b and re-evaluate the rule two ways."""
import itertools
from fractions import Fraction
import satcheck as S
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
import milp3 as M


def decode(Q, K):
    pool = IDPool()
    cl = []
    inputs = [x for x in itertools.product(range(Q), repeat=3) if len(set(x)) > 1]
    dv = {}
    for x in inputs:
        for pt in S.PATTERNS:
            for i in S.PROCS:
                v = S.view(i, x, pt)
                vs = S.values(v)
                if len(vs) >= 2 and (i, v) not in dv:
                    dv[(i, v)] = vs
                    lits = [pool.id(("d", i, v, a)) for a in vs]
                    cl.extend(CardEnc.equals(lits=lits, bound=1, vpool=pool, encoding=EncType.pairwise).clauses)
    for x in inputs:
        ys = []
        for pt in S.PATTERNS:
            y = pool.id(("y", x, pt))
            ys.append(y)
            opts = []
            for a in range(Q):
                lits = []
                ok = True
                for i in S.PROCS:
                    v = S.view(i, x, pt)
                    vs = S.values(v)
                    if len(vs) == 1:
                        if vs[0] != a:
                            ok = False
                    else:
                        if a not in vs:
                            ok = False
                        else:
                            lits.append(pool.id(("d", i, v, a)))
                if not ok:
                    continue
                z = pool.id(("z", x, pt, a))
                for l in lits:
                    cl.append([-z, l])
                opts.append(z)
            cl.append([-y] + opts)
        cl.extend(CardEnc.atleast(lits=ys, bound=64 - K, vpool=pool, encoding=EncType.seqcounter).clauses)
    s = Cadical153(bootstrap_with=cl)
    r = s.solve()
    if not r:
        return None
    mdl = set(l for l in s.get_model() if l > 0)
    rule = {}
    for (i, v), vs in dv.items():
        ch = [a for a in vs if pool.id(("d", i, v, a)) in mdl]
        assert len(ch) == 1
        rule[(i, v)] = ch[0]
    return rule


def sim(Q, rule):
    worst = 0
    for x in itertools.product(range(Q), repeat=3):
        if len(set(x)) < 2:
            continue
        bad = 0
        for pt in S.PATTERNS:
            ds = set()
            for i in S.PROCS:
                v = S.view(i, x, pt)
                vs = S.values(v)
                ds.add(vs[0] if len(vs) == 1 else rule[(i, v)])
            bad += len(ds) > 1
        worst = max(worst, bad)
    return worst


for Q, K in ((3, 47), (3, 46), (3, 44)):
    rule = decode(Q, K)
    if rule is None:
        print(f"Q={Q} K={K}: UNSAT")
        continue
    w = sim(Q, rule)
    others, links, lidx, inputs = M.setup(3, Q)
    # translate to milp3 view keys: (own, x[j] for j in others[i]) - same order as satcheck
    ev = M.exact_eval(3, Q, rule, inputs, (others, links, lidx))
    me = max(M.poly_at(c, Fraction(1, 2)) for _, c in ev)
    print(f"Q={Q} K={K}: SAT; direct-sim worst disagreements {w}/64; milp3 exact_eval at 1/2 = {me}", flush=True)
