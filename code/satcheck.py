"""Independent check of the n = 3 results with a SAT solver (shares no code with cons3/milp3/cert4).

Part A (binary, Theorem 1, all p). A rule is a bit per (process, free view). On a fixed input only 7 free
views are reachable, so per input we simulate all 128 local assignments over the 64 delivery patterns and
call an assignment GOOD if its disagreement counts by delivered links dominate the target (courteous).
A counterexample to the theorem is a rule that is bad on all six inputs; we block every good assignment
of every input and ask the solver. Expect UNSAT. Control: target raised by one at k = 3 -> expect SAT.

Part B (exact error counts at p = 1/2, where every pattern has weight 1/64). Q values; for each K ask
whether some rule has at most K disagreeing patterns on every non-constant input. Expect Q=2: K=35 UNSAT,
K=36 SAT (9/16); Q=3: K=47 UNSAT, K=48 SAT (3/4).

usage: python satcheck.py
"""
import itertools
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

PROCS = (0, 1, 2)
LINKS = [(s, r) for r in PROCS for s in PROCS if s != r]  # (sender, receiver)


def heard(i, pat):
    """senders whose message to i is delivered under pattern pat (a frozenset of links)"""
    return [s for (s, r) in LINKS if r == i and (s, r) in pat]


def view(i, x, pat):
    h = set(heard(i, pat))
    return (x[i],) + tuple((x[j] if j in h else None) for j in PROCS if j != i)


def values(v):
    return sorted({a for a in v if a is not None})


PATTERNS = [frozenset(c) for k in range(7) for c in itertools.combinations(LINKS, k)]
assert len(PATTERNS) == 64


def courteous(v):
    own, a, b = v
    got = [u for u in (a, b) if u is not None]
    if not got:
        return own
    if len(got) == 1:
        return got[0]
    vals = [own] + got
    return max(set(vals), key=vals.count)


def part_a(target):
    pool = IDPool()
    s = Cadical153()
    inputs = [x for x in itertools.product((0, 1), repeat=3) if len(set(x)) > 1]
    for x in inputs:
        free = sorted({(i, view(i, x, pt)) for pt in PATTERNS for i in PROCS if len(values(view(i, x, pt))) == 2}, key=repr)
        assert len(free) == 7, (x, len(free))
        ngood = 0
        for bits in itertools.product((0, 1), repeat=7):
            dec = dict(zip(free, bits))
            cnt = [0] * 7
            for pt in PATTERNS:
                ds = set()
                for i in PROCS:
                    v = view(i, x, pt)
                    vs = values(v)
                    ds.add(vs[0] if len(vs) == 1 else dec[(i, v)])
                if len(ds) > 1:
                    cnt[len(pt)] += 1
            if all(c >= t for c, t in zip(cnt, target)):
                ngood += 1
                # block this assignment: bit b means "decide 1" on that view
                s.add_clause([(-pool.id(f) if b else pool.id(f)) for f, b in zip(free, bits)])
        print(f"  input {''.join(map(str, x))}: {ngood} of 128 local assignments dominate the target")
    r = s.solve()
    s.delete()
    return r


def part_b(Q, K):
    pool = IDPool()
    cl = []
    inputs = [x for x in itertools.product(range(Q), repeat=3) if len(set(x)) > 1]
    dvar = {}
    for x in inputs:
        for pt in PATTERNS:
            for i in PROCS:
                v = view(i, x, pt)
                vs = values(v)
                if len(vs) >= 2 and (i, v) not in dvar:
                    dvar[(i, v)] = vs
                    lits = [pool.id(("d", i, v, a)) for a in vs]
                    cl.extend(CardEnc.equals(lits=lits, bound=1, vpool=pool, encoding=EncType.pairwise).clauses)

    def dec_lit(i, v, a):
        vs = values(v)
        if len(vs) == 1:
            return True if vs[0] == a else False
        if a not in vs:
            return False  # strong validity: an unseen value cannot be decided
        return pool.id(("d", i, v, a))

    for x in inputs:
        ys = []
        for pt in PATTERNS:
            y = pool.id(("y", x, pt))
            ys.append(y)
            opts = []
            for a in range(Q):
                lits = [dec_lit(i, view(i, x, pt), a) for i in PROCS]
                if False in lits:
                    continue
                z = pool.id(("z", x, pt, a))
                for l in lits:
                    if l is not True:
                        cl.append([-z, l])
                opts.append(z)
            cl.append([-y] + opts)
        cl.extend(CardEnc.atleast(lits=ys, bound=64 - K, vpool=pool, encoding=EncType.seqcounter).clauses)
    s = Cadical153(bootstrap_with=cl)
    r = s.solve()
    s.delete()
    return r


def main():
    court = None
    for x in itertools.product((0, 1), repeat=3):
        if len(set(x)) < 2:
            continue
        cnt = [0] * 7
        for pt in PATTERNS:
            if len({courteous(view(i, x, pt)) for i in PROCS}) > 1:
                cnt[len(pt)] += 1
        assert court is None or cnt == court
        court = cnt
    print("courteous counts by direct simulation:", court)
    print("A: theorem (expect UNSAT = False):", part_a(court))
    raised = list(court)
    raised[3] += 1
    print("A control, target c3+1 (expect SAT = True):", part_a(raised))
    for Q, K in ((2, 35), (2, 36), (3, 47), (3, 48)):
        print(f"B: Q={Q} p=1/2 error <= {K}/64 feasible:", part_b(Q, K), flush=True)


if __name__ == "__main__":
    main()
