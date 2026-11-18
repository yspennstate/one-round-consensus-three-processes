"""All-p lower bound for one-round binary consensus, n = 3, independent links.

Claim: for every deterministic algorithm A and every p in [0,1],
    (1/6) sum_x E_x(A)(p) >= E_c(p),
where x runs over the six non-constant inputs and E_c is the error of the courteous
rule (the same on every input). Hence max_x E_x(A) >= E_c and OPT = E_c.

Each E_x(A)(p) = sum_k c_k p^k q^(6-k), c_k = number of delivery patterns with k of
the six links delivered on which the decisions disagree (tables from cons3.py).
D(A) = sum_x c(A,x) - 6 c_court is a vector of 7 integers; D >= 0 componentwise is a
sufficient condition (Bernstein form). The vectors with a negative entry are
collected (distinct) and each is checked exactly on [0,1] with sympy.

modes:
  python cert3.py selfcheck N    random algorithms: table lookups vs direct evaluation
  python cert3.py cert           the certificate over all 2^30 algorithms
"""
import sys
import random
import numpy as np
import cons3 as C

COURTEOUS = 789  # bits: own 0 -> (-,1):1 (0,1):0 (1,-):1 (1,0):0 (1,1):1 ; own 1 -> (-,0):0 (0,-):0 (0,0):0 (0,1):1 (1,0):1


def decide(alg, i, vw):
    if vw in C.FIDX[i]:
        return (alg >> C.FIDX[i][vw]) & 1
    return vw[0]


def direct_counts(algs, x):
    cnt = [0] * (len(C.LINKS) + 1)
    for pat in range(1 << len(C.LINKS)):
        d = [decide(algs[i], i, C.view(i, x, pat)) for i in range(C.N)]
        if len(set(d)) > 1:
            cnt[bin(pat).count("1")] += 1
    return cnt


def build():
    tabs = []
    for x in C.INPUTS:
        R, T = C.table(x)
        tabs.append((x, T, [C.ridx(R[i]) for i in range(C.N)]))
    return tabs


def selfcheck(n):
    tabs = build()
    rnd = random.Random(20261118)
    bad = 0
    for t in range(n):
        algs = [rnd.randrange(1 << C.NB) for _ in range(C.N)] if t else [COURTEOUS] * 3
        for x, T, rid in tabs:
            a = T[rid[0][algs[0]], rid[1][algs[1]], rid[2][algs[2]]]
            b = direct_counts(algs, x)
            if list(map(int, a)) != b:
                bad += 1
    print(f"selfcheck {n} random algorithms x 6 inputs: mismatches {bad}")
    for x, T, rid in tabs:
        print("courteous", "".join(map(str, x)), list(map(int, T[rid[0][COURTEOUS], rid[1][COURTEOUS], rid[2][COURTEOUS]])))


def poly_nonneg_on_unit(d):
    """exact: sum_k d_k p^k (1-p)^(6-k) >= 0 for all p in [0,1]"""
    import sympy as sp
    p = sp.symbols("p")
    L = len(d) - 1
    P = sp.Poly(sum(int(dk) * p ** k * (1 - p) ** (L - k) for k, dk in enumerate(d)), p)
    if P.is_zero:
        return True, "zero"
    # sign changes inside (0,1) only at roots of odd multiplicity
    pts = [sp.Rational(0), sp.Rational(1)]
    _, facs = P.sqf_list()
    for f, mult in facs:
        if mult % 2 == 1:
            for r in sp.Poly(f, p).intervals():
                (a, b), m = r
                pts += [sp.Rational(a), sp.Rational(b)]
    # check sign on a dense set of rationals between isolating endpoints, plus all odd roots
    odd_inside = 0
    for f, mult in facs:
        if mult % 2 == 1:
            odd_inside += sp.Poly(f, p).count_roots(0, 1)
            # roots at the endpoints are counted by count_roots(0,1); subtract them
            for e in (0, 1):
                if sp.Poly(f, p).eval(e) == 0:
                    odd_inside -= 1
    if odd_inside == 0:
        # constant sign on (0,1): test one point that is not a root
        for t in (sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(2, 7)):
            v = P.eval(t)
            if v != 0:
                return v > 0, f"no odd root in (0,1); P({t})={v}"
        return False, "undecided"
    return False, f"{odd_inside} odd-multiplicity roots in (0,1)"


def cert():
    tabs = build()
    court = np.zeros(len(C.LINKS) + 1, dtype=np.int64)
    for x, T, rid in tabs:
        c = T[rid[0][COURTEOUS], rid[1][COURTEOUS], rid[2][COURTEOUS]]
        court += c
    print("6 x courteous counts:", list(map(int, court)))
    nb = 1 << C.NB
    i1 = np.arange(nb)
    failing = {}
    nfail = 0
    nzero = 0
    for a0 in range(nb):
        S = None
        for x, T, rid in tabs:
            V = T[rid[0][a0]].astype(np.int16)  # (s1, s2, 7)
            E = V[rid[1][:, None], rid[2][None, :], :]
            S = E if S is None else S + E
        D = S.astype(np.int32) - court.astype(np.int32)[None, None, :]
        neg = (D < 0).any(axis=-1)
        nzero += int((D == 0).all(axis=-1).sum())
        if neg.any():
            rows = D[neg]
            nfail += rows.shape[0]
            u = np.unique(rows, axis=0)
            for r in map(tuple, u):
                failing[r] = failing.get(r, 0) + 1
        if a0 % 128 == 127:
            print(f"  a0 {a0 + 1}/{nb}: algorithms with a negative Bernstein entry so far {nfail}, distinct vectors {len(failing)}")
            sys.stdout.flush()
    print(f"algorithms with D == 0 identically: {nzero}")
    print(f"algorithms with a negative Bernstein entry: {nfail}; distinct D vectors: {len(failing)}")
    ok_all = True
    nbad = 0
    for r in sorted(failing):
        ok, why = poly_nonneg_on_unit(r)
        if not ok:
            ok_all = False
            nbad += 1
            if nbad <= 20:
                print("  NEGATIVE somewhere:", r, why)
    print("CERT ALL OK" if ok_all else f"CERT FAILED on {nbad} distinct vectors")


if __name__ == "__main__":
    if sys.argv[1] == "selfcheck":
        selfcheck(int(sys.argv[2]))
    else:
        cert()
