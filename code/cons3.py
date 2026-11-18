"""Exact one-round binary consensus for three processes with independent links.

Model: processes 0,1,2 with inputs in {0,1}. In the single round every process
sends its input to the other two; each of the six directed messages arrives
independently with probability p. A deterministic algorithm maps each view
(own input, and for each other process its input or nothing) to a decision.
Validity: if all inputs are equal, only that value may be decided; since a
process that has seen no value other than its own cannot exclude that all
inputs are equal, it must decide its own input there. Views that contain the
other value are free (one bit each, ten per process).

Error of an algorithm = max over the six non-constant inputs of the probability
that not all three decisions are equal. We compute
    OPT(p) = min over all 2^30 algorithms of the error,
exactly, at rational p = a/b, in integers (values are scaled by b^6).

usage: python cons3.py a/b [a/b ...]
"""
import sys
import itertools
from fractions import Fraction
import numpy as np

N = 3
OTHERS = {i: [j for j in range(N) if j != i] for i in range(N)}
LINKS = [(j, i) for i in range(N) for j in OTHERS[i]]  # (sender, receiver)


def free_views(i):
    out = []
    for v in (0, 1):
        for sa in (None, 0, 1):
            for sb in (None, 0, 1):
                if (1 - v) in (sa, sb):
                    out.append((v, sa, sb))
    return out


FREE = {i: free_views(i) for i in range(N)}
FIDX = {i: {fv: k for k, fv in enumerate(FREE[i])} for i in range(N)}
NB = len(FREE[0])  # 10
INPUTS = [x for x in itertools.product((0, 1), repeat=N) if len(set(x)) > 1]


def view(i, x, pat):
    s = []
    for j in OTHERS[i]:
        li = LINKS.index((j, i))
        s.append(x[j] if (pat >> li) & 1 else None)
    return (x[i], s[0], s[1])


def restricted(x, i):
    """free-view indices of process i reachable under input x"""
    r = []
    for sa in (None, x[OTHERS[i][0]]):
        for sb in (None, x[OTHERS[i][1]]):
            fv = (x[i], sa, sb)
            if fv in FIDX[i]:
                r.append(FIDX[i][fv])
    return sorted(r)


def table(x):
    """counts c[a0,a1,a2,k]: number of delivery patterns with k delivered links
    on which decisions disagree, for each restricted assignment"""
    R = [restricted(x, i) for i in range(N)]
    shape = [1 << len(R[i]) for i in range(N)]
    T = np.zeros(shape + [len(LINKS) + 1], dtype=np.int64)
    pats = []
    for pat in range(1 << len(LINKS)):
        vs = [view(i, x, pat) for i in range(N)]
        pats.append((bin(pat).count("1"), vs))
    for a in itertools.product(*[range(s) for s in shape]):
        for k, vs in pats:
            d = []
            for i in range(N):
                fv = vs[i]
                if fv in FIDX[i]:
                    pos = R[i].index(FIDX[i][fv])
                    d.append((a[i] >> pos) & 1)
                else:
                    d.append(fv[0])
            if len(set(d)) > 1:
                T[a + (k,)] += 1
    return R, T


def ridx(R_i):
    """map a full 10-bit algorithm to its restricted index"""
    out = np.zeros(1 << NB, dtype=np.int64)
    for alg in range(1 << NB):
        r = 0
        for pos, b in enumerate(R_i):
            r |= ((alg >> b) & 1) << pos
        out[alg] = r
    return out


def main():
    ps = [Fraction(s) for s in sys.argv[1:]] or [Fraction(1, 2)]
    tabs = []
    for x in INPUTS:
        R, T = table(x)
        tabs.append((x, R, T, [ridx(R[i]) for i in range(N)]))
    L = len(LINKS)
    for p in ps:
        a, b = p.numerator, p.denominator
        w = np.array([a ** k * (b - a) ** (L - k) for k in range(L + 1)], dtype=object)
        assert int(b ** L) < 2 ** 62
        w = w.astype(np.int64)
        vals = [(x, (T * w).sum(axis=-1), rid) for (x, R, T, rid) in tabs]
        best = None
        nbest = 0
        arg = None
        for a0 in range(1 << NB):
            err = None
            for x, V, rid in vals:
                V2 = V[rid[0][a0]]  # shape (s1, s2)
                E = V2[rid[1][:, None], rid[2][None, :]]
                err = E if err is None else np.maximum(err, E)
            m = int(err.min())
            if best is None or m < best:
                best, nbest = m, 0
                idx = np.argwhere(err == m)[0]
                arg = (a0, int(idx[0]), int(idx[1]))
            if m == best:
                nbest += int((err == m).sum())
        val = Fraction(best, b ** L)
        print(f"p={p} OPT={val} ({float(val):.6f}) minimizers={nbest} example={arg}")
        # per-input polynomial counts of the example minimizer
        for x, R, T, rid in tabs:
            c = T[rid[0][arg[0]], rid[1][arg[1]], rid[2][arg[2]]]
            print(f"   input {''.join(map(str, x))}: counts by delivered links {list(map(int, c))}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
