"""Compare every exact q = 3 optimum in out/exact_q3*.txt with the upper bound U(p) of Proposition 2.

U(p) = min{ 1 - p^2 q^2 (1+q+q^2)  (rule s),  1 - p^2 q (1+q)  (rule t),  1 - p^2  (rule m) },  q = 1 - p.
All arithmetic is in exact rationals. Prints, for each p, the CP-SAT optimum, U(p), the rule(s) attaining U,
and whether they agree; ends with the number of distinct p and ALL OK or a list of mismatches.

usage: python code/envelope_check.py [out_dir]
"""
import glob
import os
import re
import sys
from fractions import Fraction as F

LINE = re.compile(r"^n=3 q=3 p=(\d+)/(\d+): OPTIMAL OPT=(\d+)/(\d+)")
CONTROL = re.compile(r"control t<=OPT-1: INFEASIBLE")


def terms(p):
    q = 1 - p
    return {"s": 1 - p**2 * q**2 * (1 + q + q**2), "t": 1 - p**2 * q * (1 + q), "m": 1 - p**2}


def main(out_dir):
    seen = {}
    for path in sorted(glob.glob(os.path.join(out_dir, "exact_q3*.txt"))):
        lines = open(path).read().splitlines()
        for i, ln in enumerate(lines):
            m = LINE.match(ln)
            if not m:
                continue
            p = F(int(m.group(1)), int(m.group(2)))
            opt = F(int(m.group(3)), int(m.group(4)))
            ctrl = i + 1 < len(lines) and CONTROL.search(lines[i + 1]) is not None
            if p in seen and seen[p][0] != opt:
                print(f"CONFLICT p={p}: {seen[p][0]} in {seen[p][1]} vs {opt} in {os.path.basename(path)}")
                seen[p] = (None, path, False)
                continue
            seen.setdefault(p, (opt, os.path.basename(path), ctrl))
    bad = []
    for p in sorted(seen):
        opt, src, ctrl = seen[p]
        t = terms(p)
        u = min(t.values())
        who = ",".join(k for k, v in t.items() if v == u)
        ok = opt is not None and opt == u and ctrl
        print(f"p={str(p):>7} ({float(p):.6f})  OPT={float(opt) if opt is not None else 'conflict'}  "
              f"U={float(u):.10f} by {who:<3} control={'yes' if ctrl else 'NO'}  {'OK' if ok else 'MISMATCH'}  [{src}]")
        if not ok:
            bad.append(p)
    print(f"distinct p: {len(seen)}")
    print("ALL OK" if not bad else "MISMATCH at " + ", ".join(map(str, bad)))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "out"))
