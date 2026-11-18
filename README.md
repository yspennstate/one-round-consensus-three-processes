# Optimal one-round consensus for three processes through stochastic channels

Yitzchak Shmalo, Einstein Institute of Mathematics, The Hebrew University of Jerusalem.

Code, outputs and paper source for the paper of the same title. Three synchronous processes exchange one round of
messages; each message on each of the six directed links is delivered independently with probability p, and each
process decides a value it has seen. For binary inputs the paper shows that the courteous algorithm M of
Fraigniaud, Patt-Shamir and Rajsbaum is optimal among deterministic algorithms for every p. For three input values
it gives an upper bound U(p) from three explicit rules and checks U(p) = OPT_3(p) exactly at 26 values of p.

All counts are exact integers and all reported values exact rationals.

## Layout

- `paper/main.tex` - the paper.
- `code/` - the programs. Run them from the repository root; outputs go to `out/`.
- `out/` - every output the paper uses, and the optimal rule found at each p (`rule_n3_q*_p*.json`).

## What each output is and how it was made

| Output | Command | What it shows |
| --- | --- | --- |
| `out/bin_grid.txt` | `python code/cons3.py 1/10 1/5 1/4 1/3 2/5 1/2 3/5 2/3 7/10 4/5 9/10` | binary, all 2^30 deterministic algorithms enumerated at 11 values of p; M (index 789) is the unique minimizer |
| `out/bin_ends.txt` | `python code/cons3.py 1/100 1/20 19/20 99/100` | the same near the ends of [0, 1] |
| `out/cert4.txt` | `python code/cert4.py` | Theorem 1 for every p: each algorithm has an input whose failure counts dominate those of M (Bernstein form); residue 0, `CERT4 ALL OK` |
| `out/cert4_neg.txt` | `python code/cert4.py neg` | negative control: target raised by one at k = 3; must and does fail |
| `out/selfcheck.txt` | `python code/cert3.py selfcheck 300` | tabulated counts against direct simulation, 300 random algorithms, 0 mismatches |
| `out/cert3.txt` | `python code/cert3.py cert` | the uniform-weight certificate, which fails by design (kept as a record) |
| `out/milp_q2.txt`, `out/milp_q3.txt` | `python code/milp3.py 3 2 1/5 1/2 4/5`; `python code/milp3.py 3 3 1/10 1/5 1/3 1/2 2/3 4/5 9/10` | mixed integer program (HiGHS, floating point); the rule found is re-evaluated exactly |
| `out/exact_q2.txt` | `python code/exact3.py 3 2 1/5 1/2 4/5` | CP-SAT with integer weights; the bound one below the optimum is infeasible |
| `out/exact_q3.txt` | `python code/exact3.py 3 3 1/10 1/5 1/4 1/3 2/5 9/20 1/2 2/3 4/5 9/10` | exact OPT_3 at the ten values in the table of the paper |
| `out/exact_q3_grid2.txt` | `python code/exact3.py 3 3 1/20 3/20 3/10` | further grid points |
| `out/exact_q3_ab_a.txt` .. `_d.txt` | `python code/exact3.py 3 3 6/25`, `37/150`, `25/102`, `38/155` (one per file) | points on either side of the crossing of s and t at 1 - 1/rho |
| `out/exact_q3_switch.txt` | `python code/exact3.py 3 3 7/20 3/8 19/50 77/200 39/100` | points around the crossing of t and m |
| `out/exact_q3_switch2.txt` | `python code/exact3.py 3 3 8/21 21/55 13/34 34/89` | points on either side of (3 - sqrt 5)/2 |
| `out/envelope_check.txt` | `python code/envelope_check.py out` | every exact OPT_3 above against U(p) from its closed form, with the rule attaining U; 26 points, `ALL OK` |
| `out/satcheck.txt` | `python code/satcheck.py` | a separate program (no shared code) with the SAT solver CaDiCaL through PySAT: Theorem 1 UNSAT, its control SAT; at p = 1/2, Q = 2: K = 35 UNSAT, 36 SAT; Q = 3: K = 47 UNSAT, 48 SAT |
| `out/rules3.txt` | `python code/rules3.py out/rule_n3_q3_p1_10.json out/rule_n3_q3_p1_5.json out/rule_n3_q3_p1_4.json out/rule_n3_q3_p1_3.json out/rule_n3_q3_p2_5.json out/rule_n3_q3_p1_2.json` | descriptions of the optimal three-value rules and their error on every input |
| `out/ruleAB.txt` | `python code/ruleAB.py` | the error polynomials of the rules s, t and m of the paper (A, B and m in the code) on all 24 non-constant inputs, and their crossings |
| `out/sym3_anon.txt`, `out/sym2_anon.txt` | `python code/sym3.py 3 anon 1/10 1/5 1/4 1/3 2/5 1/2`; `python code/sym3.py 2 anon 1/5 1/2` | the same optimization restricted to anonymous rules |

`code/ec_check.py` (output `out/ec_check.txt`) checks the closed form of the error of M and where the shared coin beats it; `code/satdecode.py` decodes a rule from a SAT model.

## Requirements

Two Python environments on an Apple M6 Mac mini were used:

- `cons3.py`, `cert3.py`, `cert4.py`, `milp3.py`: Python 3.13.18, NumPy 2.5.4, SymPy 1.14.0, SciPy 1.18.1.
- `exact3.py`, `satcheck.py`, `rules3.py`, `sym3.py`, `ruleAB.py`: Python 3.9.6, OR-Tools 9.15.6755, python-sat 1.9.dev23,
  SciPy 1.13.1, SymPy 1.14.0, NumPy 2.0.2.
- `envelope_check.py` needs only the Python standard library.

Each CP-SAT solve at a single p took a few seconds and under 200 MB of memory.

## Reference

The model and the two-process optimum are due to P. Fraigniaud, B. Patt-Shamir and S. Rajsbaum, Coordination through
stochastic channels, DISC 2025, LIPIcs 356, paper 32, doi:10.4230/LIPIcs.DISC.2025.32. The broadcast variant is in
their DISC 2026 paper, Consensus with stochastic broadcast, LIPIcs 397, 18:1-18:17, doi:10.4230/LIPIcs.DISC.2026.18.
