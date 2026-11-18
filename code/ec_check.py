"""Checks used in the paper: E_c from the counts equals the closed form; where the shared coin (error q) beats E_c."""
import sympy as sp

p = sp.symbols("p")
q = 1 - p
counts = [1, 4, 9, 12, 8, 2, 0]
Ec = sum(c * p ** k * q ** (6 - k) for k, c in enumerate(counts))
closed = 1 - (1 - q ** 2) * (1 - p * q) ** 2 - p ** 2 * q ** 4
print("closed form matches counts:", sp.expand(Ec - closed) == 0)
f = sp.factor(sp.expand(Ec - q))
print("E_c - q =", f)
for r in sp.Poly(sp.expand(Ec - q), p).real_roots():
    print("root", r, float(r))
print("E_c(1/2) =", Ec.subs(p, sp.Rational(1, 2)), " E_c(2/5) =", Ec.subs(p, sp.Rational(2, 5)))
