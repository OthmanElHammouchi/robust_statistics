import sympy
import sympy.stats

x, c = sympy.symbols("x c")

tukey_psi = x * (1 - x**2 / c**2) ** 2
Z = sympy.stats.Normal("z", 0, 1)

scale_fac = sympy.N(
    sympy.integrate(tukey_psi.diff(x) * sympy.stats.density(Z)(x), (x, -c, c)).subs(
        c, 4.68
    ),
    5,
) # 0.7573

print(sympy.latex(tukey_psi.diff(x) * sympy.stats.density(Z)(x), fold_short_frac = True))