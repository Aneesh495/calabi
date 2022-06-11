#!/usr/bin/env python3
"""
Python verification harness for Calabi-Yau quintic threefold combinatorics,
Picard-Fuchs Picard-Fuchs differential system, and Gromov-Witten invariants.
"""

import math
from fractions import Fraction

def factorial(n):
    return math.factorial(n)

def c_n(n):
    return factorial(5 * n) // (factorial(n) ** 5)

def harmonic(n):
    s = Fraction(0, 1)
    for k in range(1, n + 1):
        s += Fraction(1, k)
    return s

def d_n(n):
    cn = c_n(n)
    h_diff = harmonic(5 * n) - harmonic(n)
    return 5 * cn * h_diff

def verify_coefficients():
    print("=======================================================")
    print("  Python Verification: Picard-Fuchs Series Expansion   ")
    print("=======================================================")
    for n in range(0, 5):
        val_c = c_n(n)
        val_d = d_n(n) if n > 0 else 0
        print(f"  n={n}: c_{n} = {val_c:<15} d_{n} = {val_d}")

    # Ground truth checks
    assert c_n(0) == 1
    assert c_n(1) == 120
    assert c_n(2) == 113400
    assert c_n(3) == 168168000
    assert d_n(1) == 770
    assert d_n(2) == 810225
    print("\n[OK] Combinatorial coefficients match analytical values identically.")

def multiply_series(a, b, order):
    c = [0.0] * (order + 1)
    for i in range(order + 1):
        for j in range(order + 1 - i):
            c[i + j] += a[i] * b[j]
    return c

def invert_series(a, order):
    b = [0.0] * (order + 1)
    b[0] = 1.0 / a[0]
    for n in range(1, order + 1):
        s = sum(a[k] * b[n - k] for k in range(1, n + 1))
        b[n] = -s / a[0]
    return b

def substitute_series(f, g, order):
    res = [0.0] * (order + 1)
    res[0] = f[0]
    g_pow = [0.0] * (order + 1)
    g_pow[0] = 1.0
    for n in range(1, order + 1):
        g_pow = multiply_series(g_pow, g, order)
        if n < len(f):
            fn = f[n]
            for k in range(order + 1):
                res[k] += fn * g_pow[k]
    return res

def compute_instantons(order=4):
    print("\n=======================================================")
    print("  Python Verification: Mirror Map Series Reversion     ")
    print("=======================================================")

    c = [float(c_n(i)) for i in range(order + 5)]
    d = [float(d_n(i)) if i > 0 else 0.0 for i in range(order + 5)]

    # R(z) = S_1(z) / varpi_0(z)
    r = [0.0] * (order + 5)
    for n in range(1, order + 5):
        s = d[n] - sum(r[k] * c[n - k] for k in range(1, n))
        r[n] = s / c[0]

    # E(z) = exp(R(z))
    e = [0.0] * (order + 5)
    e[0] = 1.0
    for n in range(1, order + 5):
        e[n] = sum(k * r[k] * e[n - k] for k in range(1, n + 1)) / n

    q_series = [0.0] * (order + 5)
    for k in range(1, order + 5):
        q_series[k] = e[k - 1]

    # Lagrange inversion: z(q)
    z_series = [0.0] * (order + 5)
    z_series[1] = 1.0

    for n in range(2, order + 5):
        z_pow = [[0.0] * (order + 5) for _ in range(n + 1)]
        z_pow[0][0] = 1.0
        for p in range(1, n + 1):
            for i in range(p, order + 5):
                z_pow[p][i] = sum(z_series[j] * z_pow[p - 1][i - j] for j in range(1, i - p + 2))
        s = sum(q_series[k] * z_pow[k][n] for k in range(2, n + 1))
        z_series[n] = -s

    print(f"  Inverse mirror map z(q):")
    for k in range(1, 5):
        print(f"    z_{k} = {z_series[k]:.6f}")

    # dz/dq
    dz_dq = [(k + 1) * z_series[k + 1] for k in range(order + 1)]

    # q / z(q)
    z_norm = [z_series[k + 1] for k in range(order + 1)]
    inv_z_norm = invert_series(z_norm, order)

    # G = (q / z) * (dz / dq)
    g = multiply_series(inv_z_norm, dz_dq, order)
    g2 = multiply_series(g, g, order)
    g3 = multiply_series(g2, g, order)

    # 1 / varpi_0(z(q))^2
    varpi0_z = substitute_series(c, z_series, order)
    varpi0_sq = multiply_series(varpi0_z, varpi0_z, order)
    inv_varpi0_sq = invert_series(varpi0_sq, order)

    # 1 / (1 - 3125 * z(q))
    denom = [1.0] + [-3125.0 * z_series[k] for k in range(1, order + 1)]
    inv_denom = invert_series(denom, order)

    # C_{ttt} = 5 * inv_denom * inv_varpi0_sq * g3
    pre = multiply_series(inv_denom, inv_varpi0_sq, order)
    cttt = multiply_series(pre, g3, order)
    cttt = [5.0 * x for x in cttt]

    print(f"\n  Topological Yukawa Coupling C_{{ttt}}(q):")
    print(f"    C_0 = {cttt[0]:.4f} (Classical intersection number = 5)")
    for k in range(1, order + 1):
        print(f"    A_{k} = {cttt[k]:.1f}")

    # Mobius inversion
    n_counts = [0] * (order + 1)
    for k in range(1, order + 1):
        ak = cttt[k]
        div_sum = sum(n_counts[d] * (d ** 3) for d in range(1, k) if k % d == 0)
        nk = round((ak - div_sum) / (k ** 3))
        n_counts[k] = nk

    print(f"\n  Extracted Instanton Numbers n_d:")
    print(f"    n_1 (Lines)          = {n_counts[1]} (Expected: 2875)")
    print(f"    n_2 (Conics)         = {n_counts[2]} (Expected: 609250)")
    print(f"    n_3 (Twisted Cubics) = {n_counts[3]} (Expected: 317206375)")
    print(f"    n_4 (Quartics)       = {n_counts[4]} (Expected: 242467530000)")

    assert n_counts[1] == 2875
    assert n_counts[2] == 609250
    assert n_counts[3] == 317206375
    assert n_counts[4] == 242467530000
    print("\n[OK] All Gromov-Witten invariants verified with 100% precision.")

if __name__ == "__main__":
    verify_coefficients()
    compute_instantons(4)
