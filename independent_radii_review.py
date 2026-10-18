#!/usr/bin/env python3
"""Supplementary independent checks for PROOF.md; no core-checker imports.

Alexander Richard Towell, with substantial OpenAI Codex assistance.
Ricardo Pinho's axes, construction and normalization are predecessor work.
This AI review is not human peer review or a formal proof. Polynomial
identities are exact; finite geometry probes supplement the written argument.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


# Exact sparse polynomials in a fixed list of symbolic variables. This checks
# the identities as polynomial identities, rather than just at sampled points.
NAMES = ("u", "v", "q", "rho", "sigma", "a", "f", "tau")
Monomial = tuple[int, ...]
Polynomial = dict[Monomial, Q]


def constant(value: int | Q) -> Polynomial:
    return {} if value == 0 else {(0,) * len(NAMES): Q(value)}


def variable(name: str) -> Polynomial:
    powers = [0] * len(NAMES)
    powers[NAMES.index(name)] = 1
    return {tuple(powers): Q(1)}


def add(*polynomials: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for polynomial in polynomials:
        for powers, coefficient in polynomial.items():
            result[powers] = result.get(powers, Q(0)) + coefficient
    return {powers: c for powers, c in result.items() if c}


def scale(polynomial: Polynomial, coefficient: int | Q) -> Polynomial:
    return {powers: c * coefficient for powers, c in polynomial.items()
            if c * coefficient}


def multiply(*polynomials: Polynomial) -> Polynomial:
    result = constant(1)
    for polynomial in polynomials:
        product: Polynomial = {}
        for left, lc in result.items():
            for right, rc in polynomial.items():
                powers = tuple(a + b for a, b in zip(left, right))
                product[powers] = product.get(powers, Q(0)) + lc * rc
        result = {powers: c for powers, c in product.items() if c}
    return result


def derivative(polynomial: Polynomial, name: str) -> Polynomial:
    index = NAMES.index(name)
    result: Polynomial = {}
    for powers, coefficient in polynomial.items():
        if powers[index]:
            reduced = list(powers)
            reduced[index] -= 1
            result[tuple(reduced)] = coefficient * powers[index]
    return result


def polynomial_checks() -> int:
    u, v, q, rho, sigma, a, f, tau = map(variable, NAMES)
    one = constant(1)
    minus_one = constant(-1)
    d = add(one, multiply(u, u))
    w = add(u, scale(multiply(v, q), -1))
    D = add(one, multiply(u, w))
    alpha_n = add(d, multiply(v, q, add(one, scale(u, -1))))
    beta_n = add(d, scale(multiply(v, q, add(one, u)), -1))
    constraint = add(multiply(alpha_n, rho), multiply(beta_n, sigma), scale(d, -1))
    x_n = scale(multiply(q, add(sigma, multiply(u, rho))), -1)
    y_n = multiply(q, add(rho, scale(multiply(u, sigma), -1)))
    t_n = add(multiply(d, rho), multiply(v, x_n))
    H = add(one, scale(q, -1), multiply(u, w))
    P = add(scale(D, 4), scale(multiply(d, add(one, w)), -3))
    coefficient = add(scale(u, 4), constant(-3), scale(multiply(u, u), -3))
    p_tau = add(one, scale(multiply(f, tau), -3),
                multiply(add(scale(multiply(a, a), -3), scale(multiply(a, f), 4)), tau, tau),
                scale(multiply(a, a, f, tau, tau, tau), -3))
    p_tau_direct = add(scale(add(one, multiply(a, f, tau, tau)), 4),
                       scale(multiply(add(one, multiply(a, a, tau, tau)),
                                      add(one, multiply(f, tau))), -3))
    identities = {
        "alpha_plus_beta": add(alpha_n, beta_n, scale(D, -2)),
        "first_endpoint": add(multiply(alpha_n, u), beta_n,
                              scale(multiply(d, add(one, w)), -1)),
        "lower_angle": add(y_n, scale(multiply(u, x_n), -1),
                           scale(multiply(q, rho, d), -1)),
        "lower_radius": add(t_n, scale(multiply(v, x_n), -1),
                            scale(multiply(rho, d), -1)),
        "upper_angle": add(x_n, multiply(u, y_n), multiply(q, sigma, d)),
        "upper_radius_mod_constraint": add(d, scale(t_n, -1),
                                            scale(multiply(v, y_n), -1),
                                            scale(multiply(sigma, d), -1), constraint),
        "tetrahedron_lower": add(t_n, scale(y_n, -1), scale(multiply(H, rho), -1),
                                 scale(multiply(q, add(u, scale(v, -1)), sigma), -1)),
        "tetrahedron_upper_mod_constraint": add(d, scale(t_n, -1), x_n,
                                                scale(multiply(H, sigma), -1),
                                                multiply(q, add(u, scale(v, -1)), rho), constraint),
        "cutoff_polynomial": add(P, minus_one, scale(multiply(u, u), 3),
                                 scale(multiply(w, coefficient), -1)),
        "negative_w_coefficient": add(coefficient,
                                      scale(multiply(add(u, constant(Q(-2, 3))),
                                                     add(u, constant(Q(-2, 3)))), 3),
                                      constant(Q(5, 3))),
        "q_monotonicity": add(derivative(P, "q"), multiply(v, coefficient)),
        "tau_polynomial": add(p_tau, scale(p_tau_direct, -1)),
        "tau_derivative_upper_bound": add(scale(f, -3), scale(multiply(a, f, tau), 2),
                                           scale(derivative(p_tau, "tau"), -1),
                                           scale(multiply(a, add(a, scale(f, -1)), tau), -6),
                                           scale(multiply(a, a, f, tau, tau), -9)),
    }
    for name, residual in identities.items():
        require(not residual, f"Polynomial identity failed: {name}")
    return len(identities)


def phi(q: Q) -> Q:
    q2 = q * q
    return (1 - q2) * (1 - q2 / 200 + Q(27, 1250) * q2 * q2) / 4


@dataclass(frozen=True)
class Sector:
    index: int
    l: Q
    r: Q
    a: Q
    b: Q


def sectors() -> tuple[Sector, ...]:
    result = []
    for index in range(2048):
        l, r = Q(index, 1024) - 1, Q(index + 1, 1024) - 1
        b = (phi(l) - phi(r)) / (r - l)
        a = phi(l) + b * l
        result.append(Sector(index, l, r, a, b))
    return tuple(result)


def finite_hypotheses(all_sectors: tuple[Sector, ...]) -> tuple[Q, int]:
    margins = []
    for sector in all_sectors[1024:]:
        l, r, a, b = sector.l, sector.r, sector.a, sector.b
        require(0 < a < 1 and 0 <= b < 1 and 0 <= a-b <= 1-r,
                "Predecessor coefficient hypothesis failed")
        require(0 <= phi(l) < 1 and 0 <= phi(r) < 1, "Node shape hypothesis failed")
        reflected = all_sectors[2047-sector.index]
        require((reflected.l, reflected.r, reflected.a, reflected.b) == (-r, -l, a, -b),
                "Secant reflection failed")
        f = phi(l)
        E = 1 - 3*a*a + f*(4*a - 3 - 3*a*a)
        require(E >= 0, "Additional finite hypothesis failed")
        margins.append((E, sector.index))
    return min(margins)


def dot(left: tuple[Q, Q, Q], right: tuple[Q, Q, Q]) -> Q:
    return left[0]*right[0] + left[1]*right[1] + 2*left[2]*right[2]


def base_area_squared(sector: Sector, tau: Q, left_radius: Q,
                      right_radius: Q, upper: bool) -> Q:
    a, b, l, r = sector.a, sector.b, sector.l, sector.r
    if upper:
        left = (l*left_radius, Q(0), -left_radius)
        right = (r*right_radius, Q(0), -right_radius)
        axis = (-tau*a, Q(1), tau*b)
    else:
        left = (Q(0), l*left_radius, left_radius)
        right = (Q(0), r*right_radius, right_radius)
        axis = (Q(1), tau*a, tau*b)
    axis_norm = dot(axis, axis)
    gram_ll = dot(left, left) - dot(left, axis)**2 / axis_norm
    gram_rr = dot(right, right) - dot(right, axis)**2 / axis_norm
    gram_lr = dot(left, right) - dot(left, axis)*dot(right, axis)/axis_norm
    # Divide by 4 for triangle area squared, and by A_min(K)^2=2.
    return (gram_ll*gram_rr-gram_lr*gram_lr)/8


def geometric_probes(all_sectors: tuple[Sector, ...]) -> tuple[int, int]:
    nodes = tuple(Q(index, 1024)-1 for index in range(2049))
    node_values = tuple(phi(q) for q in nodes)
    geometry_count = 0
    area_count = 0
    for index in (1024, 1025, 1279, 1535, 2046, 2047):
        sector = all_sectors[index]
        a, b, l, r = sector.a, sector.b, sector.l, sector.r
        for tau in (Q(0), Q(1, 2), Q(1)):
            for upper in (False, True):
                Rl, Rr = Q(3, 7), Q(11, 13)
                expected = ((r-l)*Rl*Rr)**2/(4*(1+tau*tau*(a*a+2*b*b)))
                require(base_area_squared(sector, tau, Rl, Rr, upper) == expected,
                        "Projected-base normalization failed")
                area_count += 1
            for angular_fraction in (Q(1, 4), Q(3, 4)):
                q = l + angular_fraction*(r-l)
                u, v = tau*a, tau*b
                w, d = u-v*q, 1+u*u
                D = 1+u*w
                T, rho0, sigma0 = d/(2*D), u/(1+w), 1/(1+w)
                alpha = 1+v*q*(1-u)/d
                beta = 1-v*q*(1+u)/d
                require(0 < beta <= alpha and rho0 < T < sigma0,
                        "Segment ordering failed")
                require(sigma0 >= Q(3, 2)*T, "Large-cutoff bound failed")
                for segment_fraction in (Q(1, 5), Q(1, 2), Q(4, 5)):
                    rho = rho0 + segment_fraction*(T-rho0)
                    sigma = (1-alpha*rho)/beta
                    x = -q*(sigma+u*rho)/d
                    y = q*(rho-u*sigma)/d
                    t = rho+v*x
                    require(x < 0 < y and -x >= y and y <= t and -x <= 1-t,
                            "Test point is outside its asserted tetrahedron region")
                    lower = tuple(z*t+tau*f*x for z, f in zip(nodes, node_values))
                    upper = tuple(z*(1-t)-tau*f*y for z, f in zip(nodes, node_values))
                    lower_candidates = tuple(k for k in range(2048)
                                             if lower[k] <= y <= lower[k+1])
                    upper_candidates = tuple(k for k in range(2048)
                                             if upper[k] <= x <= upper[k+1])
                    require(lower_candidates == (index,), "Unexpected lower angular sector")
                    require(upper_candidates == (2047-index,), "Unexpected upper angular sector")
                    require((y-u*x, t-v*x, x+u*y, 1-t-v*y) ==
                            (q*rho, rho, -q*sigma, sigma), "Intercept identity failed")
                    geometry_count += 1
    return geometry_count, area_count


def failure_controls(all_sectors: tuple[Sector, ...]) -> dict[str, str]:
    # The predecessor's single diagonal point cannot imply a quadratic bound.
    A, B, T = Q(1, 4), Q(1, 2), Q(1, 2)
    require(max(A, B) >= T and A*A+B*B < 2*T*T,
            "Diagonal-only counterexample control failed")
    # The sum bound is conditional. This pair satisfies the entire segment
    # relaxation in both orientations but is NOT claimed to cover K.
    sector = all_sectors[-1]
    q = (sector.l+sector.r)/2
    u, v = sector.a, sector.b
    w = u-v*q
    T = (1+u*u)/(2*(1+u*w))
    rho0, sigma0 = u/(1+w), 1/(1+w)
    A, B = rho0/4, sigma0
    require(0 < A < rho0 and B >= T and A+B < 2*T and A*A+B*B >= 2*T*T,
            "Small-cutoff sum-bound counterexample control failed")
    return {
        "diagonal_only": "A=1/4, B=1/2, T=1/2: max>=T but A^2+B^2<2T^2",
        "unconditional_sum": "last-sector midpoint, tau=1, A=rho0/4, B=sigma0: A+B<2T",
        "interpretation": "These refute proposed implications, not the theorem; neither is a verified cover.",
    }


def main() -> None:
    identities = polynomial_checks()
    all_sectors = sectors()
    E, index = finite_hypotheses(all_sectors)
    geometry_count, area_count = geometric_probes(all_sectors)
    controls = failure_controls(all_sectors)
    result = {
        "verified": True,
        "scope": "Independent algebra and finite controls; not a full-cover or integral certificate replay",
        "source_commit": "ec5e88b8c954b8ce42ee71a993271e14dce1aded",
        "polynomial_identities": identities,
        "positive_sectors": 1024,
        "minimum_E": str(E),
        "minimum_E_sector": index,
        "rational_geometry_probes": geometry_count,
        "angular_sectors_examined_per_probe": 2048,
        "normalized_area_probes": area_count,
        "failure_controls": controls,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
