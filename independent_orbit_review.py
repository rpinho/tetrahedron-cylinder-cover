#!/usr/bin/env python3
"""Independent exact controls for the arbitrary-profile four-orbit proof.

Alexander Richard Towell, with substantial OpenAI Codex assistance.
Credits: Ricardo Pinho's fixed axes and normalization; Towell's preceding
shared/even-profile audits with AI assistance. This is a separate AI review,
not human peer review or a formal proof. No production checker is imported.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
import hashlib
from itertools import product
import json
from pathlib import Path


class ReviewFailure(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ReviewFailure(message)


Point = tuple[Q, Q, Q]
Profile = tuple[Q, ...]


def phi(q: Q) -> Q:
    q2 = q*q
    return (1-q2)*(1-q2/200+Q(27, 1250)*q2*q2)/4


@dataclass(frozen=True)
class Sector:
    index: int
    l: Q
    r: Q
    a: Q
    b: Q


@dataclass(frozen=True)
class Profiles:
    lower: Profile
    upper: Profile

    def rotate(self) -> Profiles:
        return Profiles(self.upper, self.lower[::-1])


@dataclass(frozen=True)
class Cylinder:
    vertices: frozenset[Point]
    axis: Point


@dataclass(frozen=True)
class Segment:
    sector: Sector
    tau: Q
    q: Q
    alpha: Q
    beta: Q
    rho0: Q
    sigma0: Q
    T: Q

    def point(self, rho: Q) -> tuple[Point, Q]:
        sigma = (1-self.alpha*rho)/self.beta
        u, v = self.tau*self.sector.a, self.tau*self.sector.b
        d = 1+u*u
        x = -self.q*(sigma+u*rho)/d
        y = self.q*(rho-u*sigma)/d
        return (x, y, rho+v*x), sigma

    def uncovered_interval(self, x: Q, y: Q) -> tuple[Q, Q] | None:
        # Both radial failures are strict. The open interval is nonempty
        # exactly when max(rho0,x) < min(T,(1-beta*y)/alpha).
        low = max(self.rho0, x)
        high = min(self.T, (1-self.beta*y)/self.alpha)
        return (low, high) if low < high else None

    def holds(self, x: Q, y: Q) -> bool:
        return self.uncovered_interval(x, y) is None


def make_sectors() -> tuple[Sector, ...]:
    result = []
    for j in range(2048):
        l, r = Q(j, 1024)-1, Q(j+1, 1024)-1
        b = (phi(l)-phi(r))/(r-l)
        result.append(Sector(j, l, r, phi(l)+b*l, b))
    return tuple(result)


def make_profiles() -> Profiles:
    nodes = tuple(Q(j, 1024)-1 for j in range(2049))
    profiles = Profiles(tuple(Q(2, 5)+q/17+q*q/23 for q in nodes),
                        tuple(Q(3, 4)-q/19+q*q*q/29 for q in nodes))
    require(all(value > 0 for p in (profiles.lower, profiles.upper) for value in p),
            "Control profiles must be positive")
    require(profiles.lower != profiles.lower[::-1] and
            profiles.upper != profiles.upper[::-1], "Control profiles must be non-even")
    return profiles


def profile_orbit(profiles: Profiles) -> tuple[Profiles, ...]:
    orbit = []
    current = profiles
    for _ in range(4):
        orbit.append(current)
        current = current.rotate()
    require(current == profiles, "Fourth rotation did not recover the profiles")
    return tuple(orbit)


def cutoff(profile: Profile, j: int, fraction: Q) -> Q:
    return 1/((1-fraction)/profile[j]+fraction/profile[j+1])


def rotate_point(point: Point) -> Point:
    x, y, t = point
    return -y, x, 1-t


def rotate_vector(vector: Point) -> Point:
    x, y, t = vector
    return -y, x, -t


def cylinder(sector: Sector, profile: Profile, tau: Q, upper: bool) -> Cylinder:
    j, l, r, a, b = sector.index, sector.l, sector.r, sector.a, sector.b
    if upper:
        vertices = ((Q(0), Q(0), Q(1)),
                    (l*profile[j], Q(0), 1-profile[j]),
                    (r*profile[j+1], Q(0), 1-profile[j+1]))
        axis = (-tau*a, Q(1), tau*b)
    else:
        vertices = ((Q(0), Q(0), Q(0)),
                    (Q(0), l*profile[j], profile[j]),
                    (Q(0), r*profile[j+1], profile[j+1]))
        axis = (Q(1), tau*a, tau*b)
    return Cylinder(frozenset(vertices), axis)


def geometry_transformations(sectors: tuple[Sector, ...], orbit: tuple[Profiles, ...]) -> int:
    count = 0
    for tau in (Q(0), Q(1, 2), Q(1)):
        for profiles in orbit:
            rotated = profiles.rotate()
            for sector in sectors:
                for upper in (False, True):
                    original = cylinder(sector, profiles.upper if upper else profiles.lower,
                                        tau, upper)
                    target_sector = sector if upper else sectors[2047-sector.index]
                    target_profile = rotated.lower if upper else rotated.upper
                    target = cylinder(target_sector, target_profile, tau, not upper)
                    require(frozenset(map(rotate_point, original.vertices)) == target.vertices,
                            "Node triangles do not transform to the claimed profiles")
                    expected_axis = tuple(-entry for entry in target.axis) if upper else target.axis
                    require(rotate_vector(original.axis) == expected_axis,
                            "Axis reflection/orientation mismatch")
                    count += 1
    return count


def interpolation_and_edges(sectors: tuple[Sector, ...], orbit: tuple[Profiles, ...]) -> int:
    count = 0
    initial = orbit[0]
    for sector in sectors[1024:]:
        j, reflected = sector.index, 2047-sector.index
        for fraction in (Q(1, 3), Q(2, 3)):
            for profile in (initial.lower, initial.upper):
                require(cutoff(profile[::-1], j, fraction) ==
                        cutoff(profile, reflected, 1-fraction),
                        "Reflection does not commute with inverse-radius interpolation")
            z = (cutoff(initial.lower, j, fraction),
                 cutoff(initial.upper, reflected, 1-fraction),
                 cutoff(initial.lower, reflected, 1-fraction),
                 cutoff(initial.upper, j, fraction))
            actual_edges = tuple((cutoff(p.lower, j, fraction),
                                  cutoff(p.upper, reflected, 1-fraction)) for p in orbit)
            cycle_edges = tuple((z[k], z[(k+1) % 4]) for k in range(4))
            require(sorted(actual_edges) == sorted(cycle_edges),
                    "The four profile orientations do not yield the directed cycle")
            count += 1
    return count


def make_segment(sector: Sector, tau: Q, fraction: Q = Q(1, 2)) -> Segment:
    q = sector.l+fraction*(sector.r-sector.l)
    a, b = sector.a, sector.b
    u, v = tau*a, tau*b
    w, d = u-v*q, 1+u*u
    D = 1+u*w
    alpha, beta = 1+v*q*(1-u)/d, 1-v*q*(1+u)/d
    result = Segment(sector, tau, q, alpha, beta, u/(1+w), 1/(1+w), d/(2*D))
    require(0 < beta <= alpha and result.rho0 < result.T < result.sigma0,
            "Segment hypotheses fail")
    require(result.sigma0 >= Q(3, 2)*result.T, "Endpoint compensation bound fails")
    return result


def member(point: Point, sector: Sector, profile: Profile, tau: Q, upper: bool) -> bool:
    x, y, t = point
    if upper:
        radial, angular = 1-t+tau*sector.b*y, x+tau*sector.a*y
    else:
        radial, angular = t-tau*sector.b*x, y-tau*sector.a*x
    left_weight = (sector.r*radial-angular)/((sector.r-sector.l)*profile[sector.index])
    right_weight = (angular-sector.l*radial)/((sector.r-sector.l)*profile[sector.index+1])
    return left_weight >= 0 and right_weight >= 0 and left_weight+right_weight <= 1


def segment_controls(sectors: tuple[Sector, ...], orbit: tuple[Profiles, ...],
                     segments: tuple[Segment, ...]) -> int:
    nodes = tuple(Q(j, 1024)-1 for j in range(2049))
    values = tuple(phi(q) for q in nodes)
    count = 0
    for segment in segments:
        j = segment.sector.index
        reflected = 2047-j
        for fraction in (Q(1, 4), Q(3, 4)):
            rho = segment.rho0+fraction*(segment.T-segment.rho0)
            point, sigma = segment.point(rho)
            x, y, t = point
            require(x < 0 < y and -x >= y and y <= t and -x <= 1-t,
                    "Segment point violates tetrahedron/sign hypotheses")
            lower = tuple(q*t+segment.tau*f*x for q, f in zip(nodes, values))
            upper = tuple(q*(1-t)-segment.tau*f*y for q, f in zip(nodes, values))
            require(tuple(k for k in range(2048) if lower[k] <= y <= lower[k+1]) == (j,),
                    "Unexpected lower sector")
            require(tuple(k for k in range(2048) if upper[k] <= x <= upper[k+1]) == (reflected,),
                    "Unexpected upper sector")
            for profiles in orbit:
                for is_upper, sector, profile, radius in (
                    (False, sectors[j], profiles.lower, rho),
                    (True, sectors[reflected], profiles.upper, sigma),
                ):
                    require(member(point, sector, profile, segment.tau, is_upper) ==
                            (radius <= cutoff(profile, sector.index, Q(1, 2))),
                            "Non-even triangle membership disagrees with signed radial cutoff")
            count += 1
    return count


def pairing_bound(values: tuple[Q, ...], segment: Segment) -> Q:
    require(len(values) == 4 and all(value > 0 for value in values), "Four positive radii required")
    require(all(segment.holds(values[i], values[(i+1) % 4]) for i in range(4)),
            "Incomplete directed cycle constraints")
    pairs = tuple((i, (i+1) % 4) for i, value in enumerate(values) if value < segment.T)
    used = tuple(index for pair in pairs for index in pair)
    require(len(set(used)) == len(used), "Compensation pairs overlap")
    for i, j in pairs:
        require(values[j] > segment.T and values[i]**2+values[j]**2 >= 2*segment.T**2,
                "Directed pair compensation failed")
    remaining = tuple(i for i in range(4) if i not in used)
    require(all(values[i] >= segment.T for i in remaining), "Unpaired small radius")
    total = sum((value*value for value in values), Q(0))
    require(total >= 4*segment.T**2, "Four-radius quadratic bound failed")
    return total


def matching_controls(segments: tuple[Segment, ...]) -> tuple[int, int]:
    attempted = accepted = 0
    for segment in segments:
        grid = sorted({value for value in
                       (segment.T/5, segment.rho0/2, segment.rho0,
                        (segment.rho0+segment.T)/2, segment.T,
                        Q(6, 5)*segment.T, Q(3, 2)*segment.T, 2*segment.T)
                       if value > 0})
        for values in product(grid, repeat=4):
            attempted += 1
            if all(segment.holds(values[i], values[(i+1) % 4]) for i in range(4)):
                pairing_bound(values, segment)
                accepted += 1
    require(accepted > 0, "No feasible matching controls were tested")
    return attempted, accepted


def negative_controls(segment: Segment, profiles: Profiles) -> dict[str, int | str]:
    require(segment.tau == 0, "Failure controls use the zero-tilt segment")
    T = segment.T
    omitted_edges = 0
    for missing in range(4):
        z = tuple(T/10 if i == missing else T for i in range(4))
        retained = tuple(i for i in range(4) if i != missing)
        require(all(segment.holds(z[i], z[(i+1) % 4]) for i in retained),
                "Incomplete-orbit control does not satisfy retained edges")
        require(sum(value*value for value in z) < 4*T*T, "Incomplete-orbit control too costly")
        interval = segment.uncovered_interval(z[missing], z[(missing+1) % 4])
        require(interval is not None, "The omitted edge should fail")
        rho = (interval[0]+interval[1])/2
        sigma = (1-segment.alpha*rho)/segment.beta
        require(rho > z[missing] and sigma > z[(missing+1) % 4], "Missing-edge witness failed")
        try:
            pairing_bound(z, segment)
        except ReviewFailure:
            omitted_edges += 1
        else:
            raise ReviewFailure("Incomplete orbit was accepted")
    z = (T/10, T/10, T/10, Q(19, 10)*T)
    star = ((0, 3), (1, 3), (2, 3), (3, 0))
    require(all(segment.holds(z[i], z[j]) for i, j in star), "Noninjective edges should hold")
    require(all(z[i]**2+z[3]**2 >= 2*T*T for i in range(3)), "Individual compensation should pass")
    require(sum(value*value for value in z) < 4*T*T, "Repeated compensation control too costly")
    require(len({j for _, j in star}) < len(star), "Control successor map should not be injective")
    require(profiles.rotate().upper != profiles.lower,
            "A swap without reflection must fail on non-even profiles")
    return {
        "omitted_orientation_controls": omitted_edges,
        "noninjective_compensation": "Three T/10 radii all point to one 19T/10 radius; total is 91T^2/25<4T^2",
        "incorrect_even_only_swap": "Rejected: (A,B)->(B,A) differs from the actual (B,A^-)",
        "scope": "Counterexamples to incomplete necessary-condition arguments, not verified covers",
    }


def cost_reflection_controls(sectors: tuple[Sector, ...], profiles: Profiles) -> int:
    # Work with each triangle numerator and squared weight. Equal reflected
    # denominators make summing the four unequal products exactly legitimate.
    count = 0
    for sector in sectors[1024:]:
        j, k = sector.index, 2047-sector.index
        reflected = sectors[k]
        require((reflected.l, reflected.r, reflected.a, reflected.b) ==
                (-sector.r, -sector.l, sector.a, -sector.b), "Reflected weight mismatch")
        original_products = sum((p[j]*p[j+1]+p[k]*p[k+1]
                                 for p in (profiles.lower, profiles.upper)), Q(0))
        orbit_products = sum((o.lower[j]*o.lower[j+1] for o in profile_orbit(profiles)), Q(0))
        require(original_products == orbit_products, "Four-product cost pairing lost a profile")
        require(sector.r-sector.l == reflected.r-reflected.l, "Reflected interval length differs")
        count += 1
    return count


def main() -> None:
    sectors = make_sectors()
    profiles = make_profiles()
    orbit = profile_orbit(profiles)
    transformations = geometry_transformations(sectors, orbit)
    interpolations = interpolation_and_edges(sectors, orbit)
    segments = tuple(make_segment(sectors[j], tau)
                     for j in (1024, 1535, 2047) for tau in (Q(0), Q(1, 2), Q(1)))
    geometry = segment_controls(sectors, orbit, segments)
    attempted, accepted = matching_controls(segments)
    controls = negative_controls(segments[0], profiles)
    costs = cost_reflection_controls(sectors, profiles)
    result = {
        "verified": True,
        "scope": "Independent exact orbit/interpolation/matching controls; no integral or full-cover replay",
        "accepted_source_commit": "ec5e88b8c954b8ce42ee71a993271e14dce1aded",
        "even_extension_commit": "2bff558c3c6ed29dcda80730ef2c255ee58ac477",
        "non_even_cylinder_transformations": transformations,
        "interior_angle_orbit_checks": interpolations,
        "segment_points": geometry,
        "angular_sectors_checked_per_point_per_family": 2048,
        "scalar_tuples_examined": attempted,
        "feasible_scalar_tuples_checked": accepted,
        "reflected_cost_pairs": costs,
        "negative_controls": controls,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
