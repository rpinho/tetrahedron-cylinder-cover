"""Exact arithmetic, geometric controls and failure cases for the extension.

Alexander Richard Towell, with substantial OpenAI Codex assistance.
Finite geometric controls supplement the written continuum argument.
"""
import copy
from contextlib import redirect_stdout
from fractions import Fraction as F
import io
import json
from pathlib import Path
import unittest

import independent_radii as check
import structural_independent as previous


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cylinder(sector, lower, radii, tau):
    l, r, a, b = sector.l, sector.r, sector.a, sector.b
    left, right = radii
    if lower:
        return ([(0, 0, 0), (0, l*left, left), (0, r*right, right)],
                (1, tau*a, tau*b))
    return ([(0, 0, 1), (l*left, 0, 1-left), (r*right, 0, 1-right)],
            (-tau*a, 1, tau*b))


def contains(vertices, axis, point):
    # Direct edge-cross-axis halfspaces; independent of the proof's intercept formulas.
    for i in range(3):
        a, b, c = vertices[i], vertices[(i+1) % 3], vertices[(i+2) % 3]
        normal = cross(sub(b, a), axis)
        reference = dot(normal, sub(c, a))
        if not reference:
            raise ValueError('Degenerate prism fixture')
        if reference * dot(normal, sub(point, a)) < 0:
            return False
    return True


class GeometricControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.positive = previous.make_sectors(8, previous.FINAL_FAMILY.shape)
        cls.full = tuple(previous.Sector(-s.r, -s.l, s.a, -s.b)
                         for s in reversed(cls.positive)) + cls.positive

    def test_segment_stays_inside_and_has_only_two_candidate_cylinders(self):
        # Parameterize by physical t via the two angular planes, not formula (2).
        for s in self.positive:
            for q in [(2*s.l+s.r)/3, (s.l+2*s.r)/3]:
                for tau in [F(0), F(1, 2), F(1)]:
                    u, v, w = tau*s.a, tau*s.b, tau*(s.a-s.b*q)
                    d = 1+u*u
                    T = d/(2*(1+u*w))
                    beta = 1-v*q*(1+u)/d
                    t0, t1 = w/(1+w), beta*T
                    for blend in [F(0), F(1, 4), F(1, 2), F(3, 4), F(1)]:
                        t = t0+(t1-t0)*blend
                        x = -q*(1-t+w*t)/(1+w*w)
                        y = q*(t-w*(1-t))/(1+w*w)
                        self.assertLessEqual(abs(x), 1-t)
                        self.assertLessEqual(abs(y), t)
                        rho, sigma = t-v*x, 1-t-v*y
                        self.assertEqual(y-u*x, q*rho)
                        self.assertEqual(x+u*y, -q*sigma)
                        self.assertLessEqual(rho, T)
                        self.assertGreaterEqual(sigma, T)
                        if blend:
                            lower = [k for k, z in enumerate(self.full)
                                     if contains(*cylinder(z, True, (F(100), F(100)), tau),
                                                 (x, y, t))]
                            upper = [k for k, z in enumerate(self.full)
                                     if contains(*cylinder(z, False, (F(100), F(100)), tau),
                                                 (x, y, t))]
                            self.assertEqual([self.full[k] for k in lower], [s])
                            self.assertEqual([(self.full[k].l, self.full[k].r) for k in upper],
                                             [(-s.r, -s.l)])

    def test_quarter_turn_exchanges_independent_even_profiles(self):
        def turn(p):
            return (-p[1], p[0], 1-p[2])
        def turn_axis(p):
            return (-p[1], p[0], -p[2])
        for s in self.full:
            reflected = previous.Sector(-s.r, -s.l, s.a, -s.b)
            A = lambda q: 1+q*q/7
            B = lambda q: 2+q**4/11
            for tau in [F(0), F(1, 2), F(1)]:
                vertices, axis = cylinder(s, True, (A(s.l), A(s.r)), tau)
                expected, expected_axis = cylinder(reflected, False,
                                                   (A(-s.r), A(-s.l)), tau)
                self.assertEqual(set(map(turn, vertices)), set(expected))
                self.assertEqual(turn_axis(axis), expected_axis)
                vertices, axis = cylinder(s, False, (B(s.l), B(s.r)), tau)
                expected, expected_axis = cylinder(s, True, (B(s.l), B(s.r)), tau)
                self.assertEqual(set(map(turn, vertices)), set(expected))
                self.assertEqual(turn_axis(axis), tuple(-x for x in expected_axis))

    def test_single_diagonal_point_is_not_a_cover_certificate(self):
        # tau=0: A=1/4, B=1/2 passes max(A,B)>=T=1/2 everywhere,
        # but the t=3/8 slice is uncovered. This is NOT a counterexample cover.
        self.assertGreaterEqual(max(F(1, 4), F(1, 2)), F(1, 2))
        point = (F(0), F(0), F(3, 8))
        for s in self.full:
            self.assertFalse(contains(*cylinder(s, True, (F(1, 4),)*2, F(0)), point))
            self.assertFalse(contains(*cylinder(s, False, (F(1, 2),)*2, F(0)), point))

    def test_endpoint_estimate_is_needed_for_the_scalar_argument(self):
        # A non-source sector at u=9/10,w=0 has T=181/200 and sigma0=1.
        # A=1/10,B=1 covers both test segments but violates the squared budget.
        T, A, B, rho0, sigma0 = F(181, 200), F(1, 10), F(1), F(9, 10), F(1)
        self.assertLess(A, rho0)
        self.assertGreaterEqual(B, sigma0)
        self.assertGreaterEqual(B, T)  # lower radius after exchange covers that segment
        self.assertLess(A*A+B*B, 2*T*T)
        sector = previous.Sector(F(3, 4), F(1), F(9, 10), F(9, 10))
        with self.assertRaises(previous.VerificationError):
            check.check_endpoints((sector,))

    def test_unequal_radial_profiles_have_the_full_area_factor(self):
        # Generic physical projection formula in metric diag(1,1,2), A_min^2=2.
        for s in self.full:
            for lower, radii in [(True, (F(4, 7), F(5, 8))),
                                 (False, (F(2, 3), F(7, 9)))]:
                vertices, axis = cylinder(s, lower, radii, F(3, 5))
                normal = cross(sub(vertices[1], vertices[0]), sub(vertices[2], vertices[0]))
                axis_norm_squared = axis[0]**2+axis[1]**2+2*axis[2]**2
                normalized_area_squared = dot(normal, axis)**2/(4*axis_norm_squared)
                expected = ((s.r-s.l)*radii[0]*radii[1])**2/(4*axis_norm_squared)
                self.assertEqual(normalized_area_squared, expected)
        # At zero tilt and constant unequal radii, all n sectors sum to A^2+B^2.
        A, B = F(1, 3), F(2, 3)
        area = sum((s.r-s.l)*(A*A+B*B)/2 for s in self.full)
        self.assertEqual(area, F(5, 9))


class EndpointTests(unittest.TestCase):
    def test_known_positive_margin(self):
        # Catches a sign/factor mistake in the new endpoint polynomial.
        sector = previous.Sector(F(3, 4), F(1), F(1, 2), F(1, 2))
        self.assertEqual(check.endpoint_margin(sector), F(1, 32))

    def test_endpoint_condition_rejects_an_otherwise_admissible_sector(self):
        # Base geometry alone does not imply the new endpoint estimate.
        sector = previous.Sector(F(3, 4), F(1), F(9, 10), F(9, 10))
        with self.assertRaises(previous.VerificationError):
            check.check_endpoints((sector,))

    def test_empty_sector_list_rejected(self):
        with self.assertRaises(previous.VerificationError):
            check.check_endpoints(())

    def test_complete_fixed_family_endpoint_check(self):
        sectors = previous.make_sectors(previous.FINAL_FAMILY.n,
                                        previous.FINAL_FAMILY.shape)
        result = check.check_endpoints(sectors)
        self.assertIsInstance(result, dict)
        self.assertEqual(result['sectors_checked'], len(sectors))
        self.assertGreater(F(result['minimum_endpoint_margin']), 0)
        self.assertGreaterEqual(F(result['endpoint_ratio_lower']), F(3, 2))


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(Path(__file__).with_name(
            'structural_final_family.json').read_text())

    def test_complete_replay_retains_the_normalization(self):
        # Catches a factor-of-two loss or skipping the inherited integral replay.
        result = check.verify(self.data)
        self.assertIsInstance(result, dict)
        self.assertTrue(result['verified'])
        self.assertGreaterEqual(F(result['minimum_lower_enclosure']), F(4993, 10000))
        self.assertEqual(result['q_tilt_cells_checked'], 860160)
        self.assertEqual(result['endpoint_check']['sectors_checked'], 1024)

    def test_missing_tilt_interval_rejected(self):
        data = copy.deepcopy(self.data)
        data['accepted_intervals'].pop(0)
        with self.assertRaises(previous.VerificationError):
            check.verify(data)


class CommandTests(unittest.TestCase):
    def test_bad_endpoint_control_fails(self):
        with redirect_stdout(io.StringIO()):
            self.assertEqual(check.main(['--negative-control', 'endpoint']), 1)

    def test_incomplete_tilt_control_fails(self):
        with redirect_stdout(io.StringIO()):
            self.assertEqual(check.main(['--negative-control', 'partition']), 1)


if __name__ == '__main__':
    unittest.main()
