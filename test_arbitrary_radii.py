"""Exact controls for the four-orientation extension, not a formal proof.

Alexander Richard Towell, with substantial OpenAI Codex assistance.
Predecessor axes and normalization: Ricardo Pinho.
"""
from fractions import Fraction as F
from itertools import product
import unittest

import structural_independent as previous
from test_independent_radii import cylinder


def cutoff(sector, radii, q):
    """Intersect the triangle's outer edge using barycentric coordinates."""
    weight = (q-sector.l)/(sector.r-sector.l)
    return 1/((1-weight)/radii[0]+weight/radii[1])


def segment_covered(x, y, rho0, slope):
    """Exact normalized T=1 segment criterion, including its open endpoints."""
    return x >= 1 or y >= 1+slope*(1-max(x, rho0))


class ArbitraryProfileGeometry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        positive = previous.make_sectors(8, previous.FINAL_FAMILY.shape)
        cls.sectors = tuple(previous.Sector(-s.r, -s.l, s.a, -s.b)
                            for s in reversed(positive)) + positive

    def test_four_rotations_on_non_even_node_radii(self):
        # Catches an omitted angular reflection or a wrong axis sign.
        A = lambda q: 2+q/7+q*q/11
        B = lambda q: 3-q/5+q**3/13
        self.assertNotEqual(A(F(1)), A(F(-1)))
        self.assertNotEqual(B(F(1)), B(F(-1)))
        turn = lambda p: (-p[1], p[0], 1-p[2])
        turn_axis = lambda p: (-p[1], p[0], -p[2])
        for s, tau, lower in product(self.sectors, (F(0), F(1, 2), F(1)),
                                     (True, False)):
            profile = A if lower else B
            vertices, axis = cylinder(s, lower, (profile(s.l), profile(s.r)), tau)
            original_vertices, original_axis = vertices, axis
            current_sector, current_lower = s, lower
            current_profile = profile
            for _ in range(4):
                vertices = list(map(turn, vertices))
                axis = turn_axis(axis)
                if current_lower:
                    current_sector = previous.Sector(-current_sector.r,
                        -current_sector.l, current_sector.a, -current_sector.b)
                    old_profile = current_profile
                    current_profile = lambda q, old=old_profile: old(-q)
                current_lower = not current_lower
                expected_vertices, expected_axis = cylinder(current_sector,
                    current_lower, (current_profile(current_sector.l),
                                    current_profile(current_sector.r)), tau)
                self.assertEqual(set(vertices), set(expected_vertices))
                self.assertIn(axis, (expected_axis, tuple(-x for x in expected_axis)))
            self.assertEqual(set(vertices), set(original_vertices))
            self.assertEqual(axis, original_axis)

    def test_reflection_commutes_with_reciprocal_interpolation(self):
        # Catches reversing the sector without reversing its endpoint radii.
        s = previous.Sector(F(1, 4), F(1, 2), F(1, 3), F(1, 7))
        reflected = previous.Sector(-s.r, -s.l, s.a, -s.b)
        radii, q = (F(2), F(5)), F(1, 3)
        self.assertEqual(cutoff(s, radii, q), F(5, 2))
        self.assertEqual(cutoff(reflected, tuple(reversed(radii)), -q), F(5, 2))
        self.assertNotEqual(cutoff(reflected, radii, -q), F(5, 2))

    def test_full_area_keeps_both_angular_halves(self):
        # At tau=0 the complete node formula is rational. Evenness would
        # incorrectly replace the full sum by twice the positive half.
        A, B = lambda q: 2+q/2, lambda q: 3+q/3
        term = lambda s: (s.r-s.l)*(A(s.l)*A(s.r)+B(s.l)*B(s.r))/2
        full = sum(map(term, self.sectors))
        positive = sum(term(s) for s in self.sectors if s.l >= 0)
        paired = sum(term(s)+term(previous.Sector(-s.r, -s.l, s.a, -s.b))
                     for s in self.sectors if s.l >= 0)
        self.assertEqual(full, paired)
        self.assertNotEqual(full, 2*positive)


class CycleControls(unittest.TestCase):
    def test_complete_declared_rational_domain(self):
        # Finite evidence only: checks actual segment inequalities, including
        # the small-radius branch, before comparing the four-square budget.
        values = tuple(F(i, 4) for i in range(1, 9))
        feasible = 0
        for rho0, slope in product((F(0), F(1, 4), F(1, 2)),
                                   (F(1), F(3, 2), F(2))):
            self.assertGreaterEqual(1+slope*(1-rho0), F(3, 2))
            for radii in product(values, repeat=4):
                if all(segment_covered(radii[i], radii[(i+1) % 4], rho0, slope)
                       for i in range(4)):
                    feasible += 1
                    self.assertGreaterEqual(sum(x*x for x in radii), 4)
        self.assertGreater(feasible, 0)

    def test_incomplete_cycle_is_insufficient(self):
        # Two true oriented constraints leave two small radii uncompensated.
        radii = (F(1), F(1, 10), F(1), F(1, 10))
        self.assertTrue(segment_covered(radii[0], radii[1], F(0), F(1)))
        self.assertTrue(segment_covered(radii[2], radii[3], F(0), F(1)))
        self.assertFalse(segment_covered(radii[1], radii[2], F(0), F(1)))
        self.assertLess(sum(x*x for x in radii), 4)

    def test_noninjective_compensation_is_insufficient(self):
        # Three small values can share one large successor in a general
        # directed graph; counting that large value three times is invalid.
        small, large = F(3, 4), F(5, 4)
        self.assertTrue(segment_covered(small, large, F(0), F(1)))
        self.assertGreaterEqual(small*small+large*large, 2)
        self.assertLess(3*small*small+large*large, 4)


if __name__ == '__main__':
    unittest.main()
