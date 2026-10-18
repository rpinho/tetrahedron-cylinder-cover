# A finite cylinder cover of the tetrahedron saving 3.86 × 10⁻⁴ below the half-area bound

This repository contains an explicit cover of the regular tetrahedron by **4,096 cylinders** with compact triangular
perpendicular bases, and two independently written exact checkers proving that its normalized total base area is

    0.499613830783714148 ≤ C < 0.499613830783714149,

a saving of more than **0.000386169** below 1/2. This is problem 100 of the OpenAI math release
([openai/math](https://github.com/openai/math), *Finite angular cylinder covers below the half-area bound*,
27 Sep 2026; [lean/docs/100.md](https://github.com/openai/math/blob/main/lean/docs/100.md)), whose explicit small-tilt
example saves about 5.4 × 10⁻¹⁰ with 16,000,000 cylinders. The construction here uses the same angular-sector family
as that preprint, with a finite tilt τ = 9/20 and numerically optimized rational radii.

No claim of optimality or priority is made: we found no larger rigorously published saving as of 15 Oct 2026, but the
search was not exhaustive. A target of 10⁻³ was **not** reached, and `structural_bound.py` proves C ≥ 0.4993 for the
final symmetric axis family with arbitrary shared even radial cutoffs, for every tilt 0 ≤ τ ≤ 1.

## Statement and normalization

The target is a finite cover of a convex body K in Euclidean three-space by
cylinders B_i + R u_i, with unit u_i and compact planar bases perpendicular to
u_i. The objective in this job is

    C = (sum_i area(B_i)) / A_min(K),
    A_min(K) = min_{|u|=1} area(proj_{u-perp} K).

Use coordinates (x,y,t) with physical point (x,y,sqrt(2)t), and
K = {|x| <= 1-t, |y| <= t}. Its four vertices are (±1,0,0), (0,±1,1)
in these coordinates. Physical edge length is 2. The coordinate metric is
diag(1,1,2). The minimum projection area is sqrt(2), NOT the area of an
arbitrarily selected projection and NOT the tetrahedron's surface area.

For a physical unit vector u, the projection area is
max(sqrt(2)|u_x|,|u_z|) + max(sqrt(2)|u_y|,|u_z|).
It is at least sqrt(2): setting U=max(|u_x|,|u_z|/sqrt(2)) and
V=max(|u_y|,|u_z|/sqrt(2)), (U+V)^2 >= |u|^2. Equality holds along
an x or y axis. Thus the half-area benchmark here is C=1/2.

The distinct directionwise quantity is
R = sum_i area(B_i)/area(proj_{u_i-perp} K). Always R <= C.
A cover below C=1/2 therefore also violates the directionwise half-bound.
Affine invariance applies to R; it does not transfer the same absolute
minimum-projection normalization to an arbitrary tetrahedron.

Source: angular preprint, Theorem 1.1, definitions (1.1)-(1.4), Lemma 2.1,
and lean/docs/100.md. The preprint guarantees, for 0<epsilon<=1/2000,
2 ceil(2/epsilon^2) triangular cylinders with
C=1/2-(13/6000)epsilon^2+remainder, |remainder|<=2epsilon^4.
At epsilon=1/2000 the guaranteed saving is
13/24000000000 - 1/8000000000000, approximately 5.4154e-10;
there are 16,000,000 cylinders.

The downloaded lean/ComparatorChallenges/TriangularCovering.lean encodes
these quantities and an existential asymptotic statement. Its final comparator
theorem has `sorry`; this file by itself is not a checked proof of the paper.
This is a statement comparator, not an audit of any separate formalization.
This job supplies a separate executable rational certificate.

## The construction

There are n=2048 sectors in each of two families, hence 4096 cylinders.
The exact direction parameters are

    tau = 9/20,
    q_j = -1 + 2j/2048,
    phi(q) = (1-q^2)(1-q^2/200 + 27q^4/1250)/4.

For a sector [l,r]=[q_j,q_{j+1}], set

    alpha_j = (r phi(l)-l phi(r))/(r-l),
    beta_j = (phi(l)-phi(r))/(r-l).

Let R_j be the positive rational radius at node q_j. The triangles in the
paper's rational coordinates are

    P_j = conv{(0,0,0), (0,q_j R_j,R_j),
               (0,q_{j+1} R_{j+1},R_{j+1})},
    Q_j = conv{(0,0,1), (q_j R_j,0,1-R_j),
               (q_{j+1} R_{j+1},0,1-R_{j+1})}.

The axes are (1,tau alpha_j,tau beta_j) and
(-tau alpha_j,1,tau beta_j). The cylinders are the full real spans of these
axes added to the triangles, with perpendicular bases obtained by projection.
Every rational R_j, vertex and direction is explicitly supplied in cover.json;
R_j can be read directly from the lower triangle's t coordinate. Reflection
R_j=R_{n-j} holds exactly. No parameter limit is part of this construction.

Writing z_j=1/R_j, its normalized area is exactly

    C = sum_{j=0}^{n-1} (2/n) /
        [z_j z_{j+1} sqrt(1+tau^2(alpha_j^2+2 beta_j^2))].

The saved numerical optimizer output is only the source of the proposed
rational radii. Coverage and the area bound are established by the separate
exact checker, not by optimization tolerances or by extrapolation from the
preprint's asymptotic theorem.

## Why the coverage certificate proves a cover

Write each lower cylinder's intercept triangle as the convex hull of the
origin and (0,l rho_l,rho_l), (0,r rho_r,rho_r), with -1<=l<r<=1
and positive endpoint radii. Its axis is (1,a,b). The infinite side planes
are y=l t+(a-l b)x and y=r t+(a-r b)x. The checker confirms exactly that
neighboring side planes coincide, the first plane is y=-t, and the last is
y=t. For any point of K, the finite sequence of boundary heights therefore
starts below y and ends above y. A first crossing selects an adjacent pair
bracketing y. The inequalities give l t0<=y0<=r t0 at the x=0 intercept,
and force t0>=0. No monotonicity assumption is needed.

The same argument at y=0 selects an upper angular wedge, using s=1-t.
Consequently, if a point is uncovered, it lies in a lower wedge and an upper
wedge and strictly violates both radial face inequalities. The checker builds
all three halfspace inequalities of each triangular prism directly from
cross products of its edges and axis. It does not reuse the search's
intercept-membership formulas.

For each angular pair, the auxiliary LP has variables (x,y,t,delta), the
four tetrahedron inequalities, four side-face inequalities, and

    -a_lower . (x,y,t) + delta <= -b_lower,
    -a_upper . (x,y,t) + delta <= -b_upper.

Here a_lower . p <= b_lower and a_upper . p <= b_upper are radial
membership conditions. Simultaneous strict failure gives a feasible delta>0.
The certificate specifies a support of at most four LP dual multipliers.
The checker solves the dual equations with exact fractions and confirms

    multipliers >= 0,
    sum_i multiplier_i * row_i = (0,0,0,1),
    sum_i multiplier_i * rhs_i <= 0.

Thus delta<=0 for every point in that pair. An optional exact Farkas witness
handles empty pair cells. The final certificates have no empty cells.

The affine isometry (x,y,t)->(-y,x,1-t), in the physical metric, permutes the
cover through four equivalent angular-pair cases. The checker verifies the
triangle and axis transformations exactly and checks one representative per
orbit, then checks that no orbit is missing or repeated. This reduces the
certificate count by four. All inequalities are closed; strict uncovered
points would still require delta>0, so faces, edges and vertices are included.

## Rational Euclidean data

The file cover.json uses the paper's rational (x,y,t) coordinates and metric
rather than rounded physical sqrt(2). The companion cover_euclidean.json
removes even that coordinate convention. Apply the similarity

    (X,Y,Z) = (x+y, x-y, 2t-1).

Its Euclidean metric is twice diag(1,1,2); every area scales by two.
The body has vertices (1,1,-1), (-1,-1,-1), (1,-1,1), (-1,1,1),
edge length 2 sqrt(2), and minimum projection area 2 sqrt(2).
The normalized cost is unchanged. All axes and intercept triangles are
rational in ordinary Euclidean coordinates. Each perpendicular triangle is
listed explicitly using the exact rational projection

    P_perp = P - v (P.v)/(v.v).

An axis need not be normalized in the data: its real span defines the same
cylinder. The checker verifies perpendicularity, equality of the projected
and intercept cylinders, nondegeneracy, the similarity, and the file hash. It also computes the squared
normalized area directly from the cross product of each exported perpendicular
base and checks exact agreement with the area used in the primary certificate.

## Running the checks

```sh
python verify.py cover.json certificate.json.gz --euclidean cover_euclidean.json --output verification.json   # standard library only
pip install -r requirements.txt                       # python-flint, for the independent check
cd independent_check && python independent.py --areas-only && python cells.py --workers 3 && python controls.py
```

`verify.py` checks every one of the 4,194,304 angular-pair cases with 1,048,576 exact rational LP dual certificates and a
verified fourfold symmetry, and bounds every square root rationally. `independent_check/` was written separately, without
access to `verify.py` or the dual certificate, by a referee pass whose job was to break the claim: it re-derives the problem
statement from the preprint, clips every one of the 4,194,304 pair cells exactly with FLINT rationals (no open gap remains),
recomputes every perpendicular base area exactly (C = 0.49961383078371414863…), and rejects covers with a cylinder removed
or shrunk (`controls.py`).

## Structural lower-bound audit

[REPORT.md](REPORT.md) supplies the geometric necessary-condition proof for
the final 2048-sector family, arbitrary positive shared even node radii, and
all tilts in [0,1]. An independent standard-library checker recomputes the
0.4993 lower bound using rational Newton square-root enclosures and verifies
the complete squared-tilt partition:

```sh
python structural_independent.py --output structural_independent_result.json
python -m unittest -v test_structural_independent.py
```

This is a written proof plus exact finite arithmetic, not a Lean proof or a
bound for arbitrary cylinder covers. The audit does not certify the generic
parameter settings of `structural_bound.py`, including its default tilt limit 2.
Alexander Richard Towell contributed this audit with substantial OpenAI Codex
assistance; separate AI review is not human peer review.

[INDEPENDENT_RADII.md](INDEPENDENT_RADII.md) extends the same 0.4993 barrier
to two arbitrary positive node-radius profiles, one in each family, with
no evenness requirement. A segment of test points replaces the single
diagonal point. Its constraints for four isometric images give a four-term
quadratic radius inequality. An exact endpoint check and the same complete
integral replay verify the arithmetic:

```sh
python independent_radii.py --output independent_radii_result.json
python independent_radii_review.py
python independent_orbit_review.py
python -m unittest -v test_independent_radii.py test_arbitrary_radii.py
```

The first command checks 1,024 endpoint inequalities and 860,160 rational
integration cells; its minimum enclosure is
`71956848089952005/144115188075855872 > 4993/10000`.
The separate review script checks exact identities and finite geometric
controls without importing the main checker. The continuum implication is a
written proof, not a formal proof. The separate orbit review and new tests
check the non-even transformation, complete area, and controls exposing
incomplete constraints and double-counted compensation. The axes, sector
mesh and tilt restriction remain part of the scope; no improved cover is
claimed. This extension is by Alexander Richard Towell with substantial OpenAI Codex assistance,
building on Ricardo Pinho's construction and the preceding audit. Separate
AI review is not human peer review.

Both `python -O independent_radii.py --negative-control endpoint` and
`python -O independent_radii.py --negative-control partition` must exit 1:
they reject an invalid endpoint estimate and a missing tilt interval.

## AI assistance

The search, the certificate and both checkers were developed with AI coding assistants (OpenAI Codex and Claude) under the
author's direction. Every numerical claim is checked by the code in this repository.

## License

MIT; see `LICENSE`.
