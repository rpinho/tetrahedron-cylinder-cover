# The 0.4993 barrier with arbitrary positive radial profiles

Alexander Richard Towell, with substantial OpenAI Codex assistance.
Predecessor: Ricardo Pinho's construction and the accepted shared-profile
audit in REPORT.md at ec5e88b8c954b8ce42ee71a993271e14dce1aded.
The independent-even-profile extension at
2bff558c3c6ed29dcda80730ef2c255ee58ac477 supplies the segment argument
and endpoint estimate below. The four-orientation argument removes its
evenness assumption.
This is a written proof with exact finite arithmetic, not a formal proof.

## Statement

Keep the source n=2048 equal sectors, polynomial phi, secant axes and
0<=tau<=1. Let A_j,B_j be two arbitrary positive node-radius profiles,
used respectively in the lower and upper intercept triangles.
They need not agree and neither is required to be even. If the cylinders cover
K={|x|<=1-t, |y|<=t}, their total perpendicular-base area divided by sqrt(2)
is at least 4993/10000.

All coordinate, axis and metric conventions are those of the accepted
REPORT.md. No claim is made for arbitrary axes, other sector meshes, or tilts
outside [0,1]. This does not construct an improved cover.

## Radial cutoffs and one extra finite hypothesis

On a positive-half sector [l,r], phi_j(q)=a-bq. Write A(q),B(q) for the
reciprocals of the affine interpolants of the inverse node radii. These,
not the linear interpolants of the radii, are the exact radial cutoffs.
The predecessor establishes the sector hypotheses

    0<a<1, 0<=b<1, 0<=a-b<=1-r, phi(l)>=0, phi(r)>=0.

We additionally check, for every positive-half sector with f=a-bl,

    E = 1 - 3a^2 + f(4a-3-3a^2) >= 0.                 (1)

For q in [l,r], put u=tau*a, w=tau*(a-bq), d=1+u^2, D=1+u*w.
Then 0<=w<=u<1. Define

    T=d/(2D),       sigma_0=1/(1+w),    rho_0=u/(1+w).

Condition (1) implies sigma_0 >= (3/2)T throughout the continuum of q,tau.
Indeed this inequality is equivalent to

    P=4D-3d(1+w)=1-3u^2+w(4u-3-3u^2) >= 0.

The coefficient of w is negative for all real u. Since b>=0, P is
nondecreasing in q within the sector. At q=l, with f=a-bl,

    P(tau)=1-3f*tau+(-3a^2+4af)*tau^2-3a^2*f*tau^3.

Its derivative is at most -3f+2af*tau, because f<=a and the cubic
derivative term is nonpositive. This is nonpositive for 0<=tau<=1 and
a<1. Therefore P(q,tau)>=P(l,1)=E. The checker verifies (1) by exact
rational arithmetic for all 1024 sectors. No sampling of q or tau is used
for this implication.

## A segment of test points

Fix an interior q>0 of a positive-half sector. Set v=tau*b and

    alpha=1+v*q*(1-u)/d,     beta=1-v*q*(1+u)/d.

We have alpha>=beta>0: v*q<=u gives beta>=(1-u)/d. Also
alpha+beta=2D/d=1/T. Consider the segment

    alpha*rho+beta*sigma=1,       rho_0<=rho<=T.

Its first endpoint is (rho_0,sigma_0), since
alpha*u+beta=1+w, and its second is (T,T).
Positivity of alpha,beta and u<1 on this same weighted line give
rho_0<T<sigma_0, so the segment has nonempty interior.
Throughout the segment, sigma>=rho>=0. In its open interior
rho>u*sigma and sigma>=T. Associate the spatial point

    x=-q*(sigma+u*rho)/d,
    y= q*(rho-u*sigma)/d,
    t=rho+v*x.                                             (2)

The four exact intercept identities are

    y-u*x=q*rho,       t-v*x=rho,
    x+u*y=-q*sigma,    1-t-v*y=sigma.                      (3)

Thus the lower angular coordinate is q with radius rho; the upper
coordinate is -q in the reflected sector (coefficients a,-b), with radius
sigma. For interior points x<0<y and |x|>=y.

Every point of the segment lies in K. Put
H=1-q+u*(u-v*q), which is nonnegative. Direct calculation gives

    d*(t-y)     = H*rho + q*(u-v)*sigma,
    d*(1-t+x)   = H*sigma - q*(u-v)*rho.

Here u-v=tau*(a-b)>=0, and
q*(u-v)<=a-b<=1-r<=1-q. Thus the first expression is nonnegative;
the second is at least (H-q*(u-v))*rho>=0, using sigma>=rho.
For x<=0<=y these are the two required tetrahedron inequalities.

No other angular sectors can cover an interior point of this segment.
For lower sectors use H_k=q_k*t+tau*phi(q_k)*x. At nonpositive nodes
H_k<=0<y; on positive sectors the differences are
(q_(k+1)-q_k)*(t-tau*b_k*x)>0. Identities (3) place y strictly
between the chosen sector's two boundaries.
For upper sectors use G_k=q_k*(1-t)-tau*phi(q_k)*y. At nonnegative
nodes G_k>=-tau*phi(q_k)*y>-y>=x. On negative sectors the differences
are (q_(k+1)-q_k)*(1-t+tau*b_k*y)>0, since |b_k|<1 and
1-t>=|x|>=y>0. Identities (3) place x strictly in the reflected
sector. The inherited coefficient checks imply 0<=phi(q_k)<1.

Consequently coverage forces, for every interior segment point,

    rho <= A(q)  or  sigma <= B(-q).                       (4)

Limits at the segment endpoints may be taken by continuity. The argument
does not require uniqueness of the angular sector at an endpoint.

## Four orientations and the pointwise area inequality

The physical isometry J(x,y,t)=(-y,x,1-t) preserves K. It sends a lower
sector to the reflected upper sector and an upper sector to the lower
sector with the same angular interval (axis orientation is immaterial).
Write A^-(q)=A(-q) and B^-(q)=B(-q). Without any symmetry of the radii,
the four covers have profiles

    (A,B),       (B,A^-),       (A^-,B^-),       (B^-,A).  (5)

For example, J sends (0,q A(q),A(q)) to
(-q A(q),0,1-A(q)), an upper intercept at angle -q with the same radius.
It sends (q B(q),0,1-B(q)) to (0,q B(q),B(q)). These formulas also
hold at every node, and reflection commutes with reciprocal-affine
interpolation. The corresponding axes transform to the reflected upper
axis and the negative lower axis. Thus (5) describes actual isometric
images of the cover with the same tau and fixed secant geometry. We neither assert
that the original cover is invariant nor average covers.

For fixed q, write F(x,y) for the implication
rho<=x or sigma<=y at every interior point of the segment. Applying (4)
to the four covers in (5) gives

    F(A(q),B(-q)),  F(B(q),A(q)),
    F(A(-q),B(q)),  F(B(-q),A(-q)).                       (6)

These are exactly the consecutive directed edges of the four-cycle

    z_0=A(q) -> z_1=B(-q) -> z_2=A(-q) -> z_3=B(q) -> z_0.

Whenever F(x,y) holds and x<T, the pair satisfies y>T and
x^2+y^2>=2T^2. To prove this, first suppose rho_0<=x<T and let rho
decrease to x from above in the segment. Then

    y >= (1-alpha*x)/beta
      = T+(alpha/beta)*(T-x) >= 2T-x > T.

Hence x^2+y^2>=x^2+(2T-x)^2=2T^2+2(T-x)^2.
If x<rho_0, let rho decrease to rho_0 from above instead; then
y>=sigma_0>=(3/2)T>T and x^2+y^2>=(9/4)T^2>=2T^2.

Pair every cycle vertex z_i<T with its successor z_(i+1), with indices
taken modulo four. These pairs are disjoint: the successor map is
injective, and every selected successor is >T and therefore cannot
itself be a selected small vertex. If there are k small vertices, the
pairs contribute at least 2k T^2; each of the remaining 4-2k vertices
is >=T. Consequently,

    A(q)^2+B(-q)^2+A(-q)^2+B(q)^2 >= 4T(q,tau^2)^2.       (7)

The small-radius case is essential: no unconditional sum-of-radii bound
is used there. Nor do we claim A(q)^2+B(q)^2>=2T^2 for non-even
profiles. A single diagonal point, two incomplete cycle edges, or
counting one large radius as compensation for multiple small radii
would not justify (7).

## Area and exact inherited arithmetic

Each lower perpendicular triangle has normalized area
(r-l)*A_l*A_r/(2*sqrt(1+tau^2*(a^2+2b^2))); the upper one uses B.
Integrating reciprocal affine squares therefore gives the complete cost

    C=(1/2) sum_sectors integral_l^r
      (A(q)^2+B(q)^2)/sqrt(1+tau^2*(a^2+2b^2)) dq.

The secant coefficients at reflected sectors are (a,-b), so the weight
in this integral is even whether or not the radii are even. Pairing the
integrals at q and -q and using (7), this is at least

    2 sum_positive_sectors integral_l^r
      T(q,tau^2)^2/sqrt(1+tau^2*(a^2+2b^2)) dq.

This is exactly the integral independently enclosed in the accepted
shared-profile audit. Its complete 210-interval squared-tilt partition
and 860160 rational integration cells prove a lower bound above
4993/10000. The new verifier replays that arithmetic rather than trusting
the previously saved lower values. One-sided sector limits handle all
angular endpoints; tau=0 and tau=1 are included directly throughout.

The new finite obligation is (1). The line-segment and case arguments
above are written proofs over continuous parameters. Finite identity and
failure controls supplement them and are not substitutes for them.

Credit Ricardo Pinho for the axes, explicit cover, normalization and source
certificate. The accepted shared-profile argument and independent arithmetic
are predecessor work. The earlier even-profile extension supplied the
segment and endpoint estimate. The present four-orientation argument
removes both shared-profile and evenness restrictions on node radii;
it does not change the geometric family or improve its cover.
