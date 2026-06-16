from sage.all import *

R.<t> = PolynomialRing(QQ)
S.<x> = PolynomialRing(QQ)
T.<th> = PolynomialRing(QQ)
U.<z> = PolynomialRing(QQ)

# Bernstein coefficients on [a,b].
def bernstein_coefficients(P, a, b, m=None):
    a = QQ(a); b = QQ(b)
    if m is None:
        m = P.degree()
    Q = S(0)
    for j in range(P.degree() + 1):
        Q += P[j] * (a + (b-a)*x)^j
    c = [Q.monomial_coefficient(x^j) for j in range(m+1)]
    coeffs = []
    for k in range(m+1):
        val = QQ(0)
        for j in range(k+1):
            val += c[j] * QQ(binomial(k,j)) / QQ(binomial(m,j))
        coeffs.append(val)
    return coeffs

# Local path-plus-leaf model for the fifth moment.
def q_poly_path(L):
    # vertices: a=0, u=1, then L-1 internal path vertices, plus a leaf at a
    n = L + 2
    a = 0; u = 1; leaf = n-1
    Q = matrix(T, n, n, 0)
    edges = []
    prev = u
    for k in range(L-1):
        v = 2 + k
        edges.append((prev, v))
        prev = v
    edges.append((prev, a))
    edges.append((a, leaf))
    for i,j in edges:
        Q[i,i] += 1; Q[j,j] += 1
        Q[i,j] += 1; Q[j,i] += 1
    b = matrix(T, n, 1, [1,1] + [0]*(n-2))
    M = Q + th*b*b.transpose()
    return (b.transpose() * (M^5) * b)[0,0]

qstar = 64*th^5 + 240*th^4 + 472*th^3 + 603*th^2 + 512*th + 249
for L in [2,3,4,5,6]:
    assert all((q_poly_path(L)-qstar)[i] >= 0
               for i in range((q_poly_path(L)-qstar).degree()+1))
# For L >= 6, walks of length five cannot see farther, so the L=6 value is stable.

# Expansion and one-sided integration algebra.  Use a two-variable ring here.
W.<tt,z> = PolynomialRing(QQ)
A2 = 2140 - 7*tt^4
B2 = 4414 - 4*tt^4
Pt1z = (64*(1-z)^5 + 240*(1-z)^4 + 472*(1-z)^3
        + 603*(1-z)^2 + 512*(1-z) + 249) - tt^4*(3 + 4*(1-z))
assert Pt1z == A2 - B2*z + 4099*z^2 - 2072*z^3 + 560*z^4 - 64*z^5

hmax = QQ(1573)/QQ(4090)
assert A2.subs({tt: QQ(3)}) == 1573
assert B2.subs({tt: QQ(3)}) == 4090
assert QQ(4099) - QQ(2072)*hmax - QQ(64)*hmax^3 > 3000

# Scalar envelope polynomials.
A = 2140 - 7*t^4
B = 4414 - 4*t^4
K = 5*(t^3 + 5*t^2 + 25*t + 125)
Lhat = QQ(1)/QQ(2)*A^2*B^2 + 1000*A^3

def Mconst(r):
    return QQ(binomial(2*r+2, r+1)) / QQ(r+1)

def first_branch_constant(r):
    return QQ((r-1)^(r-1)) / QQ(r^r)

Lav = QQ(1)/QQ(2) * (QQ(7)/QQ(2) - t)^2

checks = []
for name, r, a0, b0 in [
    ("LG0", 3, QQ(2), QQ(12)/QQ(5)),
    ("LG1", 5, QQ(12)/QQ(5), QQ(14)/QQ(5)),
    ("LG2", 8, QQ(14)/QQ(5), QQ(3)),
]:
    P = t^(r-1)*Lav - Mconst(r)*first_branch_constant(r)
    checks.append((name, P, a0, b0))

for name, r, a0, b0 in [
    ("LG3", 10, QQ(3), QQ(7)/QQ(2)),
    ("LG4", 30, QQ(7)/QQ(2), QQ(58)/QQ(15)),
]:
    P = t^(r-1)*Lhat - Mconst(r)*first_branch_constant(r)*B^3*K
    checks.append((name, P, a0, b0))

P = Lhat - Mconst(30) * (4-t) / QQ(4)^30 * B^3 * K
checks.append(("LG5", P, QQ(58)/QQ(15), QQ(4)))

for name, P, a0, b0 in checks:
    coeffs = bernstein_coefficients(P, a0, b0)
    assert min(coeffs) > 0, name

print("All line-graph scalar and local-moment checks passed.")
