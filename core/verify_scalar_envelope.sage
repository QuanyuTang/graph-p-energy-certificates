from sage.all import *

R.<t> = PolynomialRing(QQ)
S.<x> = PolynomialRing(QQ)

def bernstein_coefficients(P, a, b, m=None):
    a = QQ(a)
    b = QQ(b)
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

C19 = QQ(binomial(39,19)) / QQ(20)

A = 1812 - 7*t^4
B = 4001 - 4*t^4
K = 5*(t^3 + 5*t^2 + 25*t + 125)
Lhat = QQ(1)/QQ(2)*A^2*B^2 + 1000*A^3

G0 = QQ(7)/QQ(2) - 2*t
G1 = 3*t*(QQ(7)/QQ(2)-t)^2 - 5
G2 = 320*t^3*(QQ(7)/QQ(2)-t)^2 - 1701

F1 = t^18*Lhat - C19*QQ(18)^18/QQ(19)^19 * B^3*K
F2 = Lhat - C19*(4-t)/QQ(4)^19 * B^3*K

checks = [
    ("G0", G0, QQ(0), QQ(3)/QQ(2), QQ(1)/QQ(2)),
    ("G1", G1, QQ(3)/QQ(2), QQ(2), QQ(17)/QQ(2)),
    ("G2", G2, QQ(2), QQ(3), QQ(459)),
    ("F1", F1, QQ(3), QQ(72)/QQ(19), QQ(10)^20),
    ("F2a", F2, QQ(72)/QQ(19), QQ(74)/QQ(19), QQ(10)^11),
    ("F2b", F2, QQ(74)/QQ(19), QQ(75)/QQ(19), QQ(10)^10),
    ("F2c", F2, QQ(75)/QQ(19), QQ(151)/QQ(38), QQ(10)^9),
    ("F2d", F2, QQ(151)/QQ(38), QQ(4), QQ(10)^8),
]

for name, P, a, b, lower in checks:
    coeffs = bernstein_coefficients(P, a, b)
    assert min(coeffs) >= lower, name

print("All scalar Bernstein coefficient checks passed.")
