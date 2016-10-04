def weil_pairing(E, m, P, Q):
    """
    E: an elliptic curve
    m: a positive integer
    (P, Q): m-torsion points on E
    """
    if P == Q or P == 0 or Q == 0:
        return 1
    try:
        return (-1)^m * miller(E, m, P, Q) / miller(E, m, Q, P)
    except ZeroDivisionError:
        return 1

def tate_pairing(E, m, P, Q):
    """
    E: an elliptic curve
    m: a positive integer
    (P, Q): m-torsion points on E
    """
    if Q == 0:
        return 1
    if P == Q:
        R = E.random_point()
        return tate_pairing(E, m, P, P + R) / tate_pairing(E, m, P, R)
    try:
        return miller(E, m, P, Q)
    except ZeroDivisionError:
        R = E.random_point()
        return tate_pairing(E, m, P, Q + R) / tate_pairing(E, m, P, R)

def miller(E, m, P, Q):
    """
    E: an elliptic curve
    m: a positive integer
    (P, Q): m-torsion points on E
    """
    fPi, iP = 1, 0
    fPj, jP = 1, P
    while m > 0:
        if m % 2 == 1:
            fPi *= fPj * line(E, iP, jP, Q) / line(E, iP + jP, 0, Q)
            iP += jP
        fPj *= fPj * line(E, kP, kP, Q) / line(E, 2 * jP, 0, Q)
        jP += jP
        m //= 2
    return fPi

def line(E, P, Q, R):
    """
    E: an elliptic curve
    (P, Q, R): points on E
    """
    if P == 0:
        if Q == 0:
            return 1
        else:
            return R[0] - Q[0]
    elif Q == 0 or P == -Q:
        return R[0] - P[0]
    elif P == Q:
        m = (3 * P[0]^2 + E.a4()) / (2 * P[1])
        c = (-P[0]^3 + E.a4() * P[0] + 2 * E.a6()) / (2 * P[1])
    else:
        m = (P[1] - Q[1]) / (P[0] - Q[0])
        c = (Q[0] * P[1] - P[0] * Q[1]) / (Q[0] - P[0])
    return R[1] - (m * R[0] + c)
