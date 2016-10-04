load("pairings.sage")

def baby_step_giant_step(G, P, n):
    """
    G: a point
    P: a multiple of G
    n: the order of G
    """
    m = round(sqrt(n))
    log_G = {}
    jG = 0 * G
    for j in (0..m-1):
        log_G[jG] = j
        jG += G
    mG = m * G
    miG = 0 * G
    for i in (0..floor(n/m)):
        if P - miG in log_G:
            return m * i + log_G[P - miG]
        miG += mG

def pollard_rho(G, P, l, m=20):
    """
    G: a point of prime order
    P: a multiple of G
    l: the order of G
    """
    A, B = random_matrix(GF(l), 2, m)
    def f(Q, a, b):
        h = hash(Q) % m
        return Q + ZZ(A[h]) * G + ZZ(B[h]) * P, a + A[h], b + B[h]
    a_i, b_i = a_2i, b_2i = random_vector(GF(l), 2)
    Q_i = Q_2i = ZZ(a_i) * G + ZZ(b_i) * P
    for i in (1..):
        Q_i, a_i, b_i = f(Q_i, a_i, b_i)
        Q_2i, a_2i, b_2i = f(*f(Q_2i, a_2i, b_2i))
        if Q_i == Q_2i:
            if b_i != b_2i:
                return ZZ((a_i - a_2i) / (b_2i - b_i))
            else:
                return pollard_rho(G, P, l)

def pollard_lambda(G, P, l, m=20):
    """
    G: a point of prime order
    P: a multiple of G
    l: the order of G
    """
    A, B = random_matrix(GF(l), 2, m)
    def f(Q, a, b):
        h = hash(Q) % m
        return Q + ZZ(A[h]) * G + ZZ(B[h]) * P, a + A[h], b + B[h]
    T = {}
    theta_inv = floor(sqrt(l) / log(l))
    while true:
        a_i, b_i = random_vector(GF(l), 2)
        Q_i = ZZ(a_i) * G + ZZ(b_i) * P
        for i in (1..):
            Q_i, a_i, b_i = f(Q_i, a_i, b_i)
            if hash(Q_i) % theta_inv == 0:
                if Q_i in T:
                    a_j, b_j = T[Q_i]
                    if b_i != b_j:
                        return ZZ((a_i - a_j) / (b_j - b_i))
                else:
                    T[Q_i] = a_i, b_i
                    break

def pohlig_hellman_reduction(G, P, n, n_factors):
    """
    G: a point
    P: a multiple of G
    n: the order of G
    n_factors: the prime factorisation of n
    """
    residues = []
    moduli = [l^e for l, e in n_factors]
    for l, e in n_factors:
        d = 0
        G_prime = (n // l) * G
        for j in (0..e-1):
            P_prime = (n // l^(j + 1)) * (P - d * G)
            d += l^j * baby_step_giant_step(G_prime, P_prime, l)
        residues.append(d)
    return crt(residues, moduli)

def mov_reduction(G, P, l, E, n):
    """
    G: a point on E of prime order
    P: a multiple of G
    l: the order of G
    E: an elliptic curve
    n: the order of E
    """
    m = n
    while m % l == 0:
        m //= l
    g = 1
    while g == 1:
        Q = E.random_point() * m
        while Q * l != 0:
            Q *= l
        g = weil_pairing(E, l, G, Q)
    p = weil_pairing(E, l, P, Q)
    return p.log(g)

def frey_ruck_reduction(G, P, l, E):
    """
    G: a point on E of prime order
    P: a multiple of G
    l: the order of G
    E: an elliptic curve
    """
    q = E.base_field().order()
    g = 1
    while g == 1:
        Q = E.random_point()
        g = tate_pairing(E, l, G, Q)^((q - 1) // l)
    p = tate_pairing(E, l, P, Q)^((q - 1) // l)
    return p.log(g)

def anomalous_curve_reduction(G, P, E):
    """
    G: a point on E of order p
    P: a multiple of G
    E: an elliptic curve over GF(p)
    """
    p = E.base_field().order()
    a = (elliptic_log(E, P) // p).residue()
    g = (elliptic_log(E, G) // p).residue()
    return ZZ(a / g)

def elliptic_log(E, P):
    """
    E: an elliptic curve over GF(p)
    P: a point on E
    """
    p = E.base_field().order()
    Ep = E.base_extend(QQ).base_extend(Qp(p, 2))
    P_tilde = Ep.lift_x(ZZ(P[0]) + p * randint(0, p - 1))
    if P_tilde[1].residue() != P[1]:
        P_tilde *= -1
    return -(p * P_tilde)[0] / (p * P_tilde)[1]
