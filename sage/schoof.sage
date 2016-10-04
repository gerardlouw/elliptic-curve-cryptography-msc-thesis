load("endomorphisms.sage")

def schoof(E):
    """
    E: an elliptic curve
    """
    q = E.base_field().order()
    residues = [trace_of_frobenius_mod_2(E)]
    moduli = [2]
    l = 3
    while prod(moduli) <= 4 * sqrt(q):
        if q % l == 0:
            l = next_prime(l)
        residues.append(trace_of_frobenius_mod(E, l))
        moduli.append(l)
        l = next_prime(l)
    t = crt(residues, moduli)
    if t > 2 * sqrt(q):
        return t - prod(moduli)
    else:
        return t

def trace_of_frobenius_mod_2(E):
    """
    E: an elliptic curve
    """
    F = E.base_field()
    _.<x> = F[]
    f = x^3 + E.a4() * x + E.a6()
    q = F.order()
    if gcd(f, power_mod(x, q, f) - x).is_constant():
        return 1
    else:
        return 0

def trace_of_frobenius_mod(E, l):
    """
    E: an elliptic curve
    l: a prime number
    """
    psi_l = E.division_polynomial(l)
    q_l = E.base_field().order() % l
    t_l = 0
    while true:
        try:
            phi = FrobeniusEndomorphismMod(E, psi_l)
            lhs = phi^2 + q_l
            rhs = t_l * phi
            while t_l <= (l - 1) // 2:
                if lhs == rhs:
                    return t_l
                elif lhs == -rhs:
                    return -t_l
                t_l += 1
                rhs += phi
        except ZeroDivisionError as e:
            psi_l = e[0]
