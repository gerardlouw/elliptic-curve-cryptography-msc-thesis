load("discrete_logs.sage")
load("ecc.sage")
load("schoof.sage")

BITS = 32
MAX_COFACTOR = 8

if __name__ == "__main__":
    # Find an elliptic curve with a near-prime number of points
    while True:
        p = random_prime(2^BITS - 1, lbound=2^(BITS - 1))
        F = FiniteField(p)
        E = random_elliptic_curve(F)
        n = schoof(E) # E.count_points() is much faster
        if any(i.divides(n) and (n // i).is_prime() for i in (1..MAX_COFACTOR)):
            break
    
    # Find a point of prime order on the elliptic curve
    while True:
        G = E.random_point()
        l = G.order()
        if l.is_prime():
            break
    
    s, K = diffie_hellman_key_generation(G, l)
    
    print "Private key: {}".format(s)
    print "Public key: {}".format(K)
    s_ = baby_step_giant_step(G, K, l)
    print "Discrete log of public key (private key): {}".format(s_)

    m = FiniteField(l).random_element()
    print "Message: {}".format(m)
    C_1, c_2 = elgamal_encrypt_1(m, K, G, l)
    print "Encrypted message (ciphertext): {}".format((C_1, c_2))
    m_ = elgamal_decrypt_1((C_1, c_2), s, l)
    print "Decrypted ciphertext (message): {}".format(m_)

    M = E.random_point()
    print "Message: {}".format(M)
    C_1, C_2 = elgamal_encrypt_2(M, K, G, l)
    print "Encrypted message (ciphertext): {}".format((C_1, C_2))
    M_ = elgamal_decrypt_2((C_1, C_2), s)
    print "Decrypted ciphertext (message): {}".format(M_)
