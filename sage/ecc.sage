from Crypto.Hash import SHA512
from Crypto.Random import random

def random_elliptic_curve(K):
    """
    K: a field
    """
    while true:
        a_4, a_6 = random_vector(K, 2)
        if -16 * (4 * a_4^3 + 27 * a_6^2) != 0:
            return EllipticCurve([a_4, a_6])

def diffie_hellman_key_generation(G, l):
    """
    G: the base point
    l: the order of G
    """
    s = randint(0, l - 1)
    K = s * G
    return s, K

def diffie_hellman_key_agreement(s_A, K_B):
    """
    s_A: Alice's secret key
    K_B: Bob's public key
    """
    return s_A * K_B

def elgamal_encrypt_1(m, K_B, G, l):
    """
    m: the message (encoded as an element of GF(l))
    K_B: Bob's public key
    G: the base point
    l: the order of G
    """
    e = randint(1, l - 1)
    C_1 = e * G
    c_2 = (m + hash(e * K_B)) % l
    return C_1, c_2

def elgamal_decrypt_1((C_1, c_2), s_B, l):
    """
    (C_1, c_2): the ciphertext
    s_B: Bob's secret key
    l: the order of the base point
    """
    return (c_2 - hash(s_B * C_1)) % l

def elgamal_encrypt_2(M, K_B, G, l):
    """
    M: the message (encoded as a multiple of G)
    K_B: Bob's public key
    G: the base point
    l: the order of G
    """
    e = randint(1, l - 1)
    C_1 = e * G
    C_2 = M + e * K_B
    return C_1, C_2

def elgamal_decrypt_2((C_1, C_2), s_B):
    """
    (C_1, C_2): the ciphertext
    s_B: Bob's secret key
    """
    return C_2 - s_B * C_1

def schnorr_sign(m, s_A, G, l):
    """
    m: the message
    s_A: Alice's secret key
    G: the base point
    l: the order of G
    """
    e = randint(0, l - 1)
    s_1 = hash((m, e * G)) % l
    s_2 = (e - s_A * s_1) % l
    return m, s_1, s_2

def schnorr_verify((m, s_1, s_2), K_A, G, l):
    """
    m: the message
    (s_1, s_2): the signature
    K_A: Alice's public key
    G: the base point
    l: the order of G
    """
    return s_1 == hash((m, s_1 * K_A + s_2 * G)) % l

def hash(m):
    return ZZ(SHA512.new(str(m)).hexdigest(), 16)

def randint(a, b):
    return ZZ(random.randint(int(a), int(b)))
