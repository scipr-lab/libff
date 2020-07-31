#!/usr/bin/env sage -python

from sage.all import *
import sys
sys.path.append("../")
import params_generator

# Prime order of the subgroup we work in
def r(x):
    return 36*(x**4) + 36*(x**3) + 18*(x**2) + 6*x + 1

# Prime used to generate the base finite field
def q(x):
    return 36*(x**4) + 36*(x**3) + 24*(x**2) + 6*x + 1

# Compute G2 cofactor
# See: Proposition 1, Section 3.3: https://eprint.iacr.org/2015/247.pdf
def g2_h(x):
    return 36*x^4+ 36*x^3+ 30*x^2+ 6*x + 1

# Computes the order of G1, the safe subgroup of E/Fq
def g1_order(curve_order):
    decomposition = factor(curve_order)
    # Factor returns the prime decomposition and orders prime
    # factors from smaller to biggest
    biggest_factor = decomposition[-1]
    assert(biggest_factor[1] == 1)
    return biggest_factor[0]

def main():
    print("Generating parameters for alt_bn128")
    # Curve parameter
    param = 0x44e992b44a6909f1

    prime_r = r(param)
    assert(prime_r == 21888242871839275222246405745257275088548364400416034343698204186575808495617)

    prime_q = q(param)
    assert(prime_q == 21888242871839275222246405745257275088696311157297823662689037894645226208583)
    if (mod(prime_q, 6) != 1):
        raise BaseException("Unexpected: q should be = 1 (mod 6). See: https://eprint.iacr.org/2007/390.pdf")

    # Scalar field
    print('prime_r = {}'.format(prime_r))
    #params_generator.generate_libff_Fp_model_params(prime_r)
    Fr = GF(prime_r)

    # Base field
    print('prime_q = {}'.format(prime_q))
    #params_generator.generate_libff_Fp_model_params(prime_q)
    Fq = GF(prime_q)

    # E/Fq
    curve = EllipticCurve(Fq, [0, 3])
    curve_order = curve.order()

    # Cofactors
    h1 = curve_order // g1_order(curve_order)
    # G1 cofactor should be 1
    assert(h1 == 1)
    print('h1 = {}'.format(h1))
    h2 = g2_h(param)
    print('h2 = {}'.format(h2))

if __name__ == '__main__':
    main()
