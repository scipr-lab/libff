#!/usr/bin/env sage -python

from sage.all import *
import sys

# Prime order of the subgroup we work in
def r(x):
    return x**4 - x**2 + 1

# Prime used to generate the base finite field
def q(x):
    return (x-1)**2 * (x**4-x**2+1)//3 + x

# Compute G1 cofactor
def g1_h(x):
    return (x-1)**2//3

# Compute G2 cofactor
def g2_h(x):
    return (x**8-4*x**7+5*x**6-4*x**4+6*x**3-4*x**2-4*x+13)//9

# Computes the order of G1, the safe subgroup of E/Fq
def g1_order(curve_order):
    decomposition = factor(curve_order)
    # Factor returns the prime decomposition and orders prime
    # factors from smaller to biggest
    biggest_factor = decomposition[-1]
    assert(biggest_factor[1] == 1)
    return biggest_factor[0]

def main():
    print("Generating parameters for bls12_381")
    # Curve parameter
    param = -0xd201000000010000

    prime_r = r(param)
    assert(prime_r == 52435875175126190479447740508185965837690552500527637822603658699938581184513)

    prime_q = q(param)
    print prime_q
    assert(prime_q == 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787)
    if (mod(prime_q, 6) != 1):
        raise BaseException("Unexpected: q should be = 1 (mod 6).")

    # Scalar field
    print('prime_r = {}'.format(prime_r))
    #params_generator.generate_libff_Fp_model_params(prime_r)
    Fr = GF(prime_r)

    # Base field
    print('prime_q = {}'.format(prime_q))
    #params_generator.generate_libff_Fp_model_params(prime_q)
    Fq = GF(prime_q)

    # E/Fq
    curve = EllipticCurve(Fq, [0, 4])
    curve_order = curve.order()

    # Cofactors
    _h1 = curve_order // g1_order(curve_order)
    h1 = g1_h(param)
    assert(_h1 == h1)
    print('h1 = {}'.format(h1))
    h2 = g2_h(param)
    print('h2 = {}'.format(h2))

if __name__ == '__main__':
    main()
