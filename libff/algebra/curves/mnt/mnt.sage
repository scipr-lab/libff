#!/usr/bin/env sage -python

from sage.all import *
import sys
sys.path.append("../")
import params_generator

# Computes the order of G1, the safe subgroup of E/Fq
def g1_order(curve_order):
    decomposition = factor(curve_order)
    # Factor returns the prime decomposition and orders prime
    # factors from smaller to biggest
    biggest_factor = decomposition[-1]
    assert(biggest_factor[1] == 1)
    return biggest_factor[0]

def mnt4_298_params():
    print("Generating parameters for MNT4-298")
    # Scalar field
    prime_r = 475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137
    params_generator.generate_libff_Fp_model_params(prime_r)
    Fr = GF(prime_r)

    # Base field characteristic
    prime_q = 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081
    params_generator.generate_libff_Fp_model_params(prime_q)
    Fq = GF(prime_q)

    # E/Fq
    coeff_a = 2
    coeff_b = 423894536526684178289416011533888240029318103673896002803341544124054745019340795360841685
    curve = EllipticCurve(Fq, [coeff_a, coeff_b])

    # E'/Fq2
    non_residue = Fq(17)
    Fqx.<j> = PolynomialRing(Fq, 'j')
    assert(Fqx(j^2 + non_residue).is_irreducible())
    Fq2.<u> = GF(prime_q^2, modulus=j^2 - non_residue)

    twist_coeff_a = Fq2(coeff_a * non_residue)
    twist_coeff_b = Fq2(coeff_b * non_residue * u)
    twist = EllipticCurve(Fq2, [twist_coeff_a, twist_coeff_b])

    # Cofactors
    curve_order = curve.order()
    print('curve_order = {}'.format(curve_order))
    twist_order = twist.order()
    print('twist_order = {}'.format(twist_order))

    g1_h = curve_order // g1_order(curve_order)
    assert(g1_h == 1)
    assert(curve_order == prime_r)
    print('g1_h = {}'.format(g1_h))
    g2_h = twist_order // g1_order(curve_order)
    print('g2_h = {}'.format(g2_h))

def mnt6_298_params():
    print("Generating parameters for MNT6-298")
    # Scalar field
    prime_r = 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081
    params_generator.generate_libff_Fp_model_params(prime_r)
    Fr = GF(prime_r)

    # Base field characteristic
    prime_q = 475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137
    params_generator.generate_libff_Fp_model_params(prime_q)
    Fq = GF(prime_q)

    # E/Fq
    coeff_a = 11
    coeff_b = 106700080510851735677967319632585352256454251201367587890185989362936000262606668469523074
    curve = EllipticCurve(Fq, [coeff_a, coeff_b])

    # E'/Fq3
    non_residue = Fq(5)
    print('non_residue = {}'.format(non_residue))
    Fqx.<j> = PolynomialRing(Fq, 'j')
    assert(Fqx(j^3 + non_residue).is_irreducible())
    Fq3.<u> = GF(prime_q^3, modulus=j^3 - non_residue)

    twist_coeff_a = Fq3(coeff_a * u^2)
    twist_coeff_b = Fq3(coeff_b * non_residue)

    twist = EllipticCurve(Fq3, [twist_coeff_a, twist_coeff_b])

    # Cofactors
    curve_order = curve.order()
    print('curve_order = {}'.format(curve_order))
    twist_order = twist.order()
    print('twist_order = {}'.format(twist_order))

    g1_h = curve_order // g1_order(curve_order)
    assert(g1_h == 1)
    assert(curve_order == prime_r)
    print('g1_h = {}'.format(g1_h))
    g2_h = twist_order // g1_order(curve_order)
    print('g2_h = {}'.format(g2_h))

def main():
    mnt4_298_params()
    mnt6_298_params()

if __name__ == '__main__':
    main()
