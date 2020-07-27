#!/usr/bin/env sage -python

from sage.all import *

print("=== MNT4-298 ===")
# Base field characteristic
prime_q = 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081
Fq = GF(prime_q)

# G1
coeff_a = 2
coeff_b = 423894536526684178289416011533888240029318103673896002803341544124054745019340795360841685
ec_g1 = EllipticCurve(Fq, [coeff_a, coeff_b])

# G2
non_residue = Fq(17)
print('non_residue = {}'.format(non_residue))
Fqx.<j> = PolynomialRing(Fq, 'j')
assert(Fqx(j^2 + non_residue).is_irreducible())
Fq2.<u> = GF(prime_q^2, modulus=j^2 + non_residue)

twist_coeff_a = Fq2(coeff_a * non_residue)
twist_coeff_b = Fq2(coeff_b * non_residue * u)

ec_g2 = EllipticCurve(Fq2, [twist_coeff_a, twist_coeff_b])

# Cofactors
order_g1 = ec_g1.order()
order_g2 = ec_g2.order()

g2_h = order_g2 // order_g1
print('g2_h = {}'.format(g2_h))

print("=== MNT6-298 ===")
# Base field characteristic
prime_q = 475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137
Fq = GF(prime_q)

# G1
coeff_a = 11
coeff_b = 106700080510851735677967319632585352256454251201367587890185989362936000262606668469523074
ec_g1 = EllipticCurve(Fq, [coeff_a, coeff_b])

# G2
non_residue = Fq(5)
print('non_residue = {}'.format(non_residue))
Fqx.<j> = PolynomialRing(Fq, 'j')
assert(Fqx(j^3 + non_residue).is_irreducible())
Fq3.<u> = GF(prime_q^3, modulus=j^3 + non_residue)

twist_coeff_a = Fq3(coeff_a * u^2)
twist_coeff_b = Fq3(coeff_b * non_residue)

ec_g2 = EllipticCurve(Fq3, [twist_coeff_a, twist_coeff_b])

# Cofactors
order_g1 = ec_g1.order()
order_g2 = ec_g2.order()

g2_h = order_g2 // order_g1
print('g2_h = {}'.format(g2_h))
