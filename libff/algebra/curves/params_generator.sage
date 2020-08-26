#!/usr/bin/env sage -python

from sage.all import *

def generate_libff_Fp_model_params(prime):
    num_bits = ceil(log(prime, 2))
    print('num_bits = {}'.format(num_bits))

    euler = (prime-1)/2
    print('euler = {}'.format(euler))

    factorization = factor(prime-1)
    t = 0
    term_2 = factorization[0]
    counter = 0
    if term_2[0] != 2:
        raise BaseException("The prime decomposition doesn't have any factor 2."
        "The 'high 2-adicity' requirement isn't respected")
    while not(is_odd(t)):
        s = term_2[1] - counter
        t = (prime-1)/(2**s)
        counter = counter + 1
    print('s = {}'.format(s))
    is_odd(t); print('t = {}'.format(t))

    t_minus_1_over_2 = (t-1)/2
    print('t_minus_1_over_2 = {}'.format(t_minus_1_over_2))

    multiplicative_generator = primitive_root(prime)
    print('multiplicative_generator = {}'.format(multiplicative_generator))

    root_of_unity = pow(multiplicative_generator, t, prime)
    print('root_of_unity = {}'.format(root_of_unity))

    nqr = least_quadratic_nonresidue(prime)
    print('nqr = {}'.format(nqr))

    nqr_to_t = pow(nqr, t, prime)
    print('nqr_to_t = {}'.format(nqr_to_t))

    word_len_64_bits = 64
    W_64_bits = 2**(word_len_64_bits)
    k_64_bits = ceil(num_bits/word_len_64_bits)
    print('k_64_bits (nb limbs) = {}'.format(k_64_bits))
    R_64_bits = mod(W_64_bits**k_64_bits, prime); R_64_bits
    Rsquared_64_bits = R_64_bits**2
    print('Rsquared_64_bits = {}'.format(Rsquared_64_bits))
    Rcubed_64_bits = R_64_bits**3
    print('Rcubed_64_bits = {}'.format(Rcubed_64_bits))
    inv_64_bits = hex(int(mod((1/-prime), W_64_bits)))
    print('inv_64_bits = {}'.format(inv_64_bits))

    word_len_32_bits = 32
    W_32_bits = 2**(32)
    k_32_bits = ceil(num_bits/word_len_32_bits)
    print('k_32_bits (nb limbs) = {}'.format(k_32_bits))
    R_32_bits = mod(W_32_bits**k_32_bits, prime); R_32_bits
    Rsquared_32_bits = R_32_bits**2
    print('Rsquared_32_bits = {}'.format(Rsquared_32_bits))
    Rcubed_32_bits = R_32_bits**3
    print('Rcubed_32_bits = {}'.format(Rcubed_32_bits))
    inv_32_bits = hex(int(mod(1/-prime, W_32_bits)))
    print('inv_32_bits = {}'.format(inv_32_bits))