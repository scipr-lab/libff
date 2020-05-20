/** @file
 *****************************************************************************

 Implementation of interfaces for (square-and-multiply) exponentiation.

 See algorithms.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef EXPONENTIATION_TCC_
#define EXPONENTIATION_TCC_

#include <libff/common/utils.hpp>

namespace libff {

using std::size_t;

template<typename FieldT, mp_size_t m>
FieldT power(const FieldT &base, const bigint<m> &exponent)
{
    FieldT result = FieldT::one();

    bool found_one = false;

    for (long i = exponent.max_bits() - 1; i >= 0; --i)
    {
        if (found_one)
        {
            result = result * result;
        }

        if (exponent.test_bit(i))
        {
            found_one = true;
            result = result * base;
        }
    }

    return result;
}

template<typename FieldT>
FieldT power(const FieldT &base, const unsigned long exponent)
{
    return power<FieldT>(base, bigint<1>(exponent));
}

template<typename FieldT, mp_size_t n>
FieldT tonelli_shanks_sqrt(const FieldT &value)
{
    // Find s and t such that p^n - 1 = t * 2^s, with t odd.
    const mp_size_t m = n * FieldT::extension_degree();
    bigint<m> one_i = bigint<m>(1);
    bigint<m> p_to_n_minus_1 = FieldT::field_char().template power<m>(FieldT::extension_degree()) - one_i;
    size_t s = 0;
    bigint<m> t = p_to_n_minus_1;
    while ((t % 2).is_zero())
    {
        s++;
        t = t / 2;
    }

    // Find a quadratic non-residue by generating random elements and testing with Euler's criterion.
    FieldT zero = FieldT::zero();
    FieldT neg_one = -FieldT::one();
    bigint<m> euler = p_to_n_minus_1 / 2;
    FieldT nqr;
    while (true)
    {
        nqr.randomize();
        if ((nqr^euler) == neg_one)
            break;
    }

    return tonelli_shanks_sqrt(value, s, nqr^t, (t - one_i) / 2);
}

template<typename FieldT, mp_size_t n>
FieldT tonelli_shanks_sqrt(const FieldT &value, const size_t s, const FieldT &nqr_to_t, const bigint<n> &t_minus_1_over_2)
{
    if (value.is_zero())
        return FieldT::zero();

    FieldT one = FieldT::one();

    size_t v = s;
    FieldT z = nqr_to_t;
    FieldT w = value^t_minus_1_over_2;
    FieldT x = value * w;
    FieldT b = x * w; // b = (*this)^t

#if DEBUG
    // check if square with euler's criterion
    FieldT check = b;
    for (size_t i = 0; i < v-1; ++i)
    {
        check = check.squared();
    }
    if (check != one)
    {
        assert(0);
    }
#endif

    // compute square root with Tonelli--Shanks
    // (does not terminate if not a square!)

    while (b != one)
    {
        size_t m = 0;
        FieldT b2m = b;
        while (b2m != one)
        {
            /* invariant: b2m = b^(2^m) after entering this loop */
            b2m = b2m.squared();
            m += 1;
        }

        int j = v-m-1;
        w = z;
        while (j > 0)
        {
            w = w.squared();
            --j;
        } // w = z^2^(v-m-1)

        z = w.squared();
        b = b * z;
        x = x * w;
        v = m;
    }

    return x;
}

} // libff

#endif // EXPONENTIATION_TCC_
