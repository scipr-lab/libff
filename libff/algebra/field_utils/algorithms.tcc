/** @file
 *****************************************************************************

 Declaration of interfaces for (square-and-multiply) exponentiation and
 Tonelli-Shanks square root.

 See algorithms.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef ALGORITHMS_TCC_
#define ALGORITHMS_TCC_

#include "libff/common/utils.hpp"
#include "libff/common/profiling.hpp"

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
void find_tonelli_shanks_constants()
{
    // Find s and t such that p^n - 1 = t * 2^s, with t odd.
    const mp_size_t m = n * FieldT::extension_degree();
    bigint<m> one_i = bigint<m>(1); // integer one
    bigint<m> p_to_n_minus_1 = FieldT::field_char().template power<m>(FieldT::extension_degree()) - one_i;
    size_t s = 0;
    bigint<m> t = p_to_n_minus_1;
    while ((t % 2).is_zero())
    {
        s++;
        t = t / 2;
    }

    // Find a quadratic non-residue by generating random elements and testing with Euler's criterion.
    FieldT neg_one = -FieldT::one();
    bigint<m> euler = p_to_n_minus_1 / 2;
    FieldT nqr;
    while (true)
    {
        nqr.randomize();
        if ((nqr^euler) == neg_one)
            break;
    }

    FieldT::euler = euler;
    FieldT::s = s;
    FieldT::t = t;
    FieldT::t_minus_1_over_2 = (t - one_i) / 2;
    FieldT::nqr = nqr;
    FieldT::nqr_to_t = nqr^t;
}

template<typename FieldT>
FieldT tonelli_shanks_sqrt(const FieldT &value)
{
    // a few assertions to make sure s, t, and nqr are initialized
    assert(FieldT::s != 0);
    assert(!(FieldT::t % 2).is_zero()); // check that t is odd
    assert(!FieldT::nqr.is_zero());

    if (value.is_zero())
        return FieldT::zero();

    FieldT one = FieldT::one();

    size_t v = FieldT::s;
    FieldT z = FieldT::nqr_to_t;
    FieldT w = value^FieldT::t_minus_1_over_2;
    FieldT x = value * w;
    FieldT b = x * w; // b = (*this)^t

#if DEBUG
    // check if square with euler's criterion
    FieldT check = b;
    for (size_t i = 0; i < v-1; ++i)
    {
        check = check.squared();
    }
    assert(check == one);
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

#endif // ALGORITHMS_TCC_
