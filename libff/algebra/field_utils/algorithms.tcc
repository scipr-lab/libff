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

template<typename FieldT>
FieldT tonelli_shanks_sqrt(const FieldT &value)
{
    // A few assertions to make sure s, t, and nqr are initialized.
    assert(FieldT::s != 0);
    assert(!FieldT::t.is_even()); // Check that t is odd.
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
