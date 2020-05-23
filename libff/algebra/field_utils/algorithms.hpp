/** @file
 *****************************************************************************

 Declaration of interfaces for (square-and-multiply) exponentiation.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef EXPONENTIATION_HPP_
#define EXPONENTIATION_HPP_

#include <cstdint>

#include "libff/algebra/field_utils/bigint.hpp"

namespace libff {

/** Repeated squaring. */
template<typename FieldT, mp_size_t m>
FieldT power(const FieldT &base, const bigint<m> &exponent);

/** Repeated squaring. */
template<typename FieldT>
FieldT power(const FieldT &base, const unsigned long exponent);

/** Compute s, t, and quadratic non-residue as well as (p^n-1)/2, (t-1)/2, and nqr^t.
 *  Stores these constants inside static field attributes.
 *  Only needs to be called once for each field. */
template<typename FieldT, mp_size_t n>
void find_tonelli_shanks_constants();

/** Tonelli-Shanks square root with given s, t, and quadratic non-residue.
 *  Only terminates if there is a square root. */
template<typename FieldT>
FieldT tonelli_shanks_sqrt(const FieldT &value);

} // libff

#include <libff/algebra/field_utils/algorithms.tcc>

#endif // EXPONENTIATION_HPP_
