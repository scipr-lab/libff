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

#include <libff/algebra/field_utils/bigint.hpp>

namespace libff {

/** Repeated squaring. */
template<typename FieldT, mp_size_t m>
FieldT power(const FieldT &base, const bigint<m> &exponent);

/** Repeated squaring. */
template<typename FieldT>
FieldT power(const FieldT &base, const unsigned long exponent);

/** Tonelli-Shanks square root with given s, t, and quadratic non-residue.
 *  Only terminates if there is a square root. */
template<typename FieldT, mp_size_t n>
FieldT tonelli_shanks_sqrt(const FieldT &value, const std::size_t s, const FieldT &nqr_to_t, const bigint<n> &t_minus_1_over_2);

/** Tonelli-Shanks square root: compute s, t, and quadratic non-residue as part of function.
 *  Only terminates if there is a square root. */
template<typename FieldT, mp_size_t n>
FieldT tonelli_shanks_sqrt(const FieldT &value);

} // libff

#include <libff/algebra/field_utils/algorithms.tcc>

#endif // EXPONENTIATION_HPP_
