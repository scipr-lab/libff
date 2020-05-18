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

/** Tonelli-Shanks square root. */
template<typename FieldT>
FieldT tonelli_shanks_sqrt(const FieldT &value);

} // libff

#include <libff/algebra/field_utils/algorithms.tcc>

#endif // EXPONENTIATION_HPP_
