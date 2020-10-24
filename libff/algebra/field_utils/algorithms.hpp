/** @file
 *****************************************************************************

 Declaration of interfaces for (square-and-multiply) exponentiation and
 Tonelli-Shanks square root.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef ALGORITHMS_HPP_
#define ALGORITHMS_HPP_

#include <cstdint>

#include "libff/algebra/field_utils/bigint.hpp"

namespace libff {

/** Repeated squaring. */
template<typename FieldT, mp_size_t m>
FieldT power(const FieldT &base, const bigint<m> &exponent);

/** Repeated squaring. */
template<typename FieldT>
FieldT power(const FieldT &base, const unsigned long exponent);

/**
 * Tonelli-Shanks square root with given s, t, and quadratic non-residue.
 * Only terminates if there is a square root. Only works if required parameters
 * are set in the field class.
 */
template<typename FieldT>
FieldT tonelli_shanks_sqrt(const FieldT &value);

} // libff

#include <libff/algebra/field_utils/algorithms.tcc>

#endif // ALGORITHMS_HPP_
