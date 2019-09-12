/** @file
 *****************************************************************************

 Declaration of interfaces for wNAF ("width-w Non-Adjacent Form") exponentiation routines.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef WNAF_HPP_
#define WNAF_HPP_

#include <vector>

#include <libff/algebra/fields/bigint.hpp>

namespace libff {

/**
 * Find the wNAF representation of the given scalar relative to the given
 * window size, reusing the given vector to store it.
 */
template<mp_size_t n>
void update_wnaf(std::vector<long> &naf, const size_t window_size, const bigint<n> &scalar);

/**
 * Find the wNAF representation of the given scalar relative to the given window size.
 */
template<mp_size_t n>
std::vector<long> find_wnaf(const size_t window_size, const bigint<n> &scalar);

/**
 * Compute optimal window size.
 */
template<typename T>
size_t wnaf_opt_window_size(const size_t scalar_bits);

/**
 * In additive notation, use wNAF exponentiation (with the given window size)
 * to compute scalar * base, for a scalar in naf form.
 */
template<typename T>
T fixed_window_wnaf_exp(const size_t window_size, const T &base, const std::vector<long> &naf);

/**
 * In additive notation, use wNAF exponentiation (with the given window size) to compute scalar * base.
 */
template<typename T, mp_size_t n>
T fixed_window_wnaf_exp(const size_t window_size, const T &base, const bigint<n> &scalar);

/**
 * In additive notation, use wNAF exponentiation (with the window size determined by T) to compute scalar * base.
 */
template<typename T, mp_size_t n>
T opt_window_wnaf_exp(const T &base, const bigint<n> &scalar, const size_t scalar_bits);

} // libff

#include <libff/algebra/scalar_multiplication/wnaf.tcc>

#endif // WNAF_HPP_
