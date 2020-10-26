/** @file
 *****************************************************************************
 Declaration of bigint wrapper class around GMP's MPZ long integers.

 Notice that this class has no arithmetic operators. This is deliberate. All
 bigints should either be hardcoded or operated on the bit level to ensure
 high performance.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BIGINT_HPP_
#define BIGINT_HPP_
#include <cstddef>
#include <iostream>

#include <gmp.h>

#include <libff/common/serialization.hpp>

namespace libff {

template<mp_size_t n> class bigint;
template<mp_size_t n> std::ostream& operator<<(std::ostream &, const bigint<n>&);
template<mp_size_t n> std::istream& operator>>(std::istream &, bigint<n>&);

/**
 * Wrapper class around GMP's MPZ long integers. It supports arithmetic operations,
 * serialization and randomization. Serialization is fragile, see common/serialization.hpp.
 */

template<mp_size_t n>
class bigint {
public:
    static const mp_size_t N = n;

    mp_limb_t data[n] = {0};

    bigint() = default;
    bigint(const unsigned long x); /// Initalize from a small integer
    bigint(const char* s); /// Initialize from a string containing an integer in decimal notation
    bigint(const mpz_t r); /// Initialize from MPZ element

    static bigint one();

    void print() const;
    void print_hex() const;
    bool operator==(const bigint<n>& other) const;
    bool operator!=(const bigint<n>& other) const;
    bool operator<(const bigint<n>& other) const;
    void clear();
    bool is_zero() const;
    bool is_even() const;
    std::size_t max_bits() const { return n * GMP_NUMB_BITS; } /// Returns the number of bits representable by this bigint type
    std::size_t num_bits() const; /// Returns the number of bits in this specific bigint value, i.e., position of the most-significant 1

    unsigned long as_ulong() const; /// Return the last limb of the integer
    void to_mpz(mpz_t r) const;
    bool test_bit(const std::size_t bitno) const;

    bigint& randomize();

    friend std::ostream& operator<< <n>(std::ostream &out, const bigint<n> &b);
    friend std::istream& operator>> <n>(std::istream &in, bigint<n> &b);
};

} // libff
#include <libff/algebra/field_utils/bigint.tcc>
#endif
