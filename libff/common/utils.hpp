/** @file
 *****************************************************************************
 Declaration of miscellaneous math, serialization, and other common utility
 functions.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace libff {

typedef std::vector<bool> bit_vector;

template<bool B, class T = void>
struct enable_if { typedef void* type; };

template<class T>
struct enable_if<true, T> { typedef T type; };

std::size_t get_power_of_two(std::size_t n);

std::size_t round_to_next_power_of_2(const std::size_t n);
bool is_power_of_2(const std::size_t n);

/// returns ceil(log2(n)), so 1ul<<log2(n) is the smallest power of 2, that is not less than n
std::size_t log2(std::size_t n);

inline std::size_t exp2(std::size_t k) { return std::size_t(1) << k; }

std::size_t to_twos_complement(int i, std::size_t w);
int from_twos_complement(std::size_t i, std::size_t w);

std::size_t bitreverse(std::size_t n, const std::size_t l);
bit_vector int_list_to_bits(const std::initializer_list<unsigned long> &l, const std::size_t wordsize);
/* throws error if y = 0 */
long long div_ceil(long long x, long long y);

bool is_little_endian();

std::string FORMAT(const std::string &prefix, const char* format, ...);

/* A variadic template to suppress unused argument warnings */
template<typename ... Types>
void UNUSED(Types&&...) {}

#ifdef DEBUG
#define FMT libff::FORMAT
#else
#define FMT(...) (libff::UNUSED(__VA_ARGS__), "")
#endif

void serialize_bit_vector(std::ostream &out, const bit_vector &v);
void deserialize_bit_vector(std::istream &in, bit_vector &v);

/** Should not be used for fields, because the field function is named ceil_size_in_bits instead. */
template<typename CurveT>
std::size_t curve_size_in_bits(const std::vector<CurveT> &v);

/* Print a vector in the form { elem0 elem1 elem2 ... }, with a newline at the end
template<typename T>
void print_vector(std::vector<T> &vec);
template<typename T>
void print_vector(std::vector<T> vec);*/

template<typename T>
void print_vector(std::vector<T> &vec)
{
    printf("{ ");
    for (auto const& elem : vec)
    {
        std::cout << elem << " ";
    }
    printf("}\n");
}

template<typename T>
void print_vector(std::vector<T> vec)
{
    printf("{ ");
    for (auto const& elem : vec)
    {
        std::cout << elem << " ";
    }
    printf("}\n");
}

/**
 * Returns a random element of T that is not zero or one.
 * T can be a field or elliptic curve group.
 * Used for testing to generate a test example that doesn't error.
 */
template<typename T>
T random_element_non_zero_one();
/**
 * Returns a random element of T that is not zero.
 * T can be a field or elliptic curve group.
 * Used for testing to generate a test example that doesn't error.
 */
template<typename T>
T random_element_non_zero();
/**
 * Returns a random element of T that is not equal to y.
 * T can be a field or elliptic curve group.
 * Used for testing to generate a test example that doesn't error.
 */
template<typename T>
T random_element_exclude(T y);

#define ARRAY_SIZE(arr) (sizeof(arr)/sizeof(arr[0]))

} // namespace libff

#include <libff/common/utils.tcc> /* note that utils has a templatized part (utils.tcc) and non-templatized part (utils.cpp) */
#endif // UTILS_HPP_
