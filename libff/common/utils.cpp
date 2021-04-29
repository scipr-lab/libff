/** @file
 *****************************************************************************
 Implementation of misc math and serialization utility functions.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <algorithm>
#include <cassert>
#include <cstdarg>
#include <cstdint>

#include <libff/common/utils.hpp>

namespace libff {

using std::size_t;

/**
 * Round n to the next power of two.
 * If n is a power of two, return n
 */
size_t get_power_of_two(size_t n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n |= n >> 32;
    n++;

    return n;
}

/* If n is a power of 2, returns n */
size_t round_to_next_power_of_2(const size_t n)
{
    return (1ULL << log2(n));
}

bool is_power_of_2(const size_t n)
{
    return ((n != 0) && ((n & (n-1)) == 0));
}

size_t log2(size_t n)
/* returns ceil(log2(n)), so 1ul<<log2(n) is the smallest power of 2,
   that is not less than n. */
{
    size_t r = ((n & (n-1)) == 0 ? 0 : 1); // add 1 if n is not power of 2

    while (n > 1)
    {
        n >>= 1;
        r++;
    }

    return r;
}

size_t to_twos_complement(int i, size_t w)
{
    assert(i >= -(1l<<(w-1)));
    assert(i < (1l<<(w-1)));
    return (i >= 0) ? i : i + (1L<<w);
}

int from_twos_complement(size_t i, size_t w)
{
    assert(i < (1UL<<w));
    return (i < (1UL<<(w-1))) ? i : i - (1UL<<w);
}

size_t bitreverse(size_t n, const size_t l)
{
    size_t r = 0;
    for (size_t k = 0; k < l; ++k)
    {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }
    return r;
}

bit_vector int_list_to_bits(const std::initializer_list<unsigned long> &l, const size_t wordsize)
{
    bit_vector res(wordsize*l.size());
    for (size_t i = 0; i < l.size(); ++i)
    {
        for (size_t j = 0; j < wordsize; ++j)
        {
            res[i*wordsize + j] = (*(l.begin()+i) & (1UL<<(wordsize-1-j))) != 0U;
        }
    }
    return res;
}

long long div_ceil(long long x, long long y)
{
    if (y == 0)
    {
        throw std::invalid_argument("libff::div_ceil: division by zero, second argument must be non-zero");
    }
    return (x + (y-1)) / y;
}

bool is_little_endian()
{
    uint64_t a = 0x12345678;
    unsigned char *c = (unsigned char*)(&a);
    return (*c == 0x78);
}

std::string FORMAT(const std::string &prefix, const char* format, ...)
{
    const static size_t MAX_FMT = 256;
    char buf[MAX_FMT];
    va_list args;
    va_start(args, format);
    vsnprintf(buf, MAX_FMT, format, args);
    va_end(args);

    return prefix + std::string(buf);
}

void serialize_bit_vector(std::ostream &out, const bit_vector &v)
{
    out << v.size() << "\n";
    for (auto b : v)
    {
        out << b << "\n";
    }
}

void deserialize_bit_vector(std::istream &in, bit_vector &v)
{
    size_t size;
    in >> size;
    v.resize(size);
    for (size_t i = 0; i < size; ++i)
    {
        bool b;
        in >> b;
        v[i] = b;
    }
}

} // namespace libff
