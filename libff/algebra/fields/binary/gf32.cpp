#include <cstdio>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <sodium/randombytes.h>

#include "libff/algebra/field_utils/algorithms.hpp"
#include "libff/algebra/fields/binary/gf32.hpp"

#ifdef USE_ASM
#include <emmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>
#endif

namespace libff {

using std::size_t;

const uint64_t gf32::modulus_;
gf32 gf32::multiplicative_generator = gf32(2);

gf32::gf32() : value_(0)
{
}

gf32::gf32(const uint32_t value) : value_(value)
{
}

std::vector<uint64_t> gf32::to_words() const
{
    return std::vector<uint64_t>({uint64_t(this->value_)});
}

bool gf32::from_words(std::vector<uint64_t> words)
{
    this->value_ = uint32_t(words[0]);
    return true;
}

gf32& gf32::operator+=(const gf32 &other)
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    this->value_ ^= other.value_;
    return (*this);
}

gf32& gf32::operator-=(const gf32 &other)
{
#ifdef PROFILE_OP_COUNTS
    this->sub_cnt++;
#endif
    this->value_ ^= other.value_;
    return (*this);
}

// multiplication over GF(2^k) is carryless multiplication
gf32& gf32::operator*=(const gf32 &other)
{
#ifdef PROFILE_OP_COUNTS
    this->mul_cnt++;
#endif
    /* Does not require *this and other to be different, and therefore
       also works for squaring, implemented below. */

    /* Slow, but straight-forward */
    uint32_t result = 0;
    uint32_t shifted = this->value_;

    for (uint32_t i = 0; i < 32; ++i)
    {
        if ((other.value_ & (1ULL << i)) != 0U)
        {
            result ^= shifted;
        }
        if ((shifted & (1UL << 31)) != 0U)
        {
            shifted <<= 1;
            shifted ^= libff::gf32::modulus_;
        }
        else
        {
            shifted <<= 1;
        }
    }

    this->value_ = result;

    return (*this);
}

gf32& gf32::operator^=(const unsigned long pow)
{
    (*this) = *this ^ pow;
    return (*this);
}

gf32& gf32::square()
{
#ifdef PROFILE_OP_COUNTS
    this->sqr_cnt++;
    this->mul_cnt--;
#endif
    this->operator*=(*this);
    return *this;
}

gf32& gf32::invert()
{
    (*this) = inverse();
    return (*this);
}

gf32 gf32::operator+(const gf32 &other) const
{
    gf32 result(*this);
    return (result+=(other));
}

gf32 gf32::operator-(const gf32 &other) const
{
    gf32 result(*this);
    return (result-=(other));
}

gf32 gf32::operator-() const
{
    /* additive inverse matches the element itself */
    return gf32(*this);
}

gf32 gf32::operator*(const gf32 &other) const
{
    gf32 result(*this);
    return (result*=(other));
}

gf32 gf32::operator^(const unsigned long pow) const
{
    return power<gf32>(*this, pow);
}

gf32 gf32::squared() const
{
    gf32 result(*this);
    result.square();
    return result;
}

// repeatedly square pt, num_times. For use in inverse.
void square_multi(gf32* pt, int8_t num_times)
{
    for (int8_t i = 0; i < num_times; i++)
    {
        (*pt).square();
    }
}

/* calculate el^{-1} as el^{2^{32}-2}. */
gf32 gf32::inverse() const
{
#ifdef PROFILE_OP_COUNTS
    this->inv_cnt++;
    this->mul_cnt -= 9;
    this->sqr_cnt -= 31;
#endif
    assert(!this->is_zero());
    gf32 a(*this);

    gf32 result(0);
    for (size_t i = 0; i <= 4; ++i)
    {
        /* entering the loop a = el^{2^{2^i}-1} */
        gf32 b = a;
        for (size_t j = 0; j < (1UL<<i); ++j)
        {
            b.square();
        }
        /* after the loop b = a^{2^i} = el^{2^{2^i}*(2^{2^i}-1)} */
        a *= b;
        /* now a = el^{2^{2^{i+1}}-1} */

        if (i == 0)
        {
            result = b;
        }
        else
        {
            result *= b;
        }
    }
    /* now result = el^{2^32-2} */
    return result;
}

gf32 gf32::sqrt() const
{
    return (*this)^bigint<1>("2147483648"); // 2^31
}

void gf32::randomize()
{
    randombytes_buf(&this->value_, 32/8);
}

void gf32::clear()
{
    this->value_ = 0;
}

bool gf32::operator==(const gf32 &other) const
{
    return (this->value_ == other.value_);
}

bool gf32::operator!=(const gf32 &other) const
{
    return !(this->operator==(other));
}

void gf32::print() const
{
    printf("%u\n", this->value_);
}

bool gf32::is_zero() const
{
    return (this->value_ == 0);
}

gf32 gf32::zero()
{
    return gf32(0);
}

gf32 gf32::one()
{
    return gf32(1);
}

gf32 gf32::random_element()
{
    gf32 result;
    result.randomize();
    return result;
}

std::ostream& operator<<(std::ostream &out, const gf32 &el)
{
    out << el.value_;
    return out;
}

std::istream& operator>>(std::istream &in, gf32 &el)
{
    in >> el.value_;
    return in;
}

} // namespace libff
