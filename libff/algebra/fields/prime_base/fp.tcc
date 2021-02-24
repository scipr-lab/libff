/** @file
 *****************************************************************************
 Implementation of arithmetic in the finite field F[p], for prime p of fixed length.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP_TCC_
#define FP_TCC_
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

#include <libff/algebra/field_utils/field_utils.hpp>
#include <libff/algebra/field_utils/fp_aux.tcc>

namespace libff {

using std::size_t;

template<mp_size_t n, const bigint<n>& modulus>
void Fp_model<n,modulus>::mul_reduce(const bigint<n> &other)
{
    /* stupid pre-processor tricks; beware */
#if defined(__x86_64__) && defined(USE_ASM)
    if (n == 3)
    { // Use asm-optimized Comba multiplication and reduction
        mp_limb_t res[2*n];
        mp_limb_t c0, c1, c2;
        COMBA_3_BY_3_MUL(c0, c1, c2, res, this->mont_repr.data, other.data);

        mp_limb_t k;
        mp_limb_t tmp1, tmp2, tmp3;
        REDUCE_6_LIMB_PRODUCT(k, tmp1, tmp2, tmp3, inv, res, modulus.data);

        /* subtract t > mod */
        __asm__
            ("/* check for overflow */        \n\t"
             MONT_CMP(16)
             MONT_CMP(8)
             MONT_CMP(0)

             "/* subtract mod if overflow */  \n\t"
             "subtract%=:                     \n\t"
             MONT_FIRSTSUB
             MONT_NEXTSUB(8)
             MONT_NEXTSUB(16)
             "done%=:                         \n\t"
             :
             : [tmp] "r" (res+n), [M] "r" (modulus.data)
             : "cc", "memory", "%rax");
        mpn_copyi(this->mont_repr.data, res+n, n);
    }
    else if (n == 4)
    { // use asm-optimized "CIOS method"

        mp_limb_t tmp[n+1];
        mp_limb_t T0=0, T1=1, cy=2, u=3; // TODO: fix this

        __asm__ (MONT_PRECOMPUTE
                 MONT_FIRSTITER(1)
                 MONT_FIRSTITER(2)
                 MONT_FIRSTITER(3)
                 MONT_FINALIZE(3)
                 MONT_ITERFIRST(1)
                 MONT_ITERITER(1, 1)
                 MONT_ITERITER(1, 2)
                 MONT_ITERITER(1, 3)
                 MONT_FINALIZE(3)
                 MONT_ITERFIRST(2)
                 MONT_ITERITER(2, 1)
                 MONT_ITERITER(2, 2)
                 MONT_ITERITER(2, 3)
                 MONT_FINALIZE(3)
                 MONT_ITERFIRST(3)
                 MONT_ITERITER(3, 1)
                 MONT_ITERITER(3, 2)
                 MONT_ITERITER(3, 3)
                 MONT_FINALIZE(3)
                 "/* check for overflow */        \n\t"
                 MONT_CMP(24)
                 MONT_CMP(16)
                 MONT_CMP(8)
                 MONT_CMP(0)

                 "/* subtract mod if overflow */  \n\t"
                 "subtract%=:                     \n\t"
                 MONT_FIRSTSUB
                 MONT_NEXTSUB(8)
                 MONT_NEXTSUB(16)
                 MONT_NEXTSUB(24)
                 "done%=:                         \n\t"
                 :
                 : [tmp] "r" (tmp), [A] "r" (this->mont_repr.data), [B] "r" (other.data), [inv] "r" (inv), [M] "r" (modulus.data),
                   [T0] "r" (T0), [T1] "r" (T1), [cy] "r" (cy), [u] "r" (u)
                 : "cc", "memory", "%rax", "%rdx"
        );
        mpn_copyi(this->mont_repr.data, tmp, n);
    }
    else if (n == 5)
    { // use asm-optimized "CIOS method"

        mp_limb_t tmp[n+1];
        mp_limb_t T0=0, T1=1, cy=2, u=3; // TODO: fix this

        __asm__ (MONT_PRECOMPUTE
                 MONT_FIRSTITER(1)
                 MONT_FIRSTITER(2)
                 MONT_FIRSTITER(3)
                 MONT_FIRSTITER(4)
                 MONT_FINALIZE(4)
                 MONT_ITERFIRST(1)
                 MONT_ITERITER(1, 1)
                 MONT_ITERITER(1, 2)
                 MONT_ITERITER(1, 3)
                 MONT_ITERITER(1, 4)
                 MONT_FINALIZE(4)
                 MONT_ITERFIRST(2)
                 MONT_ITERITER(2, 1)
                 MONT_ITERITER(2, 2)
                 MONT_ITERITER(2, 3)
                 MONT_ITERITER(2, 4)
                 MONT_FINALIZE(4)
                 MONT_ITERFIRST(3)
                 MONT_ITERITER(3, 1)
                 MONT_ITERITER(3, 2)
                 MONT_ITERITER(3, 3)
                 MONT_ITERITER(3, 4)
                 MONT_FINALIZE(4)
                 MONT_ITERFIRST(4)
                 MONT_ITERITER(4, 1)
                 MONT_ITERITER(4, 2)
                 MONT_ITERITER(4, 3)
                 MONT_ITERITER(4, 4)
                 MONT_FINALIZE(4)
                 "/* check for overflow */        \n\t"
                 MONT_CMP(32)
                 MONT_CMP(24)
                 MONT_CMP(16)
                 MONT_CMP(8)
                 MONT_CMP(0)

                 "/* subtract mod if overflow */  \n\t"
                 "subtract%=:                     \n\t"
                 MONT_FIRSTSUB
                 MONT_NEXTSUB(8)
                 MONT_NEXTSUB(16)
                 MONT_NEXTSUB(24)
                 MONT_NEXTSUB(32)
                 "done%=:                         \n\t"
                 :
                 : [tmp] "r" (tmp), [A] "r" (this->mont_repr.data), [B] "r" (other.data), [inv] "r" (inv), [M] "r" (modulus.data),
                   [T0] "r" (T0), [T1] "r" (T1), [cy] "r" (cy), [u] "r" (u)
                 : "cc", "memory", "%rax", "%rdx"
        );
        mpn_copyi(this->mont_repr.data, tmp, n);
    }
    else
#endif
    {
        mp_limb_t res[2*n];
        mpn_mul_n(res, this->mont_repr.data, other.data, n);

        /*
          The Montgomery reduction here is based on Algorithm 14.32 in
          Handbook of Applied Cryptography
          <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.
         */
        for (size_t i = 0; i < n; ++i)
        {
            mp_limb_t k = inv * res[i];
            /* calculate res = res + k * mod * b^i */
            mp_limb_t carryout = mpn_addmul_1(res+i, modulus.data, n, k);
            carryout = mpn_add_1(res+n+i, res+n+i, n-i, carryout);
            assert(carryout == 0);
        }

        if (mpn_cmp(res+n, modulus.data, n) >= 0)
        {
            const mp_limb_t borrow = mpn_sub(res+n, res+n, n, modulus.data, n);
#ifndef NDEBUG
            assert(borrow == 0);
#else
            UNUSED(borrow);
#endif
        }

        mpn_copyi(this->mont_repr.data, res+n, n);
    }
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus>::Fp_model(const bigint<n> &b)
{
    mpn_copyi(this->mont_repr.data, Rsquared.data, n);
    mul_reduce(b);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus>::Fp_model(const long x, const bool is_unsigned)
{
    static_assert(std::numeric_limits<mp_limb_t>::max() >= static_cast<unsigned long>(std::numeric_limits<long>::max()), "long won't fit in mp_limb_t");
    if (is_unsigned || x >= 0)
    {
        this->mont_repr.data[0] = (mp_limb_t)x;
    }
    else
    {
        const mp_limb_t borrow = mpn_sub_1(this->mont_repr.data, modulus.data, n, (mp_limb_t)-x);
#ifndef NDEBUG
            assert(borrow == 0);
#else
            UNUSED(borrow);
#endif
    }

    mul_reduce(Rsquared);
}

template<mp_size_t n, const bigint<n>& modulus>
void Fp_model<n,modulus>::set_ulong(const unsigned long x)
{
    this->mont_repr.clear();
    this->mont_repr.data[0] = x;
    mul_reduce(Rsquared);
}

template<mp_size_t n, const bigint<n>& modulus>
void Fp_model<n,modulus>::clear()
{
    this->mont_repr.clear();
}

template<mp_size_t n, const bigint<n>& modulus>
void Fp_model<n,modulus>::randomize()
{
    (*this) = Fp_model<n, modulus>::random_element();
}

template<mp_size_t n, const bigint<n>& modulus>
bigint<n> Fp_model<n,modulus>::as_bigint() const
{
    bigint<n> one = bigint<n>::one();
    Fp_model<n, modulus> res(*this);
    res.mul_reduce(one);

    return (res.mont_repr);
}

template<mp_size_t n, const bigint<n>& modulus>
unsigned long Fp_model<n,modulus>::as_ulong() const
{
    return this->as_bigint().as_ulong();
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp_model<n,modulus>::operator==(const Fp_model& other) const
{
    return (this->mont_repr == other.mont_repr);
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp_model<n,modulus>::operator!=(const Fp_model& other) const
{
    return (this->mont_repr != other.mont_repr);
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp_model<n,modulus>::is_zero() const
{
    return (this->mont_repr.is_zero()); // zero maps to zero
}

template<mp_size_t n, const bigint<n>& modulus>
void Fp_model<n,modulus>::print() const
{
    Fp_model<n,modulus> tmp;
    tmp.mont_repr.data[0] = 1;
    tmp.mul_reduce(this->mont_repr);

    tmp.mont_repr.print();
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::zero()
{
    Fp_model<n,modulus> res;
    res.mont_repr.clear();
    return res;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::one()
{
    Fp_model<n,modulus> res;
    res.mont_repr.data[0] = 1;
    res.mul_reduce(Rsquared);
    return res;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::geometric_generator()
{
    Fp_model<n,modulus> res;
    res.mont_repr.data[0] = 2;
    res.mul_reduce(Rsquared);
    return res;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::arithmetic_generator()
{
    Fp_model<n,modulus> res;
    res.mont_repr.data[0] = 1;
    res.mul_reduce(Rsquared);
    return res;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus>& Fp_model<n,modulus>::operator+=(const Fp_model<n,modulus>& other)
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
#if defined(__x86_64__) && defined(USE_ASM)
    if (n == 3)
    {
        __asm__
            ("/* perform bignum addition */   \n\t"
             ADD_FIRSTADD
             ADD_NEXTADD(8)
             ADD_NEXTADD(16)
             "/* if overflow: subtract     */ \n\t"
             "/* (tricky point: if A and B are in the range we do not need to do anything special for the possible carry flag) */ \n\t"
             "jc      subtract%=              \n\t"

             "/* check for overflow */        \n\t"
             ADD_CMP(16)
             ADD_CMP(8)
             ADD_CMP(0)

             "/* subtract mod if overflow */  \n\t"
             "subtract%=:                     \n\t"
             ADD_FIRSTSUB
             ADD_NEXTSUB(8)
             ADD_NEXTSUB(16)
             "done%=:                         \n\t"
             :
             : [A] "r" (this->mont_repr.data), [B] "r" (other.mont_repr.data), [mod] "r" (modulus.data)
             : "cc", "memory", "%rax");
    }
    else if (n == 4)
    {
        __asm__
            ("/* perform bignum addition */   \n\t"
             ADD_FIRSTADD
             ADD_NEXTADD(8)
             ADD_NEXTADD(16)
             ADD_NEXTADD(24)
             "/* if overflow: subtract     */ \n\t"
             "/* (tricky point: if A and B are in the range we do not need to do anything special for the possible carry flag) */ \n\t"
             "jc      subtract%=              \n\t"

             "/* check for overflow */        \n\t"
             ADD_CMP(24)
             ADD_CMP(16)
             ADD_CMP(8)
             ADD_CMP(0)

             "/* subtract mod if overflow */  \n\t"
             "subtract%=:                     \n\t"
             ADD_FIRSTSUB
             ADD_NEXTSUB(8)
             ADD_NEXTSUB(16)
             ADD_NEXTSUB(24)
             "done%=:                         \n\t"
             :
             : [A] "r" (this->mont_repr.data), [B] "r" (other.mont_repr.data), [mod] "r" (modulus.data)
             : "cc", "memory", "%rax");
    }
    else if (n == 5)
    {
        __asm__
            ("/* perform bignum addition */   \n\t"
             ADD_FIRSTADD
             ADD_NEXTADD(8)
             ADD_NEXTADD(16)
             ADD_NEXTADD(24)
             ADD_NEXTADD(32)
             "/* if overflow: subtract     */ \n\t"
             "/* (tricky point: if A and B are in the range we do not need to do anything special for the possible carry flag) */ \n\t"
             "jc      subtract%=              \n\t"

             "/* check for overflow */        \n\t"
             ADD_CMP(32)
             ADD_CMP(24)
             ADD_CMP(16)
             ADD_CMP(8)
             ADD_CMP(0)

             "/* subtract mod if overflow */  \n\t"
             "subtract%=:                     \n\t"
             ADD_FIRSTSUB
             ADD_NEXTSUB(8)
             ADD_NEXTSUB(16)
             ADD_NEXTSUB(24)
             ADD_NEXTSUB(32)
             "done%=:                         \n\t"
             :
             : [A] "r" (this->mont_repr.data), [B] "r" (other.mont_repr.data), [mod] "r" (modulus.data)
             : "cc", "memory", "%rax");
    }
    else
#endif
    {
        mp_limb_t scratch[n+1];
        const mp_limb_t carry = mpn_add_n(scratch, this->mont_repr.data, other.mont_repr.data, n);
        scratch[n] = carry;

        if (carry || mpn_cmp(scratch, modulus.data, n) >= 0)
        {
            const mp_limb_t borrow = mpn_sub(scratch, scratch, n+1, modulus.data, n);
#ifndef NDEBUG
            assert(borrow == 0);
#else
            UNUSED(borrow);
#endif
        }

        mpn_copyi(this->mont_repr.data, scratch, n);
    }

    return *this;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus>& Fp_model<n,modulus>::operator-=(const Fp_model<n,modulus>& other)
{
#ifdef PROFILE_OP_COUNTS
    this->sub_cnt++;
#endif
#if defined(__x86_64__) && defined(USE_ASM)
    if (n == 3)
    {
        __asm__
            (SUB_FIRSTSUB
             SUB_NEXTSUB(8)
             SUB_NEXTSUB(16)

             "jnc     done%=\n\t"

             SUB_FIRSTADD
             SUB_NEXTADD(8)
             SUB_NEXTADD(16)

             "done%=:\n\t"
             :
             : [A] "r" (this->mont_repr.data), [B] "r" (other.mont_repr.data), [mod] "r" (modulus.data)
             : "cc", "memory", "%rax");
    }
    else if (n == 4)
    {
        __asm__
            (SUB_FIRSTSUB
             SUB_NEXTSUB(8)
             SUB_NEXTSUB(16)
             SUB_NEXTSUB(24)

             "jnc     done%=\n\t"

             SUB_FIRSTADD
             SUB_NEXTADD(8)
             SUB_NEXTADD(16)
             SUB_NEXTADD(24)

             "done%=:\n\t"
             :
             : [A] "r" (this->mont_repr.data), [B] "r" (other.mont_repr.data), [mod] "r" (modulus.data)
             : "cc", "memory", "%rax");
    }
    else if (n == 5)
    {
        __asm__
            (SUB_FIRSTSUB
             SUB_NEXTSUB(8)
             SUB_NEXTSUB(16)
             SUB_NEXTSUB(24)
             SUB_NEXTSUB(32)

             "jnc     done%=\n\t"

             SUB_FIRSTADD
             SUB_NEXTADD(8)
             SUB_NEXTADD(16)
             SUB_NEXTADD(24)
             SUB_NEXTADD(32)

             "done%=:\n\t"
             :
             : [A] "r" (this->mont_repr.data), [B] "r" (other.mont_repr.data), [mod] "r" (modulus.data)
             : "cc", "memory", "%rax");
    }
    else
#endif
    {
        mp_limb_t scratch[n+1];
        if (mpn_cmp(this->mont_repr.data, other.mont_repr.data, n) < 0)
        {
            const mp_limb_t carry = mpn_add_n(scratch, this->mont_repr.data, modulus.data, n);
            scratch[n] = carry;
        }
        else
        {
            mpn_copyi(scratch, this->mont_repr.data, n);
            scratch[n] = 0;
        }

        const mp_limb_t borrow = mpn_sub(scratch, scratch, n+1, other.mont_repr.data, n);
#ifndef NDEBUG
        assert(borrow == 0);
#else
        UNUSED(borrow);
#endif

        mpn_copyi(this->mont_repr.data, scratch, n);
    }
    return *this;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus>& Fp_model<n,modulus>::operator*=(const Fp_model<n,modulus>& other)
{
#ifdef PROFILE_OP_COUNTS
    this->mul_cnt++;
#endif

    mul_reduce(other.mont_repr);
    return *this;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus>& Fp_model<n,modulus>::operator^=(const unsigned long pow)
{
    (*this) = power<Fp_model<n, modulus> >(*this, pow);
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m>
Fp_model<n,modulus>& Fp_model<n,modulus>::operator^=(const bigint<m> &pow)
{
    (*this) = power<Fp_model<n, modulus>, m>(*this, pow);
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::operator+(const Fp_model<n,modulus>& other) const
{
    Fp_model<n, modulus> r(*this);
    return (r += other);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::operator-(const Fp_model<n,modulus>& other) const
{
    Fp_model<n, modulus> r(*this);
    return (r -= other);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::operator*(const Fp_model<n,modulus>& other) const
{
    Fp_model<n, modulus> r(*this);
    return (r *= other);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::operator^(const unsigned long pow) const
{
    Fp_model<n, modulus> r(*this);
    return (r ^= pow);
}

template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m>
Fp_model<n,modulus> Fp_model<n,modulus>::operator^(const bigint<m> &pow) const
{
    Fp_model<n, modulus> r(*this);
    return (r ^= pow);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::operator-() const
{
#ifdef PROFILE_OP_COUNTS
    this->sub_cnt++;
#endif

    if (this->is_zero())
    {
        return (*this);
    }
    else
    {
        Fp_model<n, modulus> r;
        mpn_sub_n(r.mont_repr.data, modulus.data, this->mont_repr.data, n);
        return r;
    }
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::squared() const
{
#ifdef PROFILE_OP_COUNTS
    this->sqr_cnt++;
    this->mul_cnt--; // zero out the upcoming mul
#endif
    /* stupid pre-processor tricks; beware */
#if defined(__x86_64__) && defined(USE_ASM)
    if (n == 3)
    { // use asm-optimized Comba squaring
        mp_limb_t res[2*n];
        mp_limb_t c0, c1, c2;
        COMBA_3_BY_3_SQR(c0, c1, c2, res, this->mont_repr.data);

        mp_limb_t k;
        mp_limb_t tmp1, tmp2, tmp3;
        REDUCE_6_LIMB_PRODUCT(k, tmp1, tmp2, tmp3, inv, res, modulus.data);

        /* subtract t > mod */
        __asm__ volatile
            ("/* check for overflow */        \n\t"
             MONT_CMP(16)
             MONT_CMP(8)
             MONT_CMP(0)

             "/* subtract mod if overflow */  \n\t"
             "subtract%=:                     \n\t"
             MONT_FIRSTSUB
             MONT_NEXTSUB(8)
             MONT_NEXTSUB(16)
             "done%=:                         \n\t"
             :
             : [tmp] "r" (res+n), [M] "r" (modulus.data)
             : "cc", "memory", "%rax");

        Fp_model<n, modulus> r;
        mpn_copyi(r.mont_repr.data, res+n, n);
        return r;
    }
    else
#endif
    {
        Fp_model<n, modulus> r(*this);
        return (r *= r);
    }
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus>& Fp_model<n,modulus>::square()
{
    (*this) = squared();
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus>& Fp_model<n,modulus>::invert()
{
#ifdef PROFILE_OP_COUNTS
    this->inv_cnt++;
#endif

    assert(!this->is_zero());

    bigint<n> g; /* gp should have room for vn = n limbs */

    mp_limb_t s[n+1]; /* sp should have room for vn+1 limbs */
    mp_size_t sn;

    bigint<n> v = modulus; // both source operands are destroyed by mpn_gcdext

    /* computes gcd(u, v) = g = u*s + v*t, so s*u will be 1 (mod v) */
    const mp_size_t gn = mpn_gcdext(g.data, s, &sn, this->mont_repr.data, n, v.data, n);
#ifndef NDEBUG
    assert(gn == 1 && g.data[0] == 1); /* inverse exists */
#else
    UNUSED(gn);
#endif

    mp_limb_t q; /* division result fits into q, as sn <= n+1 */
    /* sn < 0 indicates negative sn; will fix up later */

    if (std::abs(sn) >= n)
    {
        /* if sn could require modulus reduction, do it here */
        mpn_tdiv_qr(&q, this->mont_repr.data, 0, s, std::abs(sn), modulus.data, n);
    }
    else
    {
        /* otherwise just copy it over */
        mpn_zero(this->mont_repr.data, n);
        mpn_copyi(this->mont_repr.data, s, std::abs(sn));
    }

    /* fix up the negative sn */
    if (sn < 0)
    {
        const mp_limb_t borrow = mpn_sub_n(this->mont_repr.data, modulus.data, this->mont_repr.data, n);
#ifndef NDEBUG
        assert(borrow == 0);
#else
        UNUSED(borrow);
#endif
    }

    mul_reduce(Rcubed);
    return *this;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::inverse() const
{
    Fp_model<n, modulus> r(*this);
    return (r.invert());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::Frobenius_map(unsigned long power) const
{
    UNUSED(power); // only for API consistency
    Fp_model<n,modulus> copy = *this;
    return copy;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::random_element() /// returns random element of Fp_model
{
    /* note that as Montgomery representation is a bijection then
       selecting a random element of {xR} is the same as selecting a
       random element of {x} */
    Fp_model<n, modulus> r;
    do
    {
        r.mont_repr.randomize();

        /* clear all bits higher than MSB of modulus */
        size_t bitno = GMP_NUMB_BITS * n - 1;
        while (modulus.test_bit(bitno) == false)
        {
            const std::size_t part = bitno/GMP_NUMB_BITS;
            const std::size_t bit = bitno - (GMP_NUMB_BITS*part);

            static const mp_limb_t one = 1;
            r.mont_repr.data[part] &= ~(one<<bit);

            bitno--;
        }
    }
   /* if r.data is still >= modulus -- repeat (rejection sampling) */
    while (mpn_cmp(r.mont_repr.data, modulus.data, n) >= 0);

    return r;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp_model<n,modulus> Fp_model<n,modulus>::sqrt() const
{
    return tonelli_shanks_sqrt(*this);
}

template<mp_size_t n, const bigint<n>& modulus>
std::vector<uint64_t> Fp_model<n,modulus>::to_words() const
{
    // TODO: implement for other bit architectures
    static_assert(GMP_NUMB_BITS == 64, "Only 64-bit architectures are currently supported");

    bigint<n> repr = this->bigint_repr();
    std::vector<uint64_t> words;
	words.insert(words.begin(), std::begin(repr.data), std::end(repr.data));
    return words;
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp_model<n,modulus>::from_words(std::vector<uint64_t> words)
{
    // TODO: implement for other bit architectures
    static_assert(GMP_NUMB_BITS == 64, "Only 64-bit architectures are currently supported");

    typedef Fp_model<n, modulus> FieldT; // Without the typedef C++ doesn't compile.
    long start_bit = words.size() * 64 - FieldT::ceil_size_in_bits();
    assert(start_bit >= 0); // Check the vector is big enough.
    long start_word = start_bit / 64;
    long bit_offset = start_bit % 64;

    // Assumes mont_repr.data is just the right size to fit ceil_size_in_bits().
    std::copy(words.begin() + start_word, words.end(), this->mont_repr.data);
    // Zero out the left-most bit_offset bits.
    this->mont_repr.data[n - 1] = mp_limb_t((uint64_t(this->mont_repr.data[n - 1]) << bit_offset) >> bit_offset);
#ifndef MONTGOMERY_OUTPUT
    this->mul_reduce(Rsquared);
#endif
    return this->mont_repr < modulus;
}

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &out, const Fp_model<n, modulus> &p)
{
    out << p.bigint_repr();
    return out;
}

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &in, Fp_model<n, modulus> &p)
{
    in >> p.mont_repr;
#ifndef MONTGOMERY_OUTPUT
    p.mul_reduce(Fp_model<n,modulus>::Rsquared);
#endif
    return in;
}

template<mp_size_t n, const bigint<n>& modulus>
bigint<n> Fp_model<n,modulus>::bigint_repr() const
{
    // If the flag is defined, serialization and words output use the montgomery representation
    // instead of the human-readable value.
#ifdef MONTGOMERY_OUTPUT
    return this->mont_repr;
#else
    return this->as_bigint();
#endif
}

} // namespace libff
#endif // FP_TCC_
