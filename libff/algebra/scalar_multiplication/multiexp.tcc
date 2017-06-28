/** @file
 *****************************************************************************

 Implementation of interfaces for multi-exponentiation routines.

 See multiexp.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MULTIEXP_TCC_
#define MULTIEXP_TCC_

#include <algorithm>
#include <cassert>
#include <type_traits>

#include <libff/algebra/fields/bigint.hpp>
#include <libff/algebra/fields/fp_aux.tcc>
#include <libff/algebra/scalar_multiplication/wnaf.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

namespace libff {

template<mp_size_t n>
class ordered_exponent {
// to use std::push_heap and friends later
public:
    size_t idx;
    bigint<n> r;

    ordered_exponent(const size_t idx, const bigint<n> &r) : idx(idx), r(r) {};

    bool operator<(const ordered_exponent<n> &other) const
    {
#if defined(__x86_64__) && defined(USE_ASM)
        if (n == 3)
        {
            long res;
            __asm__
                ("// check for overflow           \n\t"
                 "mov $0, %[res]                  \n\t"
                 ADD_CMP(16)
                 ADD_CMP(8)
                 ADD_CMP(0)
                 "jmp done%=                      \n\t"
                 "subtract%=:                     \n\t"
                 "mov $1, %[res]                  \n\t"
                 "done%=:                         \n\t"
                 : [res] "=&r" (res)
                 : [A] "r" (other.r.data), [mod] "r" (this->r.data)
                 : "cc", "%rax");
            return res;
        }
        else if (n == 4)
        {
            long res;
            __asm__
                ("// check for overflow           \n\t"
                 "mov $0, %[res]                  \n\t"
                 ADD_CMP(24)
                 ADD_CMP(16)
                 ADD_CMP(8)
                 ADD_CMP(0)
                 "jmp done%=                      \n\t"
                 "subtract%=:                     \n\t"
                 "mov $1, %[res]                  \n\t"
                 "done%=:                         \n\t"
                 : [res] "=&r" (res)
                 : [A] "r" (other.r.data), [mod] "r" (this->r.data)
                 : "cc", "%rax");
            return res;
        }
        else if (n == 5)
        {
            long res;
            __asm__
                ("// check for overflow           \n\t"
                 "mov $0, %[res]                  \n\t"
                 ADD_CMP(32)
                 ADD_CMP(24)
                 ADD_CMP(16)
                 ADD_CMP(8)
                 ADD_CMP(0)
                 "jmp done%=                      \n\t"
                 "subtract%=:                     \n\t"
                 "mov $1, %[res]                  \n\t"
                 "done%=:                         \n\t"
                 : [res] "=&r" (res)
                 : [A] "r" (other.r.data), [mod] "r" (this->r.data)
                 : "cc", "%rax");
            return res;
        }
        else
#endif
        {
            return (mpn_cmp(this->r.data, other.r.data, n) < 0);
        }
    }
};

template<typename T, typename FieldT>
T naive_exp(typename std::vector<T>::const_iterator vec_start,
            typename std::vector<T>::const_iterator vec_end,
            typename std::vector<FieldT>::const_iterator scalar_start,
            typename std::vector<FieldT>::const_iterator scalar_end)
{
    T result(T::zero());

    typename std::vector<T>::const_iterator vec_it;
    typename std::vector<FieldT>::const_iterator scalar_it;

    for (vec_it = vec_start, scalar_it = scalar_start; vec_it != vec_end; ++vec_it, ++scalar_it)
    {
        bigint<FieldT::num_limbs> scalar_bigint = scalar_it->as_bigint();
        result = result + opt_window_wnaf_exp(*vec_it, scalar_bigint, scalar_bigint.num_bits());
    }
    assert(scalar_it == scalar_end);

    return result;
}

template<typename T, typename FieldT>
T naive_plain_exp(typename std::vector<T>::const_iterator vec_start,
                  typename std::vector<T>::const_iterator vec_end,
                  typename std::vector<FieldT>::const_iterator scalar_start,
                  typename std::vector<FieldT>::const_iterator scalar_end)
{
    T result(T::zero());

    typename std::vector<T>::const_iterator vec_it;
    typename std::vector<FieldT>::const_iterator scalar_it;

    for (vec_it = vec_start, scalar_it = scalar_start; vec_it != vec_end; ++vec_it, ++scalar_it)
    {
        result = result + (*scalar_it) * (*vec_it);
    }
    assert(scalar_it == scalar_end);

    return result;
}

/**
 * Simultaneous 2^w-ary method,
 * Section 2.1 of Bodo Moller, "Algorithms for multi-exponentiation", SAC '01
 */
template<typename T, typename FieldT>
T simul_2w_multi_exp(typename std::vector<T>::const_iterator bases,
                     typename std::vector<T>::const_iterator bases_end,
                     typename std::vector<FieldT>::const_iterator exponents,
                     typename std::vector<FieldT>::const_iterator exponents_end,
                     size_t chunk_length)
{
    UNUSED(exponents_end);

    size_t length = bases_end - bases;
    size_t pair_count = (length + 1) / 2;

    std::vector<T> precomp(pair_count << (2 * chunk_length));

    // precomp[(i << (2 * chunk_length)) | (a << chunk_length) | b] contains
    // ag + bh, where g and f are elements of pair i
    // (so g = bases[2*i], h = bases[2*i + 1])
    for (size_t i = 0; i < pair_count; i++)
    {
        // first, assume b = 0
        // cover a = 0, a = 1 manually, then process the rest in pairs of
        // 2j and 2j+1 to cover even/odd cases separately and thus
        // benefit from fast doubling
        size_t pair_mask = i << (2 * chunk_length);
        precomp[pair_mask | 0 << chunk_length] = T::zero();
        precomp[pair_mask | 1 << chunk_length] = bases[2 * i];
        for (size_t j = 1; j < (1u << (chunk_length - 1)); j++)
        {
            // (2j) * g = 2 * (jg)
            precomp[pair_mask | (2*j) << chunk_length] =
                precomp[pair_mask | j << chunk_length].dbl();
            // (2j + 1) * g = (2j) * g  + g
            precomp[pair_mask | (2*j + 1) << chunk_length] =
                precomp[pair_mask | (2*j) << chunk_length] + bases[2 * i];
        }

        // now consider b != 0
        if (2*i + 1 == length)
        {
            // bases[2*i + 1] does not exist
            break;
        }

        for (size_t j = 0; j < (1u << (chunk_length - 1)); j++)
        {
            // calculating (2j)g + bh,
            // also in pairs to benefit from doubling
            // note that (2j)g + 0h is already calculated, but
            // (2j)g + 1h isn't
            precomp[pair_mask | (2*j) << chunk_length | 1] =
                precomp[pair_mask | (2*j) << chunk_length] + bases[2*i + 1];
            for (size_t k = 1; k < (1u << (chunk_length - 1)); k++)
            {
                // (2j)g + (2k)h = 2 * (jg + kh)
                precomp[pair_mask | (2*j) << chunk_length | (2*k)] =
                    precomp[pair_mask | j << chunk_length | k].dbl();
                // (2j)g + (2k + 1)h = (2j)g + (2k)h + h
                precomp[pair_mask | (2*j) << chunk_length | (2*k + 1)] =
                    precomp[pair_mask | (2*j) << chunk_length | (2*k)] +
                    bases[2 * i + 1];
            }

            // calculating (2j+1)g + bh,
            // cannot benefit from doubling
            // note that (2j+1)g + 0h is already calculated
            for (size_t b = 1; b < (1u << chunk_length); b++)
            {
                // (2j + 1) * g + b * h = (2j + 1) * g + (b - 1) * h + h
                precomp[pair_mask | (2*j + 1) << chunk_length | b] =
                    precomp[pair_mask | (2*j + 1) << chunk_length | (b - 1)] +
                    bases[2 * i + 1];
            }
        }
    }

    // now time to do the actual exponentiation
    const mp_size_t exp_num_limbs =
        std::remove_reference<decltype(*exponents)>::type::num_limbs;
    std::vector<bigint<exp_num_limbs> > bn_exponents(2 * pair_count);
    size_t num_bits = 0;

    for (size_t i = 0; i < length; i++)
    {
        bn_exponents[i] = exponents[i].as_bigint();
        num_bits = std::max(num_bits, bn_exponents[i].num_bits());
    }
    // note: if length is odd, bn_exponents[length + 1] == 0,
    // which is exactly what we'll want

    // note: bigint.test_bit(x) works safely even for x bigger than
    // bigint.max_bits()


    T result = T::zero();
    // chunk is unsigned, so the loop condition can't be (chunk >= 0)
    size_t chunk_count = (num_bits + chunk_length - 1) / chunk_length;
    for (size_t chunk = chunk_count - 1; chunk < chunk_count; chunk--)
    {
        for (size_t i = 0; i < pair_count; i++)
        {
            size_t entry = i << (2 * chunk_length);
            for (size_t bit = 0; bit < chunk_length; bit++)
            {
                if (bn_exponents[2*i].test_bit(chunk_length * chunk + bit))
                {
                    entry |= 1 << (chunk_length + bit);
                }

                if (bn_exponents[2*i + 1].test_bit(chunk_length * chunk + bit))
                {
                    entry |= 1 << bit;
                }
            }

            result = result + precomp[entry];
        }

        if (chunk != 0)
        {
            for (size_t i = 0; i < chunk_length; i++)
            {
                result = result.dbl();
            }
        }
    }

    return result;
}

/*
 * A special case of Pippenger's algorithm from Page 15 of
 * https://eprint.iacr.org/2012/549.pdf
 * When compiled with USE_MIXED_ADDITION, assumes input is
 * in special form.
 */
template<typename T, typename FieldT>
T multi_exp_djb(typename std::vector<T>::const_iterator bases,
                typename std::vector<T>::const_iterator bases_end,
                typename std::vector<FieldT>::const_iterator exponents,
                typename std::vector<FieldT>::const_iterator exponents_end,
                size_t c)
{
    UNUSED(exponents_end);
    size_t length = bases_end - bases;

    if (c == 0)
    {
        // empirically, this seems to be a decent estimate of the optimal value of c
        size_t log2_length = log2(length);
        c = log2_length - (log2_length / 3 - 2);
    }

    const mp_size_t exp_num_limbs =
        std::remove_reference<decltype(*exponents)>::type::num_limbs;
    std::vector<bigint<exp_num_limbs> > bn_exponents(length);
    size_t num_bits = 0;

    for (size_t i = 0; i < length; i++)
    {
        bn_exponents[i] = exponents[i].as_bigint();
        num_bits = std::max(num_bits, bn_exponents[i].num_bits());
    }

    size_t num_groups = (num_bits + c - 1) / c;

    T result;
    bool result_nonzero = false;

    for (size_t k = num_groups - 1; k <= num_groups; k--)
    {
        if (result_nonzero)
        {
            for (size_t i = 0; i < c; i++)
            {
                result = result.dbl();
            }
        }

        std::vector<T> buckets(1 << c);
        std::vector<bool> bucket_nonzero(1 << c);

        for (size_t i = 0; i < length; i++)
        {
            size_t id = 0;
            for (size_t j = 0; j < c; j++)
            {
                if (bn_exponents[i].test_bit(k*c + j))
                {
                    id |= 1 << j;
                }
            }

            if (id == 0)
            {
                continue;
            }

            if (bucket_nonzero[id])
            {
#ifdef USE_MIXED_ADDITION
                buckets[id] = buckets[id].mixed_add(bases[i]);
#else
                buckets[id] = buckets[id] + bases[i];
#endif
            }
            else
            {
                buckets[id] = bases[i];
                bucket_nonzero[id] = true;
            }
        }

#ifdef USE_MIXED_ADDITION
        batch_to_special(buckets);
#endif

        T running_sum;
        bool running_sum_nonzero = false;

        for (size_t i = (1u << c) - 1; i > 0; i--)
        {
            if (bucket_nonzero[i])
            {
                if (running_sum_nonzero)
                {
#ifdef USE_MIXED_ADDITION
                    running_sum = running_sum.mixed_add(buckets[i]);
#else
                    running_sum = running_sum + buckets[i];
#endif
                }
                else
                {
                    running_sum = buckets[i];
                    running_sum_nonzero = true;
                }
            }

            if (running_sum_nonzero)
            {
                if (result_nonzero)
                {
                    result = result + running_sum;
                }
                else
                {
                    result = running_sum;
                    result_nonzero = true;
                }
            }
        }
    }

    return result;
}

template<typename T, typename FieldT>
T multi_exp_inner_rivest(typename std::vector<T>::const_iterator vec_start,
                         typename std::vector<T>::const_iterator vec_end,
                         typename std::vector<FieldT>::const_iterator scalar_start,
                         typename std::vector<FieldT>::const_iterator scalar_end)
{
    UNUSED(scalar_end);
    size_t length = vec_end - vec_start;

    if (length == 0)
    {
        return T::zero();
    }

    if (length == 1)
    {
        return (*scalar_start)*(*vec_start);
    }

    // first, we carefully count & group exponents into buckets
    // so that they are stored in exponents[] grouped by .num_bits()
    // s.t. the exponents with .num_bits() = k
    // are stored in exponents[bucket_start[k] .. bucket_start[k + 1] - 1]
    const mp_size_t exp_num_limbs =
        std::remove_reference<decltype(*scalar_start)>::type::num_limbs;
    std::vector<bigint<exp_num_limbs> > tmp_exponents(length);
    std::vector<size_t> bucket_start(exp_num_limbs * GMP_NUMB_BITS);
    size_t num_bits = 0;

    for (size_t i = 0; i < length; i++)
    {
        tmp_exponents[i] = scalar_start[i].as_bigint();
        size_t num_bits_here = tmp_exponents[i].num_bits();
        num_bits = std::max(num_bits, num_bits_here);
        bucket_start[num_bits_here]++;
    }

    for (size_t i = 1; i < num_bits + 2; i++)
    {
        bucket_start[i] += bucket_start[i - 1];
    }

    std::vector<T> bases(length);
    std::vector<bigint<exp_num_limbs> > exponents(length);
    for (size_t i = 0; i < length; i++)
    {
        size_t index = --bucket_start[tmp_exponents[i].num_bits()];
        exponents[index] = tmp_exponents[i];
        bases[index] = vec_start[i];
    }

    // now we do a weaker version of Bos and Coster:
    // take two terms ag + bh from the biggest bucket
    // and replace them with (a-b)g + b(g+h)
    // until all exponents but one are 0

    // we will keep num_bits updated (i.e. it's always going to be
    // the highest i for which bucket_start[i] != bucket_start[i + 1])

    T result = T::zero();

    while (bucket_start[1] < length - 1)
    {
        // take the last element of the last bucket,
        // and also the one before it (possibly from a smaller bucket)
        size_t i = bucket_start[num_bits + 1] - 1;
        size_t j = i - 1;

        // it could be that exponents[j] > exponents[i], but only
        // if they come from the same bucket.
        // that would be very inconvenient, so we test it & fix if necessary.
        if ((bucket_start[num_bits] <= j) &&
            (mpn_cmp(exponents[j].data, exponents[i].data, exp_num_limbs) > 0))
        {
            std::swap(exponents[i], exponents[j]);
            std::swap(bases[i], bases[j]);
        }

        size_t j_bits = exponents[j].num_bits();
        size_t limit = (num_bits - j_bits >= 20 ? 20 : num_bits - j_bits);

        if (j_bits < 1ul << limit)
        {
            /*
              In this case, exponentiating to the power of a is cheaper than
              subtracting b from a multiple times, so let's do it directly
              */
            result = result.add(opt_window_wnaf_exp(bases[i], exponents[i], num_bits));
            exponents[i].clear();
        }
        else
        {
            bases[j] = bases[j].add(bases[i]);
            // exponents[i] -= exponents[j];
            mpn_sub_n(exponents[i].data, exponents[i].data, exponents[j].data,
                exp_num_limbs);
        }


        // i might now belong in a different bucket
        size_t i_bucket = exponents[i].num_bits();
        size_t i_cur_bucket = num_bits;
        while (i_cur_bucket != i_bucket)
        {
            // move i to previous bucket by moving it to start of current one
            // and then advance the bucket_start pointer forwards
            if (i != bucket_start[i_cur_bucket])
            {
                std::swap(exponents[i], exponents[bucket_start[i_cur_bucket]]);
                std::swap(bases[i], bases[bucket_start[i_cur_bucket]]);
                i = bucket_start[i_cur_bucket];
            }
            bucket_start[i_cur_bucket]++;
            i_cur_bucket--;
        }

        // update num_bits in case i was the only thing in its bucket
        // and had its num_bits decreased
        while (bucket_start[num_bits] == bucket_start[num_bits + 1])
        {
            num_bits--;
        }
    }

    return result.add(opt_window_wnaf_exp(bases[bucket_start[num_bits]],
        exponents[bucket_start[num_bits]],
        num_bits));
}

/*
  The multi-exponentiation algorithm below is a variant of the Bos-Coster algorithm
  [Bos and Coster, "Addition chain heuristics", CRYPTO '89].
  The implementation uses suggestions from
  [Bernstein, Duif, Lange, Schwabe, and Yang, "High-speed high-security signatures", CHES '11].
*/
template<typename T, typename FieldT>
T multi_exp_inner(typename std::vector<T>::const_iterator vec_start,
                  typename std::vector<T>::const_iterator vec_end,
                  typename std::vector<FieldT>::const_iterator scalar_start,
                  typename std::vector<FieldT>::const_iterator scalar_end)
{
    const mp_size_t n = std::remove_reference<decltype(*scalar_start)>::type::num_limbs;

    if (vec_start == vec_end)
    {
        return T::zero();
    }

    if (vec_start + 1 == vec_end)
    {
        return (*scalar_start)*(*vec_start);
    }

    std::vector<ordered_exponent<n> > opt_q;
    const size_t vec_len = scalar_end - scalar_start;
    const size_t odd_vec_len = (vec_len % 2 == 1 ? vec_len : vec_len + 1);
    opt_q.reserve(odd_vec_len);
    std::vector<T> g;
    g.reserve(odd_vec_len);

    typename std::vector<T>::const_iterator vec_it;
    typename std::vector<FieldT>::const_iterator scalar_it;
    size_t i;
    for (i=0, vec_it = vec_start, scalar_it = scalar_start; vec_it != vec_end; ++vec_it, ++scalar_it, ++i)
    {
        g.emplace_back(*vec_it);

        opt_q.emplace_back(ordered_exponent<n>(i, scalar_it->as_bigint()));
    }
    std::make_heap(opt_q.begin(),opt_q.end());
    assert(scalar_it == scalar_end);

    if (vec_len != odd_vec_len)
    {
        g.emplace_back(T::zero());
        opt_q.emplace_back(ordered_exponent<n>(odd_vec_len - 1, bigint<n>(0ul)));
    }
    assert(g.size() % 2 == 1);
    assert(opt_q.size() == g.size());

    T opt_result = T::zero();

    while (true)
    {
        ordered_exponent<n> &a = opt_q[0];
        ordered_exponent<n> &b = (opt_q[1] < opt_q[2] ? opt_q[2] : opt_q[1]);

        const size_t abits = a.r.num_bits();

        if (b.r.is_zero())
        {
            // opt_result = opt_result + (a.r * g[a.idx]);
            opt_result = opt_result + opt_window_wnaf_exp(g[a.idx], a.r, abits);
            break;
        }

        const size_t bbits = b.r.num_bits();
        const size_t limit = (abits-bbits >= 20 ? 20 : abits-bbits);

        if (bbits < 1ul<<limit)
        {
            /*
              In this case, exponentiating to the power of a is cheaper than
              subtracting b from a multiple times, so let's do it directly
            */
            // opt_result = opt_result + (a.r * g[a.idx]);
            opt_result = opt_result + opt_window_wnaf_exp(g[a.idx], a.r, abits);
#ifdef DEBUG
            printf("Skipping the following pair (%zu bit number vs %zu bit):\n", abits, bbits);
            a.r.print();
            b.r.print();
#endif
            a.r.clear();
        }
        else
        {
            // x A + y B => (x-y) A + y (B+A)
            mpn_sub_n(a.r.data, a.r.data, b.r.data, n);
            g[b.idx] = g[b.idx] + g[a.idx];
        }

        // regardless of whether a was cleared or subtracted from we push it down, then take back up

        /* heapify A down */
        size_t a_pos = 0;
        while (2*a_pos + 2< odd_vec_len)
        {
            // this is a max-heap so to maintain a heap property we swap with the largest of the two
            if (opt_q[2*a_pos+1] < opt_q[2*a_pos+2])
            {
                std::swap(opt_q[a_pos], opt_q[2*a_pos+2]);
                a_pos = 2*a_pos+2;
            }
            else
            {
                std::swap(opt_q[a_pos], opt_q[2*a_pos+1]);
                a_pos = 2*a_pos+1;
            }
        }

        /* now heapify A up appropriate amount of times */
        while (a_pos > 0 && opt_q[(a_pos-1)/2] < opt_q[a_pos])
        {
            std::swap(opt_q[a_pos], opt_q[(a_pos-1)/2]);
            a_pos = (a_pos-1) / 2;
        }
    }

    return opt_result;
}

template<typename T, typename FieldT>
T multi_exp(typename std::vector<T>::const_iterator vec_start,
            typename std::vector<T>::const_iterator vec_end,
            typename std::vector<FieldT>::const_iterator scalar_start,
            typename std::vector<FieldT>::const_iterator scalar_end,
            const size_t chunks,
            const bool use_multiexp)
{
    const size_t total = vec_end - vec_start;
    if (total < chunks)
    {
        return naive_exp<T, FieldT>(vec_start, vec_end, scalar_start, scalar_end);
    }

    const size_t one = total/chunks;

    std::vector<T> partial(chunks, T::zero());

    if (use_multiexp)
    {
#ifdef MULTICORE
#pragma omp parallel for
#endif
        for (size_t i = 0; i < chunks; ++i)
        {
            partial[i] = multi_exp_inner<T, FieldT>(vec_start + i*one,
                                                    (i == chunks-1 ? vec_end : vec_start + (i+1)*one),
                                                    scalar_start + i*one,
                                                    (i == chunks-1 ? scalar_end : scalar_start + (i+1)*one));
        }
    }
    else
    {
#ifdef MULTICORE
#pragma omp parallel for
#endif
        for (size_t i = 0; i < chunks; ++i)
        {
            partial[i] = naive_exp<T, FieldT>(vec_start + i*one,
                                              (i == chunks-1 ? vec_end : vec_start + (i+1)*one),
                                              scalar_start + i*one,
                                              (i == chunks-1 ? scalar_end : scalar_start + (i+1)*one));
        }
    }

    T final = T::zero();

    for (size_t i = 0; i < chunks; ++i)
    {
        final = final + partial[i];
    }

    return final;
}

template<typename T, typename FieldT>
T multi_exp_with_mixed_addition(typename std::vector<T>::const_iterator vec_start,
                                typename std::vector<T>::const_iterator vec_end,
                                typename std::vector<FieldT>::const_iterator scalar_start,
                                typename std::vector<FieldT>::const_iterator scalar_end,
                                const size_t chunks,
                                const bool use_multiexp)
{
    assert(std::distance(vec_start, vec_end) == std::distance(scalar_start, scalar_end));
    enter_block("Process scalar vector");
    auto value_it = vec_start;
    auto scalar_it = scalar_start;

    const FieldT zero = FieldT::zero();
    const FieldT one = FieldT::one();
    std::vector<FieldT> p;
    std::vector<T> g;

    T acc = T::zero();

    size_t num_skip = 0;
    size_t num_add = 0;
    size_t num_other = 0;

    for (; scalar_it != scalar_end; ++scalar_it, ++value_it)
    {
        if (*scalar_it == zero)
        {
            // do nothing
            ++num_skip;
        }
        else if (*scalar_it == one)
        {
#ifdef USE_MIXED_ADDITION
            acc = acc.mixed_add(*value_it);
#else
            acc = acc + (*value_it);
#endif
            ++num_add;
        }
        else
        {
            p.emplace_back(*scalar_it);
            g.emplace_back(*value_it);
            ++num_other;
        }
    }
    print_indent(); printf("* Elements of w skipped: %zu (%0.2f%%)\n", num_skip, 100.*num_skip/(num_skip+num_add+num_other));
    print_indent(); printf("* Elements of w processed with special addition: %zu (%0.2f%%)\n", num_add, 100.*num_add/(num_skip+num_add+num_other));
    print_indent(); printf("* Elements of w remaining: %zu (%0.2f%%)\n", num_other, 100.*num_other/(num_skip+num_add+num_other));

    leave_block("Process scalar vector");

    return acc + multi_exp<T, FieldT>(g.begin(), g.end(), p.begin(), p.end(), chunks, use_multiexp);
}

template<typename T>
size_t get_exp_window_size(const size_t num_scalars)
{
    if (T::fixed_base_exp_window_table.empty())
    {
#ifdef LOWMEM
        return 14;
#else
        return 17;
#endif
    }
    size_t window = 1;
    for (long i = T::fixed_base_exp_window_table.size()-1; i >= 0; --i)
    {
#ifdef DEBUG
        if (!inhibit_profiling_info)
        {
            printf("%ld %zu %zu\n", i, num_scalars, T::fixed_base_exp_window_table[i]);
        }
#endif
        if (T::fixed_base_exp_window_table[i] != 0 && num_scalars >= T::fixed_base_exp_window_table[i])
        {
            window = i+1;
            break;
        }
    }

    if (!inhibit_profiling_info)
    {
        print_indent(); printf("Choosing window size %zu for %zu elements\n", window, num_scalars);
    }

#ifdef LOWMEM
    window = std::min((size_t)14, window);
#endif
    return window;
}

template<typename T>
window_table<T> get_window_table(const size_t scalar_size,
                                 const size_t window,
                                 const T &g)
{
    const size_t in_window = 1ul<<window;
    const size_t outerc = (scalar_size+window-1)/window;
    const size_t last_in_window = 1ul<<(scalar_size - (outerc-1)*window);
#ifdef DEBUG
    if (!inhibit_profiling_info)
    {
        print_indent(); printf("* scalar_size=%zu; window=%zu; in_window=%zu; outerc=%zu\n", scalar_size, window, in_window, outerc);
    }
#endif

    window_table<T> powers_of_g(outerc, std::vector<T>(in_window, T::zero()));

    T gouter = g;

    for (size_t outer = 0; outer < outerc; ++outer)
    {
        T ginner = T::zero();
        size_t cur_in_window = outer == outerc-1 ? last_in_window : in_window;
        for (size_t inner = 0; inner < cur_in_window; ++inner)
        {
            powers_of_g[outer][inner] = ginner;
            ginner = ginner + gouter;
        }

        for (size_t i = 0; i < window; ++i)
        {
            gouter = gouter + gouter;
        }
    }

    return powers_of_g;
}

template<typename T, typename FieldT>
T windowed_exp(const size_t scalar_size,
               const size_t window,
               const window_table<T> &powers_of_g,
               const FieldT &pow)
{
    const size_t outerc = (scalar_size+window-1)/window;
    const bigint<FieldT::num_limbs> pow_val = pow.as_bigint();

    /* exp */
    T res = powers_of_g[0][0];

    for (size_t outer = 0; outer < outerc; ++outer)
    {
        size_t inner = 0;
        for (size_t i = 0; i < window; ++i)
        {
            if (pow_val.test_bit(outer*window + i))
            {
                inner |= 1u << i;
            }
        }

        res = res + powers_of_g[outer][inner];
    }

    return res;
}

template<typename T, typename FieldT>
std::vector<T> batch_exp(const size_t scalar_size,
                         const size_t window,
                         const window_table<T> &table,
                         const std::vector<FieldT> &v)
{
    if (!inhibit_profiling_info)
    {
        print_indent();
    }
    std::vector<T> res(v.size(), table[0][0]);

#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < v.size(); ++i)
    {
        res[i] = windowed_exp(scalar_size, window, table, v[i]);

        if (!inhibit_profiling_info && (i % 10000 == 0))
        {
            printf(".");
            fflush(stdout);
        }
    }

    if (!inhibit_profiling_info)
    {
        printf(" DONE!\n");
    }

    return res;
}

template<typename T, typename FieldT>
std::vector<T> batch_exp_with_coeff(const size_t scalar_size,
                                    const size_t window,
                                    const window_table<T> &table,
                                    const FieldT &coeff,
                                    const std::vector<FieldT> &v)
{
    if (!inhibit_profiling_info)
    {
        print_indent();
    }
    std::vector<T> res(v.size(), table[0][0]);

#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < v.size(); ++i)
    {
        res[i] = windowed_exp(scalar_size, window, table, coeff * v[i]);

        if (!inhibit_profiling_info && (i % 10000 == 0))
        {
            printf(".");
            fflush(stdout);
        }
    }

    if (!inhibit_profiling_info)
    {
        printf(" DONE!\n");
    }

    return res;
}

template<typename T>
void batch_to_special(std::vector<T> &vec)
{
    enter_block("Batch-convert elements to special form");

    std::vector<T> non_zero_vec;
    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (!vec[i].is_zero())
        {
            non_zero_vec.emplace_back(vec[i]);
        }
    }

    batch_to_special_all_non_zeros<T>(non_zero_vec);
    auto it = non_zero_vec.begin();
    T zero_special = T::zero();
    zero_special.to_special();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (!vec[i].is_zero())
        {
            vec[i] = *it;
            ++it;
        }
        else
        {
            vec[i] = zero_special;
        }
    }
    leave_block("Batch-convert elements to special form");
}

} // libff

#endif // MULTIEXP_TCC_
