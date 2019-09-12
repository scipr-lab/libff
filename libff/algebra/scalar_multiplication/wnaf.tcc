/** @file
 *****************************************************************************

 Implementation of interfaces for wNAF ("weighted Non-Adjacent Form") exponentiation routines.

 See wnaf.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef WNAF_TCC_
#define WNAF_TCC_

#include <gmp.h>

namespace libff {

template<mp_size_t n>
void update_wnaf(std::vector<long> &wnaf, const size_t window_size, const bigint<n> &scalar)
{
    const size_t length = scalar.max_bits(); // upper bound
    wnaf.resize(length+1);
    bigint<n> c = scalar;
    long j = 0;
    while (!c.is_zero())
    {
        long u;
        if ((c.data[0] & 1) == 1)
        {
            u = c.data[0] % (1u << (window_size+1));
            if (u > (1 << window_size))
            {
                u = u - (1 << (window_size+1));
            }

            if (u > 0)
            {
                mpn_sub_1(c.data, c.data, n, u);
            }
            else
            {
                mpn_add_1(c.data, c.data, n, -u);
            }
        }
        else
        {
            u = 0;
        }
        wnaf[j] = u;
        ++j;

        mpn_rshift(c.data, c.data, n, 1); // c = c/2
    }

    wnaf.resize(j);
}

template<mp_size_t n>
std::vector<long> find_wnaf(const size_t window_size, const bigint<n> &scalar)
{
    std::vector<long> res;
    update_wnaf(res, window_size, scalar);
    return res;
}

template<typename T>
size_t wnaf_opt_window_size(const size_t scalar_bits)
{
    for (long i = T::wnaf_window_table.size() - 1; i >= 0; --i)
    {
        if (scalar_bits >= T::wnaf_window_table[i])
        {
            return i+1;
        }
    }

    return 0;
}

template<typename T>
T fixed_window_wnaf_exp(const size_t window_size, const T &base, const std::vector<long> &naf)
{
    std::vector<T> table(1ul<<(window_size-1));
    T tmp = base;
    T dbl = base.dbl();
    for (size_t i = 0; i < 1ul<<(window_size-1); ++i)
    {
        table[i] = tmp;
        tmp = tmp + dbl;
    }

    T res = T::zero();
    bool found_nonzero = false;
    for (long i = naf.size()-1; i >= 0; --i)
    {
        if (found_nonzero)
        {
            res = res.dbl();
        }

        if (naf[i] != 0)
        {
            found_nonzero = true;
            if (naf[i] > 0)
            {
                res = res + table[naf[i]/2];
            }
            else
            {
                res = res - table[(-naf[i])/2];
            }
        }
    }

    return res;
}

template<typename T, mp_size_t n>
T fixed_window_wnaf_exp(const size_t window_size, const T &base, const bigint<n> &scalar)
{
    std::vector<long> naf = find_wnaf(window_size, scalar);
    return fixed_window_wnaf_exp<T>(window_size, base, naf);
}

template<typename T, mp_size_t n>
T opt_window_wnaf_exp(const T &base, const bigint<n> &scalar, const size_t scalar_bits)
{
    const size_t best = wnaf_opt_window_size<T>(scalar_bits);
    if (best > 0)
    {
        return fixed_window_wnaf_exp(best, base, scalar);
    }
    else
    {
        return scalar * base;
    }
}

} // libff

#endif // WNAF_TCC_
