/**
 *****************************************************************************
 Basic tests for some of the field utils in this directory, mainly bigint
 and power.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include "libff/algebra/field_utils/bigint.hpp"
#include "libff/algebra/field_utils/algorithms.hpp"
#include "libff/algebra/fields/binary/gf64.hpp"
#include <gtest/gtest.h>

using namespace libff;

template<typename FieldT>
FieldT power_naive(const FieldT &base, const std::size_t exponent)
{
    FieldT result = FieldT::one();

    for (std::size_t i = 1; i <= exponent; ++i)
    {
        result *= base;
    }

    return result;
}


TEST(ExponentiationTest, SimpleTest) {
    typedef gf64 FieldT;

    const unsigned long max_power = 3000;
    FieldT X = FieldT::random_element();

    FieldT X_i = FieldT::one();
    for (unsigned long i = 0 ; i < max_power; ++i)
    {
        const FieldT X_i_naive = power_naive<FieldT>(X, i);
        const FieldT X_i_square_and_multiply_ul = power<FieldT>(X, i);
        const FieldT X_i_square_and_multiply_ull = power<FieldT>(X, (unsigned long long) i);

        EXPECT_EQ(X_i, X_i_naive);
        EXPECT_EQ(X_i, X_i_square_and_multiply_ul);
        EXPECT_EQ(X_i, X_i_square_and_multiply_ull);

        X_i *= X;
    }
}

TEST(FieldUtilsTest, BigintTest)
{
    bigint<3> zero = bigint<3>("0");
    bigint<3> one = bigint<3>("1");
    bigint<3> x = bigint<3>("987654567895678909876876545678909876543456");
    bigint<3> y = bigint<3>("324531232345676543272920293863628304859020");

    bigint<3> zero2;
    bigint<3> z = bigint<3>("987654567895678909876876545678909876543456");
    bigint<3> w = bigint<3>("987654567895678909876876545678909876543451");

    EXPECT_EQ(zero, zero2);
    EXPECT_EQ(one, bigint<3>::one());
    EXPECT_EQ(x, x);
    EXPECT_EQ(x, z);

    EXPECT_NE(zero, one);
    EXPECT_NE(x, one);
    EXPECT_NE(x, zero);
    EXPECT_NE(x, y);

    EXPECT_FALSE(x.is_zero());
    EXPECT_FALSE(one.is_zero());
    EXPECT_TRUE(zero.is_zero());

    EXPECT_TRUE(x.is_even());
    EXPECT_TRUE(y.is_even());
    EXPECT_TRUE(zero.is_even());

    EXPECT_FALSE(one.is_even());
    EXPECT_FALSE(w.is_even());

    EXPECT_LT(zero, one);
    EXPECT_LT(zero, x);
    EXPECT_LT(zero, y);
    EXPECT_LT(one, x);
    EXPECT_LT(y, x);
    EXPECT_LT(w, x);
    EXPECT_LT(y, w);
    EXPECT_FALSE(x < y);

    x.print();
    x.print_hex();

    x.randomize();
    EXPECT_EQ(x, x);
    x.clear();
    EXPECT_EQ(x, zero);
}
