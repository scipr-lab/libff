/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <gtest/gtest.h>
#include "libff/algebra/field_utils/bigint.hpp"

using namespace libff;

/** Returns a random element of FieldT that is not zero. */
template<mp_size_t n>
bigint<n> random_element_nonzero()
{
    bigint<n> x;
    x.randomize();
    while (x.is_zero())
        x.randomize();
    return x;
}

TEST(FieldUtilsTest, BigintTest)
{
    bigint<3> zero = bigint<3>("0");
    bigint<3> one = bigint<3>("1");
    bigint<3> x = bigint<3>("987654567895678909876876545678909876543456");
    bigint<3> y = bigint<3>("324531232345676543272920293863628304859020");
    EXPECT_EQ(x - y, bigint<3>("663123335550002366603956251815281571684436"));
    EXPECT_EQ(x / 100000, bigint<3>("9876545678956789098768765456789098765"));
    EXPECT_EQ(x % 100000, bigint<3>("43456"));

    bigint<1> z = bigint<1>("1234");
    EXPECT_EQ(z - z, bigint<1>("0"));
    EXPECT_EQ(z.power<4>(12), bigint<4>("12467572902176589255564000298710470656"));

    x.randomize();
    y = random_element_nonzero<3>();
    EXPECT_EQ(x - y, zero - (y - x));
    EXPECT_EQ(x - y - x, zero - y);
    EXPECT_EQ(x / 2 / 3 / 25, x / 5 / 30);
    EXPECT_TRUE(x % 2 == zero || x % 2 == one);
}
